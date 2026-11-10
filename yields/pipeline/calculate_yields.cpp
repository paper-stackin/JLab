YieldsOutput CalculateYields(
    const std::vector<std::vector<TH1F*>>& hist_signal,
    const std::vector<std::vector<TH1F*>>& hist_bg,
    const std::vector<std::vector<FitResult>>& fit_results,
    const std::vector<std::vector<double>>& integral_background, 
    const YieldsConfig& cfg)
{
    YieldsOutput res;
    res.w_vector.resize(cfg.q_max + 1);
    res.w_vector_error.resize(cfg.q_max + 1);
    res.integral_signal_excut.resize(cfg.q_max + 1);
    res.integral_signal_excut_error.resize(cfg.q_max + 1);
    res.integral_signal.resize(cfg.q_max + 1);
    res.integral_error.resize(cfg.q_max + 1);

    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
        double Q2_min_val = cfg.Q2_edges[q];
        double Q2_max_val = cfg.Q2_edges[q + 1];

        double del_q2 = Q2_max_val - Q2_min_val;

        double q2_mid_val = (Q2_max_val + Q2_min_val) / 2.0;

        int w_brink = cfg.WBinLimit(q);

        for (int w = 0; w < w_brink; ++w)
        {
            //W
            double w_min_val = cfg.W_min + w * cfg.W_bin_size;
            double w_max_val = cfg.W_min + (w + 1) * cfg.W_bin_size;
            double w_mid_val = (w_min_val + w_max_val) / 2.0;

            res.w_vector[q].push_back(w_mid_val);
            res.w_vector_error[q].push_back(cfg.W_bin_size);

            //Flux
            double E_prime =
                cfg.E_beam -
                (w_mid_val * w_mid_val
                 - cfg.M_proton * cfg.M_proton
                 + q2_mid_val)
                / (2.0 * cfg.M_proton);

            double nu = cfg.E_beam - E_prime;

            double sin2 =
                q2_mid_val /
                (4.0 * cfg.E_beam * E_prime);

            double cos2 = 1.0 - sin2;
            double tg2 = sin2 / cos2;

            double eps =
                1.0 /
                (
                    1.0 +
                    2.0 * (1.0 + nu * nu / q2_mid_val) * tg2
                );

            double flux =
                cfg.alpha / (4.0 * 3.1415)
                * 1.0 /
                  (cfg.M_proton * cfg.M_proton * cfg.E_beam * cfg.E_beam)
                * (
                    w_mid_val *
                    (w_mid_val * w_mid_val - cfg.M_proton * cfg.M_proton)
                  )
                / (
                    (1.0 - eps) * q2_mid_val
                  );

            double delta_flux =
                sqrt(
                    pow(
                        cfg.alpha / (4.0 * 3.1415)
                        * 1.0 /
                          (cfg.M_proton * cfg.M_proton * cfg.E_beam * cfg.E_beam)
                        * (
                            (3.0 * w_mid_val * w_mid_val
                             - cfg.M_proton * cfg.M_proton)
                            / (2.0 * q2_mid_val)
                            - 1.0
                          ),
                        2
                    ) * pow(cfg.W_bin_size, 2)
                    +
                    pow(
                        cfg.alpha / (4.0 * 3.1415)
                        * 1.0 /
                          (cfg.M_proton * cfg.M_proton * cfg.E_beam * cfg.E_beam)
                        * (
                            (w_mid_val * w_mid_val
                             - cfg.M_proton * cfg.M_proton)
                            / (2.0 * q2_mid_val)
                          ),
                        2
                    ) * pow(del_q2, 2)
                );

            // N_ex_cut from signal histogram
            double N_ex_cut = hist_signal[q][w]->Integral(
                cfg.signal_integral_bin_min,
                cfg.signal_integral_bin_max
            );

            // Yield without background subtraction

            double yield_ex_cut = N_ex_cut / del_q2 / cfg.W_bin_size / flux;
            double yield_ex_cut_error = yield_ex_cut * sqrt(
                pow(sqrt(N_ex_cut) / N_ex_cut, 2) + pow(delta_flux / flux, 2)
            );

            res.integral_signal_excut[q].push_back(yield_ex_cut);
            res.integral_signal_excut_error[q].push_back(yield_ex_cut_error);

            // Scale coefficient
            double scale_coef =
                (hist_signal[q][w] == hist_bg[q][w])
                    ? 1.0
                    : hist_signal[q][w]->Integral()
                      / hist_bg[q][w]->Integral();

            // Background
            int w_reper = cfg.w_fit_min;

            if (q >= cfg.w_reper_q45_start)
                w_reper = cfg.w_reper_q45;

            double N_bg = 0.0;

            if (w < w_reper)
            {
                // No background subtraction
                N_bg = 0.0;
            }
            else if (w < cfg.bg_fit_start)
            {
                // Interpolated background
                N_bg =
                    integral_background[q][w - w_reper]
                    * scale_coef;
            }
            else
            {
                // Background from fit
                const FitResult& result = fit_results[q][w];

                TF1* bg = new TF1(
                    Form("background_q%d_w%d", q, w),
                    background,
                    result.xmin,
                    result.xmax,
                    3
                );

                bg->SetParameters(
                    result.bg_amp,
                    result.bg_mean,
                    result.bg_sigma
                );

                N_bg =
                    bg->Integral(
                        cfg.bg_integral_xmin,
                        cfg.bg_integral_xmax
                    )
                    * cfg.bin_width_factor
                    * scale_coef;

                delete bg;
            }

            // Signal events after background subtraction
            double N = N_ex_cut - N_bg;

            // Final yield
            double yield = N / del_q2 / cfg.W_bin_size / flux;

            res.integral_signal[q].push_back(
                yield
            );

            res.integral_error[q].push_back(
                yield *
                sqrt(
                    pow(sqrt(N) / N, 2)
                    +
                    pow(delta_flux / flux, 2)
                )
            );
        }
    }

    return res;
}