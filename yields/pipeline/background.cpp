std::vector<std::vector<FitResult>> FitBackground(
    const std::vector<std::vector<TH1F*>>& hist_bg,
    const YieldsConfig& cfg)
{
    std::vector<std::vector<FitResult>> fit_results(
        cfg.n_Q2_bins,
        std::vector<FitResult>(cfg.n_W_bins)
    );

    for (int q = 0; q < cfg.n_Q2_bins; ++q)
    {
        int w_brink = cfg.WBinLimit(q);

        for (int w = cfg.bg_fit_start; w < w_brink; ++w)
        {
            double h_raw =
                hist_bg[q][w]->GetMaximum();

            double current_xmax;

            if (w < cfg.xmax_mid_start)
                current_xmax = cfg.xmax_low;
            else if (w < cfg.xmax_high_start)
                current_xmax = cfg.xmax_mid;
            else
                current_xmax = cfg.xmax_high;

            TF1* fitfunc = new TF1(
                Form("fitfunc_q%d_w%d", q, w),
                combined,
                cfg.xmin,
                current_xmax,
                7
            );

            double b_raw =
                hist_bg[q][w]->GetBinContent(cfg.bg_ampl_bin)
                * cfg.bg_ampl_factor;

            fitfunc->SetParameters(
                b_raw,
                cfg.bg_mean,
                cfg.bg_sigma,
                h_raw,
                cfg.signal_mean,
                cfg.signal_sigma,
                cfg.signal_tail
            );

            hist_bg[q][w]->Fit(fitfunc, "RQ");

            FitResult& result = fit_results[q][w];

            result.fitted = true;

            result.bg_amp = fitfunc->GetParameter(0);
            result.bg_mean = fitfunc->GetParameter(1);
            result.bg_sigma = fitfunc->GetParameter(2);

            result.bg_amp_error = fitfunc->GetParError(0);
            result.bg_mean_error = fitfunc->GetParError(1);
            result.bg_sigma_error = fitfunc->GetParError(2);

            result.signal_amp = fitfunc->GetParameter(3);
            result.signal_mean = fitfunc->GetParameter(4);
            result.signal_sigma = fitfunc->GetParameter(5);
            result.signal_tail = fitfunc->GetParameter(6);

            result.signal_amp_error = fitfunc->GetParError(3);
            result.signal_mean_error = fitfunc->GetParError(4);
            result.signal_sigma_error = fitfunc->GetParError(5);
            result.signal_tail_error = fitfunc->GetParError(6);

            result.xmin = cfg.xmin;
            result.xmax = current_xmax;

            delete fitfunc;
        }
    }

    return fit_results;
}

ExtrapolationOutput ExtrapolateBackground(
    const std::vector<std::vector<FitResult>>& fit_results,
    const YieldsConfig& cfg)
{
    ExtrapolationOutput output;
    output.integral_background.resize(cfg.n_Q2_bins);
    output.integral_background_error.resize(cfg.n_Q2_bins);
    output.w_bg_fit.resize(cfg.n_Q2_bins);
    output.results.resize(cfg.n_Q2_bins);

    for (int q = 0; q < cfg.n_Q2_bins; ++q)
    {
        int w_reper = cfg.w_fit_min;

        if (q >= cfg.w_reper_q45_start)
            w_reper = cfg.w_reper_q45;

        for (int w = w_reper; w <= cfg.w_fit_max; ++w)
        {
            int N_bg = 0;

            if (w >= cfg.bg_fit_start)
            {
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

                N_bg = bg->Integral(
                    cfg.bg_integral_xmin,
                    cfg.bg_integral_xmax
                ) * cfg.bin_width_factor;

                delete bg;
            }

            output.integral_background[q].push_back(N_bg);
            output.integral_background_error[q].push_back(sqrt(N_bg));
            output.w_bg_fit[q].push_back(w);
        }

        // Graph for fitting
        TGraphErrors* graph_fit_bg = new TGraphErrors(
            output.integral_background[q].size(),
            output.w_bg_fit[q].data(),
            output.integral_background[q].data(),
            nullptr,
            output.integral_background_error[q].data()
        );

        // Quadratic function
        TF1* fitFunction = new TF1(
            Form("fitFunction_q%d", q),
            Form(
                "[0]*(x-%d)*(x-%d)+[1]*(x-%d)",
                w_reper,
                w_reper,
                w_reper
            ),
            w_reper - 0.1,
            cfg.w_fit_max
        );

        // Initial parameters
        fitFunction->SetParameter(0, output.integral_background[q][6]);
        fitFunction->SetParameter(1, 1);

        // Fit
        graph_fit_bg->Fit(fitFunction, "RQ");

        // Save extrapolation parameters
        output.results[q].A = fitFunction->GetParameter(0);
        output.results[q].B = fitFunction->GetParameter(1);
        output.results[q].w_reper = w_reper;

        // Extrapolation
        for (int w = 0; w < cfg.bg_fit_start - w_reper; ++w)
        {
            output.integral_background[q][w] =
                fitFunction->Eval(
                    w_reper + w + 0.01
                );
        }

        delete fitFunction;
        delete graph_fit_bg;
    }

    return output;
}