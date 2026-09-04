#include "yields_config.h"

Double_t background(Double_t *x, Double_t *par) 
{

	Double_t bg_ampl = par[0];
    Double_t mean = par[1];
    Double_t sigma = par[2];

    Double_t logX = log(x[0] + 0.2); // Добавляем 0.1 для избежания отрицательных значений
    Double_t logNorm = bg_ampl / ( (x[0] + 0.2) * sigma * sqrt(2 * M_PI)) * exp(-0.5 * pow((logX - mean) / sigma, 2));

    return logNorm;
}

// Define the peak function
Double_t peak(Double_t *x, Double_t *par) 
{
    // par[0] - amplitude of signal
    // par[1] - mean
    // par[2] - sigma
    // par[3] - radiation decay
    Double_t gauss = par[0] * TMath::Gaus(x[0], par[1], par[2]);
    Double_t tail = TMath::Exp(-par[3] * (x[0] - par[1]));

    if(x[0] >= par[1])  return gauss * tail;
    else                return gauss;
}

// Define the combined function
Double_t combined(Double_t *x, Double_t *par) 
{
    return background(x, par) + peak(x, &par[3]);
}

struct ExtrapolationResult
{
    double A = 0.0;
    double B = 0.0;
    int w_reper = 0;
};

struct ExtrapolationOutput {
    std::vector<std::vector<double>> integral_background;
    std::vector<std::vector<double>> integral_background_error;
    std::vector<std::vector<double>> w_bg_fit;
    std::vector<ExtrapolationResult> extrapolation_results; // Вместо сырого массива [6] используем вектор
};

struct FitResult
{
    bool fitted = false;

    double bg_amp = 0;
    double bg_mean = 0;
    double bg_sigma = 0;

    double bg_amp_error = 0;
    double bg_mean_error = 0;
    double bg_sigma_error = 0;

    double signal_amp = 0;
    double signal_mean = 0;
    double signal_sigma = 0;
    double signal_tail = 0;

    double signal_amp_error = 0;
    double signal_mean_error = 0;
    double signal_sigma_error = 0;
    double signal_tail_error = 0;

    double xmin = 0;
    double xmax = 0;
};

struct YieldsOutput {
    std::vector<std::vector<double>> w_vector;
    std::vector<std::vector<double>> w_vector_error;
    std::vector<std::vector<double>> integral_signal_excut;
    std::vector<std::vector<double>> integral_signal_excut_error;
    std::vector<std::vector<double>> integral_signal;
    std::vector<std::vector<double>> integral_error;
};

std::vector<std::vector<TH1F*>> LoadHistograms(
    TFile* infile, 
    const char* hist_name_prefix, 
    int q_min, int q_max, 
    int w_min, int w_max) 
{
    std::vector<std::vector<TH1F*>> histograms(q_max + 1);

    for (int q = q_min; q <= q_max; ++q) 
    {
        histograms[q].resize(w_max + 1, nullptr);

        for (int w = w_min; w <= w_max; ++w) 
        {
            TString full_name = Form("%s_Q2_bin=%d_W_bin=%d", hist_name_prefix, q + 1, w + 1);
            histograms[q][w] = (TH1F*)infile->Get(full_name);
        }
    }

    return histograms;
}

ExtrapolationOutput ExtrapolateBackground(
    const std::vector<std::vector<FitResult>>& fit_results,
    const YieldsConfig& cfg)
{
    ExtrapolationOutput output;
    output.integral_background.resize(cfg.q_max + 1);
    output.integral_background_error.resize(cfg.q_max + 1);
    output.w_bg_fit.resize(cfg.q_max + 1);
    output.extrapolation_results.resize(cfg.q_max + 1);

    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
        int w_reper = cfg.w_fit_min;

        if (q == 4 || q == 5)
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
        output.extrapolation_results[q].A = fitFunction->GetParameter(0);
        output.extrapolation_results[q].B = fitFunction->GetParameter(1);
        output.extrapolation_results[q].w_reper = w_reper;

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

void DrawBackgroundExtrapolation(
    const std::vector<std::vector<double>> integral_background,
    const std::vector<std::vector<double>> integral_background_error,
    const std::vector<std::vector<double>> w_bg_fit,
    const std::vector<ExtrapolationResult>& extrapolation_results,
    const YieldsConfig& cfg,
    const TString& pdfFileName)
{
    TCanvas* canvas_fit_bg = new TCanvas(
        "canvas_fit_bg",
        "Fit bg",
        800,
        600
    );

    canvas_fit_bg->Print(pdfFileName + "[");

    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
        // ---------------------------------------------------------
        // Q2 range
        // ---------------------------------------------------------

        double Q2_min_val =
            cfg.Q2_vals[q];

        double Q2_max_val =
            cfg.Q2_vals[q + 1];

        TString namegraph = Form(
            "Background VS W_bin for Q^{2} in [%g, %g] GeV^{2}; W bin number; BG",
            Q2_min_val,
            Q2_max_val
        );

        // ---------------------------------------------------------
        // Graph
        // ---------------------------------------------------------

        TGraphErrors* graph_fit_bg = new TGraphErrors(
            integral_background[q].size(),
            w_bg_fit[q].data(),
            integral_background[q].data(),
            nullptr,
            integral_background_error[q].data()
        );

        graph_fit_bg->SetTitle(namegraph);
        graph_fit_bg->SetMinimum(0);
        graph_fit_bg->SetMarkerSize(10.0);
        graph_fit_bg->GetXaxis()->SetNdivisions(20, kTRUE);

        // ---------------------------------------------------------
        // Extrapolation function
        // ---------------------------------------------------------

        int w_reper =
            extrapolation_results[q].w_reper;

        TF1* fitFunction = new TF1(
            Form("fitFunction_draw_q%d", q),
            Form(
                "[0]*(x-%d)*(x-%d)+[1]*(x-%d)",
                w_reper,
                w_reper,
                w_reper
            ),
            w_reper - 0.1,
            cfg.w_fit_max
        );

        fitFunction->SetParameter(
            0,
            extrapolation_results[q].A
        );

        fitFunction->SetParameter(
            1,
            extrapolation_results[q].B
        );

        // ---------------------------------------------------------
        // Draw
        // ---------------------------------------------------------

        graph_fit_bg->Draw("AP");
        fitFunction->Draw("P SAME");

        graph_fit_bg->GetXaxis()->SetRangeUser(0, 20);
        graph_fit_bg->GetYaxis()->SetRangeUser(
            0,
            fitFunction->Eval(cfg.w_fit_max) * 1.1
        );

        // ---------------------------------------------------------
        // Fit information
        // ---------------------------------------------------------

        TPaveText* pave = new TPaveText(
            0.7,
            0.4,
            0.9,
            0.6,
            "NDC"
        );

        pave->SetFillColor(0);
        pave->SetTextAlign(12);
        pave->SetTextSize(0.03);

        pave->AddText(
            Form(
                "A*(x - %d)^{2}+B*(x-%d)",
                w_reper,
                w_reper
            )
        );

        pave->AddText("Fit results:");

        pave->AddText(
            Form(
                "   A = %.5f",
                extrapolation_results[q].A
            )
        );

        pave->AddText(
            Form(
                "   B = %.5f",
                extrapolation_results[q].B
            )
        );

        pave->Draw("SAME");

        // ---------------------------------------------------------
        // Save page
        // ---------------------------------------------------------

        canvas_fit_bg->Update();
        canvas_fit_bg->Print(pdfFileName);

        delete pave;
        delete fitFunction;
        delete graph_fit_bg;
    }

    canvas_fit_bg->Print(pdfFileName + "]");

    delete canvas_fit_bg;
}

std::vector<std::vector<FitResult>> FitBackground(
    const std::vector<std::vector<TH1F*>>& hist_bg,
    const YieldsConfig& cfg)
{
    int n_q = cfg.q_max - cfg.q_min + 1;

    std::vector<std::vector<FitResult>> fit_results(
        n_q,
        std::vector<FitResult>(100)
    );

    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
        int w_brink = 64;

        if (q == 4)
            w_brink = 56;
        else if (q == 5)
            w_brink = 37;

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

            FitResult& result =
                fit_results[q - cfg.q_min][w];

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

YieldsOutput CalculateYields(
    const std::vector<std::vector<TH1F*>>& hist_signal,
    const std::vector<std::vector<TH1F*>>& hist_bg,
    const std::vector<std::vector<FitResult>>& fit_results,
    const std::vector<std::vector<double>>& integral_background, // добавил ссылку &, чтобы не копировать
    const YieldsConfig& cfg)
{
    YieldsOutput res;
    res.w_vector.resize(cfg.q_max + 1);
    res.w_vector_error.resize(cfg.q_max + 1);
    res.integral_signal_excut.resize(cfg.q_max + 1);
    res.integral_signal_excut_error.resize(cfg.q_max + 1);
    res.integral_signal.resize(cfg.q_max + 1);
    res.integral_error.resize(cfg.q_max + 1);

    const double M_prot = 0.938;
    const double E_beam = 6.535;
    const double alpha = 1.0 / 137.0;
    const double del_w = 0.025;

    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
        double Q2_min_val = cfg.Q2_vals[q];
        double Q2_max_val = cfg.Q2_vals[q + 1];

        double del_q2 = Q2_max_val - Q2_min_val;

        double q2_mid_val = (Q2_max_val + Q2_min_val) / 2.0;

        int w_brink = 64;

        if (q == 4)
            w_brink = 56;
        else if (q == 5)
            w_brink = 37;

        for (int w = 0; w < w_brink; ++w)
        {
            //W
            double w_min_val = 1.4 + w * del_w;
            double w_max_val = 1.4 + (w + 1) * del_w;
            double w_mid_val = (w_min_val + w_max_val) / 2.0;

            res.w_vector[q].push_back(w_mid_val);
            res.w_vector_error[q].push_back(del_w);

            //Flux
            double E_prime =
                E_beam -
                (w_mid_val * w_mid_val
                 - M_prot * M_prot
                 + q2_mid_val)
                / (2.0 * M_prot);

            double nu = E_beam - E_prime;

            double sin2 =
                q2_mid_val /
                (4.0 * E_beam * E_prime);

            double cos2 = 1.0 - sin2;
            double tg2 = sin2 / cos2;

            double eps =
                1.0 /
                (
                    1.0 +
                    2.0 * (1.0 + nu * nu / q2_mid_val) * tg2
                );

            double flux =
                alpha / (4.0 * 3.1415)
                * 1.0 /
                  (M_prot * M_prot * E_beam * E_beam)
                * (
                    w_mid_val *
                    (w_mid_val * w_mid_val - M_prot * M_prot)
                  )
                / (
                    (1.0 - eps) * q2_mid_val
                  );

            double delta_flux =
                sqrt(
                    pow(
                        alpha / (4.0 * 3.1415)
                        * 1.0 /
                          (M_prot * M_prot * E_beam * E_beam)
                        * (
                            (3.0 * w_mid_val * w_mid_val
                             - M_prot * M_prot)
                            / (2.0 * q2_mid_val)
                            - 1.0
                          ),
                        2
                    ) * pow(del_w, 2)
                    +
                    pow(
                        alpha / (4.0 * 3.1415)
                        * 1.0 /
                          (M_prot * M_prot * E_beam * E_beam)
                        * (
                            (w_mid_val * w_mid_val
                             - M_prot * M_prot)
                            / (2.0 * q2_mid_val)
                          ),
                        2
                    ) * pow(del_q2, 2)
                );

            // N_ex_cut from signal histogram
            double N_ex_cut = hist_signal[q][w]->Integral(25, 40);

            // Yield without background subtraction

            double yield_ex_cut = N_ex_cut / del_q2 / del_w / flux;
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

            if (q == 4 || q == 5)
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
            double yield = N / del_q2 / del_w / flux;

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

void DrawMissingMassFit(
    const std::vector<std::vector<TH1F*>>& hist,
    const std::vector<std::vector<FitResult>>& fit_results,
    const char* pdfFilePrefix,
    const YieldsConfig& cfg
)
{
    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
        double Q2_min_val = cfg.Q2_vals[q];
        double Q2_max_val = cfg.Q2_vals[q + 1];

        int w_brink = 64;

        if (q == 4)
            w_brink = 56;
        else if (q == 5)
            w_brink = 37;

        TCanvas* canvas = new TCanvas(
            Form("canvas_q%d", q),
            "Canvas",
            800,
            600
        );

        TString pdfFileName =
            Form(
                "%s_in_Q2_bin_%d.pdf",
                pdfFilePrefix,
                q + 1
            );

        canvas->Print(pdfFileName + "[");

        for (int w = 0; w < w_brink; ++w)
        {
            double w_min_val = 1.4 + w * 0.025;
            double w_max_val = 1.4 + (w + 1) * 0.025;

            double h_raw = hist[q][w]->GetMaximum();

            TString title = Form(
                "LogNorm method Q^{2} [%g,%g] GeV^{2} "
                "W [%g,%g] GeV; MM_X^{2}, GeV^{2}",
                Q2_min_val,
                Q2_max_val,
                w_min_val,
                w_max_val
            );

            // ---------------------------------------------------------
            // Raw histogram
            // ---------------------------------------------------------

            hist[q][w]->SetTitle(title);

            hist[q][w]->SetStats(kFALSE);
            hist[q][w]->SetLineColor(4);
            hist[q][w]->SetLineWidth(3);

            hist[q][w]->GetYaxis()->SetRangeUser(
                0,
                h_raw * 1.05
            );

            hist[q][w]->Draw();

            // ---------------------------------------------------------
            // Draw fit
            // ---------------------------------------------------------

            const FitResult& result =
                fit_results[q][w];

            if (result.fitted)
            {
                TF1* fitfunc = new TF1(
                    Form("fitfunc_draw_q%d_w%d", q, w),
                    combined,
                    result.xmin,
                    result.xmax,
                    7
                );

                fitfunc->SetParameters(
                    result.bg_amp,
                    result.bg_mean,
                    result.bg_sigma,
                    result.signal_amp,
                    result.signal_mean,
                    result.signal_sigma,
                    result.signal_tail
                );

                fitfunc->SetLineColor(kMagenta);
                fitfunc->SetLineWidth(2);
                fitfunc->Draw("SAME");

                // Signal
                TF1* peakfunc = new TF1(
                    Form("peakfunc_draw_q%d_w%d", q, w),
                    peak,
                    result.xmin,
                    result.xmax,
                    4
                );

                peakfunc->SetParameters(
                    result.signal_amp,
                    result.signal_mean,
                    result.signal_sigma,
                    result.signal_tail
                );

                peakfunc->SetLineColor(kRed);
                peakfunc->SetLineWidth(2);
                peakfunc->Draw("SAME");

                // Background
                TF1* bg = new TF1(
                    Form("background_draw_q%d_w%d", q, w),
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

                bg->SetLineColor(kBlack);
                bg->SetLineStyle(kDashed);
                bg->SetLineWidth(2);
                bg->Draw("SAME");

                // -----------------------------------------------------
                // Vertical line at maximum
                // -----------------------------------------------------

                TLine* l_3 = new TLine(
                    hist[q][w]->GetMaximumBin() * 0.01 - 0.3,
                    0,
                    hist[q][w]->GetMaximumBin() * 0.01 - 0.3,
                    h_raw
                );

                l_3->SetLineColor(kBlack);
                l_3->SetLineStyle(kDashed);
                l_3->Draw("SAME");

                // -----------------------------------------------------
                // Fit information
                // -----------------------------------------------------

                TPaveText* pave = new TPaveText(
                    0.45,
                    0.6,
                    0.88,
                    0.9,
                    "NDC"
                );

                pave->SetFillColor(0);
                pave->SetTextAlign(12);
                pave->SetTextSize(0.03);

                pave->AddText(
                    Form(
                        "Total events: %.0f",
                        hist[q][w]->Integral()
                    )
                );

                pave->AddText("Fit results:");

                pave->AddText(
                    Form(
                        "Amplitude bg: %.5f +/- %.5f",
                        result.bg_amp,
                        result.bg_amp_error
                    )
                );

                pave->AddText(
                    Form(
                        "Mean bg: %.5f +/- %.5f",
                        result.bg_mean,
                        result.bg_mean_error
                    )
                );

                pave->AddText(
                    Form(
                        "Sigma bg: %.5f +/- %.5f",
                        result.bg_sigma,
                        result.bg_sigma_error
                    )
                );

                pave->AddText(
                    Form(
                        "Amplitude signal: %.5f +/- %.5f",
                        result.signal_amp,
                        result.signal_amp_error
                    )
                );

                pave->AddText(
                    Form(
                        "Mean signal: %.5f +/- %.5f",
                        result.signal_mean,
                        result.signal_mean_error
                    )
                );

                pave->AddText(
                    Form(
                        "Sigma signal: %.5f +/- %.5f",
                        result.signal_sigma,
                        result.signal_sigma_error
                    )
                );

                pave->AddText(
                    Form(
                        "Rad. tail: %.5f +/- %.5f",
                        result.signal_tail,
                        result.signal_tail_error
                    )
                );

                pave->Draw("SAME");

                // -----------------------------------------------------
                // Legend
                // -----------------------------------------------------

                TLegend* legend = new TLegend(
                    0.6,
                    0.4,
                    0.9,
                    0.6
                );

                legend->SetTextSize(0.02);

                legend->AddEntry(
                    hist[q][w],
                    "Pass2 Data",
                    "l"
                );

                legend->AddEntry(
                    fitfunc,
                    "Overall fit",
                    "l"
                );

                legend->AddEntry(
                    peakfunc,
                    "Signal (gaus+rad.tail)",
                    "l"
                );

                legend->AddEntry(
                    bg,
                    "Background (lognormal)",
                    "l"
                );

                legend->Draw("SAME");
            }

            // ---------------------------------------------------------
            // Integration limits
            // ---------------------------------------------------------

            TLine* l_1 = new TLine(
                -0.05,
                0,
                -0.05,
                h_raw
            );

            TLine* l_2 = new TLine(
                0.1,
                0,
                0.1,
                h_raw
            );

            l_1->SetLineColor(kBlack);
            l_2->SetLineColor(kBlack);

            l_1->Draw("SAME");
            l_2->Draw("SAME");

            canvas->Update();
            canvas->Print(pdfFileName);

            // Здесь можно удалить созданные объекты,
            // если они больше нигде не используются.
        }

        canvas->Print(pdfFileName + "]");
        delete canvas;
    }
}

void yields_pdf(void) 
{
    gErrorIgnoreLevel = kWarning;
	gROOT->SetBatch(kTRUE); 

    TFile *missing_pip_file = TFile::Open("interim/mm_pip_hists_m_pip.root");
    TFile *fullly_exclusive_file = TFile::Open("interim/mm_pip_hists_m_0.root");

    int q_min = 0; int q_max = 5;
    int w_min = 0; int w_max = 99;

    auto missing_pip_hist = LoadHistograms(
        missing_pip_file, 
        "MM", 
        q_min, q_max, w_min, w_max
    );

    auto fullly_exclusive_hist = LoadHistograms(
        fullly_exclusive_file, 
        "MMpiplus_from_MM0_topology", 
        q_min, q_max, w_min, w_max
    );

//////////////Missing Pi+/////////////////////////////////////////////////////

    std::cout << "Missing Pi+" << std::endl;

    YieldsConfig missing_pip_cfg;

    std::cout << "    FitBackground" << std::endl;

    auto missing_pip_fit_results = FitBackground(
        missing_pip_hist,
        missing_pip_cfg
    );

    std::cout << "    ExtrapolateBackground" << std::endl;

    auto [
        integral_background, 
        integral_background_error, 
        w_bg_fit, 
        missing_pip_extrapolation_results
    ] = 
    ExtrapolateBackground(
        missing_pip_fit_results,
        missing_pip_cfg
    );

    std::cout << "    CalculateYields" << std::endl;

    auto [
        w_vector,
        w_vector_error,
        intergral_signal_excut,
        intergral_signal_excut_error,
        intergral_signal,
        intergral_signal_error
    ] =
    CalculateYields(
        missing_pip_hist,
        missing_pip_hist,
        missing_pip_fit_results,
        integral_background,
        missing_pip_cfg
    );

    std::cout << "    DrawBackgroundExtrapolation" << std::endl;

    DrawBackgroundExtrapolation(
        integral_background,
        integral_background_error,
        w_bg_fit,
        missing_pip_extrapolation_results,
        missing_pip_cfg,
        "results/extrapolating_bg.pdf"
    );

    std::cout << "    DrawMissingMassFit" << std::endl;

    DrawMissingMassFit(
        missing_pip_hist, 
        missing_pip_fit_results, 
        "results/LogNormalFit",
        missing_pip_cfg
    );

//////////////Fully exclusive//////////////////////////////////////////////////////////////////////////////
    std::cout << "Fully exclusive" << std::endl;

    YieldsConfig fully_exclusive_bg_cfg;
    fully_exclusive_bg_cfg.signal_sigma = 0.05;
    fully_exclusive_bg_cfg.xmin = -0.1;
    fully_exclusive_bg_cfg.xmax_high = 0.3;
    fully_exclusive_bg_cfg.bg_ampl_bin = 50;
    fully_exclusive_bg_cfg.bg_ampl_factor = 1.0;

    std::cout << "    FitBackground" << std::endl;

    auto fully_exclusive_fit_results = FitBackground(
        fullly_exclusive_hist,
        fully_exclusive_bg_cfg
    );

    std::cout << "    ExtrapolateBackground" << std::endl;

    auto [
        background_MM0,
        background_error_MM0,
        w_bg_fit_MM0, 
        fully_exclusive_extrapolation_results
    ] =
    ExtrapolateBackground(
        fully_exclusive_fit_results,
        fully_exclusive_bg_cfg
    );

    std::cout << "    CalculateYields" << std::endl;

    auto [
        w_vector_MM0,
        w_vector_MM0_error,
        intergral_signal_excut_MM0,
        intergral_signal_excut_MM0_error,
        intergral_signal_4th_method_MM0,
        intergral_signal_4th_method_MM0_error
    ] =
    CalculateYields(
        missing_pip_hist,
        fullly_exclusive_hist,
        fully_exclusive_fit_results,
        background_MM0,
        fully_exclusive_bg_cfg
    );

    std::cout << "    DrawBackgroundExtrapolation" << std::endl;

    DrawBackgroundExtrapolation(
        background_MM0,
        background_error_MM0,
        w_bg_fit_MM0,
        fully_exclusive_extrapolation_results,
        fully_exclusive_bg_cfg,
        "results/extrapolating_bg_MM0.pdf"
    );

    std::cout << "    DrawMissingMassFit" << std::endl;

    DrawMissingMassFit(
        fullly_exclusive_hist, 
        fully_exclusive_fit_results, 
        "results/MMpiplus_from_MM0_topology_fit",
        fully_exclusive_bg_cfg
    );

//////////////GRAPHS/////////////////////////////////////////////////////////////////////////////

    std::cout << "GRAPHS" << std::endl;

    TCanvas *canvas_yields = new TCanvas("canvas_yields", "Canvas Title", 800, 600);
    TString pdfFileName_graph = Form("results/Yields_in_all_Q2_bins.pdf");  // Генерируем имя файла на основе q
    canvas_yields->Print(pdfFileName_graph + "[");

    for (int q = 0; q < 6; ++q) 
    {
	    double Q2_min_val, Q2_max_val;

	    if(q == 0)  {Q2_min_val = 0.5, Q2_max_val = 0.7;}
	    if(q == 1)  {Q2_min_val = 0.7, Q2_max_val = 1;}
	    if(q == 2)  {Q2_min_val = 1, Q2_max_val = 1.4;}
	    if(q == 3)  {Q2_min_val = 1.4, Q2_max_val = 2;}
	    if(q == 4)  {Q2_min_val = 2, Q2_max_val = 3;}
	    if(q == 5)  {Q2_min_val = 3, Q2_max_val = 5.5;}

	    char namegraph[256];
        sprintf(namegraph, "Yield VS W for Q2 in [%g,%g] GeV^2; W, GeV; Yield", Q2_min_val, Q2_max_val);

        TMultiGraph *mg = new TMultiGraph();
	
        TGraphErrors *graph_signal = new TGraphErrors(intergral_signal[q].size(), &w_vector[q][0], &intergral_signal[q][0], NULL, &intergral_signal_error[q][0]);
        graph_signal -> SetTitle(namegraph);
        graph_signal -> SetMinimum(0);
        graph_signal -> SetMarkerSize(1);
        graph_signal -> Draw("AC");

        mg -> Add(graph_signal);

        TGraphErrors *graph_signal_excut = new TGraphErrors(intergral_signal_excut[q].size(), &w_vector[q][0], &intergral_signal_excut[q][0], NULL, &intergral_signal_excut_error[q][0]);
        graph_signal_excut -> SetTitle(namegraph);
        graph_signal_excut -> SetMinimum(0);
        graph_signal_excut -> SetMarkerSize(1);
        graph_signal_excut -> SetLineColor(3);
        graph_signal_excut -> Draw("SAME");

        mg -> Add(graph_signal_excut);

        TGraphErrors *graph_signal_4 = new TGraphErrors(intergral_signal_4th_method_MM0[q].size(), &w_vector_MM0[q][0], &intergral_signal_4th_method_MM0[q][0], NULL, &intergral_signal_4th_method_MM0_error[q][0]);
        graph_signal_4 -> SetTitle(namegraph);
        graph_signal_4 -> SetMinimum(0);
        graph_signal_4 -> SetMarkerSize(1);
        graph_signal_4 -> SetLineColor(6);
        graph_signal_4 -> Draw("SAME");

        mg -> Add(graph_signal_4);
        mg -> SetTitle(namegraph);
        mg -> SetMinimum(0);
        mg -> GetXaxis() -> SetTitle("W, GeV");
        mg -> GetYaxis() -> SetTitle("Yield");

        // Нарисовать график
        mg -> Draw("AC");

        TLegend* legend = new TLegend(0.6, 0.6, 0.9, 0.9);
        legend -> SetTextSize(0.03);
        legend -> AddEntry(graph_signal, "New lognorm method", "l");
        legend -> AddEntry(graph_signal_excut, "Entries within ex.cut", "l");
        legend -> AddEntry(graph_signal_4, "Scaled bg from MM0 topology", "l");

        legend->Draw("SAME");

        canvas_yields -> Update();
        canvas_yields -> Print(pdfFileName_graph);
    }

    canvas_yields -> Print(pdfFileName_graph + "]"); // Закрываем текущий PDF-файл
    delete canvas_yields;

    gSystem -> Exit(0);
}