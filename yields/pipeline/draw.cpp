void DrawBackgroundExtrapolation(
    const ExtrapolationOutput& extrapolation_output,
    const YieldsConfig& cfg,
    const char* topology_name)
{
    TCanvas* canvas_fit_bg = new TCanvas(
        "canvas_fit_bg",
        "Fit bg",
        1280,
        720
    );

    TString pdfFileName = TString::Format("results/extrapolating_bg_%s.pdf", topology_name);

    canvas_fit_bg->Print(pdfFileName + "[");

    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
        // Q2 range
        double Q2_min_val = cfg.Q2_edges[q];
        double Q2_max_val = cfg.Q2_edges[q + 1];

        TString namegraph = Form(
            "Background for Q^{2} #in [%g, %g] GeV^{2}; W bin; Number of bg events",
            Q2_min_val,
            Q2_max_val
        );

        // Graph
        TGraphErrors* graph_fit_bg = new TGraphErrors(
            extrapolation_output.integral_background[q].size(),
            extrapolation_output.w_bg_fit[q].data(),
            extrapolation_output.integral_background[q].data(),
            nullptr,
            extrapolation_output.integral_background_error[q].data()
        );

        graph_fit_bg->SetTitle(namegraph);
        graph_fit_bg->SetMinimum(0);
        graph_fit_bg->SetMarkerSize(10.0);
        graph_fit_bg->GetXaxis()->SetNdivisions(20, kTRUE);

        // Extrapolation function
        int w_reper = extrapolation_output.results[q].w_reper;

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
            extrapolation_output.results[q].A
        );

        fitFunction->SetParameter(
            1,
            extrapolation_output.results[q].B
        );

        // Draw
        graph_fit_bg->Draw("AP");
        fitFunction->Draw("SAME");

        graph_fit_bg->GetXaxis()->SetRangeUser(0, 20);
        graph_fit_bg->GetYaxis()->SetRangeUser(
            0,
            fitFunction->Eval(cfg.w_fit_max) * 1.1
        );

        // Fit information
        TPaveText* pave = new TPaveText(
            0.7,
            0.4,
            0.71,
            0.6,
            "NDC"
        );

        pave->SetFillColor(0);
        pave->SetTextAlign(12);
        pave->SetTextSize(0.03);

        pave->AddText(
            Form(
                "A #bullet (x - %d)^{2} + B #bullet (x - %d)",
                w_reper,
                w_reper
            )
        );

        pave->AddText("Fit results:");

        pave->AddText(
            Form(
                "   A = %.2f",
                extrapolation_output.results[q].A
            )
        );

        pave->AddText(
            Form(
                "   B = %.2f",
                extrapolation_output.results[q].B
            )
        );

        pave->Draw("SAME");
        
        // Save page
        canvas_fit_bg->Update();
        canvas_fit_bg->Print(pdfFileName);

        delete pave;
        delete fitFunction;
        delete graph_fit_bg;
    }

    canvas_fit_bg->Print(pdfFileName + "]");

    delete canvas_fit_bg;
}

void DrawMissingMassFit(
    const std::vector<std::vector<TH1F*>>& hist,
    const std::vector<std::vector<FitResult>>& fit_results,
    const char* topology_name,
    const YieldsConfig& cfg
)
{
    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
        double Q2_min_val = cfg.Q2_edges[q];
        double Q2_max_val = cfg.Q2_edges[q + 1];

        int w_brink = cfg.WBinLimit(q);

        TCanvas* canvas = new TCanvas(
            Form("canvas_q%d", q),
            "Canvas",
            800,
            600
        );
        
        TString pdfFileName = TString::Format("results/%s_fit_in_Q2_bin_%d.pdf", topology_name, q + 1);

        canvas->Print(pdfFileName + "[");

        for (int w = 0; w < w_brink; ++w)
        {
            double w_min_val = cfg.W_min + w * cfg.W_bin_size;
            double w_max_val = cfg.W_min + (w + 1) * cfg.W_bin_size;

            double h_raw = hist[q][w]->GetMaximum();
            double y_axis_maximum = h_raw * 1.05;

            TString title = Form(
                "MM^{2} distribution for Q^{2} #in [%g, %g] GeV^{2}, "
                "W #in [%g, %g] GeV; MM_{#pi^{+}}^{2}, GeV^{2}",
                Q2_min_val,
                Q2_max_val,
                w_min_val,
                w_max_val
            );

            
            // Raw histogram
            hist[q][w]->SetTitle(title);

            hist[q][w]->SetStats(kFALSE);
            hist[q][w]->SetLineColor(4);
            hist[q][w]->SetLineWidth(3);
            
            hist[q][w]->GetYaxis()->SetRangeUser(
                0, y_axis_maximum
            );

            hist[q][w]->Draw();

            // Draw fit
            const FitResult& result = fit_results[q][w];

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

                // Vertical line at maximum
                double hist_bin_width =
                    (cfg.hist_max - cfg.hist_min) / cfg.n_hist_bins;

                TLine* maximum_bin_line = new TLine(
                    cfg.hist_min + hist[q][w]->GetMaximumBin() * hist_bin_width,
                    0,
                    cfg.hist_min + hist[q][w]->GetMaximumBin() * hist_bin_width,
                    y_axis_maximum
                );

                maximum_bin_line->SetLineColor(kBlack);
                maximum_bin_line->SetLineStyle(kDashed);
                maximum_bin_line->Draw("SAME");

                // Fit information
                TPaveText* pave = new TPaveText(
                    0.55, 0.55, 0.89, 0.89, "NDC"
                );

                pave->SetFillColor(0);
                pave->SetTextAlign(12);
                pave->SetTextSize(0.03);

                float events_within_ec = hist[q][w]->Integral(
                    cfg.signal_integral_bin_min,
                    cfg.signal_integral_bin_max
                );
                float bg_events_within_ec = bg->Integral(
                    cfg.bg_integral_xmin,
                    cfg.bg_integral_xmax
                ) * cfg.bin_width_factor;

                pave->AddText(
                    Form(
                        "Total events: %.0f",
                        hist[q][w]->Integral()
                    )
                );

                pave->AddText(
                    Form(
                        "Events within e.c.: %.0f",
                        events_within_ec
                    )
                );

                pave->AddText(
                    Form(
                        "Bg events within e.c.: %.0f",
                        bg_events_within_ec
                    )
                );

                pave->AddText("");

                pave->AddText(
                    Form(
                        "Amplitude bg: %.0f #pm %.0f",
                        result.bg_amp,
                        result.bg_amp_error
                    )
                );

                pave->AddText(
                    Form(
                        "Mean bg: %.4f #pm %.4f",
                        result.bg_mean,
                        result.bg_mean_error
                    )
                );

                pave->AddText(
                    Form(
                        "#sigma bg: %.4f #pm %.4f",
                        result.bg_sigma,
                        result.bg_sigma_error
                    )
                );

                pave->AddText("");

                pave->AddText(
                    Form(
                        "Amplitude signal: %.0f #pm %.0f",
                        result.signal_amp,
                        result.signal_amp_error
                    )
                );

                pave->AddText(
                    Form(
                        "Mean signal: %.5f #pm %.5f",
                        result.signal_mean,
                        result.signal_mean_error
                    )
                );

                pave->AddText(
                    Form(
                        "#sigma signal: %.5f #pm %.5f",
                        result.signal_sigma,
                        result.signal_sigma_error
                    )
                );

                pave->AddText(
                    Form(
                        "Rad. tail: %.3f #pm %.3f",
                        result.signal_tail,
                        result.signal_tail_error
                    )
                );

                pave->Draw("SAME");

                TPaveText* pave_bg_ratio = new TPaveText(
                    0.11, 0.5, 0.2, 0.6, "NDC"
                );

                pave_bg_ratio->SetFillColor(0);
                pave_bg_ratio->SetTextAlign(12);
                pave_bg_ratio->SetTextSize(0.05);
                pave_bg_ratio->SetTextColor(kRed);

                pave_bg_ratio->AddText(
                    Form(
                        "   R = %.2f ", 
                        bg_events_within_ec / events_within_ec
                    )
                );

                pave_bg_ratio->Draw("SAME");
                
                // Legend
                TLegend* legend = new TLegend(
                    0.65, 0.3, 0.88, 0.5
                );

                legend->SetTextSize(0.02);

                legend->AddEntry(
                    hist[q][w],
                    "Experimental data",
                    "l"
                );

                legend->AddEntry(
                    fitfunc,
                    "Overall fit",
                    "l"
                );

                legend->AddEntry(
                    peakfunc,
                    "Signal (gaus + rad. tail)",
                    "l"
                );

                legend->AddEntry(
                    bg,
                    "Background (lognormal)",
                    "l"
                );

                legend->Draw("SAME");
            }
            
            // Integration limits
            TLine* left_exclusivity_cut = new TLine(
                cfg.bg_integral_xmin,
                0,
                cfg.bg_integral_xmin,
                y_axis_maximum
            );

            left_exclusivity_cut->SetLineColor(kBlack);
            left_exclusivity_cut->Draw("SAME");

            TLine* right_exclusivity_cut = new TLine(
                cfg.bg_integral_xmax,
                0,
                cfg.bg_integral_xmax,
                y_axis_maximum
            );

            right_exclusivity_cut->SetLineColor(kBlack);
            right_exclusivity_cut->Draw("SAME");

            canvas->Update();
            canvas->Print(pdfFileName);

            // Здесь можно удалить созданные объекты,
            // если они больше нигде не используются.
        }

        canvas->Print(pdfFileName + "]");
        delete canvas;
    }
}

void DrawYields(
    YieldsOutput missing_pip,
    YieldsOutput fully_exclusive,
    const YieldsConfig& cfg
)
{
    std::cout << "Drawing Yields Distribution" << std::endl;

    TCanvas *canvas_yields = new TCanvas("canvas_yields", "Canvas Title", 800, 600);
    TString pdfFileName_graph = Form("results/Yields_in_all_Q2_bins.pdf");  
    canvas_yields->Print(pdfFileName_graph + "[");

    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
	    double Q2_min_val = cfg.Q2_edges[q];
	    double Q2_max_val = cfg.Q2_edges[q + 1];

	    char namegraph[256];
        sprintf(namegraph, "Yield for Q^{2} #in [%g, %g] GeV^{2}; W, GeV; Yield", Q2_min_val, Q2_max_val);

        TMultiGraph *mg = new TMultiGraph();
        mg->SetTitle(namegraph);

        TGraphErrors *yield_without_bg_estimation = new TGraphErrors(
            missing_pip.integral_signal_excut[q].size(), 
            &missing_pip.w_vector[q][0], 
            &missing_pip.integral_signal_excut[q][0], 
            NULL, 
            &missing_pip.integral_signal_excut_error[q][0]
        );
        yield_without_bg_estimation->SetMarkerSize(1);
        yield_without_bg_estimation->SetLineColor(3);
        yield_without_bg_estimation->Draw("SAME");

        mg->Add(yield_without_bg_estimation);
	
        TGraphErrors *yield_with_bg_from_mpip = new TGraphErrors(
            missing_pip.integral_signal[q].size(), 
            &missing_pip.w_vector[q][0], 
            &missing_pip.integral_signal[q][0], 
            NULL, 
            &missing_pip.integral_error[q][0] // Проверьте имя поля: integral_error или integral_signal_error
        );
        yield_with_bg_from_mpip->SetMarkerSize(1);
        yield_with_bg_from_mpip->Draw("AC");

        mg->Add(yield_with_bg_from_mpip);

        TGraphErrors *yield_with_bg_from_fe = new TGraphErrors(
            fully_exclusive.integral_signal[q].size(), 
            &fully_exclusive.w_vector[q][0], 
            &fully_exclusive.integral_signal[q][0], 
            NULL, 
            &fully_exclusive.integral_error[q][0]
        );
        yield_with_bg_from_fe->SetMarkerSize(1);
        yield_with_bg_from_fe->SetLineColor(6);
        yield_with_bg_from_fe->Draw("SAME");

        mg->Add(yield_with_bg_from_fe);

        // Нарисовать график
        mg->Draw("AC");

        TLegend* legend = new TLegend(0.55, 0.7, 0.9, 0.9);
        legend->SetTextSize(0.03);
        legend->AddEntry(yield_without_bg_estimation, "Without bg estimation", "l");
        legend->AddEntry(yield_with_bg_from_mpip, "Bg from M#pi^{+} topology", "l");
        legend->AddEntry(yield_with_bg_from_fe, "Scaled bg from FE topology", "l");

        legend->Draw("SAME");

        canvas_yields->Update();
        canvas_yields->Print(pdfFileName_graph);
    }

    canvas_yields->Print(pdfFileName_graph + "]"); // Закрываем текущий PDF-файл
    delete canvas_yields;
}