#include "PCAL_config.h"

void InterpolateLocalDips(TH1F* hist, double relative_depth = 0.25, int max_width = 3)
{
    int n_bins = hist->GetNbinsX();

    for (int i = 2; i <= n_bins - 1; ) {
        bool filled = false;

        for (int width = std::min(max_width, n_bins - i); width >= 1; --width) {
            int j = i + width - 1;
            if (j >= n_bins - 1) continue;

            double left_neighbor = hist->GetBinContent(i - 1);
            double right_neighbor = hist->GetBinContent(j + 1);
            double baseline = 0.5 * (left_neighbor + right_neighbor);

            double region_sum = 0.0;
            for (int k = i; k <= j; ++k) {
                region_sum += hist->GetBinContent(k);
            }
            double region_avg = region_sum / width;
            double depth = baseline - region_avg;

            if (left_neighbor > region_avg && right_neighbor > region_avg &&
                baseline > 0.0 &&
                (depth / baseline) > relative_depth) {
                std::vector<int> bins_to_fill;
                for (int k = i; k <= j; ++k) {
                    bins_to_fill.push_back(k);
                }
                if (i > 1) bins_to_fill.push_back(i - 1);
                if (j < n_bins) bins_to_fill.push_back(j + 1);

                for (int k : bins_to_fill) {
                    hist->SetBinContent(k, baseline);
                }

                i = j + 1;
                filled = true;
                break;
            }
        }

        if (!filled) {
            ++i;
        }
    }
}

void DrawVerticalLine(double x_pos, double y_top, Color_t color = kBlack, int style = 2)
{
    TLine* line = new TLine(x_pos, 0, x_pos, y_top);
    line->SetLineColor(color);
    line->SetLineStyle(style);
    line->Draw("same");
}

void DrawCircle(double radius, double phi1_deg, double phi2_deg, int N = 300, Color_t color = kRed, int lineWidth = 2)
{
    std::vector<double> x(N), y(N);

    double phi1 = phi1_deg * TMath::DegToRad();
    double phi2 = phi2_deg * TMath::DegToRad();

    for (int i = 0; i < N; ++i) {
        double phi = phi1 + (phi2 - phi1) * i / (N - 1);
        x[i] = radius * cos(phi);
        y[i] = radius * sin(phi);
    }

    TGraph* gCircle = new TGraph(N, x.data(), y.data());
    gCircle->SetLineColor(color);
    gCircle->SetLineWidth(lineWidth);
    gCircle->Draw("L SAME");
}

std::pair<TGraph*, TF1*> DrawGraphWithFit(const std::vector<double>& x_vals, const std::vector<double>& y_vals,
                                          const char* fit_formula, double fit_xmin, double fit_xmax,
                                          Color_t color, const char* fit_name = "fit")
{
    TF1* fit = new TF1(fit_name, fit_formula, fit_xmin, fit_xmax);

    TGraph* graph = new TGraph(x_vals.size(), x_vals.data(), y_vals.data());
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(color);

    graph->Fit(fit, "RNQ");

    fit->SetLineColor(color);
    graph->Draw("P same");
    fit->Draw("same");

    return {graph, fit};
}

void PCAL_pdf(void)
{
    // Таймер
    auto start_time = std::chrono::steady_clock::now();
    
    // Конфигурация
    PCALConfig cfg;

    // Настройки стиля
    gStyle->SetOptStat(cfg.stat_option);
    gROOT->SetBatch(cfg.use_batch ? kTRUE : kFALSE);
    gStyle->SetPalette(cfg.palette);
    gStyle->SetOptFit(cfg.fit_option);
    gErrorIgnoreLevel = kError;
    gPrintViaErrorHandler = kTRUE;

    // Файл с гистограммами
    TFile *file = TFile::Open("PCAL_hist.root");

    // Гистограммы для каждого сектора и каждого бина по x'
    TH1F *x_bins[6][cfg.N_x_bins];
    TH1F *x_bins_interp[6][cfg.N_x_bins];
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < cfg.N_x_bins; ++j)
        {
            TString name = Form("h_pcal_sec%d_bin%d", i + 1, j + 1);
            x_bins[i][j] = (TH1F*)file->Get(name);
            x_bins_interp[i][j] = nullptr;
        }
    }

    // Гистограммы для 2D распределения x' и y' для каждого сектора
    TH2F* h2_xpy[6];
    for (int i = 0; i < 6; ++i)
    {
        TString name = Form("h2_sec%d", i+1);
        h2_xpy[i] = (TH2F*)file->Get(name);
    }

    //  Интерполяция провалов в 1D гистограммах
    if (cfg.enable_dip_interpolation) {
        for (const auto& sector_cfg : cfg.dip_interpolation_sectors) {
            int sec = sector_cfg.sector;
            int sec_idx = sec - 1;
            if (sec_idx < 0 || sec_idx >= 6) continue;

            int start_bin = std::max(0, sector_cfg.start_bin - 1);
            int end_bin = (sector_cfg.end_bin > 0) ? std::min(cfg.N_x_bins - 1, sector_cfg.end_bin - 1) : cfg.N_x_bins - 1;

            for (int bin = start_bin; bin <= end_bin; ++bin) {
                TString interp_name = Form("h_pcal_sec%d_bin%d_interp", sec, bin + 1);
                x_bins_interp[sec_idx][bin] = (TH1F*)x_bins[sec_idx][bin]->Clone(interp_name);
                InterpolateLocalDips(x_bins_interp[sec_idx][bin], sector_cfg.relative_depth, cfg.dip_max_width);
            }
        }
    }

    // Графики для каждого сектора
    for (int sec = 0; sec < 6; ++sec)
    {
        TCanvas *c = new TCanvas(
            Form("c_sec%d", sec + 1),
            Form("Sector %d", sec + 1),
            cfg.canvas_width, cfg.canvas_height
        );

        TString pdf_name = Form("%s%d%s", cfg.pdf_name_prefix, sec + 1, cfg.pdf_extension);
        c->Print(pdf_name + "[");

        std::vector<double> x_vals, upper_vals, lower_vals;

        // Распределение по бинам x'
        for (int bin = 0; bin < cfg.N_x_bins; ++bin)
        {
            // Гистограммы
            x_bins[sec][bin]->SetLineColor(kBlue);
            x_bins[sec][bin]->SetLineWidth(2);
            x_bins[sec][bin]->Draw("HIST");

            if (cfg.enable_dip_interpolation && x_bins_interp[sec][bin]) {
                x_bins_interp[sec][bin]->SetLineColor(kRed);
                x_bins_interp[sec][bin]->SetLineWidth(1);
                x_bins_interp[sec][bin]->Draw("HIST SAME");
            } 

            // Квантили вертикальных линий 
            double probs[2] = {cfg.quantile_low, cfg.quantile_high}, quant[2];
            
            if (cfg.enable_dip_interpolation && x_bins_interp[sec][bin]) {
                x_bins_interp[sec][bin]->GetQuantiles(2, quant, probs);
            } else {
                x_bins[sec][bin]->GetQuantiles(2, quant, probs);
            }

            // Вертикальные линии
            double y_lower = quant[0], y_upper = quant[1];
            gPad->Update();
            double top = gPad->GetUymax();

            DrawVerticalLine(y_lower, top);
            DrawVerticalLine(y_upper, top);

            // Легенда
            TLegend* leg = new TLegend(0.6, 0.6, 0.93, 0.93);
            leg->AddEntry((TObject*)0, Form("Total events = %.0f", x_bins[sec][bin]->Integral()), "");
            if (cfg.enable_dip_interpolation && x_bins_interp[sec][bin]) {
                leg->AddEntry((TObject*)0, Form("Total events with interpolation= %.0f", x_bins_interp[sec][bin]->Integral()), "");
            }
            leg->AddEntry((TObject*)0, Form("left = %.4f cm", y_lower), "");
            leg->AddEntry((TObject*)0, Form("right = %.4f cm", y_upper), "");
            leg->Draw();

            c->Print(pdf_name);

            // Точки для построения графиков верхней и нижней границ
            double x = cfg.x_start + cfg.bin_size * (bin + 0.5);

            x_vals.push_back(x);
            lower_vals.push_back(y_lower);
            upper_vals.push_back(y_upper);
        }
        
        // 2D гистограмма сектора
        h2_xpy[sec]->Draw("COLZ");

        // Radial каты
        DrawCircle(cfg.r_min, cfg.radial_phi1_deg, cfg.radial_phi2_deg);
        DrawCircle(cfg.r_max, cfg.radial_phi1_deg, cfg.radial_phi2_deg);

        // Верхняя и нижняя границы
        std::pair<TGraph*, TF1*> lower_result = DrawGraphWithFit(x_vals, lower_vals, cfg.fit_formula, cfg.x_fit_min, cfg.x_fit_max, kRed, "f_lower");
        std::pair<TGraph*, TF1*> upper_result = DrawGraphWithFit(x_vals, upper_vals, cfg.fit_formula, cfg.x_fit_min, cfg.x_fit_max, kRed, "f_upper");
        TGraph* g_lower = lower_result.first;
        TF1* f_lower = lower_result.second;
        TGraph* g_upper = upper_result.first;
        TF1* f_upper = upper_result.second;

        // Легенда
        // h2_xpy[sec]->SetStats(0); 
        // TLegend* leg_2d = new TLegend(0.7, 0.7, 0.93, 0.93);

        // leg_2d->AddEntry((TObject*)0, Form("Total events = %.0f", h2_xpy[sec]->Integral()), "");
        // leg_2d->AddEntry((TObject*)0, Form("y = A + B * x + C * x^{2}"), "");
        // leg_2d->AddEntry((TObject*)0, Form("A_upper = %.5f", f_upper->GetParameter(0)), "");
        // leg_2d->AddEntry((TObject*)0, Form("B_upper = %.5f", f_upper->GetParameter(1)), "");
        // leg_2d->AddEntry((TObject*)0, Form("C_upper = %.5f", f_upper->GetParameter(2)), "");
        // leg_2d->AddEntry((TObject*)0, Form("A_lower = %.5f", f_lower->GetParameter(0)), "");
        // leg_2d->AddEntry((TObject*)0, Form("B_lower = %.5f", f_lower->GetParameter(1)), "");
        // leg_2d->AddEntry((TObject*)0, Form("C_lower = %.5f", f_lower->GetParameter(2)), "");
        // leg_2d->Draw();

        c->Update();
        c->Print(pdf_name);

        c->Print(pdf_name + "]");

        delete c;
    }

    // Таймер
    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    double elapsed_s = elapsed_ms / 1000.0;

    std::cout << "\nFinished in " << std::fixed << std::setprecision(3) << elapsed_s << " s" << std::endl;

    gSystem->Exit(0);
}