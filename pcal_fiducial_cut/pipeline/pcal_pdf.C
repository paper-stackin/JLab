#include "pcal_config_yaml.h"
#include "pcal_boundaries.h"

void InterpolateLocalDips(TH1F* hist, int sector, double x,
                          const std::vector<ExclusionStripSettings::ExclusionStripConfig>& strips,
                          double strip_margin)
{
    const int n_bins = hist->GetNbinsX();

    for (const auto& strip : strips) {
        if (strip.sector != sector) continue;

        // Поворачиваем прямую на угол, соответствующий сектору
        double angle = TMath::Pi() / 180.0 * (-60.0 * (sector - 1));
        double denominator = TMath::Cos(angle) - TMath::Sin(angle) * strip.slope;
        double new_slope = (TMath::Sin(angle) + TMath::Cos(angle) * strip.slope) / denominator;
        double new_c_max = (strip.c_max + strip_margin) / denominator;
        double new_c_min = (strip.c_min - strip_margin) / denominator;

        const double y_lower = new_slope * x + new_c_min;
        const double y_upper = new_slope * x + new_c_max;

        int k_first = -1;
        int k_last = -1;
        for (int k = 1; k <= n_bins; ++k) {
            const double y = hist->GetBinCenter(k);
            if (y > y_lower && y < y_upper) {
                if (k_first < 0) k_first = k;
                k_last = k;
            }
        }
        if (k_first < 0) continue;

        const int k_left = k_first - 1;
        const int k_right = k_last + 1;
        if (k_left < 1 || k_right > n_bins) continue;

        const double baseline = 0.5 * (hist->GetBinContent(k_left) + hist->GetBinContent(k_right));
        // for (int k = k_first; k <= k_last; ++k) {
        for (int k = k_left; k <= k_right; ++k) {
            hist->SetBinContent(k, baseline);
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

void pcal_pdf(const char* config_file)
{
    // Таймер
    auto start_time = std::chrono::steady_clock::now();
    
    // Конфигурация
    PCALConfig cfg = LoadConfig(config_file);

    // Настройки стиля
    gStyle->SetOptStat(cfg.style.stat_option);
    gROOT->SetBatch(cfg.style.use_batch ? kTRUE : kFALSE);
    gStyle->SetPalette(cfg.style.palette);
    gStyle->SetOptFit(cfg.style.fit_option);
    gErrorIgnoreLevel = kError;
    gPrintViaErrorHandler = kTRUE; 

    // Файл с гистограммами
    TFile *file = TFile::Open(cfg.file.intermediate_file.c_str());

    // Гистограммы для каждого сектора и каждого бина по x'
    TH1F *x_bins[6][cfg.histogram.n_x_bins];
    TH1F *x_bins_interp[6][cfg.histogram.n_x_bins];
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < cfg.histogram.n_x_bins; ++j)
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
    if (cfg.dip.enable) {
        for (const auto& sector_cfg : cfg.dip.sectors) {
            int sec = sector_cfg.sector;
            int sec_idx = sec - 1;
            if (sec_idx < 0 || sec_idx >= 6) continue;

            int start_bin = std::max(0, sector_cfg.start_bin - 1);
            int end_bin = (sector_cfg.end_bin > 0) ? std::min(cfg.histogram.n_x_bins - 1, sector_cfg.end_bin - 1) : cfg.histogram.n_x_bins - 1;

            for (int bin = start_bin; bin <= end_bin; ++bin) {
                TString interp_name = Form("h_pcal_sec%d_bin%d_interp", sec, bin + 1);
                x_bins_interp[sec_idx][bin] = (TH1F*)x_bins[sec_idx][bin]->Clone(interp_name);
                if (cfg.exclusion.enable) {
                    const double x = cfg.histogram.x_start + cfg.histogram.bin_size * (bin + 0.5);
                    InterpolateLocalDips(x_bins_interp[sec_idx][bin], sec, x,
                                         cfg.exclusion.strips, cfg.exclusion.margin);
                }
            }
        }
    }

    // Графики для каждого сектора
    for (int sec = 0; sec < 6; ++sec)
    {
        TCanvas *c = new TCanvas(
            Form("c_sec%d", sec + 1),
            Form("Sector %d", sec + 1),
            cfg.style.canvas_width, cfg.style.canvas_height
        );

        TString pdf_name = Form("%s%d.pdf", cfg.file.pdf_name_prefix.c_str(), sec + 1);
        c->Print(pdf_name + "[");

        // Точки для построения графиков верхней и нижней границ
        std::vector<double> x_vals, upper_vals, lower_vals;

        // Распределение по бинам x'
        for (int bin = 0; bin < cfg.histogram.n_x_bins; ++bin)
        {
            // Гистограммы
            x_bins[sec][bin]->SetLineColor(kBlue);
            x_bins[sec][bin]->SetLineWidth(2);
            x_bins[sec][bin]->Draw("HIST");

            if (cfg.dip.enable && x_bins_interp[sec][bin]) {
                x_bins_interp[sec][bin]->SetLineColor(kRed);
                x_bins_interp[sec][bin]->SetLineWidth(1);
                x_bins_interp[sec][bin]->Draw("HIST SAME");
            } 

            TH1F* hist_for_bounds = (cfg.dip.enable && x_bins_interp[sec][bin])
                ? x_bins_interp[sec][bin]
                : x_bins[sec][bin];

            const char* fit_name = Form("sg_fit_sec%d_bin%d", sec + 1, bin + 1);
            auto bounds = ComputeVerticalBounds(hist_for_bounds, cfg, fit_name);
            double y_lower = bounds.first;
            double y_upper = bounds.second;

            gPad->Update();
            double top = gPad->GetUymax();

            DrawVerticalLine(y_lower, top);
            DrawVerticalLine(y_upper, top);

            // Легенда
            TLegend* leg = new TLegend(0.6, 0.6, 0.93, 0.93);
            leg->AddEntry((TObject*)0, Form("Total events = %.0f", x_bins[sec][bin]->Integral()), "");
            if (cfg.dip.enable && x_bins_interp[sec][bin]) {
                leg->AddEntry((TObject*)0, Form("Total events with interpolation= %.0f", x_bins_interp[sec][bin]->Integral()), "");
            }
            leg->AddEntry((TObject*)0, Form("left = %.4f cm", y_lower), "");
            leg->AddEntry((TObject*)0, Form("right = %.4f cm", y_upper), "");
            // if (cfg.boundary.use_supergauss) {
            //     leg->AddEntry(sg_fit, Form("SuperGauss (n = %d)", cfg.boundary.supergauss_n), "l");
            // }
            leg->Draw();

            c->Print(pdf_name);

            // Точки для построения графиков верхней и нижней границ
            double x = cfg.histogram.x_start + cfg.histogram.bin_size * (bin + 0.5);

            x_vals.push_back(x);
            lower_vals.push_back(y_lower);
            upper_vals.push_back(y_upper);
        }
        
        // 2D гистограмма сектора
        h2_xpy[sec]->Draw("COLZ");

        // Radial каты
        DrawCircle(cfg.radial.r_min, cfg.radial.phi_1, cfg.radial.phi_2);
        DrawCircle(cfg.radial.r_max, cfg.radial.phi_1, cfg.radial.phi_2);

        // Верхняя и нижняя границы
        std::pair<TGraph*, TF1*> lower_result = DrawGraphWithFit(x_vals, lower_vals, cfg.boundary.fit_formula, cfg.boundary.x_fit_min, cfg.boundary.x_fit_max, kRed, "f_lower");
        std::pair<TGraph*, TF1*> upper_result = DrawGraphWithFit(x_vals, upper_vals, cfg.boundary.fit_formula, cfg.boundary.x_fit_min, cfg.boundary.x_fit_max, kRed, "f_upper");
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