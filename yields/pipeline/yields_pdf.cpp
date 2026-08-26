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

struct BackgroundFitConfig
{
    int q_min = 0;
    int q_max = 5;

    int w_fit_min = 8;
    int w_fit_max = 17;

    int w_reper_q45 = 3;
    int bg_fit_start = 13;

    double xmin = -0.2;

    double xmax_low = 0.15;
    double xmax_mid = 0.25;
    double xmax_high = 0.4;

    int xmax_mid_start = 10;
    int xmax_high_start = 17;

    double bg_mean = -1.0;
    double bg_sigma = 0.3;

    int bg_ampl_bin = 25;
    double bg_ampl_factor = 0.3;

    double signal_mean = 0.0196;
    double signal_sigma = 0.03;
    double signal_tail = 0.02;

    double bg_integral_xmin = -0.05;
    double bg_integral_xmax = 0.1;
    double bin_width_factor = 100.0;
};

void CalculateBackgroundIntegral(
    TH1F* MM_raw[12][100],
    std::vector<double>* intergral_background,
    std::vector<double>* intergral_background_error,
    std::vector<double>* w_bg_fit,
    const BackgroundFitConfig& cfg)
{
    for (int q = cfg.q_min; q <= cfg.q_max; ++q)
    {
        int w_reper = cfg.w_fit_min;

        if (q == 4 || q == 5)
            w_reper = cfg.w_reper_q45;

        for (int w = w_reper; w <= cfg.w_fit_max; ++w)
        {
            double current_xmax;

            if (w < cfg.xmax_mid_start)
                current_xmax = cfg.xmax_low;
            else if (w < cfg.xmax_high_start)
                current_xmax = cfg.xmax_mid;
            else
                current_xmax = cfg.xmax_high;

            int N_bg = 0;

            if (w >= cfg.bg_fit_start)
            {
                TF1* fitfunc = new TF1(
                    Form("fitfunc_q%d_w%d", q, w),
                    combined,
                    cfg.xmin,
                    current_xmax,
                    7
                );

                double bg_ampl =
                    MM_raw[q][w]->GetBinContent(cfg.bg_ampl_bin) * cfg.bg_ampl_factor;

                fitfunc->SetParameters(
                    bg_ampl,
                    cfg.bg_mean,
                    cfg.bg_sigma,
                    MM_raw[q][w]->GetMaximum(),
                    cfg.signal_mean,
                    cfg.signal_sigma,
                    cfg.signal_tail
                );

                MM_raw[q][w]->Fit(fitfunc, "R");

                TF1* bg = new TF1(
                    Form("background_q%d_w%d", q, w),
                    background,
                    cfg.xmin,
                    current_xmax,
                    3
                );

                bg->SetParameters(
                    fitfunc->GetParameter(0),
                    fitfunc->GetParameter(1),
                    fitfunc->GetParameter(2)
                );

                N_bg = bg->Integral(
                    cfg.bg_integral_xmin,
                    cfg.bg_integral_xmax
                ) * cfg.bin_width_factor;

                delete bg;
                delete fitfunc;
            }

            intergral_background[q].push_back(N_bg);
            intergral_background_error[q].push_back(sqrt(N_bg));
            w_bg_fit[q].push_back(w);
        }
    }
}

void ExtrapolateBackground(
    std::vector<double>* intergral_background,
    std::vector<double>* intergral_background_error,
    std::vector<int>* w_bg_fit,
    const BackgroundFitConfig& cfg,
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
        // Q2 range
        // double Q2_vals[7] = {0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.5}; 
        double Q2_min_val;
        double Q2_max_val;

        if (q == 0) {
            Q2_min_val = 0.5;
            Q2_max_val = 0.7;
        }
        else if (q == 1) {
            Q2_min_val = 0.7;
            Q2_max_val = 1.0;
        }
        else if (q == 2) {
            Q2_min_val = 1.0;
            Q2_max_val = 1.4;
        }
        else if (q == 3) {
            Q2_min_val = 1.4;
            Q2_max_val = 2.0;
        }
        else if (q == 4) {
            Q2_min_val = 2.0;
            Q2_max_val = 3.0;
        }
        else {
            Q2_min_val = 3.0;
            Q2_max_val = 5.5;
        }

        TString namegraph = Form(
            "Background VS W_bin for Q^{2} in [%g, %g] GeV^{2}; W bin number; BG",
            Q2_min_val,
            Q2_max_val
        );

        TGraphErrors* graph_fit_bg = new TGraphErrors(
            intergral_background[q].size(),
            w_bg_fit[q].data(),
            intergral_background[q].data(),
            nullptr,
            intergral_background_error[q].data()
        );

        graph_fit_bg->SetTitle(namegraph);
        graph_fit_bg->SetMinimum(0);
        graph_fit_bg->SetMarkerSize(10.0);
        graph_fit_bg->GetXaxis()->SetNdivisions(20, kTRUE);

        int w_reper = cfg.w_fit_min;

        if (q == 4 || q == 5)
            w_reper = cfg.w_reper_q45;

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
        fitFunction->SetParameter(
            0,
            intergral_background[q][6]
        );

        fitFunction->SetParameter(1, 1);

        // Fit
        graph_fit_bg->Fit(fitFunction, "R");

        graph_fit_bg->Draw("AP");
        fitFunction->Draw("P SAME");

        graph_fit_bg->GetXaxis()->SetRangeUser(0, 20);
        graph_fit_bg->GetYaxis()->SetRangeUser(
            0,
            fitFunction->Eval(cfg.w_fit_max) * 1.1
        );

        // Fit information
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
            Form("   A = %.5f", fitFunction->GetParameter(0))
        );
        pave->AddText(
            Form("   B = %.5f", fitFunction->GetParameter(1))
        );

        pave->Draw("SAME");

        // Extrapolation
        for (int w = 0; w < cfg.bg_fit_start - w_reper; ++w)
        {
            intergral_background[q][w] =
                fitFunction->Eval(
                    w_reper + w + 0.01
                );
        }

        canvas_fit_bg->Update();
        canvas_fit_bg->Print(pdfFileName);

        delete pave;
        delete fitFunction;
        delete graph_fit_bg;
    }

    canvas_fit_bg->Print(pdfFileName + "]");
    delete canvas_fit_bg;
}

void yields_pdf(void) 
{
	gROOT->SetBatch(kTRUE); 

	std::vector<double> intergral_signal[6];
	std::vector<double> intergral_signal_4th_method_MM0[6];
	std::vector<double> intergral_signal_excut[6];
	std::vector<double> intergral_background[6];
	std::vector<double> intergral_MM0_background[6];
	std::vector<double> background_MM0[6];

	std::vector<double> intergral_signal_error[6];
	std::vector<double> intergral_signal_4th_method_MM0_error[6];
	std::vector<double> intergral_signal_excut_error[6];
	std::vector<double> intergral_background_error[6];
	std::vector<double> intergral_MM0_background_error[6];
	std::vector<double> background_error_MM0[6];
	
	std::vector<double> intergral_signal_old[6];
	std::vector<double> intergral_signal_error_old[6];

	std::vector<double> w_vector[6];
	std::vector<double> w_vector_error[6];

	std::vector<double> w_vector_MM0[6];
	std::vector<double> w_vector_MM0_error[6];

	std::vector<double> w_bg_fit[6];
	std::vector<double> w_bg_fit_MM0[6];

    TFile *infile_exp = TFile::Open("interim/mm_pip_hists_m_pip.root");
    TFile *infile_ann_good = TFile::Open("interim/mm_pip_hists_m_0.root");

    TH1F *MM_raw[12][100];
    TH1F *MM_my_4th_curve[12][100];
    TH1F *MM_ann_good_4th_curve[12][100];

    for (int q = 0; q < 6; ++q)
    {
        for (int w = 0; w < 100; ++w)
        {
            TString mm_raw_name = Form("MM_Q2_bin=%d_W_bin=%d;1", q+1, w+1);
            MM_raw[q][w] = (TH1F*)infile_exp->Get(mm_raw_name);

            TString mm_ann_name = Form("MMpiplus_from_MM0_topology_Q2_bin=%d_W_bin=%d;1", q+1, w+1);
            MM_ann_good_4th_curve[q][w] = (TH1F*)infile_ann_good->Get(mm_ann_name);

            MM_my_4th_curve[q][w] = new TH1F(
                Form("MMpiplus_with_zero_topology_BG_Q2_bin=%d_W_bin=%d", q + 1, w + 1),
                "",
                100, -0.3, 0.7
            );

            MM_my_4th_curve[q][w]->SetDirectory(nullptr);
        }
    }

/////////////////////////////////////////////////////////////////////// Background fit/////////////////////////////////////////////////////

    BackgroundFitConfig missing_pip_bg_cfg;

    CalculateBackgroundIntegral(
        MM_raw,
        intergral_background,
        intergral_background_error,
        w_bg_fit,
        missing_pip_bg_cfg
    );

/////////////////////////////////////////////////////////////////////////////////GRAPH BG/////////////////////////////////////////////////////////////////////////////

    TCanvas *canvas_fit_bg = new TCanvas("canvas_fit_bg", "Fit bg", 800, 600);
    TString pdfFileName_bg = Form("results/extrapolating_bg.pdf");  
    canvas_fit_bg -> Print(pdfFileName_bg + "[");

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
        sprintf(namegraph, "Background VS W_bin for Q^{2} in [%g, %g] GeV^{2}; W bin number; BG", Q2_min_val, Q2_max_val);
	
        TGraphErrors *graph_fit_bg = new TGraphErrors(intergral_background[q].size(), &w_bg_fit[q][0], &intergral_background[q][0], NULL, &intergral_background_error[q][0]);
        graph_fit_bg -> SetTitle(namegraph);
        graph_fit_bg -> SetMinimum(0);
        graph_fit_bg -> SetMarkerSize(10.0);

        int w_reper = 8;
        if(q == 4 || q == 5)    w_reper = 3;

        TF1 *fitFunction = new TF1("fitFunction", "[0]*(x-8)*(x-8)+[1]*(x-8)", w_reper - 0.1, 17);
        if(q == 4 || q == 5)    fitFunction = new TF1("fitFunction", "[0]*(x-3)*(x-3)+[1]*(x-3)", w_reper - 0.1, 17);
        

        // Устанавливаем начальные значения параметров
        fitFunction -> SetParameter(0, intergral_background[q][6]);  
        fitFunction -> SetParameter(1, 1);   
        //fitFunction->SetParameter(2, 1);  

        graph_fit_bg -> GetXaxis() -> SetNdivisions(20, kTRUE);
        graph_fit_bg -> Fit(fitFunction, "R"); // "R" означает, что мы хотим использовать график как "range" для подгонки
        graph_fit_bg -> Draw("AP");

        fitFunction -> Draw("P SAME");

        graph_fit_bg -> GetXaxis() -> SetRangeUser(0, 20); 
        graph_fit_bg -> GetYaxis() -> SetRangeUser(0, fitFunction -> Eval(17) * 1.1);  

        TPaveText* pave = new TPaveText(0.7, 0.4, 0.9, 0.6, "NDC");
        pave -> SetFillColor(0);
        pave -> SetTextAlign(12);
        pave -> SetTextSize(0.03);

        if(q !=4 && q != 5) pave->AddText(Form("A*(x - 8)^{2}+B*(x-8)"));        
        if(q==4 || q==5)    pave->AddText(Form("A*(x - 3)^{2}+B*(x - 3)"));

        pave -> AddText(Form("Fit results:"));
        pave -> AddText(Form("   A =  %.5f ", fitFunction->GetParameter(0)));
        pave -> AddText(Form("   B =  %.5f ", fitFunction->GetParameter(1)));
        pave -> Draw("SAME");

        for (int w = 0; w < 13 - w_reper; ++w)                                 // 4,5,6,7,8,9,10,11,12  when w_reper = 4
    	intergral_background[q][w] = fitFunction -> Eval(w_reper + w + 0.01);	   // 8,9,10,11,12 when w_reper = 8       

        canvas_fit_bg -> Update();
        canvas_fit_bg -> Print(pdfFileName_bg);
    }

    canvas_fit_bg -> Print(pdfFileName_bg + "]"); // Закрываем текущий PDF-файл
    delete canvas_fit_bg;

/////////////////////////////////////////////////////////////////////// NEW LOGNORMAL Method/////////////////////////////////////////////////////

	int q_fin = 6;
	int w_brink = 64;
    
	for (int q = 0; q < q_fin; ++q) 
    {
	    // Создаем отдельный PDF-файл для каждого значения "w"
        TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 800, 600);
        TString pdfFileName = Form("results/LogNormalFit_in_Q2_bin_%d.pdf", q+1);  // Генерируем имя файла на основе q
        canvas -> Print(pdfFileName + "[");

	    double Q2_min_val;
	    double Q2_max_val;

	    double del_q2;
	    float q2_mid_val;

	    if(q == 0)  {Q2_min_val = 0.5, Q2_max_val = 0.7, del_q2 = 0.2, q2_mid_val = 0.6;}
	    if(q == 1)  {Q2_min_val = 0.7, Q2_max_val = 1, del_q2 = 0.3, q2_mid_val = 0.85;}
	    if(q == 2)  {Q2_min_val = 1, Q2_max_val = 1.4, del_q2 = 0.4, q2_mid_val = 1.2;}
	    if(q == 3)  {Q2_min_val = 1.4, Q2_max_val = 2, del_q2 = 0.6, q2_mid_val = 1.7;}
	    if(q == 4)  {Q2_min_val = 2, Q2_max_val = 3, del_q2 = 1, q2_mid_val = 2.5;}
	    if(q == 5)  {Q2_min_val = 3, Q2_max_val = 5.5, del_q2 = 2.5, q2_mid_val = 4.25;}

	    if(q == 5)  w_brink = 37;
	    if(q == 4)  w_brink = 56;

	    double del_w = 0.025;

	    for (int w = 0; w < w_brink; ++w) 
        {
	        // char namehist_raw[256];
	        // sprintf(namehist_raw, "MM_Q2_bin=%d_W_bin=%d;1", q+1, w+1);

	        // MM_raw[q][w] = (TH1F*)infile_exp -> Get(Form(namehist_raw));
	        MM_raw[q][w] -> SetStats(kFALSE);  // Don't display stats

	        float w_min_val = 1.4 + w * 0.025;
            float w_max_val = 1.4 + (w + 1) * 0.025;
            float w_mid_val = (w_min_val + w_max_val) / 2;

            float M_prot = 0.938;
            float E_beam = 6.535;
            float E_prime = E_beam - (w_mid_val * w_mid_val - M_prot * M_prot + q2_mid_val) / (2 * M_prot);
            float nu = E_beam - E_prime;
            float sin2 = q2_mid_val / (4 * E_beam * E_prime);
            float cos2 = 1 - sin2;
            float tg2 = sin2 / cos2;
            float eps = 1.0 / (1.0 + 2.0 * (1 + nu * nu / q2_mid_val) * tg2);
            float alpha = 1.0 / 137.0;

            float flux = alpha / (4 * 3.1415) * 1 / (M_prot * M_prot * E_beam * E_beam) * (w_mid_val * (w_mid_val * w_mid_val - M_prot * M_prot)) / ((1 - eps) * q2_mid_val);
            float delta_flux = sqrt(pow(alpha / (4 * 3.1415) * 1 / (M_prot * M_prot * E_beam * E_beam) * ((3 * w_mid_val * w_mid_val - M_prot * M_prot) / (2 * q2_mid_val) - 1), 2) * pow(del_w, 2) + pow(alpha / (4 * 3.1415) * 1 / (M_prot * M_prot * E_beam * E_beam) * ((w_mid_val * w_mid_val - M_prot * M_prot) / (2 * q2_mid_val)), 2) * pow(del_q2, 2));

            char namelabel_raw[256];
            sprintf(namelabel_raw, "LogNorm_method_Q2_[%g,%g]GeV^2_W_in_[%g,%g]GeV; MM_X^2, GeV^2", Q2_min_val, Q2_max_val, w_min_val, w_max_val);

            MM_raw[q][w] -> SetTitle(namelabel_raw);
            MM_raw[q][w] -> SetLineColor(4);
	        MM_raw[q][w] -> SetLineWidth(1);

	        double h_raw = MM_raw[q][w] -> GetMaximum();
	        double mean = 0.0196;
            double sigma = 0.03;
            double tail = 0.02;   
	        double b_raw = MM_raw[q][w] -> GetBinContent(25) * 0.3;

	        double xmin, xmax;
	        xmin = -0.2;

	        if(w < 10)      xmax = 0.15;
	        else if(w < 17) xmax = 0.25;
	        else            xmax = 0.4;	        

	        MM_raw[q][w] -> GetYaxis() -> SetRangeUser(0, h_raw * 1.05);
            MM_raw[q][w] -> Draw();
            MM_raw[q][w] -> SetLineWidth(3);

            w_vector[q].push_back(w_mid_val);
            w_vector_error[q].push_back(del_w);

            double N, N_ex_cut, N_old;
            int w_reper = 8;

            if(q == 4 || q == 5)    w_reper = 3;

            if(w < w_reper) 
            {
                N = MM_raw[q][w] -> Integral(25, 40);
                N_old = MM_raw[q][w] -> Integral(25, 40);

                intergral_signal[q].push_back(N / del_q2 / del_w / flux);
                intergral_signal_error[q].push_back(N / del_q2 / del_w / flux * sqrt(pow(sqrt(N) / N, 2) + pow(delta_flux / flux, 2)));

                intergral_signal_old[q].push_back(N_old / del_q2 / del_w / flux);
                intergral_signal_error_old[q].push_back(N_old / del_q2 / del_w / flux * sqrt(pow(sqrt(N_old) / N_old, 2) + pow(delta_flux / flux, 2)));

                intergral_signal_excut[q].push_back(N / del_q2 /del_w / flux);
                intergral_signal_excut_error[q].push_back(N / del_q2 / del_w / flux * sqrt(pow(sqrt(N) / N, 2) + pow(delta_flux / flux, 2)));

                TLine *l_1 = new TLine(-0.05, 0, -0.05, MM_raw[q][w]->GetMaximum());
	            TLine *l_2 = new TLine(0.1, 0, 0.1, MM_raw[q][w]->GetMaximum());
	            l_1 -> SetLineColor(kBlack);
	            l_2 -> SetLineColor(kBlack);
	            l_1 -> Draw("SAME");
	            l_2 -> Draw("SAME");
            }

            else if(w >= w_reper && w < 13 )
            {
                N = MM_raw[q][w] -> Integral(25, 40) - intergral_background[q][w-w_reper];
                N_ex_cut = MM_raw[q][w] -> Integral(25,40);
                N_old = MM_raw[q][w] -> Integral(25,40);

                intergral_signal[q].push_back(N / del_q2 / del_w / flux);
                intergral_signal_error[q].push_back(N / del_q2 / del_w / flux * sqrt(pow(sqrt(N) / N, 2) + pow(delta_flux / flux, 2)));

                intergral_signal_excut[q].push_back(N_ex_cut / del_q2 / del_w / flux);
                intergral_signal_excut_error[q].push_back(N_ex_cut / del_q2 / del_w / flux * sqrt(pow(sqrt(N_ex_cut) / N_ex_cut, 2) + pow(delta_flux / flux, 2)));

                intergral_signal_old[q].push_back(N_old / del_q2 /del_w /flux);
                intergral_signal_error_old[q].push_back(N_old / del_q2 / del_w / flux * sqrt(pow(sqrt(N_old) / N_old, 2) + pow(delta_flux / flux, 2)));

                TLine *l_1 = new TLine(-0.05, 0, -0.05, MM_raw[q][w] -> GetMaximum());
	            TLine *l_2 = new TLine(0.1, 0, 0.1, MM_raw[q][w] -> GetMaximum());
	            l_1 -> SetLineColor(kBlack);
	            l_2 -> SetLineColor(kBlack);
	            l_1 -> Draw("SAME");
	            l_2 -> Draw("SAME");
            }
 
            else 
            {
	            N_ex_cut = MM_raw[q][w] -> Integral(25, 40);
                intergral_signal_excut[q].push_back(N_ex_cut / del_q2 / del_w / flux);
                intergral_signal_excut_error[q].push_back(N_ex_cut / del_q2 / del_w / flux * sqrt(pow(sqrt(N_ex_cut) / N_ex_cut, 2) + pow(delta_flux / flux, 2)));

                TLine *l_1 = new TLine(-0.05, 0, -0.05, MM_raw[q][w] -> GetMaximum());
	            TLine *l_2 = new TLine(0.1, 0, 0.1, MM_raw[q][w] -> GetMaximum());
	            l_1 -> SetLineColor(kBlack);
	            l_2 -> SetLineColor(kBlack);
	            l_1 -> Draw("SAME");
	            l_2 -> Draw("SAME");

                TLine *l_3 = new TLine((MM_raw[q][w] -> GetMaximumBin()) * 0.01 - 0.3, 0, (MM_raw[q][w] -> GetMaximumBin()) * 0.01 - 0.3, MM_raw[q][w] -> GetMaximum());
                l_3 -> SetLineColor(kBlack);
	            l_3 -> Draw("SAME");
                l_3 -> SetLineStyle(kDashed);

                TF1 *fitfunc = new TF1("fitfunc", combined, xmin, xmax, 7);

                // Set initial parameters for the background function
                fitfunc -> SetParameter(0, b_raw);
                //fitfunc->SetParLimits(0,b_raw*0.1,b_raw*10);
                fitfunc -> SetParameter(1, -1);
                //fitfunc->SetParLimits(1,-0.5, -2);
                fitfunc -> SetParameter(2, 0.3);
                //fitfunc->SetParLimits(2,0.15,0.4);

                // Set initial parameters for the peak function
                fitfunc -> SetParameter(3, h_raw);
                fitfunc -> SetParameter(4, mean);
                fitfunc -> SetParameter(5, sigma);   
                fitfunc -> SetParameter(6, tail);
                //fitfunc->SetParLimits(6,0.5,2);
            
                // Fit the histogram
                MM_raw[q][w] -> GetYaxis() -> SetRangeUser(0, h_raw * 1.05);
                MM_raw[q][w] -> Fit(fitfunc, "R");

                // Draw the histogram and fits
                fitfunc -> SetLineColor(kMagenta);
                fitfunc -> SetLineWidth(2);
                fitfunc -> Draw("SAME");

                TF1 *peakfunc = new TF1("peakfunc", peak, xmin, xmax, 4);
                peakfunc -> SetParameters(fitfunc -> GetParameter(3), fitfunc -> GetParameter(4), fitfunc -> GetParameter(5), fitfunc -> GetParameter(6));
                peakfunc -> SetLineColor(kRed);
                peakfunc -> SetLineWidth(2);
                peakfunc -> Draw("SAME");

                TF1 *bg = new TF1("background", background, xmin, xmax, 3);
                bg -> SetParameters(fitfunc -> GetParameter(0), fitfunc -> GetParameter(1), fitfunc -> GetParameter(2));
                bg -> SetLineColor(kBlack);
                bg -> SetLineStyle(kDashed);
                bg -> SetLineWidth(2);
                bg -> Draw("SAME");

                TPaveText* pave = new TPaveText(0.45, 0.6, 0.88, 0.9, "NDC");
                pave -> SetFillColor(0);
                pave -> SetTextAlign(12);
                pave -> SetTextSize(0.03);
                pave -> AddText(Form("   Total events: %.0f ", MM_raw[q][w] -> Integral()));
                pave -> AddText(Form("Fit results:"));
                pave -> AddText(Form("   Amplitude bg: %.5f +/- %.5f", bg -> GetParameter(0), bg -> GetParError(0)));
                pave -> AddText(Form("   Mean bg: %.5f +/- %.5f", bg -> GetParameter(1), bg -> GetParError(1)));
                pave -> AddText(Form("   Sigma bg: %.5f +/- %.5f", bg -> GetParameter(2), bg -> GetParError(2)));
                pave -> AddText(Form("   Amplitude signal: %.5f +/- %.5f", fitfunc -> GetParameter(3), fitfunc -> GetParError(3)));
                pave -> AddText(Form("   Mean signal: %.5f +/- %.5f", fitfunc -> GetParameter(4), fitfunc -> GetParError(4)));
                pave -> AddText(Form("   Sigma signal: %.5f +/- %.5f", fitfunc -> GetParameter(5), fitfunc -> GetParError(5)));
                pave -> AddText(Form("   Rad. tail: %.5f +/- %.5f", fitfunc -> GetParameter(6), fitfunc -> GetParError(6)));
                pave -> Draw("SAME");

                TLegend* legend = new TLegend(0.6, 0.4, 0.9, 0.6);
                legend -> SetTextSize(0.02);
                legend -> AddEntry(MM_raw[q][w], "Pass2 Data", "l");
                legend -> AddEntry(fitfunc, "Overall fit", "l");
                legend -> AddEntry(peakfunc, "Signal (gaus+rad.tail)", "l");
                legend -> AddEntry(bg, "Background (lognormal)", "l");
                legend -> Draw("SAME");

	            N = MM_raw[q][w] -> Integral(25,40) - bg -> Integral(-0.05, 0.1) * 100;
	            N_old = peakfunc -> Integral(xmin, xmax) * 100;

                intergral_signal[q].push_back(N / del_q2 / del_w / flux);
                intergral_signal_error[q].push_back(N / del_q2 / del_w / flux * sqrt(pow(sqrt(N) / N, 2) + pow(delta_flux / flux, 2)));

                intergral_signal_old[q].push_back(N_old / del_q2 / del_w / flux);
                intergral_signal_error_old[q].push_back(N_old / del_q2 / del_w / flux * sqrt(pow(sqrt(N_old) / N_old, 2) + pow(delta_flux / flux, 2)));
            }

	        canvas -> Update();
            canvas -> Print(pdfFileName);

            // delete MM_raw[q][w];
        }

        canvas -> Print(pdfFileName + "]"); // Закрываем текущий PDF-файл
        delete canvas;
    }

/////////////////////////////////////////////////////////////////////// Background fit for MM0 case //////////////////////////////////////////////////////////////////////////////

    BackgroundFitConfig fully_exclusive_bg_cfg;

    CalculateBackgroundIntegral(
        MM_ann_good_4th_curve,
        background_MM0,
        background_error_MM0,
        w_bg_fit_MM0,
        fully_exclusive_bg_cfg
    );

//////////////////////////////////////////////////////////////////////////////////////Graph BG for MM0 case////////////////////////////////////////////////////////////////////////////////////////

    TCanvas *canvas_fit_bg_MM0 = new TCanvas("canvas_fit_bg_MM0", "Fit bg MM0", 800, 600);
    TString pdfFileName_bg_MM0 = Form("results/extrapolating_bg_MM0.pdf");  
    canvas_fit_bg_MM0->Print(pdfFileName_bg_MM0 + "[");

    for (int q = 0; q < 6; ++q) 
    {
	    double Q2_min_val, Q2_max_val;

	    if(q == 0)  {Q2_min_val = 0.5, Q2_max_val = 0.7;}
	    if(q == 1)  {Q2_min_val = 0.7, Q2_max_val = 1;}
	    if(q == 2)  {Q2_min_val = 1, Q2_max_val = 1.4;}
	    if(q == 3)  {Q2_min_val = 1.4, Q2_max_val = 2;}
	    if(q == 4)  {Q2_min_val = 2, Q2_max_val = 3;}
	    if(q == 5)  {Q2_min_val = 3, Q2_max_val = 5.5;}

        char namegraph_MM0[256];
        sprintf(namegraph_MM0, "Background from MM0 topology VS W_bin for Q2 in [%g,%g] GeV^2; W bin number; BG", Q2_min_val, Q2_max_val);
	
        TGraphErrors *graph_fit_bg_MM0 = new TGraphErrors(background_MM0[q].size(), &w_bg_fit_MM0[q][0], &background_MM0[q][0], NULL, &background_error_MM0[q][0]);
        graph_fit_bg_MM0 -> SetTitle(namegraph_MM0);
        graph_fit_bg_MM0 -> SetMinimum(0);
        graph_fit_bg_MM0 -> SetMarkerSize(10.0);

        int w_reper = 8;

        if(q==4 || q==5)    w_reper = 3;      

        TF1 *fitFunction = new TF1("fitFunction", "[0]*(x-8)*(x-8)+[1]*(x-8)", w_reper - 0.1, 17);

        if(q==4 || q==5)    fitFunction = new TF1("fitFunction", "[0]*(x-3)*(x-3)+[1]*(x-3)", w_reper-0.1, 17);
        
        // Устанавливаем начальные значения параметров
        fitFunction -> SetParameter(0, background_MM0[q][6]);  
        fitFunction -> SetParameter(1, 1);   
        //fitFunction->SetParameter(2, 1);  

        graph_fit_bg_MM0 -> GetXaxis() -> SetNdivisions(20, kTRUE);
        graph_fit_bg_MM0 -> Fit(fitFunction, "R"); // "R" означает, что мы хотим использовать график как "range" для подгонки
        graph_fit_bg_MM0 -> Draw("AP");

        fitFunction -> Draw("P SAME");

        graph_fit_bg_MM0 -> GetXaxis() -> SetRangeUser(0, 20); 
        graph_fit_bg_MM0 -> GetYaxis() -> SetRangeUser(0, fitFunction -> Eval(17) * 1.1);  

        TPaveText* pave = new TPaveText(0.7, 0.4, 0.9, 0.6, "NDC");
        pave -> SetFillColor(0);
        pave -> SetTextAlign(12);
        pave -> SetTextSize(0.03);

        if(q != 4 && q != 5)    pave -> AddText(Form("A*(x - 8)^2+B*(x-8)"));        
        if(q == 4 || q == 5)    pave->AddText(Form("A*(x - 3)^2+B*(x - 3)"));
        
        pave -> AddText(Form("Fit results:"));
        pave -> AddText(Form("   A =  %.5f ", fitFunction -> GetParameter(0)));
        pave -> AddText(Form("   B =  %.5f ", fitFunction -> GetParameter(1)));
        pave -> Draw("SAME");

        for (int w = 0; w < 13 - w_reper; ++w)                             // 4,5,6,7,8,9,10,11,12  when w_reper = 4
    	background_MM0[q][w] = fitFunction -> Eval(w_reper + w + 0.01);	   // 8,9,10,11,12 when w_reper = 8
        
        canvas_fit_bg_MM0 -> Update();
        canvas_fit_bg_MM0 -> Print(pdfFileName_bg_MM0);
    }

    canvas_fit_bg_MM0 -> Print(pdfFileName_bg_MM0 + "]"); // Закрываем текущий PDF-файл
    delete canvas_fit_bg_MM0;

// ////////////////////////////////////////////////////////////////////////////////4th method - scaling MM0 BG/////////////////////////////////////////////////////////////////////

	q_fin = 6;
	w_brink = 64;

	for (int q = 0; q < q_fin; ++q) 
    {
	    // Создаем отдельный PDF-файл для каждого значения "w"
        TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 800, 600);
        TString pdfFileName = Form("results/MMpiplus_from_MM0_topology_fit_in_Q2_bin_%d.pdf", q+1);  // Генерируем имя файла на основе q
        canvas->Print(pdfFileName + "[");

	    double Q2_min_val;
	    double Q2_max_val;

	    double del_q2;
	    float q2_mid_val;

	    if(q == 0)  {Q2_min_val = 0.5, Q2_max_val = 0.7, del_q2 = 0.2, q2_mid_val = 0.6;}
	    if(q == 1)  {Q2_min_val = 0.7, Q2_max_val = 1, del_q2 = 0.3, q2_mid_val = 0.85;}
	    if(q == 2)  {Q2_min_val = 1, Q2_max_val = 1.4, del_q2 = 0.4, q2_mid_val = 1.2;}
	    if(q == 3)  {Q2_min_val = 1.4, Q2_max_val = 2, del_q2 = 0.6, q2_mid_val = 1.7;}
	    if(q == 4)  {Q2_min_val = 2, Q2_max_val = 3, del_q2 = 1, q2_mid_val = 2.5;}
	    if(q == 5)  {Q2_min_val = 3, Q2_max_val = 5.5, del_q2 = 2.5, q2_mid_val = 4.25;}

	    if(q == 5) w_brink = 37;
	    if(q == 4) w_brink = 56;

	    double del_w = 0.025;

	    for (int w = 0; w < w_brink; ++w) 
        {      
            char namehist_raw[256];
	        sprintf(namehist_raw, "MM_Q2_bin=%d_W_bin=%d;1", q + 1, w + 1);
	
	        MM_my_4th_curve[q][w] = (TH1F*)infile_exp -> Get(Form(namehist_raw));
	        MM_my_4th_curve[q][w] -> SetStats(kFALSE);  // Don't display stats

	        MM_ann_good_4th_curve[q][w] -> SetStats(kFALSE);  // Don't display stats

	        double scale_coef = MM_my_4th_curve[q][w] -> Integral() / MM_ann_good_4th_curve[q][w] -> Integral();

	        float w_min_val = 1.4 + w * 0.025;
            float w_max_val = 1.4 + (w + 1) * 0.025;
            float w_mid_val = (w_min_val + w_max_val) / 2;

            float M_prot = 0.938;
            float E_beam = 6.535;
            float E_prime = E_beam - (w_mid_val * w_mid_val - M_prot * M_prot + q2_mid_val) / (2 * M_prot);
            float nu = E_beam - E_prime;
            float sin2 = q2_mid_val / (4 * E_beam * E_prime);
            float cos2 = 1 - sin2;
            float tg2 = sin2 / cos2;
            float eps = 1.0 / (1.0 + 2.0 * (1 + nu * nu / q2_mid_val) * tg2);
            float alpha = 1.0 / 137.0;

            float flux = alpha / (4 * 3.1415) * 1 / (M_prot * M_prot * E_beam * E_beam) * (w_mid_val * (w_mid_val * w_mid_val - M_prot * M_prot)) / ((1 - eps) * q2_mid_val);
            float delta_flux = sqrt(pow(alpha / (4 * 3.1415) * 1 / (M_prot * M_prot * E_beam * E_beam) * ((3 * w_mid_val * w_mid_val - M_prot * M_prot) / (2 * q2_mid_val) - 1), 2) * pow(del_w, 2) + pow(alpha / (4 * 3.1415) * 1 / (M_prot * M_prot * E_beam * E_beam) * ((w_mid_val * w_mid_val - M_prot * M_prot) / (2 * q2_mid_val)), 2) * pow(del_q2, 2));

            double h_raw = MM_ann_good_4th_curve[q][w] -> GetMaximum();
	        double mean = 0.0196;
            double sigma = 0.05;
            double tail = 0.02;
            double b_raw = MM_ann_good_4th_curve[q][w] -> GetBinContent(50);
            double xmin = -0.1;
            double xmax = 0.3;

	        if(w<10)        xmax = 0.15;
	        else if(w<17)   xmax = 0.25;
	        
            char namelabel_ann[256];
            sprintf(namelabel_ann, "MMpiplus_from_MM0_Q2_[%g,%g]GeV^2_W_in_[%g,%g]GeV; MM_X^2, GeV^2", Q2_min_val, Q2_max_val, w_min_val, w_max_val);

            MM_ann_good_4th_curve[q][w] -> SetTitle(namelabel_ann);
            MM_ann_good_4th_curve[q][w] -> SetLineColor(4);
	        MM_ann_good_4th_curve[q][w] -> SetLineWidth(1);

	        MM_ann_good_4th_curve[q][w] -> GetYaxis() -> SetRangeUser(0, h_raw * 1.05);
            MM_ann_good_4th_curve[q][w] -> Draw();
            MM_ann_good_4th_curve[q][w] -> SetLineWidth(3);

            w_vector_MM0[q].push_back(w_mid_val);
            w_vector_MM0_error[q].push_back(del_w);

            double N;
            int w_reper = 8;

            if(q == 4 || q == 5)    w_reper = 3;
            
            if(w < w_reper) 
            {
                N = MM_my_4th_curve[q][w] -> Integral(25, 40);

                intergral_signal_4th_method_MM0[q].push_back(N / del_q2 / del_w / flux);
                intergral_signal_4th_method_MM0_error[q].push_back(N / del_q2 / del_w / flux * sqrt(pow(sqrt(N) / N, 2) + pow(delta_flux / flux, 2)));
    
                TLine *l_1 = new TLine(-0.05, 0, -0.05, MM_ann_good_4th_curve[q][w] -> GetMaximum());
	            TLine *l_2 = new TLine(0.1, 0, 0.1, MM_ann_good_4th_curve[q][w] -> GetMaximum());
	            l_1 -> SetLineColor(kBlack);
	            l_2 -> SetLineColor(kBlack);
	            l_1 -> Draw("SAME");
	            l_2 -> Draw("SAME");
            }

            else if(w >= w_reper && w < 13) 
            {
                N = MM_my_4th_curve[q][w] -> Integral(25, 40) - background_MM0[q][w-w_reper] * scale_coef;

                intergral_signal_4th_method_MM0[q].push_back(N / del_q2 / del_w / flux);
                intergral_signal_4th_method_MM0_error[q].push_back(N / del_q2 / del_w / flux * sqrt(pow(sqrt(N) / N, 2) + pow(delta_flux / flux, 2)));

                TLine *l_1 = new TLine(-0.05, 0, -0.05, MM_ann_good_4th_curve[q][w] -> GetMaximum());
	            TLine *l_2 = new TLine(0.1, 0, 0.1, MM_ann_good_4th_curve[q][w] -> GetMaximum());
	            l_1 -> SetLineColor(kBlack);
	            l_2 -> SetLineColor(kBlack);
	            l_1 -> Draw("SAME");
	            l_2 -> Draw("SAME");
            }
  
            else 
            {
                TLine *l_1 = new TLine(-0.05, 0, -0.05, MM_ann_good_4th_curve[q][w] -> GetMaximum());
	            TLine *l_2 = new TLine(0.1, 0, 0.1, MM_ann_good_4th_curve[q][w] -> GetMaximum());
	            l_1 -> SetLineColor(kBlack);
	            l_2 -> SetLineColor(kBlack);
	            l_1 -> Draw("SAME");
	            l_2 -> Draw("SAME");

                TLine *l_3 = new TLine((MM_ann_good_4th_curve[q][w] -> GetMaximumBin()) * 0.01 - 0.3, 0, (MM_ann_good_4th_curve[q][w] -> GetMaximumBin()) * 0.01 - 0.3, MM_ann_good_4th_curve[q][w] -> GetMaximum());
                l_3 -> SetLineColor(kBlack);
	            l_3 -> Draw("SAME");
                l_3 -> SetLineStyle(kDashed);

                TF1 *fitfunc = new TF1("fitfunc", combined, xmin, xmax, 7);

                // Set initial parameters for the background function
                fitfunc -> SetParameter(0, b_raw);
                //fitfunc->SetParLimits(0,b_raw*0.1,b_raw*10);
                fitfunc -> SetParameter(1, -1);
                //fitfunc->SetParLimits(1,-0.5, -2);
                fitfunc -> SetParameter(2, 0.3);
                //fitfunc->SetParLimits(2,0.15,0.4);

                // Set initial parameters for the peak function
                fitfunc -> SetParameter(3, h_raw);
                fitfunc -> SetParameter(4, mean);
                fitfunc -> SetParameter(5, sigma);   
                fitfunc -> SetParameter(6, tail);
                //fitfunc->SetParLimits(6,0.5,2);

                // Fit the histogram
                MM_ann_good_4th_curve[q][w] -> GetYaxis() -> SetRangeUser(0, h_raw * 1.05);
                MM_ann_good_4th_curve[q][w] -> Fit(fitfunc, "R");             //// We are fitting the respective histogram from MMO topology! And using its BG to subtract from our topology!!!! 

                // Draw the histogram and fits
                fitfunc -> SetLineColor(kMagenta);
                fitfunc -> SetLineWidth(2);
                fitfunc -> Draw("SAME");
    
                TF1 *peakfunc = new TF1("peakfunc", peak, xmin, xmax, 4);
                peakfunc -> SetParameters(fitfunc -> GetParameter(3), fitfunc -> GetParameter(4), fitfunc -> GetParameter(5), fitfunc -> GetParameter(6));
                peakfunc -> SetLineColor(kRed);
                peakfunc -> SetLineWidth(2);
                peakfunc -> Draw("SAME");

                TF1 *bg = new TF1("background", background, xmin, xmax, 3);
                bg -> SetParameters(fitfunc -> GetParameter(0), fitfunc -> GetParameter(1), fitfunc -> GetParameter(2));
                bg -> SetLineColor(kBlack);
                bg -> SetLineStyle(kDashed);
                bg -> SetLineWidth(2);
                bg -> Draw("SAME");

                TPaveText* pave = new TPaveText(0.45, 0.6, 0.88, 0.9, "NDC");
                pave -> SetFillColor(0);
                pave -> SetTextAlign(12);
                pave -> SetTextSize(0.03);
                pave -> AddText(Form("   Total events: %.0f ", MM_ann_good_4th_curve[q][w] -> Integral()));
                pave -> AddText(Form("Fit results:"));
                pave -> AddText(Form("   Amplitude bg: %.5f +/- %.5f", bg -> GetParameter(0), bg -> GetParError(0)));
                pave -> AddText(Form("   Mean bg: %.5f +/- %.5f", bg -> GetParameter(1), bg -> GetParError(1)));
                pave -> AddText(Form("   Sigma bg: %.5f +/- %.5f", bg -> GetParameter(2), bg -> GetParError(2)));
                pave -> AddText(Form("   Amplitude signal: %.5f +/- %.5f", fitfunc -> GetParameter(3), fitfunc -> GetParError(3)));
                pave -> AddText(Form("   Mean signal: %.5f +/- %.5f", fitfunc -> GetParameter(4), fitfunc -> GetParError(4)));
                pave -> AddText(Form("   Sigma signal: %.5f +/- %.5f", fitfunc -> GetParameter(5), fitfunc -> GetParError(5)));
                pave -> AddText(Form("   Rad. tail: %.5f +/- %.5f", fitfunc -> GetParameter(6), fitfunc -> GetParError(6)));
                pave -> Draw("SAME");

                TLegend* legend = new TLegend(0.6, 0.4, 0.9, 0.6);
                legend -> SetTextSize(0.02);
                legend -> AddEntry(MM_ann_good_4th_curve[q][w], "Pass2 Data", "l");
                legend -> AddEntry(fitfunc, "Overall fit", "l");
                legend -> AddEntry(peakfunc, "Signal (gaus+rad.tail)", "l");
                legend -> AddEntry(bg, "Background (lognormal)", "l");
                legend -> Draw("SAME");

	            N = MM_my_4th_curve[q][w] -> Integral(25, 40) - (bg -> Integral(-0.05, 0.1) * 100) * scale_coef;

                intergral_signal_4th_method_MM0[q].push_back(N / del_q2 / del_w / flux);
                intergral_signal_4th_method_MM0_error[q].push_back(N / del_q2 / del_w / flux * sqrt(pow(sqrt(N) / N, 2) + pow(delta_flux / flux, 2)));
            }

            canvas -> Update();
            canvas -> Print(pdfFileName);

            delete MM_my_4th_curve[q][w];
        }

        canvas -> Print(pdfFileName + "]"); // Закрываем текущий PDF-файл
        delete canvas;
    }

/////////////////////////////////////////////////////////////////////////////////GRAPHS/////////////////////////////////////////////////////////////////////////////

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

        TGraphErrors *graph_signal_old = new TGraphErrors(intergral_signal_old[q].size(), &w_vector[q][0], &intergral_signal_old[q][0], NULL, &intergral_signal_error_old[q][0]);
        graph_signal_old -> SetTitle(namegraph);
        graph_signal_old -> SetMinimum(0);
        graph_signal_old -> SetMarkerSize(1);
        graph_signal_old -> SetLineColor(4);
        graph_signal_old -> Draw("SAME");

        //mg->Add(graph_signal_old);

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
        //legend->AddEntry(graph_signal_old, "Old lognorm method", "l");

        legend->Draw("SAME");

        canvas_yields -> Update();
        canvas_yields -> Print(pdfFileName_graph);
    }

    canvas_yields -> Print(pdfFileName_graph + "]"); // Закрываем текущий PDF-файл
    delete canvas_yields;

    gSystem -> Exit(0);
}