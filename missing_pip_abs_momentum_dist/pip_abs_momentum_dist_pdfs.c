void pip_abs_momentum_dist_pdfs()
{
    double n_bins = 100, first_bin_val = 0, last_bin_val = 2;
    double step = (last_bin_val - first_bin_val) / n_bins;
    float ME_cuts_08[6][64] = {0};
    float ME_cuts_09[6][64] = {0};

    TH1F *momentum_dist[6][64];
	for (int q = 0; q < 6; ++q) 
	    for (int w = 0; w < 64; ++w)
		{
			char namehist[256];
			sprintf(namehist, "momentum_dist_Q2_bin=%d_W_bin=%d", q + 1, w + 1); 
			momentum_dist[q][w] = new TH1F (namehist, namehist, n_bins, first_bin_val, last_bin_val);
		}

	TFile *infile = new TFile("pip_abs_momentum_dist_hists.root");

    int q_fin = 6;
	int w_brink = 64;
	
	for (int q = 0; q < q_fin; ++q)
    { 
		TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 800, 600);
        TString pdfFileName = Form("different_q2_bins/pip_abs_momentum_dist_q2_bin_%d.pdf", q+1);
        canvas -> Print(pdfFileName + "[");

        if(q == 5)  w_brink = 37;
		if(q == 4)  w_brink = 56;

        double Q2_min, Q2_max;

		if      (q == 0)    {Q2_min = 0.5, Q2_max = 0.7;}
		else if (q == 1)    {Q2_min = 0.7, Q2_max = 1.0;}
		else if (q == 2)    {Q2_min = 1.0, Q2_max = 1.4;}
		else if (q == 3)    {Q2_min = 1.4, Q2_max = 2.0;}
		else if (q == 4)    {Q2_min = 2.0, Q2_max = 3.0;}
		else if (q == 5)    {Q2_min = 3.0, Q2_max = 5.5;}

        for (int w = 0; w < w_brink; ++w)
		{
			float W_min = 1.4 + w * 0.025;
			float W_max = 1.4 + (w + 1) * 0.025;

            char namehist[256];
			sprintf(namehist, "momentum_dist_Q2_bin=%d_W_bin=%d", q + 1, w + 1);

            momentum_dist[q][w] = (TH1F*) infile -> Get(Form("%s", namehist));

            char namelabel[256];
			sprintf(namelabel, "#pi^{+} absolute momentum distribution in COM, Q^{2} [%g, %g] GeV^{2}, W [%g, %g] GeV; |p_{#pi^{+}}|, GeV; N", Q2_min, Q2_max, W_min, W_max);
			momentum_dist[q][w] -> SetTitle(namelabel);

            double x_08, x_09;

            Double_t *integral = momentum_dist[q][w] -> GetIntegral();
            for (int i = 0; i <= momentum_dist[q][w] -> GetNbinsX(); i++) 
            {
                //std::cout << "Bin " << i << " cumulative probability = " << integral[i] << std::endl;

                if (integral[i] <= 0.8 && integral[i + 1] >= 0.8)
                {
                    if (0.8 - integral[i] < integral[i + 1] - 0.8)
                        x_08 = i;
                    else
                        x_08 = i + 1;
                }

                if (integral[i] <= 0.9 && integral[i + 1] >= 0.9)
                {
                    if (0.9 - integral[i] < integral[i + 1] - 0.9)
                        x_09 = i;
                    else
                        x_09 = i + 1;
                }
            }

            //std::cout << "x_08 = " << x_08 << ", x_09 = " << x_09 << std::endl;

            momentum_dist[q][w] -> Draw();

            float p_08 = first_bin_val + x_08 * step;
            float p_09 = first_bin_val + x_09 * step;
            float ME_08 = std::sqrt(p_08 * p_08 + 0.13957 * 0.13957);
            float ME_09 = std::sqrt(p_09 * p_09 + 0.13957 * 0.13957);

            TLine *l_1 = new TLine(p_08, 0, p_08, momentum_dist[q][w] -> GetMaximum() * 1.0485);
	        l_1 -> SetLineColor(kRed);
            l_1 -> Draw("SAME");

            TLine *l_2 = new TLine(p_09, 0, p_09, momentum_dist[q][w] -> GetMaximum() * 1.0485);
	        l_2 -> SetLineColor(kRed);
            l_2 -> Draw("SAME");

            TPaveText* pave = new TPaveText(0.55, 0.6, 0.75, 0.7, "NDC");
            pave -> SetFillColor(0);
            pave -> SetTextAlign(12);
            pave -> SetTextSize(0.03);

            pave -> AddText(Form("80% : |p| = %.2f GeV, ME_{cut} = %.3f GeV", p_08, ME_08));
            pave -> AddText(Form("90% : |p| = %.2f GeV, ME_{cut} = %.3f GeV", p_09, ME_09));
            pave -> Draw("SAME");

            canvas -> Update();
	        canvas -> Print(pdfFileName);

            ME_cuts_08[q][w] = ME_08;
            ME_cuts_09[q][w] = ME_09;

			delete momentum_dist[q][w];
        }

		canvas -> Print(pdfFileName + "]"); 
		delete canvas;
	}

    ofstream f1, f2;
    f1.open("ME_cuts_08.txt", ios::out);
    f2.open("ME_cuts_09.txt", ios::out);

    for (int q = 0; q < 6; ++q)
    {
        for (int w = 0; w < 64; ++w)
        {
            f1 << ME_cuts_08[q][w] << " ";
            f2 << ME_cuts_09[q][w] << " ";
        }
        f1 << endl;
        f2 << endl;
    }

    f1.close();
    f2.close();
    
    gSystem -> Exit(0);
}