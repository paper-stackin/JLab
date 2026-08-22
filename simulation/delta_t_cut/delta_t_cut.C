void delta_t_cut(void)
{
	gROOT->SetBatch(kTRUE); 

	TH2F *p_delta_t_vs_p_FTOF = new TH2F("p_delta_t_FTOF", "FTOF proton p VS #Deltat; p, GeV; #Deltat, ns", 100, 0, 5, 100, -10, 10);
    TH2F *p_delta_t_vs_p_CTOF = new TH2F("p_delta_t_CTOF", "CTOF proton p VS #Deltat; p, GeV; #Deltat, ns", 100, 0, 5, 100, -10, 10);
    TH2F *pim_delta_t_vs_p_FTOF = new TH2F("pim_delta_t_FTOF", "FTOF #pi^{-} p VS #Deltat; p, GeV; #Deltat, ns", 100, 0, 5, 100, -8, 8);
    TH2F *pim_delta_t_vs_p_CTOF = new TH2F("pim_delta_t_CTOF", "CTOF #pi^{-} p VS #Deltat; p, GeV; #Deltat, ns", 100, 0, 5, 100, -5, 5);

    TH1F *p_delta_t_FTOF[20];
	for(int i = 0; i < 20; ++i) 
    {
		char namehist[256];
		sprintf(namehist, "FTOF proton #Deltat, p [%.2f, %.2f] GeV; #Deltat, ns", 0.25*i, 0.25*(i+1));
		p_delta_t_FTOF[i] = new TH1F(namehist, namehist, 100, -1, 1);
	}

    TH1F *p_delta_t_CTOF[20];
	for(int i = 0; i < 20; ++i) 
    {
		char namehist[256];
		sprintf(namehist, "CTOF proton #Deltat, p [%.2f, %.2f] GeV; #Deltat, ns", 0.25*i, 0.25*(i+1));
		p_delta_t_CTOF[i] = new TH1F(namehist, namehist, 100, -1, 1);
	}

    TH1F *pim_delta_t_FTOF[20];
	for(int i = 0; i < 20; ++i) 
    {
		char namehist[256];
		sprintf(namehist, "FTOF #pi^{-} #Deltat, p [%.2f, %.2f] GeV; #Deltat, ns", 0.25*i, 0.25*(i+1));
		pim_delta_t_FTOF[i] = new TH1F(namehist, namehist, 100, -1, 1);
	}

    TH1F *pim_delta_t_CTOF[20];
	for(int i = 0; i < 20; ++i) 
    {
		char namehist[256];
		sprintf(namehist, "CTOF #pi^{-} #Deltat, p [%.2f, %.2f] GeV; #Deltat, ns", 0.25*i, 0.25*(i+1));
		pim_delta_t_CTOF[i] = new TH1F(namehist, namehist, 100, -1, 1);
	}

	TFile *infile_exp = TFile::Open("MPPT_events_simulation_selection_with_additional_values_saved_2.root"); // открываем рут файл с импульсами частиц

	TTreeReader reader("MMpiptree", infile_exp);
    TTreeReaderValue<Float_t> px_e(reader, "px_e");
    TTreeReaderValue<Float_t> py_e(reader, "py_e");
    TTreeReaderValue<Float_t> pz_e(reader, "pz_e");
    TTreeReaderValue<Float_t> time_e(reader, "time_e");
    TTreeReaderValue<Float_t> path_e(reader, "path_e");
    TTreeReaderValue<Int_t> status_e(reader, "status_e");

    TTreeReaderValue<Float_t> px_p(reader, "px_p");
    TTreeReaderValue<Float_t> py_p(reader, "py_p");
    TTreeReaderValue<Float_t> pz_p(reader, "pz_p");
    TTreeReaderValue<Float_t> time_p(reader, "time_p");
    TTreeReaderValue<Float_t> path_p(reader, "path_p");
    TTreeReaderValue<Int_t> status_p(reader, "status_p");

    TTreeReaderValue<Float_t> px_pim(reader, "px_pim");
    TTreeReaderValue<Float_t> py_pim(reader, "py_pim");
    TTreeReaderValue<Float_t> pz_pim(reader, "pz_pim");
    TTreeReaderValue<Float_t> time_pim(reader, "time_pim");
    TTreeReaderValue<Float_t> path_pim(reader, "path_pim");
    TTreeReaderValue<Int_t> status_pim(reader, "status_pim");

    TTreeReaderValue<Float_t> weight(reader, "weight");

    TLorentzVector el_final, pi_minus, p_final;
    TLorentzVector el_initial(0, 0, 6.535, 6.535);
	TLorentzVector p_initial(0, 0, 0, 0.938);

	/////////////////////////////////////////////////////////////////
	while(reader.Next())
    {
        //cout << reader.Next() << "\n";	      
        el_final.SetXYZM(*px_e, *py_e, *pz_e, 0.);
        p_final.SetXYZM(*px_p, *py_p, *pz_p, 0.938);
        pi_minus.SetXYZM(*px_pim, *py_pim, *pz_pim, 0.139);
      
        double p_mom = p_final.P();
        double pim_mom = pi_minus.P();

        double p_beta = p_final.Beta();
        double pim_beta = pi_minus.Beta();

        if(p_beta <= 0 || pim_beta <= 0) continue;

        double c_light = 29.9792458;

        double t_e = *time_e - *path_e / c_light;
        double p_delta_t = *path_p / p_beta / c_light - *time_p + t_e;
        double pim_delta_t = *path_pim / pim_beta / c_light - *time_pim + t_e;

        //cout << p_delta_t << "\n";
        //cout << pim_delta_t << "\n";
        //cout << *status_p << "\n";
        //cout << *status_pim << "\n";
    	
        if (2000 <= *status_p && *status_p < 4000) p_delta_t_vs_p_FTOF->Fill(p_mom, p_delta_t, *weight);
        if (4000 <= *status_p && *status_p < 8000) p_delta_t_vs_p_CTOF->Fill(p_mom, p_delta_t, *weight);

        if (2000 <= *status_pim && *status_pim < 4000) pim_delta_t_vs_p_FTOF->Fill(pim_mom, pim_delta_t, *weight);
        if (4000 <= *status_pim && *status_pim < 8000) pim_delta_t_vs_p_CTOF->Fill(pim_mom, pim_delta_t, *weight);

        int p_bin = p_mom / 0.25;
        int pim_bin = pim_mom / 0.25;
        if (p_bin > 19) p_bin = 19;
        if (pim_bin > 19) pim_bin = 19;

        if (2000 <= *status_p && *status_p < 4000) p_delta_t_FTOF[p_bin]->Fill(p_delta_t, *weight);
        if (4000 <= *status_p && *status_p < 8000) p_delta_t_CTOF[p_bin]->Fill(p_delta_t, *weight);
        if (2000 <= *status_p && *status_p < 4000) pim_delta_t_FTOF[p_bin]->Fill(pim_delta_t, *weight);
        if (4000 <= *status_p && *status_p < 8000) pim_delta_t_CTOF[p_bin]->Fill(pim_delta_t, *weight);
	}

	/////////////////////////////////////////////////////////////////////

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptFit(111);

    TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 800, 600);
    TString pdfFileName = Form("delta_t_cuts_simulation.pdf");  
    canvas->Print(pdfFileName + "[");

    p_delta_t_vs_p_FTOF->Draw("COLZ");
    canvas->Update();
	canvas->Print(pdfFileName);

    p_delta_t_vs_p_CTOF->Draw("COLZ");
    canvas->Update();
	canvas->Print(pdfFileName);

    pim_delta_t_vs_p_FTOF->Draw("COLZ");
    canvas->Update();
	canvas->Print(pdfFileName);

    pim_delta_t_vs_p_CTOF->Draw("COLZ");
    canvas->Update();
	canvas->Print(pdfFileName);

    canvas->Print(pdfFileName + "]");
	delete canvas;

    double mean, sigma, leftLimit, rightLimit;
    TF1 *fit;
    TLine *l1;
    TLine *l2;

    TCanvas *canvas2 = new TCanvas("canvas", "Canvas Title", 800, 600);
    TString pdfFileName2 = Form("FTOF_proton_delta_t.pdf");  
    canvas2->Print(pdfFileName2 + "[");

    for (int i = 0; i < 20; ++i)
    {
        p_delta_t_FTOF[i] -> Fit("gaus");
        fit = p_delta_t_FTOF[i]->GetFunction("gaus");

        p_delta_t_FTOF[i] -> Draw();

        if (fit) {
            mean  = fit->GetParameter(1);
            sigma = fit->GetParameter(2);
            leftLimit  = mean - 3.5 * sigma;
            rightLimit = mean + 3.5 * sigma;

            l1 = new TLine(mean - 3.5*sigma, 0, mean - 3.5*sigma, p_delta_t_FTOF[i]->GetMaximum());
            l1->SetLineColor(kRed);
            l1->Draw(); // ROOT запомнит этот конкретный объект на канвасе

            l2 = new TLine(mean + 3.5*sigma, 0, mean + 3.5*sigma, p_delta_t_FTOF[i]->GetMaximum());
            l2->SetLineColor(kRed);
            l2->Draw();
        }

        canvas2->Update();
	    canvas2->Print(pdfFileName2);
    }

    canvas2->Print(pdfFileName2 + "]");
	delete canvas2;

    TCanvas *canvas3 = new TCanvas("canvas", "Canvas Title", 800, 600);
    TString pdfFileName3 = Form("CTOF_proton_delta_t.pdf");  
    canvas3->Print(pdfFileName3 + "[");

    for (int i = 0; i < 20; ++i)
    {
        p_delta_t_CTOF[i] -> Fit("gaus");
        fit = p_delta_t_CTOF[i]->GetFunction("gaus");

        p_delta_t_CTOF[i] -> Draw();

        if (fit) {
            mean  = fit->GetParameter(1);
            sigma = fit->GetParameter(2);
            leftLimit  = mean - 3.5 * sigma;
            rightLimit = mean + 3.5 * sigma;

            l1 = new TLine(mean - 3.5*sigma, 0, mean - 3.5*sigma, p_delta_t_CTOF[i]->GetMaximum());
            l1->SetLineColor(kRed);
            l1->Draw(); // ROOT запомнит этот конкретный объект на канвасе

            l2 = new TLine(mean + 3.5*sigma, 0, mean + 3.5*sigma, p_delta_t_CTOF[i]->GetMaximum());
            l2->SetLineColor(kRed);
            l2->Draw();
        }

        canvas3->Update();
	    canvas3->Print(pdfFileName3);
    }

    canvas3->Print(pdfFileName3 + "]");
	delete canvas3;

    TCanvas *canvas4 = new TCanvas("canvas", "Canvas Title", 800, 600);
    TString pdfFileName4 = Form("FTOF_pim_delta_t.pdf");  
    canvas4->Print(pdfFileName4 + "[");

    for (int i = 0; i < 20; ++i)
    {
        pim_delta_t_FTOF[i] -> Fit("gaus");
        fit = pim_delta_t_FTOF[i]->GetFunction("gaus");

        pim_delta_t_FTOF[i] -> Draw();

        if (fit) {
            mean  = fit->GetParameter(1);
            sigma = fit->GetParameter(2);
            leftLimit  = mean - 3.5 * sigma;
            rightLimit = mean + 3.5 * sigma;

            l1 = new TLine(mean - 3.5*sigma, 0, mean - 3.5*sigma, pim_delta_t_FTOF[i]->GetMaximum());
            l1->SetLineColor(kRed);
            l1->Draw(); // ROOT запомнит этот конкретный объект на канвасе

            l2 = new TLine(mean + 3.5*sigma, 0, mean + 3.5*sigma, pim_delta_t_FTOF[i]->GetMaximum());
            l2->SetLineColor(kRed);
            l2->Draw();
        }

        canvas4->Update();
	    canvas4->Print(pdfFileName4);
    }

    canvas4->Print(pdfFileName4 + "]");
	delete canvas4;

    TCanvas *canvas5 = new TCanvas("canvas", "Canvas Title", 800, 600);
    TString pdfFileName5 = Form("CTOF_pim_delta_t.pdf");  
    canvas5->Print(pdfFileName5 + "[");

    for (int i = 0; i < 20; ++i)
    {
        pim_delta_t_CTOF[i] -> Fit("gaus");
        fit = pim_delta_t_CTOF[i]->GetFunction("gaus");

        pim_delta_t_CTOF[i] -> Draw();

        if (fit) {
            mean  = fit->GetParameter(1);
            sigma = fit->GetParameter(2);
            leftLimit  = mean - 3.5 * sigma;
            rightLimit = mean + 3.5 * sigma;

            l1 = new TLine(mean - 3.5*sigma, 0, mean - 3.5*sigma, pim_delta_t_CTOF[i]->GetMaximum());
            l1->SetLineColor(kRed);
            l1->Draw(); // ROOT запомнит этот конкретный объект на канвасе

            l2 = new TLine(mean + 3.5*sigma, 0, mean + 3.5*sigma, pim_delta_t_CTOF[i]->GetMaximum());
            l2->SetLineColor(kRed);
            l2->Draw();
        }

        canvas5->Update();
	    canvas5->Print(pdfFileName5);
    }

    canvas5->Print(pdfFileName5 + "]");
	delete canvas5;

    gSystem -> Exit(0);
}