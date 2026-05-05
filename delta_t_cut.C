template <typename T, size_t N>
void MomentumBinsGraphs(T* (&hists)[N], const TString& pdfName, double fit_x_min, double fit_x_max)
{
    TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 800, 600);
    canvas->Print(pdfName + "[");

    for (size_t i = 0; i < N; ++i)
    {
        double xMin = hists[i]->GetXaxis()->GetXmin();
        double xMax = hists[i]->GetXaxis()->GetXmax();

        TF1* fit = new TF1("fit", "gaus", fit_x_min, fit_x_max);
        hists[i]->Fit(fit, "RQ");

        fit->SetRange(xMin, xMax);

        hists[i]->Draw();
        fit->Draw("same");

        gPad->Update();
        double yTop = gPad->GetUymax();

        if (fit) {
            double mean  = fit->GetParameter(1);
            double sigma = fit->GetParameter(2);
            double left  = mean - 3.5 * sigma;
            double right = mean + 3.5 * sigma;

            TLine *l1 = new TLine(left, 0, left, yTop);
            TLine *l2 = new TLine(right, 0, right, yTop);

            l1->SetLineColor(kRed);
            l2->SetLineColor(kRed);
            l1->Draw();
            l2->Draw();

            TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);

            leg->AddEntry((TObject*)0, Form("Mean = %.4f", mean), "");
            leg->AddEntry((TObject*)0, Form("#sigma = %.4f", sigma), "");
            leg->AddEntry((TObject*)0, Form("Mean - 3.5#sigma = %.4f", left), "");
            leg->AddEntry((TObject*)0, Form("Mean + 3.5#sigma = %.4f", right), "");
            leg->Draw();
        }

        canvas->Update();
        canvas->Print(pdfName);
    }

    canvas->Print(pdfName + "]");
    delete canvas;
}

void CreateDeltaTHists(TH1F* hists[], const TString& prefix)
{
    for (int i = 0; i < 20; ++i)
    {
        TString title = Form("%s #Deltat, p [%.2f, %.2f] GeV; #Deltat, ns",
                             prefix.Data(), 0.25*i, 0.25*(i+1));

        hists[i] = new TH1F(title, title, 100, -1, 1);
    }
}

void delta_t_cut(void)
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
	gROOT->SetBatch(kTRUE); 

	TH2F *p_delta_t_vs_p_FD = new TH2F("p_delta_t_FD", "FD proton p VS #Deltat; p, GeV; #Deltat, ns", 100, 0, 5, 100, -10, 10);
    TH2F *p_delta_t_vs_p_CD = new TH2F("p_delta_t_CD", "CD proton p VS #Deltat; p, GeV; #Deltat, ns", 100, 0, 5, 100, -10, 10);
    TH2F *pim_delta_t_vs_p_FD = new TH2F("pim_delta_t_FD", "FD #pi^{-} p VS #Deltat; p, GeV; #Deltat, ns", 100, 0, 5, 100, -8, 8);
    TH2F *pim_delta_t_vs_p_CD = new TH2F("pim_delta_t_CD", "CD #pi^{-} p VS #Deltat; p, GeV; #Deltat, ns", 100, 0, 5, 100, -5, 5);

    TH1F *p_delta_t_FD[20], *p_delta_t_CD[20], *pim_delta_t_FD[20], *pim_delta_t_CD[20];
    CreateDeltaTHists(p_delta_t_FD,  "FD proton");
    CreateDeltaTHists(p_delta_t_CD,  "CD proton");
    CreateDeltaTHists(pim_delta_t_FD,"FD #pi^{-}");
    CreateDeltaTHists(pim_delta_t_CD,"CD #pi^{-}");

	TFile *infile_exp = TFile::Open("/home/stepan/root_progs/2pion_new/clas12/data/MPPT_events_clas12_selection_with_additional_values_saved.root"); // открываем рут файл с импульсами частиц

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

    //TTreeReaderValue<Float_t> weight(reader, "weight");

    TLorentzVector el_final, pi_minus, p_final;
    TLorentzVector el_initial(0, 0, 6.535, 6.535);
	TLorentzVector p_initial(0, 0, 0, 0.938);

    int counter = 0;
	/////////////////////////////////////////////////////////////////
	while(reader.Next())
    {	      
        counter++;
        if (counter % 100000 == 0)  cout << counter / 100000 << "\n";

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
    	
        if (2000 <= *status_p && *status_p < 4000) p_delta_t_vs_p_FD->Fill(p_mom, p_delta_t);
        if (4000 <= *status_p && *status_p < 8000) p_delta_t_vs_p_CD->Fill(p_mom, p_delta_t);

        if (2000 <= *status_pim && *status_pim < 4000) pim_delta_t_vs_p_FD->Fill(pim_mom, pim_delta_t);
        if (4000 <= *status_pim && *status_pim < 8000) pim_delta_t_vs_p_CD->Fill(pim_mom, pim_delta_t);

        int p_bin = p_mom / 0.25;
        int pim_bin = pim_mom / 0.25;
        if (p_bin > 19) p_bin = 19;
        if (pim_bin > 19) pim_bin = 19;

        if (2000 <= *status_p && *status_p < 4000) p_delta_t_FD[p_bin]->Fill(p_delta_t);
        if (4000 <= *status_p && *status_p < 8000) p_delta_t_CD[p_bin]->Fill(p_delta_t);
        if (2000 <= *status_p && *status_p < 4000) pim_delta_t_FD[p_bin]->Fill(pim_delta_t);
        if (4000 <= *status_p && *status_p < 8000) pim_delta_t_CD[p_bin]->Fill(pim_delta_t);
	}

	/////////////////////////////////////////////////////////////////////

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptFit(111);

    TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 800, 600);
    TString pdfFileName = Form("delta_t_cuts.pdf");  
    canvas->Print(pdfFileName + "[");

    p_delta_t_vs_p_FD->Draw("COLZ");
    canvas->Update();
	canvas->Print(pdfFileName);

    p_delta_t_vs_p_CD->Draw("COLZ");
    canvas->Update();
	canvas->Print(pdfFileName);

    pim_delta_t_vs_p_FD->Draw("COLZ");
    canvas->Update();
	canvas->Print(pdfFileName);

    pim_delta_t_vs_p_CD->Draw("COLZ");
    canvas->Update();
	canvas->Print(pdfFileName);

    canvas->Print(pdfFileName + "]");
	delete canvas;

    MomentumBinsGraphs(p_delta_t_FD, "FD_proton_delta_t.pdf", -0.3, 0.2);
    MomentumBinsGraphs(p_delta_t_CD, "CD_proton_delta_t.pdf", -0.3, 0.2);
    MomentumBinsGraphs(pim_delta_t_FD, "FD_pim_delta_t.pdf", -0.3, 0.3);
    MomentumBinsGraphs(pim_delta_t_CD, "CD_pim_delta_t.pdf", -0.25, 0.25);

    gSystem -> Exit(0);
}