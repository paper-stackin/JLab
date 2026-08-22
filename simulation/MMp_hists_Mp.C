void MMp_hists_Mp(void)
{
	gROOT->SetBatch(kTRUE); 

	TH1F *MM_raw[12][100];
	for(int q = 0; q < 6; ++q) 
    {
        for (int w = 0; w < 100; ++w)
        {
		    char namehist[256];
		    sprintf(namehist, "MM_Q2_bin=%d_W_bin=%d", q+1, w+1);
		    MM_raw[q][w] = new TH1F(namehist, namehist, 100, -0.3, 0.7);
		}
	}

	TFile *infile_exp = TFile::Open("MPPT_events_simulation_selection_with_additional_values_saved.root"); // открываем рут файл с импульсами частиц

	TTreeReader reader("MMpiptree",infile_exp);
    TTreeReaderValue<Float_t> px_e(reader, "px_e");
    TTreeReaderValue<Float_t> py_e(reader, "py_e");
    TTreeReaderValue<Float_t> pz_e(reader, "pz_e");

    TTreeReaderValue<Float_t> px_p(reader, "px_p");
    TTreeReaderValue<Float_t> py_p(reader, "py_p");
    TTreeReaderValue<Float_t> pz_p(reader, "pz_p");

    TTreeReaderValue<Float_t> px_pim(reader, "px_pim");
    TTreeReaderValue<Float_t> py_pim(reader, "py_pim");
    TTreeReaderValue<Float_t> pz_pim(reader, "pz_pim");

    TTreeReaderValue<Float_t> weight(reader, "weight");

    TLorentzVector el_final, pi_minus, p_final;
    TLorentzVector el_initial(0, 0, 6.535, 6.535);
	TLorentzVector p_initial(0, 0, 0, 0.938);

	/////////////////////////////////////////////////////////////////
	while(reader.Next())
    {	
	    double W_current, Q2_current;
        int w,q;
       
        el_final.SetXYZM(*px_e, *py_e, *pz_e, 0.);
        p_final.SetXYZM(*px_p, *py_p, *pz_p, 0.938);
        pi_minus.SetXYZM(*px_pim, *py_pim, *pz_pim, 0.139);
    
        
        W_current = (el_initial - el_final + p_initial).M();
        Q2_current = -(el_initial - el_final).M2();

        if(W_current > 3) continue;
        if(W_current < 1.4) continue;
        if(Q2_current < 0.5) continue;

        if(Q2_current >= 0.5 && Q2_current <= 0.7) q = 0;  // бинирование по Q2, ширина бинов различная
        if(Q2_current > 0.7 && Q2_current <= 1) q = 1;
        if(Q2_current > 1 && Q2_current <= 1.4) q = 2;
        if(Q2_current > 1.4 && Q2_current <= 2) q = 3;
        if(Q2_current > 2 && Q2_current <= 3) q = 4;
        if(Q2_current > 3 ) q = 5;

        w = int((W_current - 1.4) / 0.025); // бинирование по W начинаем с 1.4 GeV, ширина бина - 25 MeV
   
        MM_raw[q][w] -> Fill((el_initial + p_initial - p_final - pi_minus - el_final).M2(), *weight); // заполняем гистограммы	
	}

	/////////////////////////////////////////////////////////////////////
    TCanvas *canvas = new TCanvas("canvas", "Canvas Title", 800, 600);
    TString pdfFileName = Form("MM_p_hists_check.pdf");  
    canvas->Print(pdfFileName + "[");

    for (int q = 0; q < 6; ++q)
    {
        for (int w = 0; w < 64; ++w)
        {
            MM_raw[q][w] -> Draw();
            canvas->Update();
	        canvas->Print(pdfFileName);
        }
    }

    canvas->Print(pdfFileName + "]");
	delete canvas;

    gSystem -> Exit(0);
}