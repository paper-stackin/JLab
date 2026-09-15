void pip_abs_momentum_dist_hists()
{
    gROOT -> SetBatch(kTRUE);

    double n_bins = 100, first_bin_val = 0, last_bin_val = 2;
    double step = (last_bin_val - first_bin_val) / n_bins;

    TH1F *momentum_dist[6][64];
	for (int q = 0; q < 6; ++q) 
	    for (int w = 0; w < 64; ++w)
		{
			char namehist[256];
			sprintf(namehist, "momentum_dist_Q2_bin=%d_W_bin=%d", q + 1, w + 1); 
			momentum_dist[q][w] = new TH1F (namehist, namehist, n_bins, first_bin_val, last_bin_val);
		}

    TFile *infile = new TFile("/home/stepan/downloads/missing_pip_simulation.root");
    TTreeReader reader("MMpiptree", infile);

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

	while (reader.Next())
    {	
        el_final.SetXYZM(*px_e, *py_e, *pz_e, 0.);
        p_final.SetXYZM(*px_p, *py_p, *pz_p, 0.938);
        pi_minus.SetXYZM(*px_pim, *py_pim, *pz_pim, 0.139);

        TLorentzVector momentum_transfer = el_initial - el_final;
        double Q2 = -momentum_transfer.M2();
        double W = (p_initial + momentum_transfer).M();

        if(W < 1.4 || W > 3 || Q2 < 0.5) continue;

        int q, w;

        if      (Q2 >= 0.5 && Q2 < 0.7)  q = 0;
        else if (Q2 >= 0.7 && Q2 < 1.0)  q = 1;
        else if (Q2 >= 1.0 && Q2 < 1.4)  q = 2;
        else if (Q2 >= 1.4 && Q2 < 2.0)  q = 3;
        else if (Q2 >= 2.0 && Q2 < 3.0)  q = 4;
        else if (Q2 >= 3.0)              q = 5;

        w = int((W - 1.4) / 0.025);

        TLorentzVector missing_pip = el_initial + p_initial - p_final - pi_minus - el_final;
        TVector3 beta = (momentum_transfer + p_initial).BoostVector();
        missing_pip.Boost(-beta);
        TVector3 p_missing_pip = missing_pip.Vect();

        momentum_dist[q][w] -> Fill(p_missing_pip.Mag(), *weight);
    }

    TFile *file = new TFile("pip_abs_momentum_dist_hists.root","RECREATE");

    for (int q = 0; q < 6; ++q) 
	    for (int w = 0; w < 64; ++w)
        {
            momentum_dist[q][w] -> Write();
            delete momentum_dist[q][w];
        }

    file -> Close();
    gSystem -> Exit(0);
}