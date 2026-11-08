void pi_plus_projected_vs_recreated_hists(bool use_cm_frame = false)
{
    auto start_time = std::chrono::steady_clock::now();

    gROOT->SetBatch(kTRUE);

    const char *input_file = "/home/stepan/downloads/missing_pip_clas12.root";
    const char *tree_name = "MMpiptree";
    const char *output_file = use_cm_frame
        ? "pi_plus_projected_vs_recreated_hists_cm.root"
        : "pi_plus_projected_vs_recreated_hists.root";
    const char *hist_name = "dE";

    std::cout << (use_cm_frame
        ? "Frame: CM (proton + virtual photon)"
        : "Frame: lab")
              << std::endl;

    const double E_beam = 6.535;
    const double M_proton = 0.93827;
    const double M_pion = 0.13957;
    const double M_electron = 0.0;

    const double Q2_edges[] = {0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.5};
    const int n_Q2_bins = 6;

    const int n_W_bins = 64;
    const double W_min = 1.4;
    const double W_max = 3.0;
    const double W_bin_size = 0.025;

    const int n_hist_bins = 100;
    const double hist_min = -0.3;
    const double hist_max = 0.5;

    const int speed = 1000000;
    const char *period = "M";

    TH1F *h_dE[n_Q2_bins][n_W_bins];
    for (int q = 0; q < n_Q2_bins; ++q)
    {
        for (int w = 0; w < n_W_bins; ++w)
        {
            char namehist[256];
            sprintf(namehist, "%s_Q2_bin=%d_W_bin=%d", hist_name, q + 1, w + 1);
            h_dE[q][w] = new TH1F(namehist, namehist, n_hist_bins, hist_min, hist_max);
            h_dE[q][w]->Sumw2();
        }
    }

    TFile *infile = TFile::Open(input_file);
    TTreeReader reader(tree_name, infile);

    TTreeReaderValue<Float_t> px_e(reader, "px_e");
    TTreeReaderValue<Float_t> py_e(reader, "py_e");
    TTreeReaderValue<Float_t> pz_e(reader, "pz_e");

    TTreeReaderValue<Float_t> px_p(reader, "px_p");
    TTreeReaderValue<Float_t> py_p(reader, "py_p");
    TTreeReaderValue<Float_t> pz_p(reader, "pz_p");

    TTreeReaderValue<Float_t> px_pim(reader, "px_pim");
    TTreeReaderValue<Float_t> py_pim(reader, "py_pim");
    TTreeReaderValue<Float_t> pz_pim(reader, "pz_pim");

    TLorentzVector beam(0, 0, E_beam, E_beam);
    TLorentzVector p_initial(0, 0, 0, M_proton);

    int totalEvents = reader.GetEntries(true) / speed;
    int counter = 0;

    while (reader.Next())
    {
        counter++;
        if (counter % speed == 0)
        {
            cout << "\rProcessed: " << counter / speed << "/" << totalEvents << period << " events" << flush;
        }

        TLorentzVector el, p_final, pim;
        el.SetXYZM(*px_e, *py_e, *pz_e, M_electron);
        p_final.SetXYZM(*px_p, *py_p, *pz_p, M_proton);
        pim.SetXYZM(*px_pim, *py_pim, *pz_pim, M_pion);

        double W_current = (beam - el + p_initial).M();
        double Q2_current = -(beam - el).M2();

        if (Q2_current < Q2_edges[0]) continue;

        int q;
        for (q = 0; q < n_Q2_bins; ++q)
        {
            if (Q2_current >= Q2_edges[q] && Q2_current <= Q2_edges[q + 1])
                break;
        }
        if (q == n_Q2_bins) continue;

        if (W_current < W_min || W_current > W_max) continue;
        int w = int((W_current - W_min) / W_bin_size);
        if (w < 0 || w >= n_W_bins) continue;

        TLorentzVector momentum_transfer = beam - el;
        TLorentzVector missing_pip = beam + p_initial - p_final - pim - el;

        if (use_cm_frame)
        {
            TVector3 beta = (momentum_transfer + p_initial).BoostVector();
            missing_pip.Boost(-beta);
        }

        double E_recreated = missing_pip.E();
        TVector3 p_missing_pip = missing_pip.Vect();
        double E_projected = sqrt(p_missing_pip.Mag2() + M_pion * M_pion);

        h_dE[q][w]->Fill(E_recreated - E_projected);
    }

    TFile *file = new TFile(output_file, "RECREATE");
    for (int q = 0; q < n_Q2_bins; ++q)
    {
        for (int w = 0; w < n_W_bins; ++w)
        {
            h_dE[q][w]->Write();
        }
    }
    file->Close();

    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    const auto hours = elapsed_seconds / 3600;
    const auto minutes = (elapsed_seconds % 3600) / 60;
    const auto seconds = elapsed_seconds % 60;

    std::cout << std::endl
              << "Finished in "
              << std::setw(2) << std::setfill('0') << hours << ":"
              << std::setw(2) << std::setfill('0') << minutes << ":"
              << std::setw(2) << std::setfill('0') << seconds << std::endl;

    gSystem->Exit(0);
}
