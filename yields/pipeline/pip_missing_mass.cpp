#include "pip_missing_mass_config.h"

void pip_missing_mass(const char* config_file)
{
    // Timer
    auto start_time = std::chrono::steady_clock::now();

    // Configuration
    PipMissingMassConfig cfg = LoadPipMissingMassConfig(config_file);

    gROOT->SetBatch(kTRUE);

    TH1F *MM_raw[cfg.n_Q2_bins][cfg.n_W_bins];

    for (int q = 0; q < cfg.n_Q2_bins; ++q)
    {
        for (int w = 0; w < cfg.n_W_bins; ++w)
        {
            char namehist[256];
            sprintf(namehist, "%s_Q2_bin=%d_W_bin=%d", cfg.hist_name.c_str(), q+1, w+1);
            MM_raw[q][w] = new TH1F(namehist, namehist, cfg.n_hist_bins, cfg.hist_min, cfg.hist_max);
        }
    }

    TFile *infile_exp = TFile::Open(cfg.input.c_str()); // открываем рут файл с импульсами частиц

	TTreeReader reader(cfg.tree_name.c_str(), infile_exp);
    TTreeReaderValue<Float_t> px_e(reader, "px_e");
    TTreeReaderValue<Float_t> py_e(reader, "py_e");
    TTreeReaderValue<Float_t> pz_e(reader, "pz_e");

    TTreeReaderValue<Float_t> px_p(reader, "px_p");
    TTreeReaderValue<Float_t> py_p(reader, "py_p");
    TTreeReaderValue<Float_t> pz_p(reader, "pz_p");

    TTreeReaderValue<Float_t> px_pim(reader, "px_pim");
    TTreeReaderValue<Float_t> py_pim(reader, "py_pim");
    TTreeReaderValue<Float_t> pz_pim(reader, "pz_pim");

    TLorentzVector beam(0, 0, cfg.E_beam, cfg.E_beam);
	TLorentzVector p_initial(0, 0, 0, cfg.M_proton);

    // Счётчик событий
    int totalEvents = reader.GetEntries(true) / cfg.speed;
    int counter = 0;

	while(reader.Next())
    {	
        counter++;
        if (counter % cfg.speed == 0) 
        {
            cout << "\rProcessed: " << counter / cfg.speed << "/" << totalEvents << cfg.period << " events" << flush;
        }

        TLorentzVector el, pim, p_final;
        el.SetXYZM(*px_e, *py_e, *pz_e, cfg.M_electron);
        p_final.SetXYZM(*px_p, *py_p, *pz_p, cfg.M_proton);
        pim.SetXYZM(*px_pim, *py_pim, *pz_pim, cfg.M_pion);
        
        double  W_current = (beam - el + p_initial).M();
        double  Q2_current = -(beam - el).M2();

        if (Q2_current < cfg.Q2_edges[0]) continue;

        int q;

        for (q = 0; q < cfg.n_Q2_bins; ++q)
        {
            if (Q2_current >= cfg.Q2_edges[q] && Q2_current <= cfg.Q2_edges[q + 1])
                break;
        }

        if (q == cfg.n_Q2_bins) continue; // не попадает ни в один Q2 бин

        if (W_current > cfg.W_max || W_current < cfg.W_min) continue;
        int w = int((W_current - cfg.W_min) / cfg.W_bin_size); // бинирование по W начинаем с 1.4 GeV, ширина бина - 25 MeV

        double  pip_rec_MM2 = (beam + p_initial - p_final - pim - el).M2();

        if (pip_rec_MM2 < -0.2 || pip_rec_MM2 > 0.5) continue;
        MM_raw[q][w] -> Fill(pip_rec_MM2); // заполняем гистограммы
	}

    TFile *file = new TFile(cfg.output.c_str(), "RECREATE");

    for (int q = 0; q < cfg.n_Q2_bins; ++q)
    {
        for (int w = 0; w < cfg.n_W_bins; ++w)
        {
            MM_raw[q][w] -> Write();
        }
    }

    file -> Close();

    // Timer
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

    gSystem -> Exit(0);
}