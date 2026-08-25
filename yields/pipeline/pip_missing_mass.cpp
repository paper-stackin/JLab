#include "pip_missing_mass_config.h"

double EvaluatePolynomial(const double* p, double x)
{
    return p[0] + x * (p[1] + x * (p[2] + x * p[3]));
}

bool PCalFiducialCut(
    int cal_layer_e,
    int cal_sect_e,
    double cal_x_e,
    double cal_y_e
)
{
    const double p_lower[6][4] = {
                    {175.997,-2.19104,0.00698198,-9.2885e-06},
                    {219.797,-2.71269,0.00891301,-1.15516e-05},
                    {226.274,-2.84365,0.00957273,-1.23558e-05},
                    {181.685,-2.2416,0.00702503,-9.08476e-06},
                    {223.327,-2.79671,0.00935171,-1.20265e-05},
                    {152.767,-1.77109,0.004802,-5.80913e-06}
                };

    const double p_upper[6][4] = {
                    {-176.431,2.16702,-0.00665704,8.51867e-06},
                    {-214.182,2.60839,-0.00830054,1.04541e-05},
                    {-227.966,2.84059,-0.00948886,1.21912e-05},
                    {-178.192,2.23739,-0.00721147,9.56152e-06},
                    {-207.242,2.58685,-0.00845863,1.083e-05},
                    {-168.544,2.08019,-0.00634868,8.22812e-06}
                };

    if (cal_layer_e != 1) return true;

    double angle = TMath::Pi() / 180.0 * (-60.0 * (cal_sect_e - 1));
    double x = cal_x_e * TMath::Cos(angle) - cal_y_e * TMath::Sin(angle);
    double y = cal_x_e * TMath::Sin(angle) + cal_y_e * TMath::Cos(angle);

    const double y_lower = EvaluatePolynomial(p_lower[cal_sect_e - 1], x);
    const double y_upper = EvaluatePolynomial(p_upper[cal_sect_e - 1], x);

    const double r = sqrt(x * x + y * y);

    return y > y_upper || y < y_lower || r < 150;
}

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
    TTreeReaderValue<Int_t> cal_sect_e(reader, "cal_sect_e");
    TTreeReaderValue<Int_t> cal_layer_e(reader, "cal_layer_e");
    TTreeReaderValue<Float_t> cal_x_e(reader, "cal_x_e");
    TTreeReaderValue<Float_t> cal_y_e(reader, "cal_y_e");

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

        if (PCalFiducialCut(*cal_layer_e, *cal_sect_e, *cal_x_e, *cal_y_e)) continue;

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