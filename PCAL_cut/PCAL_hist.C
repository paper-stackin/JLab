#include "PCAL_config.h"

void PCAL_hist(void)
{
    // Таймер
    auto start_time = std::chrono::steady_clock::now();
    
    // Конфигурация
    PCALConfig cfg;

    // Гистограммы для каждого сектора и каждого бина по x'
    TH1F *x_bins[6][cfg.N_x_bins];
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < cfg.N_x_bins; ++j)
        {
            TString name = Form("h_pcal_sec%d_bin%d", i + 1, j + 1);
            TString title = Form("%s x' = [%d, %d] cm, sec = %d%s", cfg.histogram_title_prefix, cfg.x_start + cfg.bin_size * j, cfg.x_start + cfg.bin_size * (j + 1), i + 1, cfg.histogram_title_suffix);
            x_bins[i][j] = new TH1F(name, title, cfg.N_y_bins, cfg.y_min, cfg.y_max);
        }
    }

    // Гистограммы для 2D распределения x' и y' для каждого сектора
    TH2F* h2_xpy[6];
    for (int i = 0; i < 6; ++i)
    {
        h2_xpy[i] = new TH2F(
            Form("h2_sec%d", i+1),
            Form("%s %d%s", cfg.histogram_2d_title_prefix, i+1, cfg.histogram_2d_title_suffix),
            100, cfg.x_start, cfg.x_start + cfg.bin_size * cfg.N_x_bins,
            100, cfg.y_min, cfg.y_max
        );
    }

    // Рут файл с импульсами частиц
    const char* file_name;
    if (cfg.is_simulation)  file_name = cfg.input_file_sim;
    else                   file_name = cfg.input_file_exp;
    TFile *infile_exp = TFile::Open(file_name);  

	TTreeReader reader("MMpiptree", infile_exp);
    TTreeReaderValue<Int_t> cal_sect_e(reader, "cal_sect_e");
    TTreeReaderValue<Int_t> cal_layer_e(reader, "cal_layer_e");
    TTreeReaderValue<Float_t> cal_x_e(reader, "cal_x_e");
    TTreeReaderValue<Float_t> cal_y_e(reader, "cal_y_e");

    std::unique_ptr<TTreeReaderValue<Float_t>> weight;
    if (cfg.is_simulation)  weight = std::make_unique<TTreeReaderValue<Float_t>>(reader, "weight");

    // Счётчик событий
    int totalEvents = reader.GetEntries(true) / cfg.cnt_speed;
    int counter = 0;

    // Чтение событий из дерева
	while(reader.Next())
    {	      
        // Счётчик событий
        counter++;
        if (counter % cfg.cnt_speed == 0) {
            cout << "\rProcessed: " << counter / cfg.cnt_speed << "/" << totalEvents << cfg.period << " events" << flush;
        }

        // Отбор событий по слою
        if (*cal_layer_e != 1) continue;

        // Распределение событий по секторам
        double angle = TMath::Pi() / 180.0 * (-60.0 * (*cal_sect_e - 1));
        double x = *cal_x_e * TMath::Cos(angle) - *cal_y_e * TMath::Sin(angle);
        double y = *cal_x_e * TMath::Sin(angle) + *cal_y_e * TMath::Cos(angle);

        int bin = int((x - cfg.x_start) / float(cfg.bin_size));
        if (bin < 0 || bin >= cfg.N_x_bins) continue;

        // Заполнение гистограмм
        if (cfg.is_simulation)  x_bins[*cal_sect_e - 1][bin]->Fill(y, **weight);
        else                    x_bins[*cal_sect_e - 1][bin]->Fill(y);

        if (cfg.is_simulation)  h2_xpy[*cal_sect_e - 1]->Fill(x, y, **weight);
        else                    h2_xpy[*cal_sect_e - 1]->Fill(x, y);
	}

    // Сохранение гистограмм в ROOT файл
    TFile *file = new TFile("PCAL_hist.root", "RECREATE");
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < cfg.N_x_bins; ++j)
        {   
            if (x_bins[i][j]) x_bins[i][j]->Write();
        }
        if (h2_xpy[i]) h2_xpy[i]->Write();
    }
    file->Close();

    // Таймер
    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    double elapsed_s = elapsed_ms / 1000.0;

    std::cout << "\nFinished in " << std::fixed << std::setprecision(3) << elapsed_s << " s" << std::endl;

    gSystem->Exit(0);
}