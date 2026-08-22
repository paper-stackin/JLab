#include <yaml-cpp/yaml.h>
#include "pcal_config_yaml.h"

inline bool IsInExclusionStrip(int sector, double x, double y, 
                               const std::vector<ExclusionStripSettings::ExclusionStripConfig>& strips, 
                               double margin = 0.25)
{
    for (const auto& strip : strips) {
    if (strip.sector != sector) continue;
    const double y_lower = strip.slope * x + strip.c_min - margin;
    const double y_upper = strip.slope * x + strip.c_max + margin;
    if (y > y_lower && y < y_upper) return true;
    }
    return false;
}

void pcal_hist(const char* config_file)
{
    // Таймер
    auto start_time = std::chrono::steady_clock::now();
    
    // Конфигурация
    PCALConfig cfg = LoadConfig(config_file);

    // Гистограммы для каждого сектора и каждого бина по x'
    TH1F *x_bins[6][cfg.histogram.n_x_bins];
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < cfg.histogram.n_x_bins; ++j)
        {
            TString name = Form("h_pcal_sec%d_bin%d", i + 1, j + 1);
            TString title = Form(
                "%s x' = [%d, %d] cm, sec = %d%s", 
                cfg.histogram.title_prefix.c_str(), 
                cfg.histogram.x_start + cfg.histogram.bin_size * j, 
                cfg.histogram.x_start + cfg.histogram.bin_size * (j + 1), 
                i + 1, 
                cfg.histogram.title_suffix.c_str()
            );
            x_bins[i][j] = new TH1F(name, title, cfg.histogram.n_y_bins, cfg.histogram.y_min, cfg.histogram.y_max);
        }
    }

    // Гистограммы для 2D распределения x' и y' для каждого сектора
    TH2F* h2_xpy[6];
    for (int i = 0; i < 6; ++i)
    {
        h2_xpy[i] = new TH2F(
            Form("h2_sec%d", i+1),
            Form("%s %d%s", cfg.histogram.title_2d_prefix.c_str(), i+1, cfg.histogram.title_2d_suffix.c_str()),
            100, cfg.histogram.x_start, cfg.histogram.x_start + cfg.histogram.bin_size * cfg.histogram.n_x_bins,
            100, cfg.histogram.y_min, cfg.histogram.y_max
        );
    }

    // Рут файл с импульсами частиц
    TFile *infile_exp = TFile::Open(cfg.file.input_file.c_str()); 

	TTreeReader reader("MMpiptree", infile_exp);
    TTreeReaderValue<Int_t> cal_sect_e(reader, "cal_sect_e");
    TTreeReaderValue<Int_t> cal_layer_e(reader, "cal_layer_e");
    TTreeReaderValue<Float_t> cal_x_e(reader, "cal_x_e");
    TTreeReaderValue<Float_t> cal_y_e(reader, "cal_y_e");

    std::unique_ptr<TTreeReaderValue<Float_t>> weight;
    if (cfg.use_weights)  weight = std::make_unique<TTreeReaderValue<Float_t>>(reader, "weight");

    // Счётчик событий
    int totalEvents = reader.GetEntries(true) / cfg.counter.speed;
    int counter = 0;

    // Чтение событий из дерева
	while(reader.Next())
    {	      
        // Счётчик событий
        counter++;
        if (counter % cfg.counter.speed == 0) {
            cout << "\rProcessed: " << counter / cfg.counter.speed << "/" << totalEvents << cfg.counter.period << " events" << flush;
        }

        // Отбор событий по слою
        if (*cal_layer_e != 1) continue;

        // Inefficient PCAL Regions
        if (cfg.exclusion.enable &&
            IsInExclusionStrip(*cal_sect_e, *cal_x_e, *cal_y_e, cfg.exclusion.strips, cfg.exclusion.margin)) {
            continue;
        }

        // Распределение событий по секторам
        double angle = TMath::Pi() / 180.0 * (-60.0 * (*cal_sect_e - 1));
        double x = *cal_x_e * TMath::Cos(angle) - *cal_y_e * TMath::Sin(angle);
        double y = *cal_x_e * TMath::Sin(angle) + *cal_y_e * TMath::Cos(angle);

        int bin = int((x - cfg.histogram.x_start) / float(cfg.histogram.bin_size));
        if (bin < 0 || bin >= cfg.histogram.n_x_bins) continue;

        // Заполнение гистограмм
        if (cfg.use_weights)  x_bins[*cal_sect_e - 1][bin]->Fill(y, **weight);
        else                    x_bins[*cal_sect_e - 1][bin]->Fill(y);

        if (cfg.use_weights)  h2_xpy[*cal_sect_e - 1]->Fill(x, y, **weight);
        else                    h2_xpy[*cal_sect_e - 1]->Fill(x, y);
	}

    // Сохранение гистограмм в ROOT файл
    TFile *file = new TFile(cfg.file.intermediate_file.c_str(), "RECREATE");
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < cfg.histogram.n_x_bins; ++j)
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