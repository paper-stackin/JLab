#ifndef PCAL_CONFIG_H
#define PCAL_CONFIG_H

#include <vector>
#include <utility>

struct DipSectorConfig {
    int sector = 0; // Сектор, для которого применяются настройки интерполяции провалов
    int start_bin = 1; // Начальный бин (включительно) для интерполяции провалов 
    int end_bin = -1; // Конечный бин (включительно) для интерполяции провалов
    double relative_depth = -1.0; // Относительная глубина провала, при которой будет применяться интерполяция (0.0 - 1.0)
};

struct PCALConfig {
    // Симуляция или экспериментальные данные
    bool is_simulation = false;

    // Параметры для гистограмм
    int x_start = 110;
    int N_x_bins = 26;
    int bin_size = 10;
    int N_y_bins = 150;
    int y_min = -200;
    int y_max = 200;
    const char* histogram_title_prefix = "PCAL electron";
    const char* histogram_title_suffix = "; y', cm";
    const char* histogram_2d_title_prefix = "PCAL sector";
    const char* histogram_2d_title_suffix = "; x', cm; y', cm";

    // Обнаружение провалов в 1D гистограммах
    bool enable_dip_interpolation = true;
    std::vector<DipSectorConfig> dip_interpolation_sectors = {
        // {sector, start_bin, end_bin, relative_depth}
        {1, 2, 20, 0.25},
        {2, 5, 20, 0.3},
        {4, 14, 24, 0.3},
        {6, 8, 22, 0.3}
    }; // Для каждого сектора, где нужно интерполировать провалы
    double dip_absolute_drop = 5.0; // Минимальное абсолютное падение провала (в единицах гистограммы)
    int dip_max_width = 5; // Максимальная ширина провала (в бинах), который будет интерполирован

    // Параметры для построения графиков верхней и нижней границ
    double quantile_low = 0.05;
    double quantile_high = 0.95;
    const char* fit_formula = "pol2"; // "pol1", "pol2", "pol3", "expo", "gaus", etc.
    int x_fit_min = 140;
    int x_fit_max = 240;

    // Параметры Radial cut
    int r_min = 150;
    int r_max = 360;
    double radial_phi1_deg = -45.0;
    double radial_phi2_deg = 45.0;

    // Имена файлов и выходов
    const char* input_file_sim = "check_pcal_simulation.root";
    const char* input_file_exp = "check_pcal.root";
    const char* pdf_name_prefix = "PCAL_sector_";
    const char* pdf_extension = ".pdf";

    // Настройки стиля
    bool use_batch = true; // true - без графического интерфейса, false - с графическим интерфейсом
    int palette = kRainBow; // kRainBow, kBird, kDeepSea, kDarkBodyRadiator, kLightTemperature, kBlueGreenYellow, kBlueRedYellow
    int stat_option = 0; // 0 - без статистики, 111 - с полной статистикой
    int fit_option = 0; // 0 - без информации о фитах, 111 - с полной информацией о фитах
    int canvas_width = 1500;
    int canvas_height = 1000;   

    // Частота счётчика событий
    int cnt_speed = 1e6;
    string period = "M"; 
};

#endif
