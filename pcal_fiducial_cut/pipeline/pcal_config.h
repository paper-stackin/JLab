#ifndef PCAL_CONFIG_H
#define PCAL_CONFIG_H

#include <vector>
#include <utility>

// Имена файлов входов и выходов
struct FileConfig {
    std::string main_folder = "/home/stepan/root_progs/2pion_new/pcal_fiducial_cut/";

    std::string input_file = main_folder + "input_data/check_pcal.root";

    std::string intermediate_file = main_folder + "intermediate_data/PCAL_hist.root";

    std::string pdf_name_prefix = main_folder + "output_graphs/PCAL_sector_";
};

// Параметры для гистограмм
struct HistogramConfig {
    int x_start = 110;
    int n_x_bins = 26;
    int bin_size = 10;

    int n_y_bins = 50;
    int y_min = -200;
    int y_max = 200;

    const char* title_prefix = "PCAL electron";
    const char* title_suffix = "; y', cm";

    const char* title_2d_prefix = "PCAL sector";
    const char* title_2d_suffix = "; x', cm; y', cm";
};

// Исключение событий между линиями y = m*x + c_min - margin и y = m*x + c_max + margin
struct ExclusionStripSettings {
    bool enable = true;

    double margin = 0.25;

    struct ExclusionStripConfig {
        int sector = 0;
        double slope = 0.0;   // m
        double c_min = 0.0;
        double c_max = 0.0;
    };

    std::vector<ExclusionStripConfig> strips = {
        // {sector, slope, c_min, c_max}
        {1, 0.56575, -94.4,  -92.0},
        {1, 0.56575, -103.5, -101.1},
        {1, 0.56575, -221.4, -219.0},
        {1, 0.56575, -229.4, -227.0},
        {2, 0.5900,  114.4,  120.8},
        {4, -0.5680, -236.3, -232.8},
        {6, -0.59138, -187.0, -185.0},
        {6, -0.59138, -195.5, -193.3},
    };
};

// Интерполяция провалов в 1D гистограммах
struct DipInterpolationConfig {
    bool enable = true;

    struct DipSectorConfig {
        int sector = 0; // Сектор, для которого применяются настройки интерполяции провалов
        int start_bin = 1; // Начальный бин (включительно) для интерполяции провалов 
        int end_bin = -1; // Конечный бин (включительно) для интерполяции провалов
    };

    std::vector<DipSectorConfig> sectors = {
        // {sector, start_bin, end_bin}
        {1, 2, 20},
        {2, 5, 20},
        {4, 14, 24},
        {6, 8, 22}
    };
};

// Параметры для построения графиков верхней и нижней границ
struct BoundaryConfig {
    bool use_supergauss = true;

    int supergauss_n = 6;
    double supergauss_level = 0.2;

    double quantile_low = 0.05;
    double quantile_high = 0.95;

    const char* fit_formula = "pol2";
    int x_fit_min = 140;
    int x_fit_max = 240;
};

// Параметры Radial cut
struct RadialCutConfig {
    int r_min = 150;
    int r_max = 360;
    double phi_1 = -45.0;
    double phi_2 = 45.0;
};

// Настройки стиля
struct StyleConfig {
    bool use_batch = true;

    int palette = kRainBow;

    int stat_option = 0;
    int fit_option = 0;

    int canvas_width = 1500;
    int canvas_height = 1000;
};

// Счётчик событий
struct CounterConfig {
    int speed = 1e6;
    std::string period = "M";
};

struct PCALConfig {
    // Брать ли веса
    bool use_weights = true;

    FileConfig file;

    HistogramConfig histogram;

    ExclusionStripSettings exclusion;

    DipInterpolationConfig dip;

    RadialCutConfig radial;

    BoundaryConfig boundary;

    StyleConfig style;

    CounterConfig counter; 
};

#endif
