#include "yields_config.h"
#include "yields_types.h"

#include "fit_model.cpp"
#include "load_histograms.cpp"
#include "background.cpp"
#include "calculate_yields.cpp"
#include "draw.cpp"
#include "process_topology.cpp"

void yields_pdf(
    const char* missing_pip_config_file = "missing_pip_config.yaml",
    const char* fully_exclusive_config_file = "fully_esclusive_config.yaml")
{
    gErrorIgnoreLevel = kWarning;
	gROOT->SetBatch(kTRUE); 

    YieldsConfig missing_pip_cfg = LoadYieldsConfig(missing_pip_config_file);
    YieldsConfig fully_exclusive_cfg = LoadYieldsConfig(fully_exclusive_config_file);

    TFile *missing_pip_file = TFile::Open(missing_pip_cfg.output.c_str());
    TFile *fully_exclusive_file = TFile::Open(fully_exclusive_cfg.output.c_str());

    auto missing_pip_hist = LoadHistograms(
        missing_pip_file,
        missing_pip_cfg.hist_name.c_str(),
        missing_pip_cfg.q_min, missing_pip_cfg.q_max,
        0, missing_pip_cfg.n_W_bins - 1
    );

    auto fully_exclusive_hist = LoadHistograms(
        fully_exclusive_file,
        fully_exclusive_cfg.hist_name.c_str(),
        fully_exclusive_cfg.q_min, fully_exclusive_cfg.q_max,
        0, fully_exclusive_cfg.n_W_bins - 1
    );

    // Топология 1: Missing Pi+ (Сигнал и фон — одна гистограмма)
    auto missing_pip = ProcessTopology(
        missing_pip_hist, {}, missing_pip_cfg, "missing_pip"
    );

    // Топология 2: Fully exclusive (Сигнал из первой, форма фона из второй)
    auto fully_exclusive = ProcessTopology(
        missing_pip_hist, fully_exclusive_hist, fully_exclusive_cfg, "fully_exclusive"
    );
    
    // Графики выходов
    DrawYields(missing_pip, fully_exclusive, missing_pip_cfg);

    gSystem->Exit(0);
}
