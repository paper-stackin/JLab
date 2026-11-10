#include "yields_config.h"
#include "yields_types.h"

#include "fit_model.cpp"
#include "load_histograms.cpp"
#include "background.cpp"
#include "calculate_yields.cpp"
#include "draw.cpp"
#include "process_topology.cpp"

void yields_pdf(void) 
{
    gErrorIgnoreLevel = kWarning;
	gROOT->SetBatch(kTRUE); 

    TFile *missing_pip_file = TFile::Open("interim/mm_pip_hists_m_pip.root");
    TFile *fully_exclusive_file = TFile::Open("interim/mm_pip_hists_m_0.root");

    int q_min = 0; int q_max = 5;
    int w_min = 0; int w_max = 99;

    auto missing_pip_hist = LoadHistograms(
        missing_pip_file, 
        "MM", 
        q_min, q_max, w_min, w_max
    );

    auto fully_exclusive_hist = LoadHistograms(
        fully_exclusive_file, 
        "MMpiplus_from_MM0_topology", 
        q_min, q_max, w_min, w_max
    );

    // Топология 1: Missing Pi+ (Сигнал и фон — одна гистограмма)
    YieldsConfig missing_pip_cfg;

    auto missing_pip = ProcessTopology(
        missing_pip_hist, {}, missing_pip_cfg, "missing_pip"
    );

    // Топология 2: Fully exclusive (Сигнал из первой, форма фона из второй)
    YieldsConfig fully_exclusive_cfg;
    fully_exclusive_cfg.signal_sigma = 0.05;
    fully_exclusive_cfg.xmin = -0.1;
    fully_exclusive_cfg.xmax_high = 0.3;
    fully_exclusive_cfg.bg_ampl_bin = 50;
    fully_exclusive_cfg.bg_ampl_factor = 1.0;

    auto fully_exclusive = ProcessTopology(
        missing_pip_hist, fully_exclusive_hist, fully_exclusive_cfg, "fully_exclusive"
    );
    
    // Графики выходов
    DrawYields(missing_pip, fully_exclusive);

    gSystem->Exit(0);
}