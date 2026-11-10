#include "yields_config.h"
#include "yields_types.h"

#include "fit_model.cpp"
#include "load_histograms.cpp"
#include "background.cpp"
#include "calculate_yields.cpp"
#include "draw.cpp"
#include "process_topology.cpp"

void yields_pdf()
{
    gErrorIgnoreLevel = kWarning;
	gROOT->SetBatch(kTRUE); 

    YieldsConfig missing_pip_cfg = LoadYieldsConfig("missing_pip_config.yaml");
    YieldsConfig fully_exclusive_cfg = LoadYieldsConfig("fully_esclusive_config.yaml");

    auto missing_pip_hist = LoadHistograms(missing_pip_cfg);
    auto fully_exclusive_hist = LoadHistograms(fully_exclusive_cfg);

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
