YieldsOutput ProcessTopology(
    const std::vector<std::vector<TH1F*>>& signal_hists,                     
    const std::vector<std::vector<TH1F*>>& bg_hists, 
    const YieldsConfig& cfg, 
    const char* topology_name) 
{
    const auto& actual_bg_hists = (!bg_hists.empty()) ? bg_hists : signal_hists;

    std::cout << "Processing topology: " << topology_name << std::endl;
    
    std::cout << "---> Fitting Missing Mass Background" << std::endl;
    std::vector<std::vector<FitResult>> fit_results = FitBackground(actual_bg_hists, cfg);

    std::cout << "---> Drawing Missing Mass Fit" << std::endl;
    DrawMissingMassFit(actual_bg_hists, fit_results, topology_name, cfg);

    std::cout << "---> Extrapolating Background" << std::endl;
    ExtrapolationOutput extrapolation_output = ExtrapolateBackground(fit_results, cfg);

    std::cout << "---> Drawing Background Extrapolation" << std::endl;
    DrawBackgroundExtrapolation(extrapolation_output, cfg, topology_name);

    std::cout << "---> Calculating Yields" << std::endl;
    return CalculateYields(signal_hists, actual_bg_hists, fit_results, extrapolation_output.integral_background, cfg);
}