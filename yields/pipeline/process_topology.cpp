YieldsOutput ProcessTopology(
    const std::vector<std::vector<TH1F*>>& signal_hists,                     
    const std::vector<std::vector<TH1F*>>& bg_hists, 
    const YieldsConfig& cfg)
{
    std::cout << "Processing topology: " << cfg.topology_name.c_str() << std::endl;
    
    std::cout << "---> Fitting Missing Mass Background" << std::endl;
    std::vector<std::vector<FitResult>> fit_results = FitBackground(bg_hists, cfg);

    std::cout << "---> Drawing Missing Mass Fit" << std::endl;
    DrawMissingMassFit(bg_hists, fit_results, cfg);

    std::cout << "---> Extrapolating Background" << std::endl;
    ExtrapolationOutput extrapolation_output = ExtrapolateBackground(fit_results, cfg);

    std::cout << "---> Drawing Background Extrapolation" << std::endl;
    DrawBackgroundExtrapolation(extrapolation_output, cfg);

    std::cout << "---> Calculating Yields" << std::endl;
    return CalculateYields(signal_hists, bg_hists, fit_results, extrapolation_output, cfg);
}