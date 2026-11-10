YieldsOutput ProcessTopology(
    const std::vector<std::vector<TH1F*>>& signal_hists,                     
    const std::vector<std::vector<TH1F*>>& bg_hists, 
    const YieldsConfig& cfg, 
    const char* diag_prefix) 
{
    const auto& actual_bg_hists = (!bg_hists.empty()) ? bg_hists : signal_hists;

    std::cout << "Processing topology: " << diag_prefix << std::endl;
    
    std::cout << "---> Fitting Background" << std::endl;
    auto fit_results = FitBackground(actual_bg_hists, cfg);

    std::cout << "---> Extrapolating Background" << std::endl;
    auto [integral_background, integral_background_error, w_bg_fit, extrapolation_results] = 
        ExtrapolateBackground(fit_results, cfg);

    TString bg_pdf_path = TString::Format("results/extrapolating_bg_%s.pdf", diag_prefix);
    TString fit_pdf_prefix = TString::Format("results/%s_fit", diag_prefix);

    std::cout << "---> Drawing Background Extrapolation" << std::endl;
    DrawBackgroundExtrapolation(
        integral_background, integral_background_error, w_bg_fit, 
        extrapolation_results, cfg, bg_pdf_path
    );

    std::cout << "---> Drawing Missing Mass Fit" << std::endl;
    DrawMissingMassFit(
        actual_bg_hists, fit_results, fit_pdf_prefix.Data(), cfg // Передаем .Data() безопасной строки
    );

    std::cout << "---> Calculating Yields" << std::endl;
    return CalculateYields(signal_hists, actual_bg_hists, fit_results, integral_background, cfg);
}