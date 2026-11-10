std::vector<std::vector<TH1F*>> LoadHistograms(const YieldsConfig& cfg) 
{
    TFile *infile = TFile::Open(cfg.output.c_str());
    std::vector<std::vector<TH1F*>> histograms(cfg.n_Q2_bins);

    for (int q = 0; q < cfg.n_Q2_bins; ++q) 
    {
        histograms[q].resize(cfg.n_W_bins, nullptr);

        for (int w = 0; w < cfg.n_W_bins; ++w) 
        {
            TString full_name = Form("%s_Q2_bin=%d_W_bin=%d", cfg.hist_name.c_str(), q + 1, w + 1);
            
            TH1F* hist = (TH1F*)infile->Get(full_name);
            if (hist) {
                hist->SetDirectory(nullptr); 
            }
            histograms[q][w] = hist;
        }
    }

    infile->Close();
    delete infile; 

    return histograms;
}