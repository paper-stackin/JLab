std::vector<std::vector<TH1F*>> LoadHistograms(
    TFile* infile, 
    const char* hist_name_prefix, 
    int q_min, int q_max, 
    int w_min, int w_max) 
{
    std::vector<std::vector<TH1F*>> histograms(q_max + 1);

    for (int q = q_min; q <= q_max; ++q) 
    {
        histograms[q].resize(w_max + 1, nullptr);

        for (int w = w_min; w <= w_max; ++w) 
        {
            TString full_name = Form("%s_Q2_bin=%d_W_bin=%d", hist_name_prefix, q + 1, w + 1);
            histograms[q][w] = (TH1F*)infile->Get(full_name);
        }
    }

    return histograms;
}