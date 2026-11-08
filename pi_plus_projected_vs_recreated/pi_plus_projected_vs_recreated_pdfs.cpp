void pi_plus_projected_vs_recreated_pdfs(bool use_cm_frame = false)
{
    auto start_time = std::chrono::steady_clock::now();

    gROOT->SetBatch(kTRUE);
    gErrorIgnoreLevel = kWarning;
    gStyle->SetOptStat(0);

    const char *input_file = use_cm_frame
        ? "pi_plus_projected_vs_recreated_hists_cm.root"
        : "pi_plus_projected_vs_recreated_hists.root";
    const char *hist_name = "dE";
    const char *frame_label = use_cm_frame ? "CM" : "Lab";

    const double Q2_edges[] = {0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.5};
    const int n_Q2_bins = 6;

    const int n_W_bins = 64;
    const double W_min = 1.4;
    const double W_bin_size = 0.025;

    TFile *infile = TFile::Open(input_file);
    if (!infile || infile->IsZombie())
    {
        std::cerr << "Cannot open " << input_file << std::endl;
        gSystem->Exit(1);
    }

    for (int q = 0; q < n_Q2_bins; ++q)
    {
        TCanvas *canvas = new TCanvas(Form("canvas_q%d", q + 1), "Canvas", 800, 600);
        TString pdfFileName = use_cm_frame
            ? TString::Format("dE_cm_Q2_bin_%d.pdf", q + 1)
            : TString::Format("dE_Q2_bin_%d.pdf", q + 1);
        canvas->Print(pdfFileName + "[");

        for (int w = 0; w < n_W_bins; ++w)
        {
            TString full_name = Form("%s_Q2_bin=%d_W_bin=%d", hist_name, q + 1, w + 1);
            TH1F *hist = (TH1F *)infile->Get(full_name);
            if (!hist)
            {
                std::cerr << "Missing histogram " << full_name << std::endl;
                continue;
            }

            double W_lo = W_min + w * W_bin_size;
            double W_hi = W_lo + W_bin_size;

            TString title = Form(
                "E_{rec} - E_{proj} for Q^{2} #in [%g, %g] GeV^{2}, "
                "W #in [%g, %g] GeV in %s frame; "
                "E_{#pi^{+}}^{rec} - E_{#pi^{+}}^{proj}, GeV; N",
                Q2_edges[q], Q2_edges[q + 1], W_lo, W_hi, frame_label);

            hist->SetTitle(title);
            hist->SetLineColor(kBlue);
            hist->SetLineWidth(2);
            hist->Draw("HIST");

            TLegend *leg = new TLegend(0.58, 0.68, 0.88, 0.88);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextSize(0.03);
            leg->AddEntry((TObject *)0, Form("Entries = %.0f", hist->GetEntries()), "");
            // leg->AddEntry((TObject *)0, Form("Mean = %.4f", hist->GetMean()), "");
            // leg->AddEntry((TObject *)0, Form("RMS = %.4f", hist->GetRMS()), "");
            leg->Draw();

            canvas->Update();
            canvas->Print(pdfFileName);
        }

        canvas->Print(pdfFileName + "]");
        delete canvas;
    }

    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();
    const auto hours = elapsed_seconds / 3600;
    const auto minutes = (elapsed_seconds % 3600) / 60;
    const auto seconds = elapsed_seconds % 60;

    std::cout << std::endl
              << "Finished in "
              << std::setw(2) << std::setfill('0') << hours << ":"
              << std::setw(2) << std::setfill('0') << minutes << ":"
              << std::setw(2) << std::setfill('0') << seconds << std::endl;

    gSystem->Exit(0);
}
