void InitRadialHistograms(TH1F *hR[6], const char *prefix)
{
    for (int i = 0; i < 6; i++)
    {
        hR[i] = new TH1F(Form("%s_%d", prefix, i+1),
                         Form("Sector %d", i+1),
                         100, 100, 400);
    }
}

void FillRadialHistogram(TFile *infile, TH1F *hR[6], bool use_weight)
{
    TTreeReader reader("MMpiptree", infile);
    TTreeReaderValue<Int_t> cal_sect_e(reader, "cal_sect_e");
    TTreeReaderValue<Int_t> cal_layer_e(reader, "cal_layer_e");
    TTreeReaderValue<Float_t> x(reader, "cal_x_e");
    TTreeReaderValue<Float_t> y(reader, "cal_y_e");
    std::unique_ptr<TTreeReaderValue<Float_t>> weight;
    if (use_weight)  weight = std::make_unique<TTreeReaderValue<Float_t>>(reader, "weight");

    int speed = 1e6;
    int totalEvents = reader.GetEntries(true) / speed;
    int counter = 0;

    while (reader.Next())
    {
        // Счётчик событий
        counter++;
        if (counter % speed == 0) {
            cout << "\rProcessed: " << counter / speed << "/" << totalEvents << "M events" << flush;
        }

        if (*cal_layer_e != 1) continue;

        double rx = *x;
        double ry = *y;
        double r = sqrt(rx * rx + ry * ry);
        if (r < 150) continue;

        int sector = *cal_sect_e;
        if (use_weight) hR[sector - 1]->Fill(r, **weight);
        else            hR[sector - 1]->Fill(r);
    }
}

void StyleRadialHistogram(TH1F *hist, Color_t color, int sector)
{
    if (hist->Integral() > 0) hist->Scale(1.0 / hist->Integral());

    hist->SetStats(0);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);

    hist->SetTitle(Form("Sector %d", sector));
    hist->GetXaxis()->SetTitle("r, cm");
    hist->GetYaxis()->SetTitle("N");
}

void radial_cuts_det(void)
{
    // Таймер
    auto start_time = std::chrono::steady_clock::now();

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
	gROOT->SetBatch(kTRUE); 
    gErrorIgnoreLevel = kError;
    gPrintViaErrorHandler = kTRUE; 

    TH1F *hR_exp[6];
    TH1F *hR_mc[6];

    InitRadialHistograms(hR_exp, "hR_exp");
    InitRadialHistograms(hR_mc, "hR_mc");

    // Открываем рут файл с импульсами частиц
    std::string mc_file_name = "input_data/check_pcal_simulation.root";
    std::string exp_file_name = "input_data/check_pcal.root";

    TFile *infile_exp = TFile::Open(exp_file_name.c_str());
    TFile *infile_mc = TFile::Open(mc_file_name.c_str());

    FillRadialHistogram(infile_exp, hR_exp, false);
    FillRadialHistogram(infile_mc, hR_mc, true);

    TCanvas *c = new TCanvas("c", "", 800, 600);
    TString pdf_name = Form("output_graphs/radial_distributions.pdf");

    c->Print(pdf_name + "[");

    for (int i = 0; i < 6; i++)
    {
        StyleRadialHistogram(hR_exp[i], kBlack, i + 1);
        StyleRadialHistogram(hR_mc[i], kRed, i + 1);

        hR_mc[i]->Draw("HIST");
        hR_exp[i]->Draw("HIST SAME");

        TLegend *leg = new TLegend(0.65, 0.75, 0.88, 0.88);
        leg->AddEntry(hR_exp[i], "Experiment", "l");
        leg->AddEntry(hR_mc[i], "Simulation", "l");
        leg->Draw();

        c->Print(pdf_name);
        delete leg;
    }

    c->Print(pdf_name + "]");

    // Таймер
    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
    double elapsed_s = elapsed_ms / 1000.0;

    std::cout << "\nFinished in " << std::fixed << std::setprecision(3) << elapsed_s << " s" << std::endl;
    
    gSystem -> Exit(0);
}