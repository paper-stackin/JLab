void PCAL_cut(void)
{
    // gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
	gROOT->SetBatch(kTRUE); 

    TH1F *x_bins[6][30];
    for (int i = 0; i < 6; ++i)
    {
        for (int j = 0; j < 30; ++j)
        {
            TString name = Form("h_pcal_sec%d_bin%d", i + 1, j + 1);
            TString title = Form("PCAL electron x' = [%d, %d] cm, sec = %d; y', cm", 20 * j, 20 * (j + 1), i + 1);
            x_bins[i][j] = new TH1F(name, title, 100, -75, 75);
        }
    }

    TH2F* h2_xpy[6];
    for (int i = 0; i < 6; ++i)
    {
        h2_xpy[i] = new TH2F(
            Form("h2_sec%d", i+1),
            Form("PCAL sector %d; x', cm; y', cm", i+1),
            100, 0, 600,
            100, -200, 200
        );
    }

	TFile *infile_exp = TFile::Open("/home/stepan/root_progs/2pion_new/clas12/data/MPPT_events_clas12_selection_with_additional_values_saved_new.root"); // открываем рут файл с импульсами частиц

	TTreeReader reader("MMpiptree", infile_exp);

    TTreeReaderValue<Int_t> cal_sect_e(reader, "cal_sect_e");
    TTreeReaderValue<Int_t> cal_layer_e(reader, "cal_layer_e");
    TTreeReaderValue<Float_t> cal_x_e(reader, "cal_x_e");
    TTreeReaderValue<Float_t> cal_y_e(reader, "cal_y_e");

    int counter = 0;
    int layer_1_counter = 0;

	while(reader.Next())
    {	      
        counter++;
        if (counter % 100000 == 0)  cout << counter / 100000 << "\n";

        if (*cal_layer_e != 1) continue;
        layer_1_counter++;

        // cout << *cal_layer_e << "\n";

        double angle = TMath::Pi() / 180.0 * (-60.0 * (*cal_sect_e - 1));
        
        double x = *cal_x_e * TMath::Cos(angle) - *cal_y_e * TMath::Sin(angle);
        double y = *cal_x_e * TMath::Sin(angle) + *cal_y_e * TMath::Cos(angle);
        
        int bin = int(x / 20.0);
        if (bin < 0 || bin >= 30) continue;
        x_bins[*cal_sect_e - 1][bin]->Fill(y);

        h2_xpy[*cal_sect_e - 1]->Fill(x, y);
	}

    gStyle->SetPalette(kRainBow);
    gStyle->SetOptFit(111);

    for (int sec = 0; sec < 6; ++sec)
    {
        TCanvas *c = new TCanvas(
            Form("c_sec%d", sec + 1),
            Form("Sector %d", sec + 1),
            800,
            600
        );

        TString pdf_name = Form("PCAL_sector_%d.pdf", sec + 1);

        c->Print(pdf_name + "[");

        for (int bin = 0; bin < 30; ++bin)
        {
            x_bins[sec][bin]->Draw();
            c->Print(pdf_name);
        }

        h2_xpy[sec]->Draw("COLZ");
        c->Print(pdf_name);

        c->Print(pdf_name + "]");

        delete c;
    }

    cout << counter << " " << layer_1_counter << " " << layer_1_counter / counter;

    gSystem -> Exit(0);
}