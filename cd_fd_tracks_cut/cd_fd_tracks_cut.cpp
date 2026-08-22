#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <cmath>
#include <iostream>

bool IsFD(int status)
{
    status = std::abs(status);
    return ((status / 2000) % 2) == 1;
}

bool IsCD(int status)
{
    status = std::abs(status);
    return ((status / 4000) % 2) == 1;
}

// Приводим phi к диапазону [-180, 180]
double DeltaPhi(double phi1, double phi2)
{
    double dphi = phi1 - phi2;

    while (dphi > 180.0)
        dphi -= 360.0;

    while (dphi < -180.0)
        dphi += 360.0;

    return dphi;
}


void cd_fd_tracks_cut()
{
    // ============================================================
    // Input
    // ============================================================

    const char* inputFile = "/home/stepan/downloads/missing_pip.root";
    const char* treeName  = "MMpiptree";
    const char* outputPDF = "removal_cd_fd.pdf";

    TFile* infile = TFile::Open(inputFile);

    if (!infile || infile->IsZombie())
    {
        std::cerr << "Error: cannot open file "
                  << inputFile << std::endl;
        return;
    }

    // ============================================================
    // TTreeReader
    // ============================================================

    TTreeReader reader(treeName, infile);

    // ------------------------------------------------------------
    // Proton
    // ------------------------------------------------------------

    TTreeReaderValue<Float_t> px_p(reader, "px_p");
    TTreeReaderValue<Float_t> py_p(reader, "py_p");
    TTreeReaderValue<Float_t> pz_p(reader, "pz_p");
    TTreeReaderValue<Int_t>   status_p(reader, "status_p");

    // ------------------------------------------------------------
    // Pi-
    // ------------------------------------------------------------

    TTreeReaderValue<Float_t> px_pim(reader, "px_pim");
    TTreeReaderValue<Float_t> py_pim(reader, "py_pim");
    TTreeReaderValue<Float_t> pz_pim(reader, "pz_pim");
    TTreeReaderValue<Int_t>   status_pim(reader, "status_pim");

    // ============================================================
    // Histograms
    // ============================================================

    // ------------------------------------------------------------
    // Proton
    // ------------------------------------------------------------

    TH1F* h_dp_p = new TH1F(
        "h_dp_p",
        "Proton;#Delta p = p_{FD} - p_{CD} [GeV];Events",
        200, -1.0, 1.0
    );

    TH1F* h_dtheta_p = new TH1F(
        "h_dtheta_p",
        "Proton;#Delta#theta = #theta_{FD} - #theta_{CD} [deg];Events",
        200, -30.0, 30.0
    );

    TH1F* h_dphi_p = new TH1F(
        "h_dphi_p",
        "Proton;#Delta#phi = #phi_{FD} - #phi_{CD} [deg];Events",
        200, -60.0, 60.0
    );

    // ------------------------------------------------------------
    // Pi-
    // ------------------------------------------------------------

    TH1F* h_dp_pim = new TH1F(
        "h_dp_pim",
        "#pi^{-};#Delta p = p_{FD} - p_{CD} [GeV];Events",
        200, -1.0, 1.0
    );

    TH1F* h_dtheta_pim = new TH1F(
        "h_dtheta_pim",
        "#pi^{-};#Delta#theta = #theta_{FD} - #theta_{CD} [deg];Events",
        200, -30.0, 30.0
    );

    TH1F* h_dphi_pim = new TH1F(
        "h_dphi_pim",
        "#pi^{-};#Delta#phi = #phi_{FD} - #phi_{CD} [deg];Events",
        200, -60.0, 60.0
    );

    // ============================================================
    // Event loop
    // ============================================================

    Long64_t counter = 0;

    while (reader.Next())
    {
        counter++;

        if (counter % 1000000 == 0)
        {
            std::cout << "\rProcessed: "
                      << counter
                      << std::flush;
        }

        // ========================================================
        // Proton
        // ========================================================

        bool pFD = IsFD(*status_p);
        bool pCD = IsCD(*status_p);

        /*
         * В текущем TTree у нас одна запись proton.
         *
         * Поэтому здесь есть принципиальная проблема:
         *
         *   pFD == true
         *   pCD == true
         *
         * означает, что status говорит о наличии обоих
         * detector contributions у ОДНОЙ REC::Particle.
         *
         * Это НЕ две независимые реконструкции:
         *
         *   p_FD
         *   p_CD
         *
         * Поэтому Delta p/theta/phi для CD-FD пары
         * из этого TTree вычислить нельзя.
         */

        if (pFD && pCD)
        {
            // Здесь нет отдельных p_FD и p_CD.
            //
            // Поэтому не заполняем гистограммы.
        }


        // ========================================================
        // Pi-
        // ========================================================

        bool pimFD = IsFD(*status_pim);
        bool pimCD = IsCD(*status_pim);

        /*
         * Аналогично для pi-.
         *
         * Одна REC::Particle может иметь одновременно
         * FD и CD contributions, но это не означает,
         * что в TTree находятся две реконструкции частицы.
         */
    }

    std::cout << "\nProcessed: "
              << counter
              << " events" << std::endl;


    // ============================================================
    // Drawing
    // ============================================================

    gStyle->SetOptStat(0);

    TCanvas* c = new TCanvas(
        "c",
        "CD-FD track removal",
        1000,
        800
    );

    c->Print(Form("%s[", outputPDF));

    // Proton
    h_dp_p->Draw();
    c->Print(outputPDF);

    h_dtheta_p->Draw();
    c->Print(outputPDF);

    h_dphi_p->Draw();
    c->Print(outputPDF);

    // Pi-
    h_dp_pim->Draw();
    c->Print(outputPDF);

    h_dtheta_pim->Draw();
    c->Print(outputPDF);

    h_dphi_pim->Draw();
    c->Print(outputPDF);

    c->Print(Form("%s]", outputPDF));

    std::cout << "Saved: "
              << outputPDF
              << std::endl;

    // ============================================================
    // Cleanup
    // ============================================================

    delete c;

    delete h_dp_p;
    delete h_dtheta_p;
    delete h_dphi_p;

    delete h_dp_pim;
    delete h_dtheta_pim;
    delete h_dphi_pim;

    infile->Close();
    delete infile;
}