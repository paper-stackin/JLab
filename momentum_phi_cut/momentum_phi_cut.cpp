#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <cmath>
#include <iostream>

bool IsFD(int status)
{
    status = std::abs(status);

    // FD contribution
    return ((status / 2000) % 2) == 1;
}

bool IsCD(int status)
{
    status = std::abs(status);

    // CD contribution
    return ((status / 4000) % 2) == 1;
}

void momentum_phi_cut()
{
    // ============================================================
    // Input
    // ============================================================

    const char* inputFile = "/home/stepan/downloads/missing_pip.root";
    const char* treeName  = "MMpiptree";
    const char* outputPDF = "pt_vs_phi.pdf";

    TFile* infile = TFile::Open(inputFile);

    if (!infile || infile->IsZombie())
    {
        std::cerr << "Error: cannot open file "
                  << inputFile << std::endl;
        return;
    }

    TTreeReader reader(treeName, infile);

    // ============================================================
    // Proton
    // ============================================================

    TTreeReaderValue<Float_t> px_p(reader, "px_p");
    TTreeReaderValue<Float_t> py_p(reader, "py_p");
    TTreeReaderValue<Float_t> pz_p(reader, "pz_p");
    TTreeReaderValue<Int_t>   status_p(reader, "status_p");

    // ============================================================
    // Pi-
    // ============================================================

    TTreeReaderValue<Float_t> px_pim(reader, "px_pim");
    TTreeReaderValue<Float_t> py_pim(reader, "py_pim");
    TTreeReaderValue<Float_t> pz_pim(reader, "pz_pim");
    TTreeReaderValue<Int_t>   status_pim(reader, "status_pim");

    // ============================================================
    // Histograms
    // ============================================================

    // phi: -180 ... 180 degrees
    // pT:  0 ... 2 GeV
    //
    // При необходимости потом легко изменить диапазон.
    
    TH2F* h_pt_phi_p_FD = new TH2F(
        "h_pt_phi_p_FD",
        "Proton FD;#phi [deg];p_{T} [GeV]",
        180, -180, 180,
        200, 0, 2
    );

    TH2F* h_pt_phi_p_CD = new TH2F(
        "h_pt_phi_p_CD",
        "Proton CD;#phi [deg];p_{T} [GeV]",
        180, -180, 180,
        200, 0, 2
    );

    TH2F* h_pt_phi_pim_FD = new TH2F(
        "h_pt_phi_pim_FD",
        "#pi^{-} FD;#phi [deg];p_{T} [GeV]",
        180, -180, 180,
        200, 0, 2
    );

    TH2F* h_pt_phi_pim_CD = new TH2F(
        "h_pt_phi_pim_CD",
        "#pi^{-} CD;#phi [deg];p_{T} [GeV]",
        180, -180, 180,
        200, 0, 2
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

        // --------------------------------------------------------
        // Proton
        // --------------------------------------------------------

        double pt_p = std::sqrt(
            (*px_p) * (*px_p) +
            (*py_p) * (*py_p)
        );

        double phi_p = std::atan2(
            *py_p,
            *px_p
        ) * 180.0 / M_PI;

        if (IsFD(*status_p))
        {
            h_pt_phi_p_FD->Fill(phi_p, pt_p);
        }

        if (IsCD(*status_p))
        {
            h_pt_phi_p_CD->Fill(phi_p, pt_p);
        }

        // --------------------------------------------------------
        // Pi-
        // --------------------------------------------------------

        double pt_pim = std::sqrt(
            (*px_pim) * (*px_pim) +
            (*py_pim) * (*py_pim)
        );

        double phi_pim = std::atan2(
            *py_pim,
            *px_pim
        ) * 180.0 / M_PI;

        if (IsFD(*status_pim))
        {
            h_pt_phi_pim_FD->Fill(phi_pim, pt_pim);
        }

        if (IsCD(*status_pim))
        {
            h_pt_phi_pim_CD->Fill(phi_pim, pt_pim);
        }
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
        "pT vs phi",
        1000,
        800
    );

    c->Print(Form("%s[", outputPDF));

    // ------------------------------------------------------------
    // Proton FD
    // ------------------------------------------------------------

    h_pt_phi_p_FD->Draw("COLZ");
    c->Print(outputPDF);

    // ------------------------------------------------------------
    // Proton CD
    // ------------------------------------------------------------

    h_pt_phi_p_CD->Draw("COLZ");
    c->Print(outputPDF);

    // ------------------------------------------------------------
    // Pi- FD
    // ------------------------------------------------------------

    h_pt_phi_pim_FD->Draw("COLZ");
    c->Print(outputPDF);

    // ------------------------------------------------------------
    // Pi- CD
    // ------------------------------------------------------------

    h_pt_phi_pim_CD->Draw("COLZ");
    c->Print(outputPDF);

    c->Print(Form("%s]", outputPDF));

    std::cout << "Saved: " << outputPDF << std::endl;

    // ============================================================
    // Cleanup
    // ============================================================

    delete c;

    delete h_pt_phi_p_FD;
    delete h_pt_phi_p_CD;
    delete h_pt_phi_pim_FD;
    delete h_pt_phi_pim_CD;

    infile->Close();
    delete infile;

    gSystem->Exit(0);
}