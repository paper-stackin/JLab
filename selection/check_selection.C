#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TStyle.h>

#include <iostream>
#include <cmath>

void check_selection()
{
    const char* input_file = "MPPT_events_simulation.root";
    const double Ebeam = 6.535; // GeV

    TFile* file = TFile::Open(input_file, "READ");

    if (!file || file->IsZombie())
    {
        std::cerr << "Error: cannot open " << input_file << std::endl;
        return;
    }

    TTree* tree = nullptr;
    file->GetObject("MMpiptree", tree);

    if (!tree)
    {
        std::cerr << "Error: TTree MMpiptree not found" << std::endl;
        file->Close();
        return;
    }

    // -------------------------
    // Branches
    // -------------------------

    float px_e, py_e, pz_e;
    float px_p, py_p, pz_p;
    float px_pim, py_pim, pz_pim;
    float weight;

    tree->SetBranchAddress("px_e",   &px_e);
    tree->SetBranchAddress("py_e",   &py_e);
    tree->SetBranchAddress("pz_e",   &pz_e);

    tree->SetBranchAddress("px_p",   &px_p);
    tree->SetBranchAddress("py_p",   &py_p);
    tree->SetBranchAddress("pz_p",   &pz_p);

    tree->SetBranchAddress("px_pim", &px_pim);
    tree->SetBranchAddress("py_pim", &py_pim);
    tree->SetBranchAddress("pz_pim", &pz_pim);

    tree->SetBranchAddress("weight", &weight);

    // -------------------------
    // Masses, GeV
    // -------------------------

    const double m_e   = 0.000511;
    const double m_p   = 0.938272;
    const double m_pim = 0.139570;

    // -------------------------
    // Histograms
    // -------------------------

    TH1D* h_p_e = new TH1D(
        "h_p_e",
        "Electron momentum; p_{e'} [GeV]; Events",
        100, 0, 7
    );

    TH1D* h_p_p = new TH1D(
        "h_p_p",
        "Proton momentum; p_{p'} [GeV]; Events",
        100, 0, 5
    );

    TH1D* h_p_pim = new TH1D(
        "h_p_pim",
        "#pi^{-} momentum; p_{#pi^{-}} [GeV]; Events",
        100, 0, 5
    );

    TH1D* h_theta_e = new TH1D(
        "h_theta_e",
        "Electron polar angle; #theta_{e'} [deg]; Events",
        90, 0, 90
    );

    TH1D* h_theta_p = new TH1D(
        "h_theta_p",
        "Proton polar angle; #theta_{p'} [deg]; Events",
        90, 0, 180
    );

    TH1D* h_theta_pim = new TH1D(
        "h_theta_pim",
        "#pi^{-} polar angle; #theta_{#pi^{-}} [deg]; Events",
        90, 0, 180
    );

    TH1D* h_phi_e = new TH1D(
        "h_phi_e",
        "Electron azimuthal angle; #phi_{e'} [deg]; Events",
        72, -180, 180
    );

    TH1D* h_phi_p = new TH1D(
        "h_phi_p",
        "Proton azimuthal angle; #phi_{p'} [deg]; Events",
        72, -180, 180
    );

    TH1D* h_phi_pim = new TH1D(
        "h_phi_pim",
        "#pi^{-} azimuthal angle; #phi_{#pi^{-}} [deg]; Events",
        72, -180, 180
    );

    TH1D* h_W = new TH1D(
        "h_W",
        "Invariant mass W; W [GeV]; Events",
        100, 1.0, 3.0
    );

    TH1D* h_Q2 = new TH1D(
        "h_Q2",
        "Four-momentum transfer; Q^{2} [GeV^{2}]; Events",
        100, 0, 6
    );

    TH1D* h_M_p_pim = new TH1D(
        "h_M_p_pim",
        "p'#pi^{-} invariant mass; M(p'#pi^{-}) [GeV]; Events",
        100, 1.0, 3.0
    );

    TH1D* h_M_ep_pim = new TH1D(
        "h_M_ep_pim",
        "e'p'#pi^{-} invariant mass; M(e'p'#pi^{-}) [GeV]; Events",
        100, 0, 4
    );

    // -------------------------
    // Event loop
    // -------------------------

    Long64_t nentries = tree->GetEntries();

    std::cout << "Entries: " << nentries << std::endl;

    TLorentzVector beam(0, 0, Ebeam, Ebeam);
    TLorentzVector target(0, 0, 0, m_p);

    for (Long64_t i = 0; i < nentries; i++)
    {
        tree->GetEntry(i);

        TLorentzVector electron;
        TLorentzVector proton;
        TLorentzVector pim;

        electron.SetXYZM(px_e, py_e, pz_e, m_e);
        proton.SetXYZM(px_p, py_p, pz_p, m_p);
        pim.SetXYZM(px_pim, py_pim, pz_pim, m_pim);

        // -------------------------
        // Basic kinematics
        // -------------------------

        double p_e   = electron.P();
        double p_p   = proton.P();
        double p_pim = pim.P();

        double theta_e   = electron.Theta() * 180.0 / M_PI;
        double theta_p   = proton.Theta()   * 180.0 / M_PI;
        double theta_pim = pim.Theta()      * 180.0 / M_PI;

        double phi_e   = electron.Phi() * 180.0 / M_PI;
        double phi_p   = proton.Phi()   * 180.0 / M_PI;
        double phi_pim = pim.Phi()      * 180.0 / M_PI;

        h_p_e->Fill(p_e, weight);
        h_p_p->Fill(p_p, weight);
        h_p_pim->Fill(p_pim, weight);

        h_theta_e->Fill(theta_e, weight);
        h_theta_p->Fill(theta_p, weight);
        h_theta_pim->Fill(theta_pim, weight);

        h_phi_e->Fill(phi_e, weight);
        h_phi_p->Fill(phi_p, weight);
        h_phi_pim->Fill(phi_pim, weight);

        // -------------------------
        // Q2 and W
        // -------------------------

        TLorentzVector virtualPhoton = beam - electron;

        double Q2 = -virtualPhoton.M2();
        double W  = (target + virtualPhoton).M();

        if (Q2 > 0)
            h_Q2->Fill(Q2, weight);

        if (W > 0)
            h_W->Fill(W, weight);

        // -------------------------
        // Invariant masses
        // -------------------------

        double M_p_pim = (proton + pim).M();
        double M_ep_pim = (electron + proton + pim).M();

        h_M_p_pim->Fill(M_p_pim, weight);
        h_M_ep_pim->Fill(M_ep_pim, weight);
    }

    // -------------------------
    // Draw
    // -------------------------

    gStyle->SetOptStat(1110);

    TCanvas* c1 = new TCanvas("c1", "Basic distributions OG", 1200, 800);
    c1->Divide(3, 3);

    c1->cd(1);
    h_p_e->Draw();

    c1->cd(2);
    h_p_p->Draw();

    c1->cd(3);
    h_p_pim->Draw();

    c1->cd(4);
    h_theta_e->Draw();

    c1->cd(5);
    h_theta_p->Draw();

    c1->cd(6);
    h_theta_pim->Draw();

    c1->cd(7);
    h_phi_e->Draw();

    c1->cd(8);
    h_phi_p->Draw();

    c1->cd(9);
    h_phi_pim->Draw();

    c1->SaveAs("basic_distributions.pdf");

    TCanvas* c2 = new TCanvas("c2", "Invariant distributions OG", 1200, 800);
    c2->Divide(2, 2);

    c2->cd(1);
    h_W->Draw();

    c2->cd(2);
    h_Q2->Draw();

    c2->cd(3);
    h_M_p_pim->Draw();

    c2->cd(4);
    h_M_ep_pim->Draw();

    c2->SaveAs("invariant_distributions.pdf");

    // -------------------------
    // Save histograms
    // -------------------------

    TFile* output = new TFile("distributions.root", "RECREATE");

    h_p_e->Write();
    h_p_p->Write();
    h_p_pim->Write();

    h_theta_e->Write();
    h_theta_p->Write();
    h_theta_pim->Write();

    h_phi_e->Write();
    h_phi_p->Write();
    h_phi_pim->Write();

    h_W->Write();
    h_Q2->Write();

    h_M_p_pim->Write();
    h_M_ep_pim->Write();

    output->Close();
    file->Close();

    std::cout << "Done." << std::endl;
}