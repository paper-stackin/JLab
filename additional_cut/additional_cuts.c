bool PCalFiducialCut(
    int cal_layer_e,
    int cal_sect_e,
    double cal_x_e,
    double cal_y_e
)
{
    const double p_lower[6][4] = {
                    {175.997,-2.19104,0.00698198,-9.2885e-06},
                    {219.797,-2.71269,0.00891301,-1.15516e-05},
                    {226.274,-2.84365,0.00957273,-1.23558e-05},
                    {181.685,-2.2416,0.00702503,-9.08476e-06},
                    {223.327,-2.79671,0.00935171,-1.20265e-05},
                    {152.767,-1.77109,0.004802,-5.80913e-06}
                };

    const double p_upper[6][4] = {
                    {-176.431,2.16702,-0.00665704,8.51867e-06},
                    {-214.182,2.60839,-0.00830054,1.04541e-05},
                    {-227.966,2.84059,-0.00948886,1.21912e-05},
                    {-178.192,2.23739,-0.00721147,9.56152e-06},
                    {-207.242,2.58685,-0.00845863,1.083e-05},
                    {-168.544,2.08019,-0.00634868,8.22812e-06}
                };

    if (cal_layer_e != 1) return true;

    double angle = TMath::Pi() / 180.0 * (-60.0 * (cal_sect_e - 1));
    double x = cal_x_e * TMath::Cos(angle) - cal_y_e * TMath::Sin(angle);
    double y = cal_x_e * TMath::Sin(angle) + cal_y_e * TMath::Cos(angle);

    const double y_lower = EvaluatePolynomial(p_lower[cal_sect_e - 1], x);
    const double y_upper = EvaluatePolynomial(p_upper[cal_sect_e - 1], x);

    const double r = sqrt(x * x + y * y);

    return y > y_upper || y < y_lower || r < 150;
}

struct Kinematics
{
    double W;
    double Q2;
    double ME;
    double MM_sq;

    int q;
    int w;
};

Kinematics CalculateKinematics(
    const TLorentzVector& el_initial,
    const TLorentzVector& el_final,
    const TLorentzVector& p_initial,
    const TLorentzVector& p_final,
    const TLorentzVector& pi_minus)
{
    Kinematics result;

    TLorentzVector virtual_photon = el_initial - el_final;
    TLorentzVector missing_pip =
        virtual_photon + p_initial - p_final - pi_minus;

    result.W = (virtual_photon + p_initial).M();
    result.Q2 = -virtual_photon.M2();

    TVector3 beta = (virtual_photon + p_initial).BoostVector();

    missing_pip.Boost(-beta);

    result.ME = missing_pip.E();
    result.MM_sq = missing_pip.M2();

    result.q = -1;

    if (result.Q2 >= 0.5 && result.Q2 <= 0.7)
        result.q = 0;
    else if (result.Q2 > 0.7 && result.Q2 <= 1.0)
        result.q = 1;
    else if (result.Q2 > 1.0 && result.Q2 <= 1.4)
        result.q = 2;
    else if (result.Q2 > 1.4 && result.Q2 <= 2.0)
        result.q = 3;
    else if (result.Q2 > 2.0 && result.Q2 <= 3.0)
        result.q = 4;
    else if (result.Q2 > 3.0)
        result.q = 5;

    result.w = static_cast<int>((result.W - 1.4) / 0.025);

    return result;
}

bool MECut(
    const TLorentzVector& el_initial,
    const TLorentzVector& el_final,
    const TLorentzVector& p_initial,
    const TLorentzVector& p_final,
    const TLorentzVector& pi_minus,
    const float ME_cuts[6][64]
)
{
    Kinematics kin = CalculateKinematics(
        el_initial,
        el_final,
        p_initial,
        p_final,
        pi_minus
    );

    if (kin.q < 0 || kin.q >= 6)
        return false;

    if (kin.w < 0 || kin.w >= 64)
        return false;

    if (kin.ME < 0.175 || kin.ME > ME_cuts[kin.q][kin.w])
        return false;

    return true;
}

void additional_cuts()
{
    TFile* file_in = TFile::Open("home/stepan/downloads/missing_pip_clas12.root", "READ");
    TTree* tree_in = (TTree*)file_in->Get("MMpiptree");
    
    // Read branches
    Float_t px_e, py_e, pz_e; 
    Int_t cal_sect_e, cal_layer_e; 
    Float_t cal_x_e, cal_y_e; 
    Float_t time_e, path_e; 
    
    Float_t px_p, py_p, pz_p; 
    Float_t time_p, path_p; 
    Int_t status_p;

    Float_t px_pim, py_pim, pz_pim; 
    Float_t time_pim, path_pim; 
    Int_t status_pim;

    tree_in->SetBranchAddress("px_e", &px_e); 
    tree_in->SetBranchAddress("py_e", &py_e); 
    tree_in->SetBranchAddress("pz_e", &pz_e); 
    tree_in->SetBranchAddress("cal_sect_e", &cal_sect_e); 
    tree_in->SetBranchAddress("cal_layer_e", &cal_layer_e); 
    tree_in->SetBranchAddress("cal_x_e", &cal_x_e); 
    tree_in->SetBranchAddress("cal_y_e", &cal_y_e); 
    tree_in->SetBranchAddress("time_e", &time_e); 
    tree_in->SetBranchAddress("path_e", &path_e);

    tree_in->SetBranchAddress("px_p", &px_p); 
    tree_in->SetBranchAddress("py_p", &py_p); 
    tree_in->SetBranchAddress("pz_p", &pz_p); 
    tree_in->SetBranchAddress("time_p", &time_p); 
    tree_in->SetBranchAddress("path_p", &path_p); 
    tree_in->SetBranchAddress("status_p", &status_p); 

    tree_in->SetBranchAddress("px_pim", &px_pim); 
    tree_in->SetBranchAddress("py_pim", &py_pim); 
    tree_in->SetBranchAddress("pz_pim", &pz_pim); 
    tree_in->SetBranchAddress("time_pim", &time_pim); 
    tree_in->SetBranchAddress("path_pim", &path_pim); 
    tree_in->SetBranchAddress("status_pim", &status_pim);
    
    // Output
    TFile* file_out = TFile::Open("output.root", "RECREATE");
    TTree* tree_out = tree_in->CloneTree(0);

    // Event loop
    TLorentzVector el_initial(0, 0, 6.535, 6.535);
    TLorentzVector p_initial(0, 0, 0, 0.938);

    Long64_t n_entries = tree_in->GetEntries();

    for (Long64_t i = 0; i < n_entries; i++)
    {
        tree_in->GetEntry(i);

        el_final.SetXYZM(px_e, py_e, pz_e, 0.);
        p_final.SetXYZM(px_p, py_p, pz_p, 0.938);
        pi_minus.SetXYZM(px_pim, py_pim, pz_pim, 0.139);

        // Additional cuts
        if(
            PCalFiducialCut(
                cal_layer_e,
                cal_sect_e,
                cal_x_e,
                cal_y_e
            )
        )
        {
            continue;
        }

        float ME_cuts[6][64];

        std::ifstream f("ME_cuts_08.txt");

        for (int i = 0; i < 6; ++i)
        {
            std::string line;
            std::getline(f, line);

            std::istringstream iss(line);

            for (int j = 0; j < 64; ++j)
                iss >> ME_cuts[i][j];
        }

        f.close();

        if (!MECut(
                el_initial,
                el_final,
                p_initial,
                p_final,
                pi_minus,
                ME_cuts))
        {
            continue;
        }

        tree_out->Fill();
    }


    
    // Save
    

    file_out->cd();
    tree_out->Write();

    file_out->Close();
    file_in->Close();

    std::cout << "Input events:  " << n_entries << std::endl;
    std::cout << "Output events: " << tree_out->GetEntries() << std::endl;
}