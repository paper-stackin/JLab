#include <memory>

#include "constants.h"
#include "electron_cuts.h"
#include "hadron_cuts.h"
#include "config.h"

// Конфигурация
Constants mc;
SelectionConfig cfg;

TLorentzVector MakeParticle(int index, double mass, const hipo::bank& PART)
{
    TLorentzVector p;
    p.SetXYZM(
        PART.getFloat("px", index),
        PART.getFloat("py", index),
        PART.getFloat("pz", index),
        mass
    );
    return p;
}

std::pair<double, double> Rotate(double x, double y, double angle)
{
    return {
        x * std::cos(angle) - y * std::sin(angle),
        x * std::sin(angle) + y * std::cos(angle)
    };
}

void GetDCInfo(const hipo::bank& Track,
                       const hipo::bank& Traj,
                       int part_index,
                       int i,
                       int& DC_sect,
                       double DC_x[3],
                       double DC_y[3])
{
    const int dc_layers[3] = {6, 18, 36};

    for (int k = 0; k < Track.getRows(); k++)
    {
        bool in_DC = (Track.getInt("pindex", k) == part_index &&
                         Track.getInt("detector", k) == 6);

        if (!in_DC) continue;

        int layer = Traj.getInt("layer", i);
        float x = Traj.getFloat("x", i);
        float y = Traj.getFloat("y", i);

        if (layer == dc_layers[0])
            DC_sect = Track.getInt("sector", k);

        for (int j = 0; j < 3; j++)
        {
            if (layer == dc_layers[j])
            {
                DC_x[j] = x;
                DC_y[j] = y;
                break;
            }
        }
    }
}

void selection()
{
    float px_e, py_e, pz_e; // Final electron
    float px_p, py_p, pz_p; // Final proton
    float px_pim, py_pim, pz_pim;   // Pi minus
    float weight;   // Weight (used in simulation)

    // Initializing TTree
    TTree tree(cfg.file.tree_name.c_str(), cfg.file.tree_title.c_str());

    // Final electron
    tree.Branch("px_e", &px_e, "px_e");
    tree.Branch("py_e", &py_e, "py_e");
    tree.Branch("pz_e", &pz_e, "pz_e");
    // Final proton
    tree.Branch("px_p", &px_p, "px_p");
    tree.Branch("py_p", &py_p, "py_p");
    tree.Branch("pz_p", &pz_p, "pz_p");
    // Pi minus
    tree.Branch("px_pim", &px_pim, "px_pim");
    tree.Branch("py_pim", &py_pim, "py_pim");
    tree.Branch("pz_pim", &pz_pim, "pz_pim");

    if(cfg.use_weights) {
        tree.Branch("weight", &weight, "weight");
    } 
    // tree.Branch("code", &code, "code/I"); /I - for integer

    // Forming data file name
    std::vector<std::string> data;
    for (int i = cfg.file.input_start; i <= cfg.file.input_end; i++)
    {
        data.push_back(cfg.file.input_prefix + std::__cxx11::to_string(i) + ".hipo");
    } 

    // Reading data
    for(int k = 0; k < data.size(); k++)
    {
        hipo::reader reader;        
        std::ifstream infile(data[k]);
        if (!infile.good()) continue;
        reader.open(data[k].c_str());

        std::cout << std::endl << std::endl;

        hipo::dictionary factory;
        reader.readDictionary(factory);

        // Reading banks
        hipo::event event;
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank Track(factory.getSchema("REC::Track"));
        hipo::bank Traj(factory.getSchema("REC::Traj"));
        hipo::bank Scint(factory.getSchema("REC::Scintillator"));
        hipo::bank ECAL(factory.getSchema("REC::Calorimeter"));
        std::unique_ptr<hipo::bank> mc_bank;
        if (cfg.use_weights) {
            mc_bank = std::make_unique<hipo::bank>(factory.getSchema("MC::Event"));
        }

        int totalEvents = reader.getEntries() / cfg.counter.speed;
        int counter = 0;

        //Analyzing events
        while(reader.next())
        {
            counter++;
            if (counter % cfg.counter.speed == 0) {
                cout << "\rFile: " << k + 1 << "/" << data.size() << 
                "        Processed: " << counter / cfg.counter.speed << "/" << totalEvents << cfg.counter.period << " events" << flush;
            }

            reader.read(event);
            event.getStructure(PART);
            event.getStructure(Track);
            event.getStructure(Traj);
            event.getStructure(Scint);
            event.getStructure(ECAL);
            if (cfg.use_weights) {
                event.getStructure(*mc_bank);
            }

            int part_rows = PART.getRows();
    
            if (part_rows == 0) continue;
            
            int el_cnt = 0, p_cnt = 0, pip_cnt = 0, pim_cnt = 0, unid_cnt = 0;
            int charged_cnt = 0;
            int el_index, p_index, pim_index, pip_index; 

            for(int i = 0; i < part_rows; i++) 
            {
                int part_id = PART.getInt("pid", i);

                if (part_id == mc.e_id) {el_index = i; el_cnt++;}
                else if (part_id == mc.p_id) {p_index = i; p_cnt++;}
                else if (part_id == mc.pip_id) {pip_index = i; pip_cnt++;}
                else if (part_id == mc.pim_id) {pim_index = i; pim_cnt++;}
                else if (part_id == mc.unid_id) {unid_cnt++;}   // Unidentified particles counter
                if (PART.getInt("charge", i)) charged_cnt++;        
            }

            if (TopologyCut(cfg.topology, el_cnt, p_cnt, pip_cnt, pim_cnt, unid_cnt, charged_cnt)) continue; // Topology cut

            if (!ElectronIsFirst(el_index)) continue; // First registered particle should be electron

            if (!ElectronInFD(PART, el_index)) continue; // Electron should be in FD
            if (!HadronInFDOrCD(PART, p_index)) continue; // Proton should be in FD or CD
            if (!HadronInFDOrCD(PART, pim_index)) continue; // Pi minus should be in FD or CD

            // Setting registered particles parameteres
            TLorentzVector beam(0, 0, cfg.E_beam, cfg.E_beam);
            TLorentzVector p_in(0, 0, 0, mc.M_p);
            TLorentzVector  el = MakeParticle(el_index, mc.M_el, PART);
            TLorentzVector  pim = MakeParticle(pim_index, mc.M_pi, PART);
            TLorentzVector  p_fin = MakeParticle(p_index, mc.M_p, PART);

            if (DoublePionWThresholdCut(beam, p_in, el)) continue; // W > 1.2 GeV, 2 pion threshold

            double el_mom = el.Rho();

            if (ElectronChi2Cut(PART, el_index)) continue; // Chi-Squared PID cut
            if (ElectronMomentumCut(el_mom)) continue; // Minimum momentum cut

            if (ElectronVertexCut(PART, el_index)) continue;  // Vertex position cut
            
            if (ElectronTOFCut(PART, Scint, el_index)) continue;  // TOF cut

            double  E_PCAL = 0., EC_in = 0., EC_out =0.;
            int     cal_el_sect = 0;
            float    lu_el, lv_el, lw_el;
            
            for (int i = 0; i < ECAL.getRows(); i++)
            {
                if (ECAL.getInt("pindex", i) != el_index)   continue;

                cal_el_sect = ECAL.getInt("sector", i);
                int ecal_layer = ECAL.getInt("layer", i);
                float ecal_energy = ECAL.getFloat("energy", i);

                if (ecal_layer == 1) 
                {
                    E_PCAL = ecal_energy;
                    lu_el = ECAL.getFloat("lu", i);
                    lv_el = ECAL.getFloat("lv", i);
                    lw_el = ECAL.getFloat("lw", i);
                }
                else if (ecal_layer == 4)
                {
                    EC_in = ecal_energy;
                }
                else if (ecal_layer == 7)
                {
                    EC_out = ecal_energy;
                }
            }

            if (ElectronUVWCut(lu_el, lv_el, lw_el)) continue; // UVW cut for the electron

            if (PiMinusContaminationCut(EC_in, E_PCAL, el_mom)) continue;  // Pi minus contamination cut
            
            // Sampling fraction cut
            double Etot_e = E_PCAL + EC_in + EC_out;

            if (SamplingFractionCut(Etot_e, el_mom, cal_el_sect)) continue;
               
            // Chi-Squared PID cuts
            if (HadronChi2Cut(PART, p_index)) continue;
            if (HadronChi2Cut(PART, pim_index)) continue;
                
            int pim_status = PART.getInt("status", pim_index);
            bool pim_in_FD = (pim_status > 2000 && pim_status < 4000); // pi-minus in FD
                
            int p_status = PART.getInt("status", p_index);
            bool p_in_FD = (p_status > 2000 && p_status < 4000); // proton in FD

            // Hadron momentum cut
            if (HadronMomentumCut(pim.Rho(), pim_in_FD)) continue;
            if (HadronMomentumCut(p_fin.Rho(), p_in_FD)) continue;

            // Beta Vs p cut
            if (HadronBetaVsPCut(pim.Beta(), pim_in_FD)) continue;
            if (HadronBetaVsPCut(p_fin.Beta(), p_in_FD)) continue;

            // Z-vertex for the hadrons
            if (HadronVertexCut(PART, p_index)) continue;
            if (HadronVertexCut(PART, pim_index)) continue;

            // DC fiducial cut
            // Drift chamber location of final electron
            std::array<double, 3> DC_el_x = {0.0, 0.0, 0.0};
            std::array<double, 3> DC_el_y = {0.0, 0.0, 0.0};
            std::array<double, 3> DC_el_x_new = {0.0, 0.0, 0.0};
            std::array<double, 3> DC_el_y_new = {0.0, 0.0, 0.0};
            int DC_el_sect = 0;
            // Drift chamber location of piminus
            std::array<double, 3> DC_pim_x = {0.0, 0.0, 0.0};
            std::array<double, 3> DC_pim_y = {0.0, 0.0, 0.0};
            std::array<double, 3> DC_pim_x_new = {0.0, 0.0, 0.0};
            std::array<double, 3> DC_pim_y_new = {0.0, 0.0, 0.0};
            int DC_p_sect = 0;
            // Drift chamber location of final proton
            std::array<double, 3> DC_p_x = {0.0, 0.0, 0.0};
            std::array<double, 3> DC_p_y = {0.0, 0.0, 0.0};
            std::array<double, 3> DC_p_x_new = {0.0, 0.0, 0.0};
            std::array<double, 3> DC_p_y_new = {0.0, 0.0, 0.0};
            int DC_pim_sect = 0;

            for(int i = 0; i < Traj.getRows(); i++) 
            {
                int part_index = Traj.getInt("pindex", i);

                if (Traj.getInt("detector", i) != 6) continue;

                // Electron
                if(part_index == el_index)                  GetDCInfo(Track, Traj, part_index, i, DC_el_sect, DC_el_x.data(), DC_el_y.data());
                // Proton 
                if(p_in_FD && part_index == p_index)        GetDCInfo(Track, Traj, part_index, i, DC_p_sect, DC_p_x.data(), DC_p_y.data());
                // Pi- 
                if(pim_in_FD && part_index == pim_index)    GetDCInfo(Track, Traj, part_index, i, DC_pim_sect, DC_pim_x.data(), DC_pim_y.data());
            }
            
            // DC fiducial cut for electron
            double el_angle = (-60.0 * (DC_el_sect - 1)) * TMath::Pi() / 180.0;

            for (int r = 0; r < 3; ++r)
            {
                std::pair<double, double> rotated = Rotate(DC_el_x[r], DC_el_y[r], el_angle);
                DC_el_x_new[r] = rotated.first;
                DC_el_y_new[r] = rotated.second;
            }

            if (ElectronDCFiducialCut(DC_el_x_new, DC_el_y_new)) continue;
                                
            // DC cuts for the proton in R1, R2, R3
            double p_angle = (-60.0 * (DC_p_sect - 1)) * TMath::Pi() / 180.0;

            for (int r = 0; r < 3; ++r)
            {
                std::pair<double, double> rotated = Rotate(DC_p_x[r], DC_p_y[r], p_angle);
                DC_p_x_new[r] = rotated.first;
                DC_p_y_new[r] = rotated.second;
            }

            if (HadronDCFiducialCut(DC_p_x_new, DC_p_y_new, DC_p_sect, 1)) continue;
                            
            // DC for the pi- in R1, R2, R3
            double pim_angle = (-60.0 * (DC_pim_sect - 1)) * TMath::Pi() / 180.0;

            for (int r = 0; r < 3; ++r)
            {
                std::pair<double, double> rotated = Rotate(DC_pim_x[r], DC_pim_y[r], pim_angle);
                DC_pim_x_new[r] = rotated.first;
                DC_pim_y_new[r] = rotated.second;
            }

            if (HadronDCFiducialCut(DC_pim_x_new, DC_pim_y_new, DC_pim_sect, -1)) continue;
                            
            px_e = PART.getFloat("px", el_index);
            py_e = PART.getFloat("py", el_index);
            pz_e = PART.getFloat("pz", el_index);

            px_p = PART.getFloat("px", p_index);
            py_p = PART.getFloat("py", p_index);
            pz_p = PART.getFloat("pz", p_index);

            px_pim = PART.getFloat("px", pim_index);
            py_pim = PART.getFloat("py", pim_index);
            pz_pim = PART.getFloat("pz", pim_index);

            if (cfg.use_weights) {
                weight = mc_bank->getFloat("weight", 0);
            } 

            tree.Fill();
        }
    }

    TFile* hfile =  new TFile(cfg.file.output.c_str(), "RECREATE");
    tree.Write();
    hfile->Close();

    TFile* file = TFile::Open(cfg.file.output.c_str(), "READ");
    TTree* tree2 = nullptr;
    file->GetObject(cfg.file.tree_name.c_str(), tree2);
    Long64_t nentries = tree2->GetEntries();

    std::cout << std::endl << std::endl;
    std::cout << "Entries: " << nentries << std::endl;
    file->Close();

    gSystem->Exit(0);
}