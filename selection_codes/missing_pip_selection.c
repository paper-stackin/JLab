#include <cstdlib>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TH1.h>
#include <TLine.h>
#include <TH2.h>
#include <TChain.h>
#include <TCanvas.h>
#include "clas12reader.h"
#include <string>
#include <TMath.h>

using namespace clas12;

void SetLorentzVector(TLorentzVector &p4,clas12::region_part_ptr rp)
{
  p4.SetXYZM(rp->par()->getPx(),rp->par()->getPy(),rp->par()->getPz(),p4.M());
}

void missing_pip_selection()
{
    float px_e, py_e, pz_e, px_p, py_p, pz_p, px_pim, py_pim, pz_pim, weight;
    int ppindex, pimpindex, pippindex;

    float Etot_e, mm2, DCchi2_e, DCndf_e, DCechi2;
    bool lu1, lv1, lw1;
    float DCexR1, DCeyR1, DCexR2, DCeyR2, DCexR3, DCeyR3;
    float DCexR1_new, DCeyR1_new, DCexR2_new, DCeyR2_new, DCexR3_new, DCeyR3_new;

    float DCpimxR1, DCpimyR1, DCpimxR2, DCpimyR2, DCpimxR3, DCpimyR3, DCpimchi2;
    float DCpimxR1_new, DCpimyR1_new;
    float DCpimxR2_new, DCpimyR2_new;
    float DCpimxR3_new, DCpimyR3_new;

    float DCpxR1, DCpyR1, DCpxR2, DCpyR2, DCpxR3, DCpyR3, DCpchi2;
    float DCpxR1_new, DCpyR1_new;
    float DCpxR2_new, DCpyR2_new;
    float DCpxR3_new, DCpyR3_new;

    int DCpndf, DCpimndf, code;

    auto Missing_pi_plus_hist = new TH1F("Missing_pi_plus_hist","Missing_pi_plus_hist",100,-0.3,0.7);
    auto W_hist = new TH1F("W_hist","W_hist",100,0,4);
    auto Q2_hist = new TH1F("Q2_hist","Q2_hist",100,0,6);
    auto Q2_vs_W_hist = new TH2F("Q2_vs_W_hist","Q2_vs_W_hist",100,0,4,100,0,6);
    auto MM_pi_plus_hist = new TH1F("MM_pi_plus_hist","MM_pi_plus_hist",100,-0.3,0.7);

    TH1F   *MM[12][70];
    for (int q=0;q<6;++q) 
    {
        for (int w=0;w<64;++w)
        {
            char namehist[256];
            sprintf(namehist,"MM_Q2_bin=%d_W_bin=%d",q+1,w+1);
            MM[q][w] = new TH1F (namehist,namehist,100,-0.3,0.7);
        }
    }

    TTree tree("MMpiptree","60nA 6.535 RGK runs MM0 topology");
    //tree.Branch("event",&particle,"px:py:pz");
    tree.Branch("px_e",&px_e,"px_e");
    tree.Branch("py_e",&py_e,"py_e");
    tree.Branch("pz_e",&pz_e,"pz_e");

    tree.Branch("px_p",&px_p,"px_p");
    tree.Branch("py_p",&py_p,"py_p");
    tree.Branch("pz_p",&pz_p,"pz_p");

    tree.Branch("px_pim",&px_pim,"px_pim");
    tree.Branch("py_pim",&py_pim,"py_pim");
    tree.Branch("pz_pim",&pz_pim,"pz_pim");

    //tree.Branch("weight",&weight,"weight");

    tree.Branch("code",&code,"code/I");

    double Mp = 0.93827, Mpi = 0.13957;
    TLorentzVector beam(0, 0, 6.535, 6.535);
    TLorentzVector pin(0,0,0,Mp);
    TLorentzVector el, pim, pfin, prec, piprec, pimrec, Sum;

    int counter=0;
  
    std::cout << " reading file example program (HIPO) "  << __cplusplus << std::endl;


    std::vector<string> data;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////   PASS 2

    string str = "/cache/clas12/rg-k/production/recon/fall2018/torus+1/6535MeV/pass2/v0/dst/train/skim30/skim30_00";
    string fin = ".hipo";
    for (int i = 5860; i <= 6000; i++) data.push_back(str + std::__cxx11::to_string(i) + fin);

    /////////////////////////////////////////////////////////////////////////////////////////////////////////////

    for(int r=0;r<data.size();r++)
    {
        hipo::reader reader;
        std::ifstream infile(data[r]);
        if(! infile.good()) continue;
        reader.open(data[r].c_str());

        hipo::dictionary factory;
        reader.readDictionary(factory);

        hipo::event event;
        hipo::bank PART(factory.getSchema("REC::Particle"));
        hipo::bank Track(factory.getSchema("REC::Track"));
        hipo::bank Traj(factory.getSchema("REC::Traj"));
        hipo::bank Scint(factory.getSchema("REC::Scintillator"));
        hipo::bank ECAL(factory.getSchema("REC::Calorimeter"));
        //hipo::bank MC(factory.getSchema("MC::Event"));

        while(reader.next()==true)
        {
            counter++;
            if (counter % 100000 == 0) cout<<counter<<"\n";

            reader.read(event);
            event.getStructure(PART);
            event.getStructure(Track);
            event.getStructure(Traj);
            event.getStructure(Scint);
            event.getStructure(ECAL);
            //event.getStructure(MC);

            int nrows = PART.getRows(), fd = 0, cd = 0, charge = 0;
    
            if (nrows == 0/* || Traj.getRows() == 0 || Track.getRows() == 0*/) continue;
        
            std::vector<int> pid;
            int pcount = 0, pipcount = 0, pimcount = 0, unid = 0;

            bool elstatus = (2000 <= abs(PART.getInt("status",0))) && (abs(PART.getInt("status",0)) < 4000);
    
            for(int i=0; i<nrows; i++) 
            {
                pid.push_back(PART.getInt("pid",i));
                if (PART.getInt("pid",i) == 2212) {ppindex =i; pcount++;}
                if (PART.getInt("pid",i) == 211) {pippindex = i; pipcount++;}
                if (PART.getInt("pid",i) == -211) {pimpindex = i; pimcount++;}
                if (PART.getInt("pid",i) == 0) { unid++;}
                if(PART.getInt("charge",i)) charge ++;        
            }

            for(int i=0; i<nrows; i++) 
            {
                if (2000 <= abs(PART.getInt("status",i)) && abs(PART.getInt("status",i)) < 4000) fd++;
                if (4000 <= abs(PART.getInt("status",i)) && abs(PART.getInt("status",i)) < 8000) cd++;
            }

            if(pid[0]==11 && elstatus && fd+cd == nrows && pcount == 1 && pipcount == 0 && unid == 0 && pimcount == 1 && charge == 3) // missing pi plus
            {

                /*
                cout<<"Event "<<endl;
                for (int i = 0; i < pid.size(); ++i)
                {
                    cout<<pid[i]<<endl; 
                }
                */

                el.SetXYZM(PART.getFloat("px",0),PART.getFloat("py",0),PART.getFloat("pz",0),0);
                pim.SetXYZM(PART.getFloat("px",pimpindex),PART.getFloat("py",pimpindex),PART.getFloat("pz",pimpindex),Mpi);
                pfin.SetXYZM(PART.getFloat("px",ppindex),PART.getFloat("py",ppindex),PART.getFloat("pz",ppindex),Mp);

                if (PART.getFloat("chi2pid",0) >= 5) continue;
                if (el.Rho() < 1.0 || el.Rho() > 6.535) continue; //final electron energy > 1 GeV
                if((beam+pin-el).M()< 1.2) continue;// W > 1.2 GeV, 2 pion threshold

                double vz_el = PART.getFloat("vz",0 );// z-vertex
                if(vz_el < -10 || vz_el > 2) continue;// z-vertex cut
            
                double hit_el = 0.; // time of hit
                for (int i = 0; i < Scint.getRows(); i++)
                {
                    if (Scint.getInt("pindex",i) == 0 && Scint.getFloat("time",i) > hit_el) hit_el = Scint.getFloat("time",i);
                }
                double vt_el = PART.getFloat("vt",0);// vertex time
                double tof_e = hit_el - vt_el;
                if (tof_e <= 21 || tof_e >= 26 ) continue;// TOF cut

                double E_PCAL = 0., EC_in = 0., EC_out =0.;
                int CALesect = 0;
                for (int i = 0; i < ECAL.getRows(); i++)
                {
                    if (ECAL.getInt("pindex",i) == 0)
                    {
                        CALesect = ECAL.getInt("sector",i);
                        if (ECAL.getInt("layer",i) == 1) 
                        {
                            E_PCAL = ECAL.getFloat("energy",i);
                            lu1 = ECAL.getFloat("lu",i) > 40 && ECAL.getFloat("lu",i) < 400;
                            lv1 = ECAL.getFloat("lv",i) > 15;
                            lw1 = ECAL.getFloat("lw",i) > 15;
                        }
                        if (ECAL.getInt("layer",i) == 4)
                        {
                            EC_in = ECAL.getFloat("energy",i);
                        }
                        if (ECAL.getInt("layer",i) == 7)
                        {
                            EC_out = ECAL.getFloat("energy",i);
                        }
                    }
                }
                if (!(lu1 && lv1 && lw1)) continue;
                if ( (EC_in/el.Rho()) < (-0.84 * E_PCAL/el.Rho() + 0.17) ) continue; // pi minus contamination
            
                //sampling fraction cut:
                Etot_e = E_PCAL + EC_in + EC_out;
                bool SFcut = true;
            
                switch (CALesect)
                {
                    case 1:
                    {
                        if(Etot_e/el.Rho() < (0.145 + 0.0216*el.Rho()-0.00181*el.Rho()*el.Rho())) SFcut = false;
                    }
                    break;
                    case 2:
                    {
                        if(Etot_e/el.Rho() < (0.134 + 0.023*el.Rho()-0.00168*el.Rho()*el.Rho())) SFcut = false;
                    }
                    break;
                    case 3:
                    {
                        if(Etot_e/el.Rho() < (0.145 + 0.0211*el.Rho()-0.00166*el.Rho()*el.Rho())) SFcut = false;
                    }
                    break;
                    case 4:
                    {
                        if(Etot_e/el.Rho() < (0.152 + 0.0161*el.Rho()-0.00127*el.Rho()*el.Rho())) SFcut = false;
                    }
                    break;
                    case 5:
                    {
                        if(Etot_e/el.Rho() < (0.141 + 0.021*el.Rho()-0.00174*el.Rho()*el.Rho())) SFcut = false;
                    }
                    break;
                    case 6:
                    {
                        if(Etot_e/el.Rho() < (0.141 + 0.02152*el.Rho()-0.0017*el.Rho()*el.Rho())) SFcut = false;
                    }
                    break;
                }
                if (!SFcut || !CALesect) continue;
               
                if (PART.getFloat("chi2pid",ppindex) >= 5) continue;
                if (PART.getFloat("chi2pid",pimpindex) >= 5) continue;

                code = 0;
                
                int pimstatus = PART.getInt("status",pimpindex);
                if (pimstatus > 2000 && pimstatus < 4000) code += 10; // pi-minus in FD
                
                int pstatus = PART.getInt("status",ppindex);
                if (pstatus > 2000 && pstatus < 4000) code += 100; // proton in FD

                // Hadron momentum cut
                if ((code /10 % 10) && pim.Rho() < 0.4) continue; //FD case pim
                if (!(code /10 % 10) && pim.Rho() < 0.2) continue; //СD case pim
                if ((code /100) && pfin.Rho() < 0.4) continue; //FD case pfin
                if (!(code /100) && pfin.Rho() < 0.2) continue; //СD case pfin

                // Beta Vs p cut
                if ((code /10 % 10) && ((pim.Beta() < 0.4) || (pim.Beta() > 1.1))) continue; //FD case pim
                if (!(code /10 % 10) && ((pim.Beta() < 0.2) || (pim.Beta() > 1.1))) continue; //CD case pim
                if ((code /100) && ((pfin.Beta() < 0.4) || (pfin.Beta() > 1.1))) continue; //FD case pfin
                if (!(code /100) && ((pfin.Beta() < 0.2) || (pfin.Beta() > 1.1))) continue; //CD case pfin

                // Z-vertex for the hadrons
                double vz_p = PART.getFloat("vz", ppindex);
                if(vz_p < -10 || vz_p > 2) continue;
                double vz_pim = PART.getFloat("vz", pimpindex);
                if(vz_pim < -10 || vz_pim > 2) continue;
                
                px_e = PART.getFloat("px",0);
                py_e = PART.getFloat("py",0);
                pz_e = PART.getFloat("pz",0);

                px_p = PART.getFloat("px",ppindex);
                py_p = PART.getFloat("py",ppindex);
                pz_p = PART.getFloat("pz",ppindex);

                px_pim = PART.getFloat("px",pimpindex);
                py_pim = PART.getFloat("py",pimpindex);
                pz_pim = PART.getFloat("pz",pimpindex);

                //weight = MC.getFloat("weight",0);

                //DCfiducial cut
                int DCesect = 0, DCpsect = 0, DCpimsect = 0;
                for(int i=0; i<Traj.getRows(); i++) 
                {
                    if(Traj.getInt("pindex",i) == 0 && Traj.getInt("detector",i) == 6)
                    {    
                        for (int k=0; k<Track.getRows(); k++)
                        {
                            if(Track.getInt("pindex", k) == 0 && Track.getInt("detector",k) == 6)//DC
                            {
                                if (Traj.getInt("layer",i) == 6)
                                {
                                    DCexR1 = Traj.getFloat("x",i);
                                    DCeyR1 = Traj.getFloat("y",i);
                                    DCesect = Track.getInt("sector",k);
                                    DCndf_e = Track.getInt("NDF",k);
                                    DCchi2_e = Track.getFloat("chi2",k);
                                }
                                if (Traj.getInt("layer",i) == 18)
                                {
                                    DCexR2 = Traj.getFloat("x",i);
                                    DCeyR2 = Traj.getFloat("y",i);
                                }
                                if (Traj.getInt("layer",i) == 36)
                                {
                                    DCexR3 = Traj.getFloat("x",i);
                                    DCeyR3 = Traj.getFloat("y",i);
                                }
                            }
                        }
                    }
                
                    //proton 
                    if((code /100) && Traj.getInt("pindex",i) == ppindex && Traj.getInt("detector",i) == 6)
                    {
                        for (int k=0; k<Track.getRows(); k++)
                        {
                            if(Track.getInt("pindex", k) == ppindex && Track.getInt("detector",k) == 6)//DC
                            {
                                if (Traj.getInt("layer",i) == 6)
                                {
                                    DCpxR1 = Traj.getFloat("x",i);
                                    DCpyR1 = Traj.getFloat("y",i);
                                    DCpsect = Track.getInt("sector",k);
                                    DCpchi2 = Track.getFloat("chi2",k);
                                    DCpndf = Track.getInt("NDF",k);
                                }
                                if (Traj.getInt("layer",i) == 18)
                                {
                                    DCpxR2 = Traj.getFloat("x",i);
                                    DCpyR2 = Traj.getFloat("y",i);
                                }
                                if (Traj.getInt("layer",i) == 36)
                                {
                                    DCpxR3 = Traj.getFloat("x",i);
                                    DCpyR3 = Traj.getFloat("y",i);
                                }
                            }
                        }
                    }
                       
                    //pi- 
                    if((code /10 % 10) && Traj.getInt("pindex",i) == pimpindex && Traj.getInt("detector",i) == 6)
                    {
                        for (int k=0; k<Track.getRows(); k++)
                        {
                            if(Track.getInt("pindex", k) == pimpindex && Track.getInt("detector",k) == 6)//DC
                            {
                                if (Traj.getInt("layer",i) == 6)
                                {
                                    DCpimxR1 = Traj.getFloat("x",i);
                                    DCpimyR1 = Traj.getFloat("y",i);
                                    DCpimsect = Track.getInt("sector",k);
                                    DCpimchi2 = Track.getFloat("chi2",k);
                                    DCpimndf = Track.getInt("NDF",k);
                                }
                                if (Traj.getInt("layer",i) == 18)
                                {
                                    DCpimxR2 = Traj.getFloat("x",i);
                                    DCpimyR2 = Traj.getFloat("y",i);
                                }
                                if (Traj.getInt("layer",i) == 36)
                                {
                                    DCpimxR3 = Traj.getFloat("x",i);
                                    DCpimyR3 = Traj.getFloat("y",i);
                                }
                            }
                        }
                    }
                }
            
                //DC fiducial cut for electrons
                switch (DCesect)
                {
                    case 1:
                    {
                        DCexR1_new = DCexR1;
                        DCeyR1_new = DCeyR1;
                
                        DCexR2_new = DCexR2;
                        DCeyR2_new = DCeyR2;
                
                        DCexR3_new = DCexR3;
                        DCeyR3_new = DCeyR3;
                    }
                    break;
                    case 2:
                    {
                        DCexR1_new = DCexR1 * cos(-60 * TMath::Pi() / 180) - DCeyR1 * sin(-60 * TMath::Pi() / 180);
                        DCeyR1_new = DCexR1 * sin(-60 * TMath::Pi() / 180) + DCeyR1 * cos(-60 * TMath::Pi() / 180);

                        DCexR2_new = DCexR2 * cos(-60 * TMath::Pi() / 180) - DCeyR2 * sin(-60 * TMath::Pi() / 180);
                        DCeyR2_new = DCexR2 * sin(-60 * TMath::Pi() / 180) + DCeyR2 * cos(-60 * TMath::Pi() / 180);

                        DCexR3_new = DCexR3 * cos(-60 * TMath::Pi() / 180) - DCeyR3 * sin(-60 * TMath::Pi() / 180);
                        DCeyR3_new = DCexR3 * sin(-60 * TMath::Pi() / 180) + DCeyR3 * cos(-60 * TMath::Pi() / 180);
                    }
                    break;
                    case 3:
                    {
                        DCexR1_new = DCexR1 * cos(-120 * TMath::Pi() / 180) - DCeyR1 * sin(-120 * TMath::Pi() / 180);
                        DCeyR1_new = DCexR1 * sin(-120 * TMath::Pi() / 180) + DCeyR1 * cos(-120 * TMath::Pi() / 180);

                        DCexR2_new = DCexR2 * cos(-120 * TMath::Pi() / 180) - DCeyR2 * sin(-120 * TMath::Pi() / 180);
                        DCeyR2_new = DCexR2 * sin(-120 * TMath::Pi() / 180) + DCeyR2 * cos(-120 * TMath::Pi() / 180);

                        DCexR3_new = DCexR3 * cos(-120 * TMath::Pi() / 180) - DCeyR3 * sin(-120 * TMath::Pi() / 180);
                        DCeyR3_new = DCexR3 * sin(-120 * TMath::Pi() / 180) + DCeyR3 * cos(-120 * TMath::Pi() / 180);
                    }
                    break;
                    case 4:
                    {
                        DCexR1_new = DCexR1 * cos(-180 * TMath::Pi() / 180) - DCeyR1 * sin(-180 * TMath::Pi() / 180);
                        DCeyR1_new = DCexR1 * sin(-180 * TMath::Pi() / 180) + DCeyR1 * cos(-180 * TMath::Pi() / 180);

                        DCexR2_new = DCexR2 * cos(-180 * TMath::Pi() / 180) - DCeyR2 * sin(-180 * TMath::Pi() / 180);
                        DCeyR2_new = DCexR2 * sin(-180 * TMath::Pi() / 180) + DCeyR2 * cos(-180 * TMath::Pi() / 180);

                        DCexR3_new = DCexR3 * cos(-180 * TMath::Pi() / 180) - DCeyR3 * sin(-180 * TMath::Pi() / 180);
                        DCeyR3_new = DCexR3 * sin(-180 * TMath::Pi() / 180) + DCeyR3 * cos(-180 * TMath::Pi() / 180);
                    }
                    break;
                    case 5:
                    {
                        DCexR1_new = DCexR1 * cos(120 * TMath::Pi() / 180) - DCeyR1 * sin(120 * TMath::Pi() / 180);
                        DCeyR1_new = DCexR1 * sin(120 * TMath::Pi() / 180) + DCeyR1 * cos(120 * TMath::Pi() / 180);

                        DCexR2_new = DCexR2 * cos(120 * TMath::Pi() / 180) - DCeyR2 * sin(120 * TMath::Pi() / 180);
                        DCeyR2_new = DCexR2 * sin(120 * TMath::Pi() / 180) + DCeyR2 * cos(120 * TMath::Pi() / 180);

                        DCexR3_new = DCexR3 * cos(120 * TMath::Pi() / 180) - DCeyR3 * sin(120 * TMath::Pi() / 180);
                        DCeyR3_new = DCexR3 * sin(120 * TMath::Pi() / 180) + DCeyR3 * cos(120 * TMath::Pi() / 180);
                    }
                    break;
                    case 6:
                    {
                        DCexR1_new = DCexR1 * cos(60 * TMath::Pi() / 180) - DCeyR1 * sin(60 * TMath::Pi() / 180);
                        DCeyR1_new = DCexR1 * sin(60 * TMath::Pi() / 180) + DCeyR1 * cos(60 * TMath::Pi() / 180);

                        DCexR2_new = DCexR2 * cos(60 * TMath::Pi() / 180) - DCeyR2 * sin(60 * TMath::Pi() / 180);
                        DCeyR2_new = DCexR2 * sin(60 * TMath::Pi() / 180) + DCeyR2 * cos(60 * TMath::Pi() / 180);

                        DCexR3_new = DCexR3 * cos(60 * TMath::Pi() / 180) - DCeyR3 * sin(60 * TMath::Pi() / 180);
                        DCeyR3_new = DCexR3 * sin(60 * TMath::Pi() / 180) + DCeyR3 * cos(60 * TMath::Pi() / 180);
                    }
                    break;
                }
        
                if (DCeyR1_new >= (0.556*DCexR1_new - 6.878) ||  DCeyR1_new <= (-0.56*DCexR1_new + 7.482) || DCexR1_new <= 24.052) continue;
                if (DCeyR2_new >= (0.578*DCexR2_new - 13,898) ||  DCeyR2_new <= (-0.577*DCexR2_new + 14.851) || DCexR2_new <= 39.705) continue;
                if (DCeyR3_new >= (0.591*DCexR3_new - 27.459) ||  DCeyR3_new <= (-0.588*DCexR3_new + 26.912) || DCexR3_new <= 77.755) continue;
                                 
                ///DC cuts for the proton in R1, R2, R3
                switch (DCpsect)
                {
                    case 1:
                    {
                        DCpxR1_new = DCpxR1;
                        DCpyR1_new = DCpyR1;
                    
                        DCpxR2_new = DCpxR2;
                        DCpyR2_new = DCpyR2;
                    
                        DCpxR3_new = DCpxR3;
                        DCpyR3_new = DCpyR3;
                    }
                    break;
                    case 2:
                    {
                        DCpxR1_new = DCpxR1 * cos(-60 * TMath::Pi() / 180) - DCpyR1 * sin(-60 * TMath::Pi() / 180);
                        DCpyR1_new = DCpxR1 * sin(-60 * TMath::Pi() / 180) + DCpyR1 * cos(-60 * TMath::Pi() / 180);

                        DCpxR2_new = DCpxR2 * cos(-60 * TMath::Pi() / 180) - DCpyR2 * sin(-60 * TMath::Pi() / 180);
                        DCpyR2_new = DCpxR2 * sin(-60 * TMath::Pi() / 180) + DCpyR2 * cos(-60 * TMath::Pi() / 180);

                        DCpxR3_new = DCpxR3 * cos(-60 * TMath::Pi() / 180) - DCpyR3 * sin(-60 * TMath::Pi() / 180);
                        DCpyR3_new = DCpxR3 * sin(-60 * TMath::Pi() / 180) + DCpyR3 * cos(-60 * TMath::Pi() / 180);
                    }
                    break;
                    case 3:
                    {
                        DCpxR1_new = DCpxR1 * cos(-120 * TMath::Pi() / 180) - DCpyR1 * sin(-120 * TMath::Pi() / 180);
                        DCpyR1_new = DCpxR1 * sin(-120 * TMath::Pi() / 180) + DCpyR1 * cos(-120 * TMath::Pi() / 180);

                        DCpxR2_new = DCpxR2 * cos(-120 * TMath::Pi() / 180) - DCpyR2 * sin(-120 * TMath::Pi() / 180);
                        DCpyR2_new = DCpxR2 * sin(-120 * TMath::Pi() / 180) + DCpyR2 * cos(-120 * TMath::Pi() / 180);
                                                            
                        DCpxR3_new = DCpxR3 * cos(-120 * TMath::Pi() / 180) - DCpyR3 * sin(-120 * TMath::Pi() / 180);
                        DCpyR3_new = DCpxR3 * sin(-120 * TMath::Pi() / 180) + DCpyR3 * cos(-120 * TMath::Pi() / 180);
                    }
                    break;
                    case 4:
                    {
                        DCpxR1_new = DCpxR1 * cos(-180 * TMath::Pi() / 180) - DCpyR1 * sin(-180 * TMath::Pi() / 180);
                        DCpyR1_new = DCpxR1 * sin(-180 * TMath::Pi() / 180) + DCpyR1 * cos(-180 * TMath::Pi() / 180);

                        DCpxR2_new = DCpxR2 * cos(-180 * TMath::Pi() / 180) - DCpyR2 * sin(-180 * TMath::Pi() / 180);
                        DCpyR2_new = DCpxR2 * sin(-180 * TMath::Pi() / 180) + DCpyR2 * cos(-180 * TMath::Pi() / 180);
                                                            
                        DCpxR3_new = DCpxR3 * cos(-180 * TMath::Pi() / 180) - DCpyR3 * sin(-180 * TMath::Pi() / 180);
                        DCpyR3_new = DCpxR3 * sin(-180 * TMath::Pi() / 180) + DCpyR3 * cos(-180 * TMath::Pi() / 180);
                    }
                    break;
                    case 5:
                    {
                        DCpxR1_new = DCpxR1 * cos(120 * TMath::Pi() / 180) - DCpyR1 * sin(120 * TMath::Pi() / 180);
                        DCpyR1_new = DCpxR1 * sin(120 * TMath::Pi() / 180) + DCpyR1 * cos(120 * TMath::Pi() / 180);

                        DCpxR2_new = DCpxR2 * cos(120 * TMath::Pi() / 180) - DCpyR2 * sin(120 * TMath::Pi() / 180);
                        DCpyR2_new = DCpxR2 * sin(120 * TMath::Pi() / 180) + DCpyR2 * cos(120 * TMath::Pi() / 180);
                                                            
                        DCpxR3_new = DCpxR3 * cos(120 * TMath::Pi() / 180) - DCpyR3 * sin(120 * TMath::Pi() / 180);
                        DCpyR3_new = DCpxR3 * sin(120 * TMath::Pi() / 180) + DCpyR3 * cos(120 * TMath::Pi() / 180);
                    }
                    break;
                    case 6:
                    {
                        DCpxR1_new = DCpxR1 * cos(60 * TMath::Pi() / 180) - DCpyR1 * sin(60 * TMath::Pi() / 180);
                        DCpyR1_new = DCpxR1 * sin(60 * TMath::Pi() / 180) + DCpyR1 * cos(60 * TMath::Pi() / 180);

                        DCpxR2_new = DCpxR2 * cos(60 * TMath::Pi() / 180) - DCpyR2 * sin(60 * TMath::Pi() / 180);
                        DCpyR2_new = DCpxR2 * sin(60 * TMath::Pi() / 180) + DCpyR2 * cos(60 * TMath::Pi() / 180);
                                                            
                        DCpxR3_new = DCpxR3 * cos(60 * TMath::Pi() / 180) - DCpyR3 * sin(60 * TMath::Pi() / 180);
                        DCpyR3_new = DCpxR3 * sin(60 * TMath::Pi() / 180) + DCpyR3 * cos(60 * TMath::Pi() / 180);
                    }
                    break;
                }
                            
                if (DCpsect && (DCpyR1_new > (0.610*DCpxR1_new - 12.720) || DCpyR1_new < (-0.604*DCpxR1_new + 12.159) || (DCpxR1_new < 38.02))) continue;
                if (DCpsect && (DCpyR2_new > (0.573*DCpxR2_new - 13,949) ||  DCpyR2_new < (-0.569*DCpxR2_new + 13.891) || (DCpxR2_new < 54.88))) continue;
                if (DCpsect && (DCpyR3_new > (0.527*DCpxR3_new - 11,998) ||  DCpyR3_new < (-0.530*DCpxR3_new + 13.372) || (DCpxR3_new < 49.))) continue;

                ///DC for the pi- in R1, R2, R3
                switch (DCpimsect)
                {
                    case 1:
                    {
                        DCpimxR1_new = DCpimxR1;
                        DCpimyR1_new = DCpimyR1;
                    
                        DCpimxR2_new = DCpimxR2;
                        DCpimyR2_new = DCpimyR2;
                    
                        DCpimxR3_new = DCpimxR3;
                        DCpimyR3_new = DCpimyR3;
                    }
                    break;
                    case 2:
                    {
                        DCpimxR1_new = DCpimxR1 * cos(-60 * TMath::Pi() / 180) - DCpimyR1 * sin(-60 * TMath::Pi() / 180);
                        DCpimyR1_new = DCpimxR1 * sin(-60 * TMath::Pi() / 180) + DCpimyR1 * cos(-60 * TMath::Pi() / 180);

                        DCpimxR2_new = DCpimxR2 * cos(-60 * TMath::Pi() / 180) - DCpimyR2 * sin(-60 * TMath::Pi() / 180);
                        DCpimyR2_new = DCpimxR2 * sin(-60 * TMath::Pi() / 180) + DCpimyR2 * cos(-60 * TMath::Pi() / 180);

                        DCpimxR3_new = DCpimxR3 * cos(-60 * TMath::Pi() / 180) - DCpimyR3 * sin(-60 * TMath::Pi() / 180);
                        DCpimyR3_new = DCpimxR3 * sin(-60 * TMath::Pi() / 180) + DCpimyR3 * cos(-60 * TMath::Pi() / 180);
                    }
                    break;
                    case 3:
                    {
                        DCpimxR1_new = DCpimxR1 * cos(-120 * TMath::Pi() / 180) - DCpimyR1 * sin(-120 * TMath::Pi() / 180);
                        DCpimyR1_new = DCpimxR1 * sin(-120 * TMath::Pi() / 180) + DCpimyR1 * cos(-120 * TMath::Pi() / 180);

                        DCpimxR2_new = DCpimxR2 * cos(-120 * TMath::Pi() / 180) - DCpimyR2 * sin(-120 * TMath::Pi() / 180);
                        DCpimyR2_new = DCpimxR2 * sin(-120 * TMath::Pi() / 180) + DCpimyR2 * cos(-120 * TMath::Pi() / 180);
                                                            
                        DCpimxR3_new = DCpimxR3 * cos(-120 * TMath::Pi() / 180) - DCpimyR3 * sin(-120 * TMath::Pi() / 180);
                        DCpimyR3_new = DCpimxR3 * sin(-120 * TMath::Pi() / 180) + DCpimyR3 * cos(-120 * TMath::Pi() / 180);
                    }
                    break;
                    case 4:
                    {
                        DCpimxR1_new = DCpimxR1 * cos(-180 * TMath::Pi() / 180) - DCpimyR1 * sin(-180 * TMath::Pi() / 180);
                        DCpimyR1_new = DCpimxR1 * sin(-180 * TMath::Pi() / 180) + DCpimyR1 * cos(-180 * TMath::Pi() / 180);

                        DCpimxR2_new = DCpimxR2 * cos(-180 * TMath::Pi() / 180) - DCpimyR2 * sin(-180 * TMath::Pi() / 180);
                        DCpimyR2_new = DCpimxR2 * sin(-180 * TMath::Pi() / 180) + DCpimyR2 * cos(-180 * TMath::Pi() / 180);
                                                            
                        DCpimxR3_new = DCpimxR3 * cos(-180 * TMath::Pi() / 180) - DCpimyR3 * sin(-180 * TMath::Pi() / 180);
                        DCpimyR3_new = DCpimxR3 * sin(-180 * TMath::Pi() / 180) + DCpimyR3 * cos(-180 * TMath::Pi() / 180);
                    }
                    break;
                    case 5:
                    {
                        DCpimxR1_new = DCpimxR1 * cos(120 * TMath::Pi() / 180) - DCpimyR1 * sin(120 * TMath::Pi() / 180);
                        DCpimyR1_new = DCpimxR1 * sin(120 * TMath::Pi() / 180) + DCpimyR1 * cos(120 * TMath::Pi() / 180);

                        DCpimxR2_new = DCpimxR2 * cos(120 * TMath::Pi() / 180) - DCpimyR2 * sin(120 * TMath::Pi() / 180);
                        DCpimyR2_new = DCpimxR2 * sin(120 * TMath::Pi() / 180) + DCpimyR2 * cos(120 * TMath::Pi() / 180);
                                                            
                        DCpimxR3_new = DCpimxR3 * cos(120 * TMath::Pi() / 180) - DCpimyR3 * sin(120 * TMath::Pi() / 180);
                        DCpimyR3_new = DCpimxR3 * sin(120 * TMath::Pi() / 180) + DCpimyR3 * cos(120 * TMath::Pi() / 180);
                    }
                    break;
                    case 6:
                    {
                        DCpimxR1_new = DCpimxR1 * cos(60 * TMath::Pi() / 180) - DCpimyR1 * sin(60 * TMath::Pi() / 180);
                        DCpimyR1_new = DCpimxR1 * sin(60 * TMath::Pi() / 180) + DCpimyR1 * cos(60 * TMath::Pi() / 180);

                        DCpimxR2_new = DCpimxR2 * cos(60 * TMath::Pi() / 180) - DCpimyR2 * sin(60 * TMath::Pi() / 180);
                        DCpimyR2_new = DCpimxR2 * sin(60 * TMath::Pi() / 180) + DCpimyR2 * cos(60 * TMath::Pi() / 180);
                                                            
                        DCpimxR3_new = DCpimxR3 * cos(60 * TMath::Pi() / 180) - DCpimyR3 * sin(60 * TMath::Pi() / 180);
                        DCpimyR3_new = DCpimxR3 * sin(60 * TMath::Pi() / 180) + DCpimyR3 * cos(60 * TMath::Pi() / 180);
                    }
                    break;
                }
                            
                if (DCpimsect && (DCpimyR1_new >= (0.556*DCpimxR1_new - 6.878) ||  DCpimyR1_new <= (-0.56*DCpimxR1_new + 7.482) || DCpimxR1_new <= 24.052)) continue;
                if (DCpimsect && (DCpimyR2_new >= (0.578*DCpimxR2_new - 13,898) ||  DCpimyR2_new <= (-0.577*DCpimxR2_new + 14.851) || DCpimxR2_new <= 39.705)) continue;
                if (DCpimsect && (DCpimyR3_new >= (0.591*DCpimxR3_new - 27.459) ||  DCpimyR3_new <= (-0.588*DCpimxR3_new + 26.912) || DCpimxR3_new <= 77.755)) continue;

                tree.Fill();
            }
        }
    }

    TFile* hfile =  new TFile("MMPiPlusTree_pass2_CHECK.root","RECREATE");
    tree.Write();

    W_hist->Write();
    Q2_hist->Write();
    Q2_vs_W_hist->Write();
    MM_pi_plus_hist->Write();

    hfile->Close();
}