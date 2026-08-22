#ifndef ELECTRON_CUTS_H
#define ELECTRON_CUTS_H

#include <array>

#include "config.h"

extern SelectionConfig cfg;

inline bool ElectronInFD(const hipo::bank& PART, int index = 0)
{
    int el_status = abs(PART.getInt("status", index));
    return (2000 <= el_status) && (el_status < 4000);
}

inline bool ElectronIsFirst(int index)
{
    return (index == 0);
}

inline bool ElectronMomentumCut(double el_momentum)
{
    return (el_momentum < 1.0 || el_momentum > cfg.E_beam);
}

inline bool ElectronVertexCut(const hipo::bank& PART, int index = 0)
{
    double vz_el = PART.getFloat("vz", index);
    return (vz_el < -10 || vz_el > 2);
}

bool SamplingFractionCut(double E_tot_e, double el_mom, int sector)
{
    const double p[6][3] = {
                    {0.145, 0.0216, 0.00181},
                    {0.134, 0.0230, 0.00168},
                    {0.145, 0.0211, 0.00166},
                    {0.152, 0.0161, 0.00127},
                    {0.141, 0.0210, 0.00174},
                    {0.141, 0.02152, 0.00170}
                };

    if (sector < 1 || sector > 6)
        return false;

    const double& p0 = p[sector - 1][0];
    const double& p1 = p[sector - 1][1];
    const double& p2 = p[sector - 1][2];

    bool SFcut = (E_tot_e / el_mom) < (p0 + p1 * el_mom - p2 * el_mom * el_mom);
    return SFcut || !sector;
}

bool ElectronDCFiducialCut(const std::array<double, 3>& DC_el_x_new,
                           const std::array<double, 3>& DC_el_y_new)
{
    const double p[3][5] = {
        {0.556, -6.878, -0.56, 7.482, 24.052},
        {0.578, -13.898, -0.577, 14.851, 39.705},
        {0.591, -27.459, -0.588, 26.912, 77.755}
    };

    for (int r = 0; r < 3; ++r)
    {
        if (DC_el_y_new[r] >= (p[r][0] * DC_el_x_new[r] + p[r][1]) ||
            DC_el_y_new[r] <= (p[r][2] * DC_el_x_new[r] + p[r][3]) ||
            DC_el_x_new[r] <= p[r][4])
        {
            return true;
        }
    }

    return false;
}

inline bool ElectronTOFCut(const hipo::bank& PART, const hipo::bank& Scint, int index = 0)
{
    double hit_el = 0.0; // Time of hit
    for (int i = 0; i < Scint.getRows(); i++)
    {
        float scint_time = Scint.getFloat("time", i);
        if (Scint.getInt("pindex", i) == index && scint_time > hit_el)
            hit_el = scint_time;
    }
    double vt_el = PART.getFloat("vt", index);  // Vertex time
    double tof_e = hit_el - vt_el;
    return (tof_e < 21 || tof_e > 26);
}

inline bool PiMinusContaminationCut(double EC_in, double E_PCAL, double el_mom)
{
    return (EC_in / el_mom) < (-0.84 * E_PCAL / el_mom + 0.17);
}

inline bool ElectronChi2Cut(const hipo::bank& PART, int index = 0)
{
    return (PART.getFloat("chi2pid", index) >= 5);
}

inline bool ElectronUVWCut(double lu, double lv, double lw)
{
    return (lu <= 40 || lu >= 400 || lv <= 15 || lw <= 15);
}

inline bool DoublePionWThresholdCut(const TLorentzVector& beam, const TLorentzVector& p_in, const TLorentzVector& el)
{
    double W = (beam + p_in - el).M();
    return (W < 1.2);
}

#endif
