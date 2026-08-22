#ifndef HADRON_CUTS_H
#define HADRON_CUTS_H

#include <array>

inline bool HadronInFDOrCD(const hipo::bank& PART, int index)
{
    int status = abs(PART.getInt("status", index));
    return (status >= 2000 && status < 8000);
}

inline bool HadronMomentumCut(double momentum, bool in_FD)
{
    if (in_FD)
        return (momentum < 0.4);
    else    // CD case
        return (momentum < 0.2);
}

inline bool HadronVertexCut(const hipo::bank& PART, int index)
{
    double vz = PART.getFloat("vz", index);
    return (vz < -10 || vz > 2);
}

inline bool HadronChi2Cut(const hipo::bank& PART, int index)
{
    return (PART.getFloat("chi2pid", index) >= 5);
}

bool HadronDCFiducialCut(const std::array<double, 3>& DC_x_new,
                         const std::array<double, 3>& DC_y_new,
                         int DC_sect, int charge)
{
    const double p_negative[3][5] = {
        {0.556, -6.878, -0.56, 7.482, 24.052},
        {0.578, -13.898, -0.577, 14.851, 39.705},
        {0.591, -27.459, -0.588, 26.912, 77.755}
    };

    const double p_positive[3][5] = {
        {0.610, -12.720, -0.604, 12.159, 38.02},
        {0.573, -13.949, -0.569, 13.891, 54.88},
        {0.527, -11.998, -0.530, 13.372, 49.}
    };

    const double (*p)[5] = (charge < 0) ? p_negative : p_positive;

    for (int r = 0; r < 3; ++r)
    {
        if (DC_sect &&
            (DC_y_new[r] >= (p[r][0] * DC_x_new[r] + p[r][1]) ||
            DC_y_new[r] <= (p[r][2] * DC_x_new[r] + p[r][3]) ||
            DC_x_new[r] <= p[r][4]))
        {
            return true;
        }
    }

    return false;
}

inline bool HadronBetaVsPCut(double beta, bool in_FD)
{
    if (in_FD)
        return (beta < 0.4 || beta > 1.1);
    else    // CD case
        return (beta < 0.2 || beta > 1.1);
}

inline bool TopologyCut(std::string topology, int el_cnt, int p_cnt, int pip_cnt, int pim_cnt, int unid_cnt, int charged_cnt) 
{ 
    if (unid_cnt > 0) return true; 
    if (el_cnt != 1) return true; 
    
    if (topology == "full") 
        return !(p_cnt == 1 && pip_cnt == 1 && pim_cnt == 1 && charged_cnt == 4); 
    else if (topology == "missing_pim") 
        return !(p_cnt == 1 && pip_cnt == 1 && pim_cnt == 0 && charged_cnt == 3); 
    else if (topology == "missing_pip") 
        return !(p_cnt == 1 && pip_cnt == 0 && pim_cnt == 1 && charged_cnt == 3); 
    else if (topology == "missing_p") 
        return !(p_cnt == 0 && pip_cnt == 1 && pim_cnt == 1 && charged_cnt == 3); 
    else 
        return true; 
}

#endif
