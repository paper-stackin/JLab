#ifndef PCAL_BOUNDARIES_H
#define PCAL_BOUNDARIES_H

#include "pcal_config.h"
#include <cmath>
#include <utility>

inline std::pair<double, double> ComputeQuantileBounds(TH1F* hist, double q_low, double q_high)
{
    double probs[2] = {q_low, q_high};
    double quant[2];
    hist->GetQuantiles(2, quant, probs);
    return {quant[0], quant[1]};
}

inline double SuperGauss(double* x, double* par)
{
    const double A = par[0];
    const double mu = par[1];
    const double sigma = par[2];
    const double n = par[3];

    const double arg = (x[0] - mu) / sigma;
    return A * std::exp(-std::pow(std::fabs(arg), n));
}

inline int EnsureEvenSuperGaussN(int n)
{
    if (n < 2) return 2;
    return (n % 2 == 0) ? n : n - 1;
}

inline std::pair<double, double> ComputeSuperGaussBounds(TH1F* hist, double level, int n, const char* fit_name)
{
    const int n_even = EnsureEvenSuperGaussN(n);

    TF1* fit = new TF1(fit_name, SuperGauss,
            hist->GetXaxis()->GetXmin(),
            hist->GetXaxis()->GetXmax(), 4);

    fit->SetParameters(
        hist->GetMaximum(),
        hist->GetBinCenter(hist->GetMaximumBin()),
        hist->GetRMS(),
        static_cast<double>(n_even)
    );
    fit->FixParameter(3, static_cast<double>(n_even));

    hist->Fit(fit, "RQ0");

    const double mu = fit->GetParameter(1);
    const double sigma = fit->GetParameter(2);

    const double dx = sigma * std::pow(-std::log(level), 1.0 / n_even);
    return {mu - dx, mu + dx};
}

inline std::pair<double, double> ComputeVerticalBounds(TH1F* hist, const PCALConfig& cfg, const char* fit_name = "sg_fit")
{
    if (cfg.boundary.use_supergauss) {
        auto bounds = ComputeSuperGaussBounds(hist, cfg.boundary.supergauss_level, cfg.boundary.supergauss_n, fit_name);
        TF1* fit = hist->GetFunction(fit_name);
        if (fit) {
            fit->SetLineColor(kRed);
            fit->SetLineWidth(2);
            fit->Draw("same");
        }
        return bounds;
    }
    return ComputeQuantileBounds(hist, cfg.boundary.quantile_low, cfg.boundary.quantile_high);
}

#endif
