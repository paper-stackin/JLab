#ifndef YIELDS_TYPES_H
#define YIELDS_TYPES_H

#include <vector>

struct ExtrapolationResult
{
    double A = 0.0;
    double B = 0.0;
    int w_reper = 0;
};

struct ExtrapolationOutput {
    std::vector<std::vector<double>> integral_background;
    std::vector<std::vector<double>> integral_background_error;
    std::vector<std::vector<double>> w_bg_fit;
    std::vector<ExtrapolationResult> results; 
};

struct FitResult
{
    bool fitted = false;

    double bg_amp = 0;
    double bg_mean = 0;
    double bg_sigma = 0;

    double bg_amp_error = 0;
    double bg_mean_error = 0;
    double bg_sigma_error = 0;

    double signal_amp = 0;
    double signal_mean = 0;
    double signal_sigma = 0;
    double signal_tail = 0;

    double signal_amp_error = 0;
    double signal_mean_error = 0;
    double signal_sigma_error = 0;
    double signal_tail_error = 0;

    double xmin = 0;
    double xmax = 0;
};

struct YieldsOutput {
    std::vector<std::vector<double>> w_vector;
    std::vector<std::vector<double>> w_vector_error;
    std::vector<std::vector<double>> integral_signal_excut;
    std::vector<std::vector<double>> integral_signal_excut_error;
    std::vector<std::vector<double>> integral_signal;
    std::vector<std::vector<double>> integral_error;
};

inline int WBinLimit(int q) 
{
    if (q == 4)
        return 56;
    else if (q == 5)
        return 37;
    else
        return 64;
}

#endif