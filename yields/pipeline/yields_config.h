#ifndef YIELDS_CONFIG_H
#define YIELDS_CONFIG_H

struct YieldsConfig
{
    int q_min = 0;
    int q_max = 5;
    double Q2_vals[7] = {0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.5};

    int w_fit_min = 8;
    int w_fit_max = 17;

    int w_reper_q45 = 3;
    int bg_fit_start = 13;

    double xmin = -0.2;

    double xmax_low = 0.15;
    double xmax_mid = 0.25;
    double xmax_high = 0.4;

    int xmax_mid_start = 10;
    int xmax_high_start = 17;

    double bg_mean = -1.0;
    double bg_sigma = 0.3;

    int bg_ampl_bin = 25;
    double bg_ampl_factor = 0.3;

    double signal_mean = 0.0196;
    double signal_sigma = 0.03;
    double signal_tail = 0.02;

    double bg_integral_xmin = -0.05;
    double bg_integral_xmax = 0.1;
    double bin_width_factor = 100.0;
};

#endif