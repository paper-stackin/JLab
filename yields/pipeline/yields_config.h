#ifndef YIELDS_CONFIG_H
#define YIELDS_CONFIG_H

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

struct YieldsConfig
{
    // Input / output for the missing-mass histograms
    std::string input = "/home/stepan/root_progs/2pion/root_trees/MMPiPlusTree_pass2_CHECK.root";
    std::string tree_name = "MMpiptree";
    std::string output = "mm_pip_hists_m_pip.root";
    std::string hist_name = "MM";
    std::string topology_name = "missing_pip";

    int n_hist_bins = 100;
    double hist_min = -0.3;
    double hist_max = 0.7;

    // Q2 / W binning
    std::vector<double> Q2_edges = {0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.5};
    int n_Q2_bins = static_cast<int>(Q2_edges.size()) - 1;
    int q_min = 0;
    int q_max = n_Q2_bins - 1;

    int n_W_bins = 64;
    double W_min = 1.4;
    double W_max = 3.0;
    double W_bin_size = 0.025;

    std::vector<int> W_bin_limits = {64, 64, 64, 64, 56, 37};

    // Background fit
    int w_fit_min = 8;
    int w_fit_max = 17;

    int w_reper_q45 = 3;
    int w_reper_q45_start = 4;
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

    int signal_integral_bin_min = 25;
    int signal_integral_bin_max = 40;

    // Physics constants
    double E_beam = 6.535;
    double M_proton = 0.93827;
    double M_pion = 0.13957;
    double M_electron = 0.0;
    double alpha = 1.0 / 137.0;

    // Event loop counter
    int speed = 1000000;
    std::string period = "M";

    int WBinLimit(int q) const
    {
        if (q >= 0 && q < static_cast<int>(W_bin_limits.size()))
            return W_bin_limits[q];
        return n_W_bins;
    }
};

inline YieldsConfig LoadYieldsConfig(const std::string& filename)
{
    YAML::Node root = YAML::LoadFile(filename);
    YieldsConfig cfg;

    if (root["input"])
        cfg.input = root["input"].as<std::string>();
    if (root["tree_name"])
        cfg.tree_name = root["tree_name"].as<std::string>();
    if (root["output"])
        cfg.output = root["output"].as<std::string>();
    if (root["hist_name"])
        cfg.hist_name = root["hist_name"].as<std::string>();
    if (root["topology_name"])
        cfg.topology_name = root["topology_name"].as<std::string>();

    if (root["n_hist_bins"])
        cfg.n_hist_bins = root["n_hist_bins"].as<int>();
    if (root["hist_min"])
        cfg.hist_min = root["hist_min"].as<double>();
    if (root["hist_max"])
        cfg.hist_max = root["hist_max"].as<double>();

    if (root["Q2_edges"])
        cfg.Q2_edges = root["Q2_edges"].as<std::vector<double>>();
    cfg.n_Q2_bins = static_cast<int>(cfg.Q2_edges.size()) - 1;
    cfg.q_min = 0;
    cfg.q_max = cfg.n_Q2_bins - 1;

    if (root["n_W_bins"])
        cfg.n_W_bins = root["n_W_bins"].as<int>();
    if (root["W_min"])
        cfg.W_min = root["W_min"].as<double>();
    if (root["W_max"])
        cfg.W_max = root["W_max"].as<double>();
    if (root["W_bin_size"])
        cfg.W_bin_size = root["W_bin_size"].as<double>();
    else if (cfg.n_W_bins > 0)
        cfg.W_bin_size = (cfg.W_max - cfg.W_min) / cfg.n_W_bins;

    if (root["E_beam"])
        cfg.E_beam = root["E_beam"].as<double>();
    if (root["M_proton"])
        cfg.M_proton = root["M_proton"].as<double>();
    if (root["M_pion"])
        cfg.M_pion = root["M_pion"].as<double>();
    if (root["M_electron"])
        cfg.M_electron = root["M_electron"].as<double>();

    if (root["speed"])
        cfg.speed = root["speed"].as<int>();
    if (root["period"])
        cfg.period = root["period"].as<std::string>();

    if (YAML::Node yields = root["yields"])
    {
        if (yields["w_fit_min"])
            cfg.w_fit_min = yields["w_fit_min"].as<int>();
        if (yields["w_fit_max"])
            cfg.w_fit_max = yields["w_fit_max"].as<int>();

        if (yields["w_reper_q45"])
            cfg.w_reper_q45 = yields["w_reper_q45"].as<int>();
        if (yields["w_reper_q45_start"])
            cfg.w_reper_q45_start = yields["w_reper_q45_start"].as<int>();
        if (yields["bg_fit_start"])
            cfg.bg_fit_start = yields["bg_fit_start"].as<int>();

        if (yields["W_bin_limits"])
            cfg.W_bin_limits = yields["W_bin_limits"].as<std::vector<int>>();

        if (yields["xmin"])
            cfg.xmin = yields["xmin"].as<double>();
        if (yields["xmax_low"])
            cfg.xmax_low = yields["xmax_low"].as<double>();
        if (yields["xmax_mid"])
            cfg.xmax_mid = yields["xmax_mid"].as<double>();
        if (yields["xmax_high"])
            cfg.xmax_high = yields["xmax_high"].as<double>();

        if (yields["xmax_mid_start"])
            cfg.xmax_mid_start = yields["xmax_mid_start"].as<int>();
        if (yields["xmax_high_start"])
            cfg.xmax_high_start = yields["xmax_high_start"].as<int>();

        if (yields["bg_mean"])
            cfg.bg_mean = yields["bg_mean"].as<double>();
        if (yields["bg_sigma"])
            cfg.bg_sigma = yields["bg_sigma"].as<double>();

        if (yields["bg_ampl_bin"])
            cfg.bg_ampl_bin = yields["bg_ampl_bin"].as<int>();
        if (yields["bg_ampl_factor"])
            cfg.bg_ampl_factor = yields["bg_ampl_factor"].as<double>();

        if (yields["signal_mean"])
            cfg.signal_mean = yields["signal_mean"].as<double>();
        if (yields["signal_sigma"])
            cfg.signal_sigma = yields["signal_sigma"].as<double>();
        if (yields["signal_tail"])
            cfg.signal_tail = yields["signal_tail"].as<double>();

        if (yields["bg_integral_xmin"])
            cfg.bg_integral_xmin = yields["bg_integral_xmin"].as<double>();
        if (yields["bg_integral_xmax"])
            cfg.bg_integral_xmax = yields["bg_integral_xmax"].as<double>();
        if (yields["bin_width_factor"])
            cfg.bin_width_factor = yields["bin_width_factor"].as<double>();

        if (yields["signal_integral_bin_min"])
            cfg.signal_integral_bin_min = yields["signal_integral_bin_min"].as<int>();
        if (yields["signal_integral_bin_max"])
            cfg.signal_integral_bin_max = yields["signal_integral_bin_max"].as<int>();

        if (yields["alpha"])
            cfg.alpha = yields["alpha"].as<double>();
    }

    return cfg;
}

#endif
