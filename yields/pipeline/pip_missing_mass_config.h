#pragma once

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

struct PipMissingMassConfig {
    std::string input = "/home/stepan/root_progs/2pion/root_trees/MMPiPlusTree_pass2_CHECK.root";
    std::string tree_name = "MMpiptree";
    std::string output = "mm_pip_hists_m_pip.root";
    std::string hist_name = "MM";

    int n_hist_bins = 100;
    double hist_min = -0.3;
    double hist_max = 0.7;

    std::vector<double> Q2_edges = {0.5, 0.7, 1.0, 1.4, 2.0, 3.0, 5.5};
    int n_Q2_bins = static_cast<int>(Q2_edges.size() - 1);

    int n_W_bins = 64;
    double W_min = 1.4;
    double W_max = 3.0;
    double W_bin_size = 0.025;

    double E_beam = 6.535;
    double M_proton = 0.93827;
    double M_pion = 0.13957;
    double M_electron = 0.0;

    int speed = 1000000;
    std::string period = "M";
};

inline PipMissingMassConfig LoadPipMissingMassConfig(const std::string& filename)
{
    YAML::Node root = YAML::LoadFile(filename);
    PipMissingMassConfig cfg;

    if (root["input"])
        cfg.input = root["input"].as<std::string>();
    if (root["tree_name"])
        cfg.tree_name = root["tree_name"].as<std::string>();
    if (root["output"])
        cfg.output = root["output"].as<std::string>();
    if (root["hist_name"])
        cfg.hist_name = root["hist_name"].as<std::string>();

    if (root["n_hist_bins"])
        cfg.n_hist_bins = root["n_hist_bins"].as<int>();
    if (root["hist_min"])
        cfg.hist_min = root["hist_min"].as<double>();
    if (root["hist_max"])
        cfg.hist_max = root["hist_max"].as<double>();

    if (root["Q2_edges"])
        cfg.Q2_edges = root["Q2_edges"].as<std::vector<double>>();
    cfg.n_Q2_bins = static_cast<int>(cfg.Q2_edges.size() - 1);

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

    return cfg;
}
