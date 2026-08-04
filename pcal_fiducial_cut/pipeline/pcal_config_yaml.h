#ifndef PCAL_CONFIG_YAML_H
#define PCAL_CONFIG_YAML_H

#include <yaml-cpp/yaml.h>
#include "pcal_config.h"

namespace YAML {
    template<>
    struct convert<FileConfig>
    {
        static bool decode(const Node& node, FileConfig& cfg)
        {
            cfg.main_folder = node["main_folder"].as<std::string>();
            cfg.input_file = cfg.main_folder + node["input_file"].as<std::string>();
            cfg.intermediate_file = cfg.main_folder + node["intermediate_file"].as<std::string>();
            cfg.pdf_name_prefix = cfg.main_folder + node["pdf_name_prefix"].as<std::string>();

            return true;
        }
    };

    template<>
    struct convert<HistogramConfig>
    {
        static bool decode(const Node& node, HistogramConfig& cfg)
        {
            cfg.x_start   = node["x_start"].as<int>();
            cfg.n_x_bins  = node["n_x_bins"].as<int>();
            cfg.bin_size  = node["bin_size"].as<int>();

            cfg.n_y_bins  = node["n_y_bins"].as<int>();
            cfg.y_min     = node["y_min"].as<int>();
            cfg.y_max     = node["y_max"].as<int>();

            return true;
        }
    };

    template<>
    struct convert<BoundaryConfig>
    {
        static bool decode(const Node& node, BoundaryConfig& cfg)
        {
            cfg.use_supergauss  = node["use_supergauss"].as<bool>();
            cfg.supergauss_n    = node["supergauss_n"].as<int>();
            cfg.supergauss_level= node["supergauss_level"].as<double>();

            cfg.quantile_low    = node["quantile_low"].as<double>();
            cfg.quantile_high   = node["quantile_high"].as<double>();

            cfg.x_fit_min       = node["x_fit_min"].as<int>();
            cfg.x_fit_max       = node["x_fit_max"].as<int>();

            return true;
        }
    };

    template<>
    struct convert<RadialCutConfig>
    {
        static bool decode(const Node& node, RadialCutConfig& cfg)
        {
            cfg.r_min = node["r_min"].as<int>();
            cfg.r_max = node["r_max"].as<int>();

            cfg.phi_1 = node["phi_1"].as<double>();
            cfg.phi_2 = node["phi_2"].as<double>();

            return true;
        }
    };

    template<>
    struct convert<CounterConfig>
    {
        static bool decode(const Node& node, CounterConfig& cfg)
        {
            cfg.speed  = node["speed"].as<int>();
            cfg.period = node["period"].as<std::string>();

            return true;
        }
    };

    template<>
    struct convert<StyleConfig>
    {
        static bool decode(const Node& node, StyleConfig& cfg)
        {
            cfg.use_batch = node["use_batch"].as<bool>();

            cfg.palette = node["palette"].as<int>();

            cfg.stat_option = node["stat_option"].as<int>();
            cfg.fit_option  = node["fit_option"].as<int>();

            cfg.canvas_width  = node["canvas_width"].as<int>();
            cfg.canvas_height = node["canvas_height"].as<int>();

            return true;
        }
    };

    template<>
    struct convert<DipInterpolationConfig::DipSectorConfig>
    {
        static bool decode(const Node& node,
                        DipInterpolationConfig::DipSectorConfig& cfg)
        {
            cfg.sector    = node["sector"].as<int>();
            cfg.start_bin = node["start_bin"].as<int>();
            cfg.end_bin   = node["end_bin"].as<int>();

            return true;
        }
    };

    template<>
    struct convert<DipInterpolationConfig>
    {
        static bool decode(const Node& node,
                        DipInterpolationConfig& cfg)
        {
            cfg.enable = node["enable"].as<bool>();

            cfg.sectors =
                node["sectors"].as<std::vector<DipInterpolationConfig::DipSectorConfig>>();

            return true;
        }
    };

    template<>
    struct convert<ExclusionStripSettings::ExclusionStripConfig>
    {
        static bool decode(const Node& node,
                        ExclusionStripSettings::ExclusionStripConfig& cfg)
        {
            cfg.sector = node["sector"].as<int>();

            cfg.slope = node["slope"].as<double>();

            cfg.c_min = node["c_min"].as<double>();
            cfg.c_max = node["c_max"].as<double>();

            return true;
        }
    };

    template<>
    struct convert<ExclusionStripSettings>
    {
        static bool decode(const Node& node,
                        ExclusionStripSettings& cfg)
        {
            cfg.enable = node["enable"].as<bool>();

            cfg.margin = node["margin"].as<double>();

            cfg.strips =
                node["strips"].as<std::vector<ExclusionStripSettings::ExclusionStripConfig>>();

            return true;
        }
    };
}

PCALConfig LoadConfig(const std::string& filename)
{
    YAML::Node root = YAML::LoadFile(filename);

    PCALConfig cfg;  // сначала дефолтные значения из pcal_config.h

    if (root["use_weights"])
        cfg.use_weights = root["use_weights"].as<bool>();

    if (root["file"])
        cfg.file = root["file"].as<FileConfig>();

    if (root["histogram"])
        cfg.histogram = root["histogram"].as<HistogramConfig>();

    if (root["exclusion"])
        cfg.exclusion = root["exclusion"].as<ExclusionStripSettings>();

    if (root["dip"])
        cfg.dip = root["dip"].as<DipInterpolationConfig>();

    if (root["radial"])
        cfg.radial = root["radial"].as<RadialCutConfig>();

    if (root["boundary"])
        cfg.boundary = root["boundary"].as<BoundaryConfig>();

    if (root["style"])
        cfg.style = root["style"].as<StyleConfig>();

    if (root["counter"])
        cfg.counter = root["counter"].as<CounterConfig>();

    return cfg;
}

#endif