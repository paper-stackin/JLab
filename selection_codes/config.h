#ifndef CONFIG_H
#define CONFIG_H

#include <string>

struct FileConfig {
    std::string input_prefix = "savkin/10717-";
    int input_start = 0;
    int input_end = 19;

    // std::string input_prefix = "/cache/clas12/rg-k/production/recon/fall2018/torus+1/6535MeV/pass2/v0/dst/train/skim30/skim30_00";
    // int input_start = 5860;
    // int input_end = 6000;
    
    std::string output = "MPPT_events_simulation.root";
    std::string tree_name = "MMpiptree";
    std::string tree_title = "Double pion chanel";
};

struct CounterConfig {
    int speed = 1e3;
    std::string period = "k";
};

struct SelectionConfig{
    std::string topology = "missing_pip";
    double E_beam = 6.535;
    bool use_weights = true;

    FileConfig file;
    
    CounterConfig counter;
};

#endif
