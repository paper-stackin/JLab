#include <chrono>
#include <cstdio>
#include <fcntl.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <unistd.h>
#include <unordered_set>
#include <vector>

class SilenceStdout {
public:
    SilenceStdout()
    {
        fflush(stdout);
        saved_fd = dup(fileno(stdout));
        int devnull = open("/dev/null", O_WRONLY);
        dup2(devnull, fileno(stdout));
        close(devnull);
    }

    ~SilenceStdout()
    {
        fflush(stdout);
        dup2(saved_fd, fileno(stdout));
        close(saved_fd);
    }

private:
    int saved_fd;
};

void generated_events_selection()
{
    std::ios::sync_with_stdio(false);

    // const std::string input_prefix = "/volatile/clas12/osg/savkin/11847/output-11847-";
    // const int input_start = 0;
    // const int input_end = 9999;
    std::string input_prefix = "savkin/10717-";
    int input_start = 0;
    int input_end = 19;
    const std::unordered_set<int> files_to_skip = {6492};

    const std::string output = "generated_events.root";
    const std::string tree_name = "MMpiptree";
    const std::string tree_title = "Generated double pion events";

    const int counter_speed = 1000;

    auto start_time = std::chrono::steady_clock::now();

    float px_e, py_e, pz_e;
    float px_p, py_p, pz_p;
    float px_pim, py_pim, pz_pim;
    float px_pip, py_pip, pz_pip;

    TFile* hfile = new TFile(output.c_str(), "RECREATE");
    TTree tree(tree_name.c_str(), tree_title.c_str());
    tree.SetAutoSave(1e6);

    tree.Branch("px_e", &px_e, "px_e/F");
    tree.Branch("py_e", &py_e, "py_e/F");
    tree.Branch("pz_e", &pz_e, "pz_e/F");

    tree.Branch("px_p", &px_p, "px_p/F");
    tree.Branch("py_p", &py_p, "py_p/F");
    tree.Branch("pz_p", &pz_p, "pz_p/F");

    tree.Branch("px_pim", &px_pim, "px_pim/F");
    tree.Branch("py_pim", &py_pim, "py_pim/F");
    tree.Branch("pz_pim", &pz_pim, "pz_pim/F");

    tree.Branch("px_pip", &px_pip, "px_pip/F");
    tree.Branch("py_pip", &py_pip, "py_pip/F");
    tree.Branch("pz_pip", &pz_pip, "pz_pip/F");

    std::vector<std::string> data;
    data.reserve(input_end - input_start + 1);

    for (int i = input_start; i <= input_end; i++)
    {
        if (files_to_skip.count(i)) continue;
        data.push_back(input_prefix + std::to_string(i) + ".hipo");
    }

    for (int k = 0; k < static_cast<int>(data.size()); k++)
    {
        std::ifstream infile(data[k]);
        if (!infile.good()) continue;

        hipo::reader reader;
        hipo::dictionary factory;
        {
            SilenceStdout quiet;
            reader.open(data[k].c_str());
            reader.readDictionary(factory);
        }

        hipo::event event;
        hipo::bank LUND(factory.getSchema("MC::Lund"));

        int totalEvents = reader.getEntries() / counter_speed;
        int counter = 0;

        while (reader.next())
        {
            counter++;
            if (counter % counter_speed == 0)
            {
                std::cout << "\rFile: " << k + 1 << "/" << data.size()
                          << "        Processed: " << counter / counter_speed
                          << "/" << totalEvents << "k events" << std::flush;
            }

            reader.read(event);
            event.getStructure(LUND);

            const int n_lund = LUND.getRows();
            if (n_lund == 0) continue;

            int el_index = -1;
            int p_index = -1;
            int pip_index = -1;
            int pim_index = -1;

            for (int i = 0; i < n_lund; i++)
            {
                const int lund_index = LUND.getInt("index", i);
                if (lund_index == 1) el_index = i;
                else if (lund_index == 2) p_index = i;
                else if (lund_index == 3) pip_index = i;
                else if (lund_index == 4) pim_index = i;
            }

            if (el_index < 0 || p_index < 0 || pip_index < 0 || pim_index < 0) continue;

            px_e = LUND.getFloat("px", el_index);
            py_e = LUND.getFloat("py", el_index);
            pz_e = LUND.getFloat("pz", el_index);

            px_p = LUND.getFloat("px", p_index);
            py_p = LUND.getFloat("py", p_index);
            pz_p = LUND.getFloat("pz", p_index);

            px_pip = LUND.getFloat("px", pip_index);
            py_pip = LUND.getFloat("py", pip_index);
            pz_pip = LUND.getFloat("pz", pip_index);

            px_pim = LUND.getFloat("px", pim_index);
            py_pim = LUND.getFloat("py", pim_index);
            pz_pim = LUND.getFloat("pz", pim_index);

            tree.Fill();
        }
    }

    tree.Write();
    hfile->Close();

    TFile* file = TFile::Open(output.c_str(), "READ");
    TTree* tree2 = nullptr;
    file->GetObject(tree_name.c_str(), tree2);
    Long64_t nentries = tree2->GetEntries();

    std::cout << std::endl << std::endl;
    std::cout << "Entries: " << nentries << std::endl;
    file->Close();

    auto end_time = std::chrono::steady_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time).count();

    const auto hours = elapsed_seconds / 3600;
    const auto minutes = (elapsed_seconds % 3600) / 60;
    const auto seconds = elapsed_seconds % 60;

    std::cout << "Finished in "
              << std::setw(2) << std::setfill('0') << hours << ":"
              << std::setw(2) << std::setfill('0') << minutes << ":"
              << std::setw(2) << std::setfill('0') << seconds << std::endl;

    gSystem->Exit(0);
}
