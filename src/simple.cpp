#ifndef __MDN_SIMPLE__
#define __MDN_SIMPLE__

#include "nlohmann/json.hpp"
#include "adios2.h"

#include "mpi.h"

#include "constants.cpp"
#include "utils.cpp"

#include <fstream>

namespace mdn::simple
{
    namespace fs = std::filesystem;
    int main(const utils::MC &sts, const std::map<fs::path, int> &storages, const int Natoms, const std::array<double, 3> &dims, const double Volume)
    {
        int total_steps = 0;
        for (const auto &[key, val] : storages)
            total_steps += val;

        sts.logger->info("Distributing storages: {} steps over {} workers", total_steps, sts.size);
        std::map<int, std::map<fs::path, std::pair<int, int>>> distribution;
        const int average_steps = ceil(total_steps / sts.size);
        const int remaining_steps = total_steps % sts.size;
        int this_worker_steps = 0;
        int storage_remaining_steps = 0;
        for (const auto &[storage, storage_steps] : storages)
        {
            
        }
        return 0;
    }

}
#endif // !__MDN_SIMPLE__