#ifndef __MDN_ROOT__
#define __MDN_ROOT__

#include "nlohmann/json.hpp"
#include "adios2.h"

#include "mpi.h"

#include "constants.cpp"
#include "utils.cpp"
#include "simple.cpp"

#include <fstream>
// #include <vector>
// #include <utility>

namespace mdn::root
{
    using json = nlohmann::json;
    namespace fs = std::filesystem;

    int get_total_steps(adios2::IO &io, const fs::path &storage)
    {
        adios2::Engine reader = io.Open(storage, adios2::Mode::Read);
        // reader.BeginStep();
        // adios2::Variable<std::string> varGreeting =
        //     io.InquireVariable<std::string>("Greeting");
        // std::string greeting;
        // reader.Get(varGreeting, greeting);
        // reader.EndStep();
        int steps = reader.Steps();
        reader.Close();
        // return greeting;
        return steps;
    }

    const std::map<fs::path, int> storage_rsolve(const utils::MC &sts, const json &_storages)
    {
        adios2::ADIOS adios = adios2::ADIOS();
        adios2::IO io = adios.DeclareIO("Steps_reader");
        sts.logger->debug("ADIOS2.IO initialized");
        std::map<fs::path, int> storages;
        int steps = 0;

        for (const auto &storage : _storages.items())
        {
            steps = get_total_steps(io, storage.key());
            storages.emplace(fs::absolute(storage.key()), steps);
            sts.logger->trace("Storage '{}' contains {} steps", storage.key(), steps);
        }

        return storages;
    }

    const int bearbeit(const utils::MC &sts, const fs::path &sto, std::array<double, 3> *dims)
    {
        adios2::ADIOS adios = adios2::ADIOS();
        adios2::IO io = adios.DeclareIO("PReader");
        adios2::Engine reader = io.Open(sto, adios2::Mode::Read);
        reader.BeginStep();
        adios2::Variable<int> _Natoms = io.InquireVariable<int>(std::string(lcf::natoms));
        adios2::Variable<double> _boxxhi = io.InquireVariable<double>(std::string(lcf::boxxhi));
        adios2::Variable<double> _boxyhi = io.InquireVariable<double>(std::string(lcf::boxyhi));
        adios2::Variable<double> _boxzhi = io.InquireVariable<double>(std::string(lcf::boxzhi));
        int Natoms;
        // double boxxhi;
        // double boxyhi;
        // double boxzhi;
        reader.Get(_Natoms, Natoms);
        reader.Get(_boxxhi, dims->data());
        reader.Get(_boxyhi, dims->data() + 1);
        reader.Get(_boxzhi, dims->data() + 2);
        reader.EndStep();
        reader.Close();
        return Natoms;
    }

    int main(const utils::MC &sts)
    {
        sts.logger->debug("Reading datafile");
        fs::path datafile_path = sts.cwd / files::data;
        if (!fs::exists(datafile_path))
        {
            sts.logger->error("File: '" + datafile_path.string() + "' does not exists");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        std::ifstream datafile_in(datafile_path);
        json data = json::parse(datafile_in);
        sts.logger->debug("Successfully read datafile");

        sts.logger->debug("Resolving storages");
        const std::map<fs::path, int> storages = storage_rsolve(sts, data.at(fields::storages));
        // for (auto const &[key, val] : storages)
        // {
        //     sts.logger->debug(std::string("Storage: ") + key.string() + std::string(" has ") + std::to_string(val) + std::string(" steps."));
        // }

        fs::path sto_check = storages.begin()->first;
        std::array<double, 3> dims;
        const int Natoms = bearbeit(sts, sto_check, &dims);
        sts.logger->debug("Number of atoms: {}", Natoms);
        sts.logger->debug("Dimensions: [{}, {}, {}]", dims[0], dims[1], dims[2]);
        data[fields::Natoms] = Natoms;
        data[fields::Dims] = dims;
        const double Volume = dims[0] * dims[1] * dims[2];
        data[fields::Volume] = Volume;
        sts.logger->debug("Updating datafile");
        std::ofstream datafile_out(datafile_path);
        datafile_out << std::setw(4) << data << std::endl;

        fs::path data_processing_folder = sts.cwd / data.at(fields::data_processing_folder);
        sts.logger->debug("Creating post-processing directory: " + data_processing_folder.string());
        if (!fs::exists(data_processing_folder))
        {
            std::error_code ec;
            if (!(fs::create_directory(data_processing_folder, ec) && ec))
            {
                sts.logger->error("Cannot create directory: " + data_processing_folder.string());
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }

        int exit_code = 1;

        if (sts.args.mode == 0)
        {
            exit_code = simple::main(sts, storages, Natoms, dims, Volume);
        }

        return exit_code;
    }
}

#endif // !__MDN_ROOT__