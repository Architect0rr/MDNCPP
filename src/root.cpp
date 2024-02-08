#ifndef __MDN_ROOT__
#define __MDN_ROOT__

#include "nlohmann/json.hpp"
#include "adios2.h"

#include "mpi.h"

#include "constants.cpp"
#include "utils.cpp"
#include "simple.cpp"
#include "barriers.cpp"
#include "logic.cpp"

#include <fstream>

namespace mdn::root
{
    using json = nlohmann::json;
    namespace fs = std::filesystem;

    int get_total_steps(adios2::IO &io, const fs::path &storage)
    {
        adios2::Engine reader = io.Open(storage, adios2::Mode::Read);
        int steps = reader.Steps();
        reader.Close();
        return steps;
    }

    const std::map<fs::path, int> storage_rsolve(utils::MC &sts, const json &_storages)
    {
        adios2::ADIOS adios = adios2::ADIOS();
        adios2::IO io = adios.DeclareIO("Steps_reader");
        sts.logger.debug("ADIOS2.IO initialized");
        std::map<fs::path, int> storages;
        int steps = 0;
        try{
            for (const auto &storage : _storages.items())
            {
                sts.logger.trace("Storage '{}', resolved path: {}", storage.key(), fs::absolute(storage.key()).string());
                steps = get_total_steps(io, storage.key());
                storages.emplace(fs::absolute(storage.key()), steps);
                sts.logger.trace("Storage '{}' contains {} steps", storage.key(), steps);
            }
        }catch(std::exception& e){
            sts.logger.error("Some std::exception was thrown while resolving storages, probably aborting");
            sts.logger.error(e.what());
        }catch (...){
            sts.logger.error("Some unknown exception was thrown while resolving storages, probably aborting");
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
        reader.Get(_Natoms, Natoms);
        reader.Get(_boxxhi, dims->data());
        reader.Get(_boxyhi, dims->data() + 1);
        reader.Get(_boxzhi, dims->data() + 2);
        reader.EndStep();
        reader.Close();
        return Natoms;
    }

    RETURN_CODES distribute(utils::MC &sts, std::map<fs::path, std::pair<int, int>> *mydistrib, int *NNatoms)
    {
        sts.logger.debug("Reading datafile");
        fs::path datafile_path = sts.cwd / files::data;
        if (!fs::exists(datafile_path))
        {
            sts.logger.error("File: '" + datafile_path.string() + "' does not exists");
            return RETURN_CODES::ERROR;
        }
        std::ifstream datafile_in(datafile_path);
        json data = json::parse(datafile_in);
        sts.logger.debug("Successfully read datafile");

        sts.logger.debug("Resolving storages");
        const std::map<fs::path, int> storages = storage_rsolve(sts, data.at(fields::storages));
        // {
        //     sts.logger.debug(std::string("Storage: ") + key.string() + std::string(" has ") + std::to_string(val) + std::string(" steps."));
        // }

        fs::path sto_check = storages.begin()->first;
        std::array<double, 3> dims;
        int Natoms = bearbeit(sts, sto_check, &dims);
        sts.logger.debug("Number of atoms: {}", Natoms);
        sts.logger.debug("Dimensions: [{}, {}, {}]", dims[0], dims[1], dims[2]);
        data[fields::Natoms] = Natoms;
        data[fields::Dims] = dims;
        const double Volume = dims[0] * dims[1] * dims[2];
        data[fields::Volume] = Volume;
        sts.logger.debug("Updating datafile");
        std::ofstream datafile_out(datafile_path);
        datafile_out << std::setw(4) << data << std::endl;

        fs::path data_processing_folder = sts.cwd / data.at(fields::data_processing_folder);
        sts.logger.debug("Creating post-processing directory: " + data_processing_folder.string());
        if (!fs::exists(data_processing_folder))
        {
            std::error_code ec;
            if (!(fs::create_directory(data_processing_folder, ec) && ec))
            {
                sts.logger.error("Cannot create directory: " + data_processing_folder.string());
                return RETURN_CODES::ERROR;
            }
        }

        sts.logger.debug("Distributing storages");

        std::map<int, std::map<fs::path, std::pair<int, int>>> distribution = simple::dns(sts, storages);

        *mydistrib = distribution.at(0);
        *NNatoms = Natoms;

        return RETURN_CODES::OK;
    }

    RETURN_CODES sanity(const utils::MC &sts)
    {
        int rnd = std::rand();
        std::cout << "      Bcasting: " << rnd << std::endl;
        MPI_Bcast(&rnd, 1, MPI_INT, cs::mpi_root, MPI_COMM_WORLD);
        // int responce = 0;
        int *responces = new int[sts.size];
        std::cout << "      Gathering..." << std::endl;
        MPI_Gather(&rnd, 1, MPI_INT, responces, 1, MPI_INT, cs::mpi_root, MPI_COMM_WORLD);

        if (responces[0] != rnd)
        {
            std::cerr << "Root sanity doesn't passed" << std::endl;
            return RETURN_CODES::ERROR;
        }

        // sts.logger.info(std::string("Random int: " + std::to_string(rnd)));
        // sts.logger.info(std::string("Received arrray:"));
        for (size_t i = 0; i < sts.size; ++i)
        {
            // sts.logger.info(std::to_string(responces[i]));
            if (responces[i] - i != rnd)
            {
                std::cerr << "Root sanity doesn't passed" << std::endl;
                return RETURN_CODES::ERROR;
            }
        }
        std::cout << "      Root sanity OK" << std::endl;
        std::cout << "      Releasing sanity barrier" << std::endl;
        SANITY_BARRIER(MPI_COMM_WORLD);
        std::cout << "      Sanity barrier released" << std::endl;

        return RETURN_CODES::OK;
    }

    RETURN_CODES setup(utils::MC &sts)
    {
        std::cout << "  Sanity:" << std::endl;
        if (sanity(sts) != RETURN_CODES::OK)
            return RETURN_CODES::ERROR;

        std::error_code ec;
        fs::path log_folder = sts.cwd / "logs";
        if (!fs::exists(log_folder))
        {
            if (!fs::create_directory(log_folder, ec))
            {
                std::cerr << "Cannot create non-existent directory: " << log_folder.string() << std::endl;
                std::cerr << "Error code: " << ec.value() << std::endl;
                std::cerr << "Message: " << ec.message() << std::endl;
                return RETURN_CODES::ERROR;
            }
            else
            {
                std::cout << "  Created logging directory: " << log_folder << std::endl;
            }
        }
        std::cout << "  Releasing setup barrier" << std::endl;
        SETUP_BARRIER(MPI_COMM_WORLD);
        std::cout << "  Setup barrier released" << std::endl;

        sts.setup_logger(log_folder, sts.rank);

        std::cout << "  Root setup OK" << std::endl;
        return RETURN_CODES::OK;
    }

    RETURN_CODES entry(utils::MC &sts)
    {
        std::cout << "Setup..." << std::endl;
        if (setup(sts) != RETURN_CODES::OK)
            return RETURN_CODES::ERROR;

        sts.logger.info("Initialized");

        std::map<fs::path, std::pair<int, int>> *distribution;
        int Natoms{};

        if (distribute(sts, distribution, &Natoms) != RETURN_CODES::OK)
            return RETURN_CODES::ERROR;

        MPI_Bcast(&Natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);

        fs::path ss = sts.args.outfile / (sts.args.outfile.filename().string() + "." + std::to_string(sts.rank));

        logic::run(sts, ss, *distribution, Natoms);

        return RETURN_CODES::OK;
    }
}

#endif // !__MDN_ROOT__