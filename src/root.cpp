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
#include "MDN.hpp"

#include <fstream>
#include <memory>

namespace mdn {

    using json = nlohmann::json;
    namespace fs = std::filesystem;

    const int MDN_root::get_total_steps(adios2::IO &io, const fs::path &storage)
    {
        adios2::Engine reader = io.Open(storage, adios2::Mode::Read);
        int steps = reader.Steps();
        reader.Close();
        return steps;
    }

    const std::map<fs::path, int> MDN_root::storage_rsolve(const json &_storages)
    {
        adios2::ADIOS adios = adios2::ADIOS(MPI_COMM_SELF);
        adios2::IO io = adios.DeclareIO("Steps_reader");
        logger.debug("ADIOS2.IO initialized");
        std::map<fs::path, int> storages;
        int steps = 0;
        fs::path rsolved_storage;
        for (const auto &storage : _storages.items())
        {
            TRY
            rsolved_storage = fs::absolute(storage.key());
            logger.trace("Storage '{}', resolved path: {}", storage.key(), rsolved_storage.string());
            steps = get_total_steps(io, rsolved_storage);
            storages.emplace(rsolved_storage, steps);
            logger.trace("Storage '{}' contains {} steps", rsolved_storage.string(), steps);
            CATCH("Error while resolving storage: " + storage.key())
        }

        return storages;
    }

    const uint64_t MDN_root::bearbeit(const fs::path &storage, double *dims)
    {
        adios2::ADIOS adios = adios2::ADIOS();
        adios2::IO io = adios.DeclareIO("PReader");
        adios2::Engine reader;
        TRY
            reader = io.Open(storage, adios2::Mode::Read);
        CATCH("Error while opening storage" + storage.string())

        reader.BeginStep();

        adios2::Variable<uint64_t> varNatoms;
        adios2::Variable<double>   varboxxhi;
        adios2::Variable<double>   varboxyhi;
        adios2::Variable<double>   varboxzhi;

        TRY
            varNatoms = io.InquireVariable<uint64_t>(std::string(lcf::natoms));
            varboxxhi = io.InquireVariable<double>(std::string(lcf::boxxhi));
            varboxyhi = io.InquireVariable<double>(std::string(lcf::boxyhi));
            varboxzhi = io.InquireVariable<double>(std::string(lcf::boxzhi));
        CATCH("Error while inquiring variables")

        uint64_t Natoms = 0;

        TRY
            reader.Get(varNatoms, Natoms);
            reader.Get(varboxxhi, dims + 0);
            reader.Get(varboxyhi, dims + 1);
            reader.Get(varboxzhi, dims + 2);
        CATCH("Error while getting variables")

        reader.EndStep();
        reader.Close();
        return Natoms;
    }

    uint64_t MDN_root::dns(std::map<fs::path, std::pair<int, int>> &_distrib)
    {
        logger.debug("Reading datafile");
        fs::path datafile_path = cwd / files::data;
        if (!fs::exists(datafile_path))
        {
            logger.error("Datafile (file: '" + datafile_path.string() + "') does not exists");
            throw std::runtime_error("Datafile (file: '" + datafile_path.string() + "') does not exists");
        }
        std::ifstream datafile_in(datafile_path);
        json data;
        TRY
            data = json::parse(datafile_in);
        CATCH("Error while parsing json from datafile")
        logger.debug("Successfully read datafile");

        logger.debug("Resolving storages");
        std::map<fs::path, int> storages;
        TRY
            storages = storage_rsolve(data.at(fields::storages));
        CATCH("Error while resolving storages")
        // {
        //     logger.debug(std::string("Storage: ") + key.string() + std::string(" has ") + std::to_string(val) + std::string(" steps."));
        // }

        fs::path sto_check = storages.begin()->first;
        logger.debug("Simulation parameters will be obtained from storage: {}", sto_check.string());

        std::unique_ptr<double[]> dims(new double[3]);
        uint64_t Natoms{};
        TRY
            Natoms = bearbeit(sto_check, dims.get());
        CATCH("Error while getting simulation parameters from one of storages")

        const std::array<double, 3> _dims = {*dims.get(), *(dims.get() + 1), *(dims.get() + 2)};
        const double Volume = dims[0] * dims[1] * dims[2];
        logger.debug("Number of atoms: {}", Natoms);
        logger.debug("Dimensions: [{}, {}, {}]", dims[0], dims[1], dims[2]);
        logger.debug("System volume: {}", Volume);

        TRY
            data[fields::Natoms] = Natoms;
            data[fields::Dims] = _dims;
            data[fields::Volume] = Volume;
            logger.debug("Updating datafile");
            std::ofstream datafile_out(datafile_path);
            datafile_out << std::setw(4) << data << std::endl;
        CATCH("Error while writing datafile")

        fs::path data_processing_folder = cwd / folders::post_process_folder;
        logger.debug("Creating post-processing directory: " + data_processing_folder.string());

        std::error_code ec;
        if (!fs::exists(data_processing_folder, ec)){
            if (!ec){
                if (!(fs::create_directory(data_processing_folder, ec) && ec)){
                    logger.error("Can not create directory '{}' due to error with code: {}. Message", data_processing_folder.string(), ec.value());
                    logger.error(ec.message());
                    throw std::runtime_error("Can not create post-processing directory: " + data_processing_folder.string());
                }
            }else{
                logger.error("Can not check directory '{}' existense due to error with code: {}. Message", data_processing_folder.string(), ec.value());
                logger.error(ec.message());
                throw std::runtime_error("Can not check post-processing directory existense: " + data_processing_folder.string());
            }
        }

        logger.debug("Distributing storages");
        TRY
            distribute(storages, _distrib);
        CATCH("Error while distributing storages")

        logger.debug("Bcasting Natoms");
        MPI_Bcast(&Natoms, 1, MPI_UINT64_T, cs::mpi_root, wcomm);
        logger.debug("Sent Natoms");
        logger.info("Distribution complete");

        return Natoms;
    }

    RETURN_CODES MDN_root::sanity()
    {
        int rnd = std::rand();
        std::cout << "      Bcasting: " << rnd << std::endl;
        MPI_Bcast(&rnd, 1, MPI_INT, cs::mpi_root, wcomm);
        // int responce = 0;
        int *responces = new int[size];
        std::cout << "      Gathering..." << std::endl;
        MPI_Gather(&rnd, 1, MPI_INT, responces, 1, MPI_INT, cs::mpi_root, wcomm);

        if (responces[0] != rnd)
        {
            std::cerr << "Root sanity doesn't passed" << std::endl;
            throw std::runtime_error("Root sanity doesn't passed");
        }

        // logger.info(std::string("Random int: " + std::to_string(rnd)));
        // logger.info(std::string("Received arrray:"));
        for (int i = 0; i < size; ++i)
        {
            // logger.info(std::to_string(responces[i]));
            if (responces[i] - i != rnd)
            {
                std::cerr << "Root sanity doesn't passed" << std::endl;
                throw std::runtime_error("Root sanity doesn't passed");
            }
        }
        std::cout << "      Root sanity OK" << std::endl;
        std::cout << "      Releasing sanity barrier" << std::endl;
        SANITY_BARRIER(wcomm);
        std::cout << "      Sanity barrier released" << std::endl;

        return RETURN_CODES::OK;
    }

    RETURN_CODES MDN_root::setup()
    {
        std::cout << "Sanity:" << std::endl;
        TRY
            sanity();
        CATCH_NOLOGGER("Error while checking sanity")

        std::error_code ec;
        fs::path log_folder = cwd / folders::log;
        if (!fs::exists(log_folder))
        {
            if (!fs::create_directories(log_folder, ec))
            {
                std::cerr << "Cannot create non-existent directories: " << log_folder.string() << std::endl;
                std::cerr << "Error code: " << ec.value() << std::endl;
                std::cerr << "Message: " << ec.message() << std::endl;
                throw std::runtime_error("Cannot create non-existent directory: " + log_folder.string());
            }
            else
            {
                std::cout << "Created logging directory: " << log_folder << std::endl;
            }
        }
        else
        {
            std::cout << "Logging folder exists" << std::endl;
        }

        std::cout << "  Releasing setup barrier" << std::endl;
        SETUP_BARRIER(wcomm);
        std::cout << "  Setup barrier released" << std::endl;

        TRY
            setup_logger(rank);
        CATCH("Error while setting up logger")

        std::cout << "  Root setup OK" << std::endl;
        return RETURN_CODES::OK;
    }

    RETURN_CODES MDN_root::entry()
    {
        std::cout << "Setup..." << std::endl;
        TRY
        setup();
        CATCH_NOLOGGER("Error wahile setting up")
        logger.info("Initialized");

        std::map<fs::path, std::pair<int, int>> distribution;
        uint64_t Natoms{};
        TRY
            Natoms = dns(distribution);
        CATCH("Error while getting distribution")

        fs::path ss = args.outfile / (args.outfile.filename().string() + "." + std::to_string(rank));
        logger.info("Output data will be written to {}", ss.string());

        logger.info("Starting calculations...");
        // TRY
        //     run(ss, distribution, Natoms);
        // CATCH("Error while doing caculations")

        return RETURN_CODES::OK;
    }
} // namespace mdn

#endif // !__MDN_ROOT__
