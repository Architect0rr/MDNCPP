// STL libraries
#include <iostream>
#include <filesystem>

// Third-party librariesq
#include "argparse/argparse.hpp"
#include "spdlog/spdlog.h"

// Convenient libraries
#include "mpi.h"

// This project headers
#include "MDN.hpp"
#include "utils.hpp"
#include "setup.hpp"
#include "logic.hpp"

namespace mdn{

    RETURN_CODES MDN::entry(){
        json data = parse_json(cwd / files::data, "Datafile");
        uint64_t Natoms = data.at(fields::Natoms).get<uint64_t>();
        if (rank == cs::mpi_root) logger.debug("Number of atoms: {}", Natoms);

        json dist = parse_json(cwd / files::distribution, "Dictribution");
        std::map<int, std::map<fs::path, std::pair<int, int>>> distribution = dist.get<std::map<int, std::map<fs::path, std::pair<int, int>>>>();

        fs::path outfile = args.outfile.parent_path() / (args.outfile.filename().string() + "." + std::to_string(rank));
        logger.info("Output data will be written to {}", outfile.string());

        uint64_t max_cluster_size = 0;
        logger.info("Starting calculations...");
        TRY
            run(outfile, distribution.at(rank), Natoms, max_cluster_size);
        CATCH("Error while doing caculations")
        logger.info("Calculations ended, gathering result info");

        uint64_t *max_sizes = NULL;
        MPI_Gather(&max_cluster_size, 1, MPI_UINT64_T, max_sizes, 1, MPI_UINT64_T, cs::mpi_root, wcomm);
        logger.debug("Gathered max cluster size");

        MPI_Send(outfile.string().c_str(), outfile.string().length(), MPI_CHAR, cs::mpi_root, MPI_TAGS::STOR, wcomm);
        logger.debug("Gathered storages");

        logger.info("Exiting entry point. NO RETURN");
        return RETURN_CODES::OK;
    }

    RETURN_CODES MDN_root::entry(){
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

        json data = parse_json(cwd / files::data, "Datafile");
        uint64_t Natoms = data.at(fields::Natoms).get<uint64_t>();
        if (rank == cs::mpi_root) logger.debug("Number of atoms: {}", Natoms);

        json dist = parse_json(cwd / files::distribution, "Dictribution");
        std::map<int, std::map<fs::path, std::pair<int, int>>> distribution = dist.get<std::map<int, std::map<fs::path, std::pair<int, int>>>>();

        fs::path outfile = args.outfile.parent_path() / (args.outfile.filename().string() + "." + std::to_string(rank));
        logger.info("Output data will be written to {}", outfile.string());

        uint64_t max_cluster_size = 0;
        logger.info("Starting calculations...");
        TRY
            run(outfile, distribution.at(rank), Natoms, max_cluster_size);
        CATCH("Error while doing caculations")
        logger.info("Calculations ended, gathering result info");

        uint64_t *max_sizes = new uint64_t[size];
        MPI_Gather(&max_cluster_size, 1, MPI_UINT64_T, max_sizes, 1, MPI_UINT64_T, cs::mpi_root, wcomm);
        max_cluster_size = *std::max_element(max_sizes, max_sizes + size);
        delete max_sizes;
        logger.debug("Gathered max cluster size");

        std::map<int, std::string> storages = gather_storages();
        storages.emplace(rank, outfile.string().c_str());
        logger.debug("Gathered storages");

        TRY
            json ddata = parse_json(cwd / files::data, "Datafile");
            ddata[fields::maxclsize] = max_cluster_size;
            ddata[fields::outfiles] = storages;
            logger.debug("Updating datafile");
            std::ofstream datafile_out(cwd / files::data);
            datafile_out << std::setw(4) << ddata << std::endl;
        CATCH("Error while writing datafile")

        logger.info("Exiting entry point. NO RETURN");
        return RETURN_CODES::OK;
    }

    RETURN_CODES MDN::start(const int argc, const char ***argv){
        TRY
            parse_args(argc, argv);
        CATCH_NOLOGGER("Error while parsing args")
        if (rank == cs::mpi_root) std::cout << "Setup..." << std::endl;
        TRY
            setup();
        CATCH_NOLOGGER("Error while setting up")
        logger.info("-!-Initialized-!-");
        TRY
            entry();
        CATCH_NOLOGGER("Error in entry point")
        logger.info("Worker shutdown");
        spdlog::shutdown();
        return RETURN_CODES::OK;
    }

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size >= 2) {
        try {
            if (rank == mdn::cs::mpi_root) {
                mdn::MDN_root mmd(rank, size, MPI_COMM_WORLD, "root");
                mmd.start(argc, const_cast<const char***>(&argv));
            } else {
                mdn::MDN mmd(rank, size, MPI_COMM_WORLD, "nonroot");
                mmd.start(argc, const_cast<const char***>(&argv));
            }
        } catch (const std::exception &e) {
            std::cerr << "Some low-level std::exception. Aborting..." << std::endl;
            std::cerr << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
        } catch (...) {
            std::cerr << "Some low-level unknown exception. Aborting..." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
        }
    } else {
        std::cerr << "It is required at least 2 processes to be runned, shutting down" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
    }

    MPI_Finalize();
    return 0;
}
