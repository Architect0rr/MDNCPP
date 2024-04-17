// STL libraries
#include <iostream>
#include <filesystem>

// Third-party libraries
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

    namespace fs = std::filesystem;

    void MDN::pre_process(){
        #ifdef __MDN_PROFILING__
            logger.debug("Profiling enabled");
        #else
            logger.debug("Profiling disabled");
        #endif // __MDN_PROFILING__
        #ifdef __CALC_ENTHROPY__
            logger.debug("Enthropy caltulation:  on");
        #else
            logger.debug("Enthropy caltulation:  off");
        #endif // __CALC_ENTHROPY__
        #ifdef __KE_PE_PRESENT__
            logger.debug("KE and PE are present: yes");
        #else
            logger.debug("KE and PE are present: no");
        #endif // __KE_PE_PRESENT__
        #ifdef __CALC_PE__
            logger.debug("PE caltulation:        on");
        #else
            logger.debug("PE caltulation:        off");
        #endif // __CALC_PE__
        #ifdef __MDN_TRACE_OUT__
            logger.debug("Trace out:             on");
        #else
            logger.debug("Trace out:             off");
        #endif // __MDN_TRACE_OUT__

        json data = parse_json(cwd / files::data, "Datafile");
        if (data.contains(fields::Natoms)) _Natoms = data.at(fields::Natoms).get<uint64_t>();
        else throw std::runtime_error("Datafile does not contains '" + std::string(fields::Natoms) + "' field");

        json dist = parse_json(cwd / files::distribution, "Distribution");
        std::map<int, std::map<fs::path, std::pair<int, int>>> distribution;
        TRY
            distribution = dist.get<std::map<int, std::map<fs::path, std::pair<int, int>>>>();
        CATCH("Error while parsing distribution")
        storages = distribution.at(rank);

        if (data.contains(fields::Dhash)){
            const unsigned long Dhash = data.at(fields::Dhash).get<const unsigned long>();
            const unsigned long hash = std::hash<json>{}(dist);
            if (Dhash == hash) cont = true;
        }

        wfile = cwd / "stat";
        int amode = MPI_MODE_RDWR;
        if (!fs::exists(wfile)) {
            amode = MPI_MODE_CREATE | amode;
            cont = false;
        }
        off = rank * (sizeof(int) + sizeof(uint64_t));
    }

    void MDN_root::pre_process(){
        MDN::pre_process();

        json data = parse_json(cwd / files::data, "Datafile");
        json dist = parse_json(cwd / files::distribution, "Distribution");
        std::vector<std::string> outfiles;
        for (int i = 0; i < size; ++i) outfiles.push_back((args.outfile_base.parent_path() / (args.outfile_base.filename().string() + "." + std::to_string(i))).string());
        data[fields::outfiles] = outfiles;
        const unsigned long hash = std::hash<json>{}(dist);
        data[fields::Dhash] = hash;
        logger.debug("Updating datafile");
        TRY
            std::ofstream datafile_out(cwd / files::data);
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
    }

    void MDN::post_process(){
        uint64_t *max_sizes = new uint64_t[size];
        MPI_Gather(&max_cluster_size, 1, MPI_UINT64_T, max_sizes, 1, MPI_UINT64_T, cs::mpi_root, wcomm);
        max_cluster_size = *std::max_element(max_sizes, max_sizes + size);
        delete[] max_sizes;
        logger.debug("Gathered max cluster size");

        // MPI_Send(args.outfile.string().c_str(), args.outfile.string().length(), MPI_CHAR, cs::mpi_root, MPI_TAGS::STOR, wcomm);
        // logger.debug("Gathered storages");

        logger.info("Exiting entry point. NO RETURN");
    }

    void MDN_root::post_process(){
        MDN::post_process();

        // std::map<int, std::string> storages = gather_storages();
        // logger.debug("Gathered storages");

        json ddata = parse_json(cwd / files::data, "Datafile");
        ddata[fields::maxclsize] = max_cluster_size;
        ddata[fields::outfiles] = storages;

        logger.debug("Updating datafile");
        TRY
            std::ofstream datafile_out(cwd / files::data);
            datafile_out << std::setw(4) << ddata << std::endl;
        CATCH("Error while writing datafile")

        logger.info("Exiting entry point. NO RETURN");
    }

    void MDN::start(const int argc, const char ***argv){
        TRY
            if (parse_args(argc, argv) == RETURN_CODES::EXIT) return;
        CATCH_NOLOGGER("Error while parsing args")
        if (rank == cs::mpi_root) std::cout << "Setup..." << std::endl;
        TRY
            setup();
            logger.info("-!-Initialized-!-");
        CATCH_NOLOGGER("Error while setting up")
        TRY
            pre_process();
        CATCH("Error in entry point")
        TRY
            logger.info("Starting calculations...");
            run();
            logger.info("Calculations ended");
        CATCH("Error while doing caculations")
        TRY
            post_process();
        CATCH("Error in entry point")
        logger.info("-!-Worker shutdown-!-");
        spdlog::shutdown();
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
