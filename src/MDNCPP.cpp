#ifndef __MDNCPP_VERSION__
#define __MDNCPP_VERSION__ "0.1.0"
#endif // !__MDNCPP_VERSION__

// Third-party libraries
#include "argparse/argparse.hpp"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"

// Convenient libraries
#include "mpi.h"

// This project libraries
#include "utils.cpp"
#include "constants.cpp"
#include "root.cpp"

// STL libraries
#include <iostream>
#include <filesystem>

namespace mdn
{
    namespace fs = std::filesystem;
    spdlog::logger setup_logger(const fs::path &log_folder, const int rank)
    {
        fs::path log_file = log_folder / (std::to_string(rank) + ".log");

        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::trace);

        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true);
        file_sink->set_level(spdlog::level::trace);

        spdlog::logger logger = spdlog::logger(std::string("Rank") + std::to_string(rank), {console_sink, file_sink});
        logger.set_level(spdlog::level::trace);

        return logger;
    }

    int root_entry(utils::MC &sts)
    {
        // sanity
        int rnd = std::rand();
        MPI_Bcast(&rnd, 1, MPI_INT, cs::mpi_root, MPI_COMM_WORLD);
        int responce = 0;
        int *responces = new int[sts.size];
        MPI_Gather(&responce, 1, MPI_INT, responces, 1, MPI_INT, cs::mpi_root, MPI_COMM_WORLD);

        if (responces[0] != 0)
        {
            std::cerr << "Root sanity doesn't passed" << std::endl;
            return 1;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        // sts.logger.info(std::string("Random int: " + std::to_string(rnd)));
        // sts.logger.info(std::string("Received arrray:"));
        for (size_t i = 1; i < sts.size; ++i)
        {
            // sts.logger.info(std::to_string(responces[i]));
            if (responces[i] - i != rnd)
            {
                std::cerr << "Root sanity doesn't passed" << std::endl;
                return 1;
            }
        }

        std::error_code ec;
        fs::path cwd = fs::current_path(ec);

        if (ec)
        {
            std::cerr << "Cannot get cwd due to error with code: " << ec.value() << std::endl;
            std::cerr << "Message: " << ec.message() << std::endl;
            return 1;
        }
        else
        {
            std::cout << "Working in directiry: " << cwd.string() << std::endl;
        }

        fs::path log_folder = cwd / "logs";
        if (!fs::exists(log_folder))
        {
            if (!fs::create_directory(log_folder, ec))
            {
                std::cerr << "Cannot create non-existent directory: " << log_folder.string() << std::endl;
                std::cerr << "Error code: " << ec.value() << std::endl;
                std::cerr << "Message: " << ec.message() << std::endl;
                return 1;
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);

        spdlog::logger logger = setup_logger(log_folder, sts.rank);

        sts.set_cwd(cwd);
        sts.add_logger(&logger);

        sts.logger->info("Initialized");

        // size_t received = 0;
        // int rcv = 0;
        // while (received < sts.size - 1)
        // {
        //     rcv = 0;
        //     MPI_Iprobe(received, 10, MPI_COMM_WORLD, &rcv, &status);
        //     if (rcv)
        //     {
        //         MPI_Recv(&(responces[received]), 1, MPI_INT, received, 10, MPI_COMM_WORLD, &status);
        //         sts.logger.info(std::to_string(responces[received]));
        //         received += 1;
        //     }
        // }
        root::main(sts);

        return 0;
    }

    int nonroot_entry(utils::MC &sts)
    {
        // sanity
        int rnd;
        MPI_Bcast(&rnd, 1, MPI_INT, cs::mpi_root, MPI_COMM_WORLD);
        // sts.logger.info(std::string("Got number: " + std::to_string(rnd)));
        rnd += sts.rank;
        // sts.logger.info(std::string("Sending back: " + std::to_string(rnd)));
        int *responces = NULL;
        MPI_Gather(&rnd, 1, MPI_INT, responces, 0, MPI_INT, cs::mpi_root, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);

        std::error_code ec;
        fs::path cwd = fs::current_path(ec);

        if (ec)
            return 1;

        fs::path log_folder = cwd / "logs";
        MPI_Barrier(MPI_COMM_WORLD);

        spdlog::logger logger = setup_logger(log_folder, sts.rank);

        sts.set_cwd(cwd);
        sts.add_logger(&logger);

        sts.logger->info("Initialized");

        return 0;
    }

    int main(const int argc, char ***argv, const int rank, const int size)
    {
        argparse::ArgumentParser args = argparse::ArgumentParser("MDNCPP", __MDNCPP_VERSION__, argparse::default_arguments::help & argparse::default_arguments::version, false);

        args.add_argument("-d", "--debug").help("Turn on debug").flag();
        args.add_argument("-t", "--trace").help("Turn on trace (enables debug)").flag();
        args.add_argument("-m", "--mode").help("Mode to run").scan<'i', int>().default_value(0);
        args.add_argument("-v", "--verbose").help("Print log to stdout").flag();
        args.add_argument("-h", "--help").flag().help("Show help (this) message");
        args.add_argument("--version").flag().help("Show version");
        args.add_description(std::string("Generate cluster distribution matrix from ADIOS2 LAMMPS data. Version: ") + std::string(__MDNCPP_VERSION__));
        args.add_epilog("Author: Perevoshchikov Egor, 2023");

        try
        {
            args.parse_args(argc, *argv);
        }
        catch (const std::exception &err)
        {
            if (rank == 0)
            {
                std::cerr << err.what() << std::endl;
                std::cerr << args << std::endl;
            }
            return 1;
        }

        if (argc <= 1 || args["--help"] == true)
        {
            if (rank == 0)
                std::cout << args.help().str();
            return 0;
        }
        if (args["--version"] == true)
        {
            if (rank == 0)
                std::cout << __MDNCPP_VERSION__ << std::endl;
            return 0;
        }

        utils::Arguments argss;
        if (args["--trace"] == true || args["--debug"] == true)
        {
            argss.debug = true;
        }
        if (args["--trace"] == true)
        {
            argss.trace = true;
        }
        if (args["--verbose"] == true)
        {
            argss.verbose = true;
        }
        argss.mode = args.get<int>("--mode");

        utils::MC sts(rank, size, argss);

        int exit_code = 1;
        if (size >= 2)
        {
            if (rank == 0)
            {
                exit_code = root_entry(sts);
            }
            else
            {
                exit_code = nonroot_entry(sts);
            }
        }
        else
        {
            std::cerr << "It is required at least 2 processes to be runned, shutting down" << std::endl;
        }

        return exit_code;
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int exit_code = mdn::main(argc, &argv, rank, size);

    MPI_Finalize();
    return exit_code;
}
