#ifndef __MDNCPP_VERSION__
#define __MDNCPP_VERSION__ "0.1.0"
#endif // !__MDNCPP_VERSION__

// Third-party libraries
#include "argparse/argparse.hpp"
#include "spdlog/spdlog.h"
// #include "spdlog/sinks/stdout_color_sinks.h"
// #include "spdlog/sinks/basic_file_sink.h"

// Convenient libraries
#include "mpi.h"

// This project libraries
#include "utils.cpp"
// #include "constants.cpp"
#include "root.cpp"
#include "nonroot.cpp"

// STL libraries
#include <iostream>
#include <filesystem>

namespace mdn
{
    namespace fs = std::filesystem;
    RETURN_CODES parse_args(const int argc, char ***argv, const int rank, const int size, utils::Arguments *argss, fs::path *pwd)
    {
        argparse::ArgumentParser args = argparse::ArgumentParser("MDNCPP", __MDNCPP_VERSION__, argparse::default_arguments::help & argparse::default_arguments::version, false);

        args.add_argument("-d", "--debug").help("Turn on debug").flag();
        args.add_argument("-t", "--trace").help("Turn on trace (enables debug)").flag();
        args.add_argument("-m", "--mode").help("Mode to run").scan<'i', int>().default_value(0);
        args.add_argument("-o", "--output").help("Output file basename").default_value(std::string_view("data.bp"));
        args.add_argument("-w", "--working_directory").help("Current working directory").default_value(std::string_view("./"));
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
                std::cout << args.help().str();
            }
            throw;
        }

        std::error_code ec;
        fs::path cwd = fs::current_path(ec);

        if (ec)
        {
            std::cerr << "Cannot get cwd due to error with code: " << ec.value() << std::endl;
            std::cerr << "Message: " << ec.message() << std::endl;
            throw std::runtime_error("Cannot get current working directory due to error");
        }
        else
        {
            if (args.is_used("--working_directory"))
            {
                fs::path abs_spec_path = fs::absolute(args.get<std::string>("--working_directory"), ec);
                if (ec)
                {
                    std::cerr << "Cannot resolve path '" << args.get<std::string>("--working_directory") << "' due to error with code: " << ec.value() << std::endl;
                    std::cerr << "Message: " << ec.message() << std::endl;
                    throw std::runtime_error("Cannot resolve absolute path to specified working directory due to error");
                }
                if (abs_spec_path != cwd)
                {
                    fs::current_path(abs_spec_path, ec);
                    if (ec)
                    {
                        std::cerr << "Cannot change working directory to '" << abs_spec_path.string() << "' due to error with code: " << ec.value() << std::endl;
                        std::cerr << "Message: " << ec.message() << std::endl;
                        throw std::runtime_error("Cannot get current working directory due to error");
                    }
                    cwd = abs_spec_path;
                }
            }
            std::cout << "  Working in directiry: " << cwd.string() << std::endl;
        }

        *pwd = cwd;

        if (argc <= 1 || args["--help"] == true)
        {
            if (rank == 0)
                std::cout << args.help().str();
            return RETURN_CODES::EXIT;
        }
        if (args["--version"] == true)
        {
            if (rank == 0)
                std::cout << __MDNCPP_VERSION__ << std::endl;
            return RETURN_CODES::EXIT;
        }

        if (args["--trace"] == true || args["--debug"] == true)
            argss->debug = true;

        if (args["--trace"] == true)
            argss->trace = true;

        if (args["--verbose"] == true)
            argss->verbose = true;

        argss->mode = args.get<int>("--mode");

        if (args.is_used("--output"))
        {
            fs::path outf = fs::absolute(args.get<std::string>("--output"), ec);
            if (ec)
            {
                std::cerr << "Cannot resolve path '" << args.get<std::string>("--output") << "' due to error with code: " << ec.value() << std::endl;
                std::cerr << "Message: " << ec.message() << std::endl;
                throw std::runtime_error("Cannot resolve specified output file basename due to error");
            }
            argss->outfile = outf;
        }
        else
        {
            argss->outfile = cwd / "data.bp";
            if (fs::exists(argss->outfile)){
                std::cerr << "Default output file '" << argss->outfile.string() << "' already exists. Error code: " << ec.value() << std::endl;
                std::cerr << "Message: " << ec.message() << std::endl;
                throw std::runtime_error("Default output file already exists");
            }
        }

        return RETURN_CODES::OK;
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    mdn::RETURN_CODES exit_code = mdn::RETURN_CODES::ERROR;
    mdn::utils::Arguments args;
    std::filesystem::path cwd;
    mdn::RETURN_CODES ecc;
    try
    {
        ecc = mdn::parse_args(argc, &argv, rank, size, &args, &cwd);
    }
    catch (const std::exception &err)
    {
        std::cerr << "Rank: " << rank << "Execution was aborted due to low-level exception." << std::endl;
        std::cerr << err.what() << std::endl;
        MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
    }

    if (ecc == mdn::RETURN_CODES::OK)
    {
        mdn::utils::MC sts(rank, size, args, cwd);
        if (size >= 2)
        {
            if (rank == mdn::cs::mpi_root)
            {
                try
                {
                    exit_code = mdn::root::entry(sts);
                }
                catch (const std::exception &err)
                {
                    std::cerr << "Rank: " << rank << "(root) Execution was aborted due to low-level exception." << std::endl;
                    std::cerr << err.what() << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
                }
            }
            else
            {
                try
                {
                    exit_code = mdn::nonroot::entry(sts);
                }
                catch (const std::exception &err)
                {
                    std::cerr << "Rank: " << rank << "(nonroot) Execution was aborted due to low-level exception." << std::endl;
                    std::cerr << err.what() << std::endl;
                    MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
                }
            }

            if (exit_code != mdn::RETURN_CODES::OK)
            {
                std::cerr << "Rank: " << rank << " Execution was aborted due to non-zero return code." << std::endl;
                MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
            }
        }
        else
        {
            std::cerr << "It is required at least 2 processes to be runned, shutting down" << std::endl;
        }
    }
    else if (ecc == mdn::RETURN_CODES::EXIT)
        exit_code = mdn::RETURN_CODES::OK;

    MPI_Finalize();
    return exit_code;
}
