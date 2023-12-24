#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include "argparse/argparse.hpp"
#include "mpi.h"

#include <filesystem>
#include <system_error>
#include <streambuf>
#include <string>
#include <ostream>
#include <iostream>

namespace fs = std::filesystem;

class existing_string_buf : public std::streambuf
{
public:
    // Somehow store a pointer to to_append.
    explicit existing_string_buf(std::string &to_append) : m_to_append(&to_append) {}

    virtual int_type overflow(int_type c)
    {
        if (c != EOF)
        {
            m_to_append->push_back(c);
        }
        return c;
    }

    virtual std::streamsize xsputn(const char *s, std::streamsize n)
    {
        m_to_append->insert(m_to_append->end(), s, s + n);
        return n;
    }

private:
    std::string *m_to_append;
};

struct MC
{
    bool success = true;
    std::error_code ec;
    spdlog::logger logger;
    const int rank = 0;
    const int size = 0;
    const fs::path cwd = "/";
    argparse::ArgumentParser args = argparse::ArgumentParser("MDNCPP", "0.1.0", argparse::default_arguments::all, false);
    MC(const int argc, const char **argv, const int rank, const int size) : logger("Empty"), rank(rank), size(size), cwd(fs::current_path(ec))
    {

        if (ec)
        {
            std::cerr << "Rank: " << rank << ", Cannot get cwd due to error with code: " << ec.value() << std::endl;
            std::cerr << "Rank: " << rank << ", Message: " << ec.message() << std::endl;
        }
        else
        {
            std::cout << "Rank: " << rank << ", Working in directiry: " << cwd.string() << std::endl;
        }

        fs::path log_folder = cwd / "logs";
        if (rank == 0)
        {
            if (!fs::exists(log_folder))
            {
                if (!fs::create_directory(log_folder, ec))
                {
                    std::cerr << "Cannot create nonexisting directory: " << log_folder.string() << std::endl;
                    std::cerr << "Error code: " << ec.value() << std::endl;
                    std::cerr << "Message: " << ec.message() << std::endl;
                }
            }
        }

        fs::path log_file = log_folder / (rank + ".log");

        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::trace);

        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true);
        file_sink->set_level(spdlog::level::trace);

        logger = spdlog::logger(std::string("Rank") + std::to_string(rank), {console_sink, file_sink});
        logger.set_level(spdlog::level::debug);

        try
        {
            args.parse_args(argc, argv);
        }
        catch (const std::exception &err)
        {
            logger.error(err.what());
            std::string s;
            existing_string_buf b(s);
            std::ostream o(&b);

            o << args;
            logger.error(s);
            success = false;
        }
    }
};

int runner(MC &sts)
{

    return 0;
}

int main(int argc, char **argv)
{
    std::cout << "Started" << std::endl;
    std::cout << "Initializing MPI" << std::endl;
    MPI_Init(&argc, &argv);
    std::cout << "MPI initialized" << std::endl;

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MC sts = MC(argc, const_cast<const char **>(argv), rank, size);

    int exit_code = 0;
    if (sts.success)
    {
        sts.logger.info("Init successful");
        exit_code = runner(sts);
    }
    else
    {
        exit_code = 0;
    }

    std::cout << "Rank: " << rank << ", Finalizing MPI" << std::endl;
    MPI_Finalize();
    std::cout << "Rank: " << rank << ", MPI finalized" << std::endl;
    return exit_code;
}