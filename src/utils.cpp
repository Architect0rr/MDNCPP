// #ifndef __MDNCPP_VERSION__
// #define __MDNCPP_VERSION__ "0.1.0"
// #endif // !__MDNCPP_VERSION__

#ifndef __MDN_utils_CPP__
#define __MDN_utils_CPP__

#include <filesystem>

#include "argparse/argparse.hpp"

#include "spdlog/spdlog.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"

#include "MDN.hpp"

// namespace mdn {
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


    void MDN::setup_logger(const fs::path &log_folder, const int rank)
    {
        fs::path log_file = log_folder / (std::to_string(rank) + ".log");

        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::trace);

        auto file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true);
        file_sink->set_level(spdlog::level::trace);

        logger = spdlog::logger(std::string("Rank") + std::to_string(rank), {console_sink, file_sink});
        logger.set_level(spdlog::level::trace);
    }

    RETURN_CODES MDN::parse_args(const int argc, const char ***argv)
    {
        argparse::ArgumentParser _args = argparse::ArgumentParser("MDNCPP", __MDNCPP_VERSION__, argparse::default_arguments::help & argparse::default_arguments::version, false);

        _args.add_argument("-d", "--debug").help("Turn on debug").flag();
        _args.add_argument("-t", "--trace").help("Turn on trace (enables debug)").flag();
        _args.add_argument("-m", "--mode").help("Mode to run").scan<'i', int>().default_value(0);
        _args.add_argument("-o", "--output").help("Output file basename").default_value(std::string_view("data.bp"));
        _args.add_argument("-w", "--working_directory").help("Current working directory").default_value(std::string_view("./"));
        _args.add_argument("-v", "--verbose").help("Print log to stdout").flag();
        _args.add_argument("-h", "--help").flag().help("Show help (this) message");
        _args.add_argument("--version").flag().help("Show version");
        _args.add_description(std::string("Generate cluster distribution matrix from ADIOS2 LAMMPS data. Version: ") + std::string(__MDNCPP_VERSION__));
        _args.add_epilog("Author: Perevoshchikov Egor, 2023");

        try
        {
            _args.parse_args(argc, *argv);
        }
        catch (const std::exception &err)
        {
            if (rank == 0)
            {
                std::cerr << err.what() << std::endl;
                std::cerr << _args << std::endl;
                std::cout << _args.help().str();
            }
            throw;
        }

        std::error_code ec;
        cwd = fs::current_path(ec);

        if (ec)
        {
            std::cerr << "Cannot get cwd due to error with code: " << ec.value() << std::endl;
            std::cerr << "Message: " << ec.message() << std::endl;
            throw std::runtime_error("Cannot get current working directory due to error");
        }
        else
        {
            if (_args.is_used("--working_directory"))
            {
                fs::path abs_spec_path = fs::absolute(_args.get<std::string>("--working_directory"), ec);
                if (ec)
                {
                    std::cerr << "Cannot resolve path '" << _args.get<std::string>("--working_directory") << "' due to error with code: " << ec.value() << std::endl;
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

        // *pwd = cwd;

        if (argc <= 1 || _args["--help"] == true)
        {
            if (rank == 0)
                std::cout << _args.help().str();
            return RETURN_CODES::EXIT;
        }
        if (_args["--version"] == true)
        {
            if (rank == 0)
                std::cout << __MDNCPP_VERSION__ << std::endl;
            return RETURN_CODES::EXIT;
        }

        if (_args["--trace"] == true || _args["--debug"] == true)
            args.debug = true;

        if (_args["--trace"] == true)
            args.trace = true;

        if (_args["--verbose"] == true)
            args.verbose = true;

        args.mode = _args.get<int>("--mode");

        if (_args.is_used("--output"))
        {
            fs::path outf = fs::absolute(_args.get<std::string>("--output"), ec);
            if (ec)
            {
                std::cerr << "Cannot resolve path '" << _args.get<std::string>("--output") << "' due to error with code: " << ec.value() << std::endl;
                std::cerr << "Message: " << ec.message() << std::endl;
                throw std::runtime_error("Cannot resolve specified output file basename due to error");
            }
            args.outfile = outf;
        }
        else
        {
            args.outfile = cwd / "data.bp";
            if (fs::exists(args.outfile))
            {
                std::cerr << "Default output file '" << args.outfile.string() << "' already exists. Error code: " << ec.value() << std::endl;
                std::cerr << "Message: " << ec.message() << std::endl;
                throw std::runtime_error("Default output file already exists");
            }
        }

        return RETURN_CODES::OK;
    }
// } // namespace

#endif // !__MDN_utils_CPP__