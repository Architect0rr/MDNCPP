#ifndef __MDN_utils_HPP__
#define __MDN_utils_HPP__

#include <filesystem>
#include <chrono>
#include <fstream>
#include <span>

#include "argparse/argparse.hpp"

#include "spdlog/spdlog.h"

#include "MDN.hpp"
#include "constants.hpp"

namespace mdn{
    namespace fs = std::filesystem;

    class existing_string_buf : public std::streambuf {
    public:

        explicit existing_string_buf(std::string &to_append) : m_to_append(&to_append) {}

        virtual int_type overflow(int_type c) {
            if (c != EOF) {
                m_to_append->push_back(c);
            }
            return c;
        }

        virtual std::streamsize xsputn(const char *s, std::streamsize n) {
            m_to_append->insert(m_to_append->end(), s, s + n);
            return n;
        }

    private:
        std::string *m_to_append;
    };

    void dp2s(std::map<int, std::map<fs::path, std::pair<int, int>>> &distribution, std::map<int, std::map<std::string, std::pair<int, int>>> &_distribution){
        for (const auto &[worker, stors] : distribution){
            for (const auto &[stor, bounds] : stors){
                _distribution[worker].insert(std::make_pair(stor.string(), bounds));
            }
        }
    }

    void ds2p(std::map<int, std::map<fs::path, std::pair<int, int>>> &_distribution, std::map<int, std::map<std::string, std::pair<int, int>>> &distribution){
        for (const auto &[worker, stors] : distribution){
            for (const auto &[stor, bounds] : stors){
                _distribution[worker].insert(std::make_pair(fs::path(stor), bounds));
            }
        }
    }

    void MDN::setup_logger(const int rank)
    {
        fs::path log_folder = cwd / folders::log;
        fs::path log_file = log_folder / (std::to_string(rank) + ".log");

        console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(spdlog::level::trace);

        file_sink = std::make_shared<spdlog::sinks::basic_file_sink_mt>(log_file, true);
        file_sink->set_level(spdlog::level::trace);

        logger = spdlog::logger(std::string("Rank") + std::to_string(rank), {console_sink, file_sink});
        // logger = spdlog::logger(std::string("Rank") + std::to_string(rank), {file_sink});
        logger.set_level(spdlog::level::trace);
        // spdlog::flush_every(std::chrono::seconds(10));
    }

    RETURN_CODES MDN::parse_args(const int argc, const char ***argv) {
        argparse::ArgumentParser _args = argparse::ArgumentParser("MDNCPP", __MDNCPP_VERSION__, argparse::default_arguments::help & argparse::default_arguments::version, false);

        _args.add_argument("-d", "--debug").help("Turn on debug").flag();
        _args.add_argument("-t", "--trace").help("Turn on trace (enables debug)").flag();
        _args.add_argument("-a", "--adios_conf").help("Path to ADIOS2 xml config file").default_value(std::string_view("./adios2_BP4_config.xml"));
        _args.add_argument("-D", "--distribution").help("Path to distribution file").default_value(std::string_view("./dist.json"));
        _args.add_argument("-m", "--mode").help("Mode to run").scan<'i', int>().default_value(0);
        _args.add_argument("-o", "--output").help("Output file basename").default_value(std::string_view("data.bp"));
        _args.add_argument("-w", "--working_directory").help("Set current working directory").default_value(std::string_view("./"));
        _args.add_argument("-v", "--verbose").help("Print log to stdout").flag();
        _args.add_argument("-c", "--cache").help("Cache some internal variables (does not affect heavy computations, may speedup startup)").flag();
        _args.add_argument("-h", "--help").flag().help("Show help (this) message");
        _args.add_argument("--version").flag().help("Show version");
        _args.add_description(std::string("Generate cluster distribution matrix from ADIOS2 LAMMPS data. Version: ") + std::string(__MDNCPP_VERSION__));
        _args.add_epilog("Author: Perevoshchikov Egor, 2023-2024");

        if (rank == cs::mpi_root) std::cout << "Rank " << cs::mpi_root << "(root): " << "Parsing args:" << std::endl;

        try
        {
            _args.parse_args(argc, *argv);
        }
        catch (const std::exception &err)
        {
            if (rank == 0) {
                std::cerr << "  └─An error occured while parsing args: " << std::endl;
                std::cerr << err.what() << std::endl;
                std::cerr << _args << std::endl;
                std::cout << _args.help().str();
            }
            throw;
        }

        std::error_code ec;
        cwd = fs::current_path(ec);

        if (ec) {
            if (rank == cs::mpi_root) {
                std::cerr << "  ├─Cannot get cwd due to error with code: " << ec.value() << std::endl;
                std::cerr << "  └─Message: " << ec.message() << std::endl;
            }
            throw std::runtime_error("Cannot get current working directory due to error");
        } else {
            if (_args.is_used("--working_directory")) {
                fs::path abs_spec_path = fs::absolute(_args.get<std::string>("--working_directory"), ec);

                if (ec) {
                    if (rank == cs::mpi_root) {
                        std::cerr << "  ├─Cannot resolve path '" << _args.get<std::string>("--working_directory") << "' due to error with code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: " << ec.message() << std::endl;
                    }
                    throw std::runtime_error("Cannot resolve absolute path to specified working directory due to error");
                }
                if (abs_spec_path != cwd) {
                    fs::current_path(abs_spec_path, ec);
                    if (ec) {
                        if (rank == cs::mpi_root) {
                            std::cerr << "  ├─Cannot change working directory to '" << abs_spec_path.string() << "' due to error with code: " << ec.value() << std::endl;
                            std::cerr << "  └─Message: " << ec.message() << std::endl;
                        }
                        throw std::runtime_error("Cannot get current working directory due to error");
                    }
                    cwd = abs_spec_path;
                }
            }
            if (rank == cs::mpi_root) std::cout << "  ├─Working in directiry: " << cwd.string() << std::endl;
        }

        if (argc <= 1 || _args["--help"] == true) {
            if (rank == 0)
                std::cout << _args.help().str();
            return RETURN_CODES::EXIT;
        }
        if (_args["--version"] == true) {
            if (rank == 0)
                std::cout << __MDNCPP_VERSION__ << std::endl;
            return RETURN_CODES::EXIT;
        }

        if (_args["--trace"] == true || _args["--debug"] == true){
            args.debug = true;
            if (rank == cs::mpi_root) std::cout << "  ├─Enabling debug" << std::endl;
        }

        if (_args["--trace"] == true){
            args.trace = true;
            if (rank == cs::mpi_root) std::cout << "  ├─Enabling trace" << std::endl;
        }

        if (_args["--verbose"] == true){
            args.verbose = true;
            if (rank == cs::mpi_root) std::cout << "  ├─Enabling verbose" << std::endl;
        }

        args.mode = _args.get<int>("--mode");
        if (rank == cs::mpi_root) std::cout << "  ├─Working mode: " << args.mode << std::endl;

        args.distribution = fs::absolute(_args.get<std::string>("--distribution"), ec);
        if (ec) {
            if (rank == cs::mpi_root) {
                std::cerr << "  ├─Cannot resolve distribution file path '" << _args.get<std::string>("--output") << "' due to error with code: " << ec.value() << std::endl;
                std::cerr << "  └─Message: " << ec.message() << std::endl;
            }
            throw std::runtime_error("Cannot resolve specified output file basename due to error");
        }
        if (!fs::exists(args.distribution, ec)) {
            if (ec) {
                if (rank == cs::mpi_root) {
                    std::cerr << "  ├─Can not check distribution file '" << args.outfile.string() << "' existense. Error code: " << ec.value() << std::endl;
                    std::cerr << "  └─Message: " << ec.message() << std::endl;
                }
                throw std::runtime_error("Can not check distribution file existense.");
            } else {
                if (rank == cs::mpi_root) {
                    std::cerr << "  ├─Distribution file '" << args.outfile.string() << "' does not exists. Error code: " << ec.value() << std::endl;
                    std::cerr << "  └─Message: " << ec.message() << std::endl;
                }
                throw std::runtime_error("Distribution file does not exists.");
            }
        }
        if (rank == cs::mpi_root) std::cout << "  ├─Absolute distribution file path: " << args.distribution.string() << std::endl;

        if (_args.is_used("--adios_conf")) {
            args.adios_conf = fs::absolute(_args.get<std::string>("--adios_conf"), ec);

            if (ec) {
                if (rank == cs::mpi_root) {
                    std::cerr << "  ├─Cannot resolve output file path '" << _args.get<std::string>("--output") << "' due to error with code: " << ec.value() << std::endl;
                    std::cerr << "  └─Message: " << ec.message() << std::endl;
                }
                throw std::runtime_error("Cannot resolve specified output file basename due to error");
            }
            if (!fs::exists(args.adios_conf, ec)) {
                if (ec) {
                    if (rank == cs::mpi_root) {
                        std::cerr << "  ├─Can not check distribution file '" << args.outfile.string() << "' existense. Error code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: " << ec.message() << std::endl;
                    }
                    throw std::runtime_error("Can not check distribution file existense.");
                } else {
                    if (rank == cs::mpi_root) {
                        std::cerr << "  ├─Distribution file '" << args.outfile.string() << "' does not exists. Error code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: " << ec.message() << std::endl;
                    }
                    throw std::runtime_error("Distribution file does not exists.");
                }
            }
            if (rank == cs::mpi_root) std::cout << "  ├─ADIOS2 configuration: " << args.adios_conf.string() << std::endl;
        } else {
            if (rank == cs::mpi_root) std::cout << "  ├─ADIOS2 configuration: builtin" << std::endl;
        }

        if (_args.is_used("--output")) {
            args.outfile_base = fs::absolute(_args.get<std::string>("--output"), ec);

            if (ec) {
                if (rank == cs::mpi_root) {
                    std::cerr << "  ├─Cannot resolve output file path '" << _args.get<std::string>("--output") << "' due to error with code: " << ec.value() << std::endl;
                    std::cerr << "  └─Message: " << ec.message() << std::endl;
                }
                throw std::runtime_error("Cannot resolve specified output file basename due to error");
            }
        } else {
            args.outfile_base = cwd / folders::post_process_folder / files::output_datafile_basename;
            if (fs::exists(args.outfile_base, ec)) {
                if (ec) {
                    if (rank == cs::mpi_root) {
                        std::cerr << "  ├─Can not check output file '" << args.outfile_base.string() << "' existense. Error code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: " << ec.message() << std::endl;
                    }
                    throw std::runtime_error("Can not check output file existense.");
                } else {
                    if (rank == cs::mpi_root) {
                        std::cerr << "  ├─Default output file '" << args.outfile_base.string() << "' already exists. Error code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: " << ec.message() << std::endl;
                    }
                    throw std::runtime_error("Default output file already exists.");
                }
            }
        }
        if (rank == cs::mpi_root) std::cout << "  ├─Template output file path: " << args.outfile_base.string() << std::endl;
        args.outfile = args.outfile_base.parent_path() / (args.outfile_base.filename().string() + "." + std::to_string(rank));

        fs::path ofpp = args.outfile.parent_path();
        if (rank == cs::mpi_root){
            if (!fs::exists(ofpp, ec)) {
                if (!ec) {
                    if(!fs::create_directories(ofpp, ec)) {
                        std::cerr << "  ├─Cannot create non-existent post-processing directory(ies): " << ofpp.string() << std::endl;
                        std::cerr << "  ├─Error code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: "    << ec.message() << std::endl;
                        throw std::runtime_error("Cannot create non-existent post-processing directory: " + ofpp.string());
                    }
                    std::cout << "  ├─Created post-processing directory(ies): " << ofpp.string() << std::endl;
                } else {
                    std::cerr << "  ├─Cannot check post-processing directory(ies) existense: " << ofpp.string() << std::endl;
                    std::cerr << "  ├─Error code: " << ec.value() << std::endl;
                    std::cerr << "  └─Message: " << ec.message() << std::endl;
                    throw std::runtime_error("Cannot check post-processing directory existense: " + ofpp.string());
                }
            }
        }

        if (_args["--cache"] == true){
            args.cache = true;
            if (rank == cs::mpi_root) std::cout << "  ├─Cache enabled" << std::endl;
        }

        if (rank == cs::mpi_root) std::cout << "  └─End of parsing args" << std::endl;

        return RETURN_CODES::OK;
    } // parse_args

    json MDN::parse_json(const fs::path& path, const std::string& name){
        logger.debug("Reading {}", name);

        if (!fs::exists(path)){
            logger.error("{} (file: '{}') does not exists", name, path.string().c_str());
            throw std::runtime_error(name + " (file: '" + path.string() + "') does not exists");
        }
        std::ifstream infilestream(path);
        json data;
        TRY
            data = json::parse(infilestream);
        CATCH("Error while parsing json from datafile")
        logger.debug("Successfully read {}", name);

        return data;
    }

    struct timer{
        uint64_t counter = 0;
        std::chrono::system_clock::time_point t1 = std::chrono::high_resolution_clock::now();
        std::chrono::system_clock::time_point t2 = std::chrono::high_resolution_clock::now();
        long double s = 0;
        void start(){
            t1 = std::chrono::high_resolution_clock::now();
            t2 = t1;
            counter = 0;
        }
        void cutoff(){
            counter += 1;
        }
        void b(){
            t1 = std::chrono::high_resolution_clock::now();
            t2 = t1;
        }
        void e(){
            t2 = std::chrono::high_resolution_clock::now();
            counter += 1;
            s += std::chrono::duration<long double, std::milli>(t2 - t1).count();
        }
        long double sum(){
            return s / counter;
        }
        void check(){
            t2 = std::chrono::high_resolution_clock::now();
        }
        long double current(){
            check();
            if (counter > 0)
                return std::chrono::duration<long double, std::milli>(t2 - t1).count() / counter;
            else
                return std::chrono::duration<long double, std::milli>(0).count();
        }
        long double count(){
            return std::chrono::duration<long double, std::milli>(t2 - t1).count();
        }
    };

} // namespace

#endif // !__MDN_utils_HPP__