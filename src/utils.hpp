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

    RETURN_CODES MDN::create_d_if_not(const fs::path& folder, const std::string& prefix = ""){
        // uint64_t token = std::chrono::system_clock::now().time_since_epoch().count();
        std::error_code ec;
        if (!fs::exists(folder, ec)){
            if (!ec){
                logger.info("{} directory '{}' does not exists, creating", prefix, folder.string());
                if (!fs::create_directories(folder, ec)){
                    logger.warn("Can not create {} directory '{}' due to error with code: {}. Message:", prefix, folder.string(), ec.value());
                    logger.warn(ec.message());
                    return RETURN_CODES::ERROR;
                }else{
                    return RETURN_CODES::OK;
                }
            }else{
                logger.warn("Can not check {} directory '{}' existense due to error with code: {}. Message:", prefix, folder.string(), ec.value());
                logger.warn(ec.message());
                return RETURN_CODES::ERROR;
            }
        } else return RETURN_CODES::OK;
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

    RETURN_CODES MDN::parse_args(const int argc, const char ***argv)
    {
        argparse::ArgumentParser _args = argparse::ArgumentParser("MDNCPP", __MDNCPP_VERSION__, argparse::default_arguments::help & argparse::default_arguments::version, false);

        _args.add_argument("-d", "--debug").help("Turn on debug").flag();
        _args.add_argument("-t", "--trace").help("Turn on trace (enables debug)").flag();
        _args.add_argument("-m", "--mode").help("Mode to run").scan<'i', int>().default_value(0);
        _args.add_argument("-o", "--output").help("Output file basename").default_value(std::string_view("data.bp"));
        _args.add_argument("-w", "--working_directory").help("Current working directory").default_value(std::string_view("./"));
        _args.add_argument("-v", "--verbose").help("Print log to stdout").flag();
        _args.add_argument("-c", "--cache").help("Cache some internal variables (does not affect heavy computations, may speedup startup)").flag();
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
            args.outfile = cwd / folders::post_process_folder / files::output_datafile_basename;
            if (fs::exists(args.outfile, ec)){
                if (ec){
                    std::cerr << "Can not check output file '" << args.outfile.string() << "' existense. Error code: " << ec.value() << std::endl;
                    std::cerr << "Message: " << ec.message() << std::endl;
                    throw std::runtime_error("Can not check output file existense.");
                }else{
                    std::cerr << "Default output file '" << args.outfile.string() << "' already exists. Error code: " << ec.value() << std::endl;
                    std::cerr << "Message: " << ec.message() << std::endl;
                    throw std::runtime_error("Default output file already exists.");
                }
            }
        }
        fs::path ofpp = args.outfile.parent_path();
        if (rank == cs::mpi_root){
            if (!fs::exists(ofpp, ec)){
                if (!ec){
                    if(!fs::create_directories(ofpp, ec)){
                        std::cerr << "Cannot create non-existent post-processing directory(ies): " << ofpp.string() << std::endl;
                        std::cerr << "Error code: " << ec.value() << std::endl;
                        std::cerr << "Message: "    << ec.message() << std::endl;
                        throw std::runtime_error("Cannot create non-existent post-processing directory: " + ofpp.string());
                        std::cout << "Created post-processing directory(ies): " << ofpp << std::endl;
                    }
                }else{
                    std::cerr << "Cannot check post-processing directory(ies) existense: " << ofpp.string() << std::endl;
                    std::cerr << "Error code: " << ec.value() << std::endl;
                    std::cerr << "Message: " << ec.message() << std::endl;
                    throw std::runtime_error("Cannot check post-processing directory existense: " + ofpp.string());
                }
            }
        }

        if (_args["--cache"] == true)
            args.cache = true;

        return RETURN_CODES::OK;
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

    struct csvWriter{
        std::ofstream str;
        csvWriter(const fs::path& file){str.open(file, std::ios_base::out);}
        ~csvWriter() {str.close();}

        template<typename T>
        csvWriter& operator<<(const std::span<T>& arr){
        for (size_t i = 0; i < arr.size() - 1; ++i)
            str << arr[i] << ',';
            str << arr[arr.size() - 1];
            return *this;
        }

        template<typename T>
        csvWriter& operator<<(const T& el){
            str << el << ',';
            return *this;
        }

        void flush(){
            str.flush();
        }
    };

} // namespace

#endif // !__MDN_utils_HPP__