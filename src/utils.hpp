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

    std::map<int, std::map<fs::path, std::pair<int, int>>> parse_dist(const json& j) {
        std::map<int, std::map<fs::path, std::pair<int, int>>> result;
        std::map<fs::path, std::pair<int, int>> innerMap;

        for (const auto& [key, value] : j.items()) {
            for (const auto& [innerKey, innerValue] : value.items())
                innerMap[fs::path(innerKey)] = std::make_pair(static_cast<int>(innerValue[0]), static_cast<int>(innerValue[1]));

            result[std::stoi(key)] = innerMap;
            innerMap.clear();
        }
        return result;
    }

    RETURN_CODES MDN::parse_args(const int argc, const char ***argv) {
        argparse::ArgumentParser _args = argparse::ArgumentParser("MDNCPP", __MDNCPP_VERSION__, argparse::default_arguments::help & argparse::default_arguments::version, false);

        _args.add_argument("-d", "--debug").help("Turn on debug").flag();
        _args.add_argument("-t", "--trace").help("Turn on trace (enables debug)").flag();
        _args.add_argument("-a", "--adios_conf").help("Path to ADIOS2 xml config file").default_value("./adios2_BP4_config.xml");
        _args.add_argument("-D", "--distribution").help("Path to distribution file").default_value("./dist.json");
        _args.add_argument("-m", "--mode").help("Mode to run").scan<'i', int>().default_value(0);
        _args.add_argument("-p", "--cut").help("Cut output matrices to this cluster size").scan<'i', uint64_t>().default_value(0);
        _args.add_argument("-w", "--working_directory").help("Set current working directory").default_value("./");
        _args.add_argument("-v", "--verbose").help("Print log to stdout").flag();
        _args.add_argument("-h", "--help").help("Show help (this) message").flag();
        _args.add_argument("--version").help("Show version").flag();
        _args.add_description(std::string("Generate cluster distribution matrix from ADIOS2 LAMMPS data. Version: ") + std::string(__MDNCPP_VERSION__));
        _args.add_epilog("Author: Perevoshchikov Egor, 2023-2024");

        if (mpi_rank == cs::mpi_root) std::cout << "Rank " << cs::mpi_root << "(root): " << "Parsing args:" << std::endl;

        try
        {
            _args.parse_args(argc, *argv);
        }
        catch (const std::exception &err)
        {
            if (mpi_rank == 0) {
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
            if (mpi_rank == cs::mpi_root) {
                std::cerr << "  ├─Cannot get cwd due to error with code: " << ec.value() << std::endl;
                std::cerr << "  └─Message: " << ec.message() << std::endl;
            }
            throw std::runtime_error("Cannot get current working directory due to error");
        } else {
            if (_args.is_used("--working_directory")) {
                fs::path abs_spec_path = fs::absolute(_args.get<std::string>("--working_directory"), ec);

                if (ec) {
                    if (mpi_rank == cs::mpi_root) {
                        std::cerr << "  ├─Cannot resolve path '" << _args.get<std::string>("--working_directory") << "' due to error with code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: " << ec.message() << std::endl;
                    }
                    throw std::runtime_error("Cannot resolve absolute path to specified working directory due to error");
                }
                if (abs_spec_path != cwd) {
                    fs::current_path(abs_spec_path, ec);
                    if (ec) {
                        if (mpi_rank == cs::mpi_root) {
                            std::cerr << "  ├─Cannot change working directory to '" << abs_spec_path.string() << "' due to error with code: " << ec.value() << std::endl;
                            std::cerr << "  └─Message: " << ec.message() << std::endl;
                        }
                        throw std::runtime_error("Cannot get current working directory due to error");
                    }
                    cwd = abs_spec_path;
                }
            }
            if (mpi_rank == cs::mpi_root) std::cout << "  ├─Working in directiry: " << cwd.string() << std::endl;
        }

        if (argc <= 1 || _args["--help"] == true) {
            if (mpi_rank == 0)
                std::cout << _args.help().str();
            return RETURN_CODES::EXIT;
        }
        if (_args["--version"] == true) {
            if (mpi_rank == 0)
                std::cout << __MDNCPP_VERSION__ << std::endl;
            return RETURN_CODES::EXIT;
        }

        if (_args["--trace"] == true || _args["--debug"] == true){
            args.debug = true;
            if (mpi_rank == cs::mpi_root) std::cout << "  ├─Enabling debug" << std::endl;
        }

        if (_args["--trace"] == true){
            args.trace = true;
            if (mpi_rank == cs::mpi_root) std::cout << "  ├─Enabling trace" << std::endl;
        }

        if (_args["--verbose"] == true){
            args.verbose = true;
            if (mpi_rank == cs::mpi_root) std::cout << "  ├─Enabling verbose" << std::endl;
        }

        args.mode = _args.get<int>("--mode");
        if (mpi_rank == cs::mpi_root) std::cout << "  ├─Working mode: " << args.mode << std::endl;

        args.cut = _args.get<uint64_t>("--cut");

        args.distribution = fs::absolute(_args.get<std::string>("--distribution"), ec);
        if (ec) {
            if (mpi_rank == cs::mpi_root) {
                std::cerr << "  ├─Cannot resolve distribution file path '" << _args.get<std::string>("--output") << "' due to error with code: " << ec.value() << std::endl;
                std::cerr << "  └─Message: " << ec.message() << std::endl;
            }
            throw std::runtime_error("Cannot resolve specified output file basename due to error");
        }
        if (!fs::exists(args.distribution, ec)) {
            if (ec) {
                if (mpi_rank == cs::mpi_root) {
                    std::cerr << "  ├─Can not check distribution file '" << args.distribution.string() << "' existense. Error code: " << ec.value() << std::endl;
                    std::cerr << "  └─Message: " << ec.message() << std::endl;
                }
                throw std::runtime_error("Can not check distribution file existense.");
            } else {
                if (mpi_rank == cs::mpi_root) {
                    std::cerr << "  ├─Distribution file '" << args.distribution.string() << "' does not exists. Error code: " << ec.value() << std::endl;
                    std::cerr << "  └─Message: " << ec.message() << std::endl;
                }
                throw std::runtime_error("Distribution file does not exists.");
            }
        }
        if (mpi_rank == cs::mpi_root) std::cout << "  ├─Absolute distribution file path: " << args.distribution.string() << std::endl;

        if (_args.is_used("--adios_conf")) {
            args.adios_conf = fs::absolute(_args.get<std::string>("--adios_conf"), ec);

            if (ec) {
                if (mpi_rank == cs::mpi_root) {
                    std::cerr << "  ├─Cannot resolve output file path '" << _args.get<std::string>("--output") << "' due to error with code: " << ec.value() << std::endl;
                    std::cerr << "  └─Message: " << ec.message() << std::endl;
                }
                throw std::runtime_error("Cannot resolve specified output file basename due to error");
            }
            if (!fs::exists(args.adios_conf, ec)) {
                if (ec) {
                    if (mpi_rank == cs::mpi_root) {
                        std::cerr << "  ├─Can not check distribution file '" << args.outfile.string() << "' existense. Error code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: " << ec.message() << std::endl;
                    }
                    throw std::runtime_error("Can not check distribution file existense.");
                } else {
                    if (mpi_rank == cs::mpi_root) {
                        std::cerr << "  ├─Distribution file '" << args.outfile.string() << "' does not exists. Error code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: " << ec.message() << std::endl;
                    }
                    throw std::runtime_error("Distribution file does not exists.");
                }
            }
            if (mpi_rank == cs::mpi_root) std::cout << "  ├─ADIOS2 configuration: " << args.adios_conf.string() << std::endl;
        } else {
            if (mpi_rank == cs::mpi_root) std::cout << "  ├─ADIOS2 configuration: builtin" << std::endl;
        }

        args.outfile_base = cwd / folders::post_process_folder;
        if (mpi_rank == cs::mpi_root){
            if (!fs::exists(args.outfile_base, ec)) {
                if (!ec) {
                    if(!fs::create_directories(args.outfile_base, ec)) {
                        std::cerr << "  ├─Cannot create non-existent post-processing directory(ies): " << args.outfile_base.string() << std::endl;
                        std::cerr << "  ├─Error code: " << ec.value() << std::endl;
                        std::cerr << "  └─Message: "    << ec.message() << std::endl;
                        throw std::runtime_error("Cannot create non-existent post-processing directory: " + args.outfile_base.string());
                    }
                    std::cout << "  ├─Created post-processing directory(ies): " << args.outfile_base.string() << std::endl;
                } else {
                    std::cerr << "  ├─Cannot check post-processing directory(ies) existense: " << args.outfile_base.string() << std::endl;
                    std::cerr << "  ├─Error code: " << ec.value() << std::endl;
                    std::cerr << "  └─Message: " << ec.message() << std::endl;
                    throw std::runtime_error("Cannot check post-processing directory existense: " + args.outfile_base.string());
                }
            }
        }

        args.outfile = args.outfile_base / ("data.bp." + std::to_string(mpi_rank));
        args.dist_csv = args.outfile_base / (std::string(files::dist_csv) + "." + std::to_string(mpi_rank));
        args.temp_csv = args.outfile_base / (std::string(files::temp_csv) + "." + std::to_string(mpi_rank));
        args.data_csv = args.outfile_base / (std::string(files::comp_data) + "." + std::to_string(mpi_rank));

        if (mpi_rank == cs::mpi_root) std::cout << "  └─End of parsing args" << std::endl;

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
            if (counter == 0) return 0;
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
                return std::chrono::duration<long double, std::milli>(t2 - t1).count();
        }
        std::chrono::duration<long double, std::milli> current_dur(){
            check();
            if (counter > 0)
                return std::chrono::duration<long double, std::milli>(t2 - t1) / counter;
            else
                return std::chrono::duration<long double, std::milli>(t2 - t1);
        }
        long double count(){
            return std::chrono::duration<long double, std::milli>(t2 - t1).count();
        }
        // std::chrono::duration<uint64_t, std::milli> count_dur(){
        //     return std::chrono::duration<uint64_t, std::milli>(t2 - t1).count();
        // }
    };

    std::string format_duration(std::chrono::duration<long double, std::milli> ms) {
        std::chrono::seconds secs = std::chrono::duration_cast<std::chrono::seconds>(ms);
        ms -= std::chrono::duration_cast<std::chrono::milliseconds>(secs);
        std::chrono::minutes mins = std::chrono::duration_cast<std::chrono::minutes>(secs);
        secs -= std::chrono::duration_cast<std::chrono::seconds>(mins);
        std::chrono::hours hour = std::chrono::duration_cast<std::chrono::hours>(mins);
        mins -= std::chrono::duration_cast<std::chrono::minutes>(hour);

        return std::to_string(hour.count()) + ":" + std::to_string(mins.count()) + ":" + std::to_string(secs.count()) + ":" + std::to_string(ms.count());
    }

    void trace(spdlog::logger& logger, const char *file, int line, const std::string& mess){
        logger.trace(mess + ": " + file + ":" + std::to_string(line));
        logger.flush();
    }

} // namespace

#endif // !__MDN_utils_HPP__