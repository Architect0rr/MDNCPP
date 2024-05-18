#ifndef __MDN_HPP__
#define __MDN_HPP__

#ifndef __MDN_TRY_CATCH_MACRO__
#define __MDN_TRY_CATCH_MACRO__
#define TRY \
    try {
#define CATCH_NOLOGGER(mess)                                                                         \
    }catch (const std::exception &e){                                                                \
        std::cerr << "Rank: " << mpi_rank << " (" << prefix << ") "                                  \
                  << "Some std::exception was thrown, probably aborting. Context:" << std::endl;     \
        std::cerr << mess << std::endl;                                                              \
        std::cerr << e.what() << std::endl;                                                          \
        MPI_Abort(wcomm, static_cast<int>(RETURN_CODES::ERROR));                                     \
    }catch (...){                                                                                    \
        std::cerr << "Rank: " << mpi_rank << " (" << prefix << ") "                                  \
                  << "Some unknown exception was thrown, probably aborting. Context:." << std::endl; \
        std::cerr << mess << std::endl;                                                              \
        MPI_Abort(wcomm, static_cast<int>(RETURN_CODES::ERROR));                                     \
    }
#define CATCH(mess)                                                                     \
    }catch (std::exception & e){                                                        \
        logger.error("Some std::exception was thrown, probably aborting. Context:");    \
        logger.error(mess);                                                             \
        logger.error(e.what());                                                         \
        throw;                                                                          \
    }catch (...){                                                                       \
        logger.error("Some unknown exception was thrown, probably aborting. Context:"); \
        logger.error(mess);                                                             \
        throw;                                                                          \
    }
#define CATCH_NOTHROW(mess)                                          \
    }catch (std::exception & e){                                     \
        logger.error("Some std::exception was thrown. Context:");    \
        logger.error(mess);                                          \
        logger.error(e.what());                                      \
    }catch (...){                                                    \
        logger.error("Some unknown exception was thrown. Context:"); \
        logger.error(mess);                                          \
    }

#define CATCH_NOTHROW_RET(mess, obj)                                 \
    }catch (std::exception & e){                                     \
        logger.error("Some std::exception was thrown. Context:");    \
        logger.error(mess);                                          \
        logger.error(e.what());                                      \
        return obj;                                                  \
    }catch (...){                                                    \
        logger.error("Some unknown exception was thrown. Context:"); \
        logger.error(mess);                                          \
        return obj;                                                  \
    }

#endif // !__MDN_TRY_CATCH_MACRO__


#include <map>
#include <span>
#include <memory>
#include <filesystem>

#include "adios2.h"
#include "nlohmann/json.hpp"

#include "config.hpp"
#include "memory.hpp"
#include "enums.hpp"

namespace mdn{
    extern spdlog::logger logger;
    namespace fs = std::filesystem;
    using json = nlohmann::json;

    struct Arguments{
        bool debug = false;
        bool trace = false;
        int mode = 0;
        uint64_t cut = 0;
        bool verbose = false;
        fs::path outfile = "/";
        fs::path outfile_base = "/";
        fs::path distribution = "/";
        fs::path adios_conf = "";
        fs::path dist_csv = "";
        fs::path temp_csv = "";
        fs::path data_csv = "";
        // bool cache = false;
    };

    class MDN{
    public:
        MDN() = delete;
        MDN(const int rank = 0, const int size = 0, MPI_Comm comm = MPI_COMM_WORLD, const char *prefix = "unknown") noexcept : mpi_rank(rank), mpi_size(size), wcomm(comm), prefix(prefix) {}
        ~MDN() = default;

        RETURN_CODES parse_args(const int argc, const char ***argv);
        void start(const int argc, const char ***argv);

    protected:
        // spdlog::logger logger = spdlog::logger("default");
        // std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> console_sink;
        // std::shared_ptr<spdlog::sinks::basic_file_sink_mt> file_sink;

        const int mpi_rank = 0;
        const int mpi_size = 0;
        const MPI_Comm wcomm = MPI_COMM_WORLD;

        fs::path cwd = "/";
        fs::path wfile = "/";
        std::string_view prefix;
        Arguments args;

        uint64_t _Natoms = 0;
        double _Volume = 0;
        uint64_t max_cluster_size = 0;
        std::map<fs::path, std::pair<int, int>> storages;
        double time_step;
        uint64_t distance;

        MPI_File fh;
        MPI_Status status;
        MPI_Offset off;
        MPI_Request req;
        bool cont = false;
        int amode;

        Memory memory;

        // uint64_t done_steps_primary = 0;
        // std::unique_ptr<uint64_t[]> done_steps;
        uint64_t total_steps = 0;
        int done_stages = 0;

        void setup_logger(const int);
        void run();
        void sstage();

        virtual void sanity();
        virtual void setup();
        virtual void pre_process();
        virtual void post_process();

        // utility functions
        json parse_json(const fs::path&, const std::string&);
        std::tuple<uint64_t, uint64_t, double, double, double, double, double, double, double, double>
        get_row(const uint64_t&, const std::vector<uint64_t>& sizes, const std::span<const uint64_t>&, const double&, const uint64_t&, const double&, const uint64_t&);
    };

    class MDN_root : public MDN{
    public:
        using MDN::MDN;

    private:
        virtual void sanity() override;
        virtual void setup() override;

        virtual void pre_process() override;
        virtual void post_process() override;

        std::map<int, std::string> gather_storages();
    };

} // namespace mdn

#endif // !__MDN_HPP__
