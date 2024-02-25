#ifndef __MDN_HPP__
#define __MDN_HPP__

#ifndef __MDN_TRY_CATCH_MACRO__
#define __MDN_TRY_CATCH_MACRO__
#define TRY \
    try {
#define CATCH_NOLOGGER(mess)                                                                         \
    }catch (const std::exception &e){                                                                \
        std::cerr << "Rank: " << rank << " (" << prefix << ") "                                      \
                  << "Some std::exception was thrown, probably aborting. Context:" << std::endl;     \
        std::cerr << mess << std::endl;                                                              \
        std::cerr << e.what() << std::endl;                                                          \
        MPI_Abort(wcomm, RETURN_CODES::ERROR);                                                       \
    }catch (...){                                                                                    \
        std::cerr << "Rank: " << rank << " (" << prefix << ") "                                      \
                  << "Some unknown exception was thrown, probably aborting. Context:." << std::endl; \
        std::cerr << mess << std::endl;                                                              \
        MPI_Abort(wcomm, RETURN_CODES::ERROR);                                                       \
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

#include <filesystem>
#include <map>

#include "utils.hpp"
#include "adios2.h"
#include "nlohmann/json.hpp"

#include "spdlog/spdlog.h"

#include "config.hpp"

namespace mdn{

    namespace fs = std::filesystem;
    using json = nlohmann::json;

    struct Arguments{
        bool debug = false;
        bool trace = false;
        int mode = 0;
        bool verbose = false;
        fs::path outfile = "/";
        bool cache = false;
    };

    class MDN{
    public:
        MDN() = delete;
        MDN(const int rank = 0, const int size = 0, MPI_Comm comm = MPI_COMM_WORLD, const char *prefix = "unknown") noexcept : rank(rank), size(size), wcomm(comm), prefix(prefix) {}
        ~MDN() = default;

        RETURN_CODES parse_args(const int argc, const char ***argv);
        RETURN_CODES start(const int argc, const char ***argv);

    protected:
        spdlog::logger logger = spdlog::logger("default");
        std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> console_sink;
        std::shared_ptr<spdlog::sinks::basic_file_sink_mt> file_sink;
        fs::path cwd = "/";
        Arguments args;
        const int rank = 0;
        const int size = 0;
        const MPI_Comm wcomm = MPI_COMM_WORLD;
        std::string_view prefix;

        void setup_logger(const int);
        RETURN_CODES run(const fs::path &, const std::map<fs::path, std::pair<int, int>> &, const uint64_t);

        // funcitons different for (non) root workers
        virtual RETURN_CODES entry() = 0;
        virtual RETURN_CODES sanity() = 0;
        virtual RETURN_CODES setup() = 0;
        virtual uint64_t dns(std::map<fs::path, std::pair<int, int>> &) = 0;

        // utility functions
        RETURN_CODES create_d_if_not(const fs::path&, const std::string&);
        RETURN_CODES check_existense(const fs::path& folder, const std::string&);
    };

    class MDN_root : public MDN{
    public:
        using MDN::MDN;

    private:
        RETURN_CODES entry() override;
        RETURN_CODES sanity() override;
        RETURN_CODES setup() override;
        uint64_t dns(std::map<fs::path, std::pair<int, int>> &) override;

        const uint64_t bearbeit(const fs::path &, double *);
        const std::map<fs::path, int> storage_rsolve(const json &);
        const int get_total_steps(adios2::IO &, const fs::path &);
        RETURN_CODES distribute(const std::map<fs::path, int> &, std::map<fs::path, std::pair<int, int>> &);
        std::map<int, std::map<fs::path, std::pair<int, int>>> make_distribution(const std::map<fs::path, int> &);
        const bool save_distribution(const std::map<int, std::map<std::string, std::pair<int, int>>>&);
        const bool load_distribution(std::map<int, std::map<std::string, std::pair<int, int>>>&);
        const bool save_storages(const std::map<fs::path, int>&, const unsigned long);
        const bool load_storages(std::map<fs::path, int>&, const unsigned long);
    };

    class MDN_nonroot : public MDN{
    public:
        using MDN::MDN;

    private:
        RETURN_CODES entry() override;
        RETURN_CODES sanity() override;
        RETURN_CODES setup() override;
        uint64_t dns(std::map<fs::path, std::pair<int, int>> &) override;
    };

} // namespace mdn

#endif // !__MDN_HPP__
