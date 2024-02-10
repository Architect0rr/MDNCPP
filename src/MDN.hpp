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

#endif // !__MDN_TRY_CATCH_MACRO__

#include <filesystem>
#include <map>

#include "utils.hpp"
#include "adios2.h"
#include "nlohmann/json.hpp"

#include "spdlog/spdlog.h"

#include "config.hpp"

// namespace mdn {

namespace fs = std::filesystem;
using json = nlohmann::json;

struct Arguments
{
    bool debug = false;
    bool trace = false;
    int mode = 0;
    bool verbose = false;
    fs::path outfile = "/";
};

class MDN
{
public:
    MDN() = delete;
    MDN(const int rank = 0, const int size = 0, MPI_Comm comm = MPI_COMM_WORLD, const char *prefix = "unknown") noexcept : rank(rank), size(size), wcomm(comm), prefix(prefix) {}
    ~MDN() = default;

    RETURN_CODES parse_args(const int argc, const char ***argv);
    RETURN_CODES start(const int argc, const char ***argv);

protected:
    spdlog::logger logger = spdlog::logger("default");
    fs::path cwd = "/";
    Arguments args;
    const int rank = 0;
    const int size = 0;
    const MPI_Comm wcomm = MPI_COMM_WORLD;
    std::string_view prefix;

    virtual RETURN_CODES entry() = 0;
    virtual RETURN_CODES sanity() = 0;
    virtual RETURN_CODES setup() = 0;
    virtual uint64_t dns(std::map<fs::path, std::pair<int, int>> *) = 0;
    void setup_logger(const fs::path &, const int);
    RETURN_CODES run(const fs::path &, const std::map<fs::path, std::pair<int, int>> &, const uint64_t);
};

class MDN_root : public MDN
{
public:
    using MDN::MDN;

private:
    RETURN_CODES entry() override;
    RETURN_CODES sanity() override;
    RETURN_CODES setup() override;
    uint64_t dns(std::map<fs::path, std::pair<int, int>> *) override;

    const uint64_t bearbeit(const fs::path &, double *);
    const std::map<fs::path, int> storage_rsolve(const json &);
    const int get_total_steps(adios2::IO &, const fs::path &);
    RETURN_CODES distribute(const std::map<fs::path, int> &, std::map<fs::path, std::pair<int, int>> *);
    std::map<int, std::map<fs::path, std::pair<int, int>>> make_distribution(const std::map<fs::path, int> &);
};

class MDN_nonroot : public MDN
{
public:
    using MDN::MDN;

private:
    RETURN_CODES entry() override;
    RETURN_CODES sanity() override;
    RETURN_CODES setup() override;
    uint64_t dns(std::map<fs::path, std::pair<int, int>> *) override;
};

// } // namespace mdn

#endif // !__MDN_HPP__
