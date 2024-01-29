#ifndef __MDN_NONROOT__
#define __MDN_NONROOT__

// #include <yas/serialize.hpp>
// #include <yas/std_types.hpp>
#include <yas/mem_streams.hpp>
#include <yas/binary_iarchive.hpp>
#include <yas/binary_oarchive.hpp>

#include "mpi.h"

#include "utils.cpp"
#include "constants.cpp"
#include "barriers.cpp"
#include "logic.cpp"

namespace mdn::nonroot
{
    namespace fs = std::filesystem;
    RETURN_CODES sanity(const utils::MC &sts)
    {
        int rnd;
        MPI_Bcast(&rnd, 1, MPI_INT, cs::mpi_root, MPI_COMM_WORLD);
        // sts.logger.info(std::string("Got number: " + std::to_string(rnd)));
        rnd += sts.rank;
        // sts.logger.info(std::string("Sending back: " + std::to_string(rnd)));
        int *responces = NULL;
        MPI_Gather(&rnd, 1, MPI_INT, responces, 0, MPI_INT, cs::mpi_root, MPI_COMM_WORLD);
        SANITY_BARRIER(MPI_COMM_WORLD);

        return RETURN_CODES::OK;
    }

    RETURN_CODES setup(utils::MC &sts)
    {
        if (sanity(sts) != RETURN_CODES::OK)
            return RETURN_CODES::ERROR;

        SETUP_BARRIER(MPI_COMM_WORLD);

        fs::path log_folder = sts.cwd / "logs";

        sts.setup_logger(log_folder, sts.rank);

        return RETURN_CODES::OK;
    }

    std::map<int, std::map<fs::path, std::pair<int, int>>> dns(utils::MC &sts)
    {
        // DISTRIBUTION_BARRIER(MPI_COMM_WORLD);
        // MPI_Status status;
        // int iop{};
        // std::size_t count{};
        // MPI_Iprobe(0, 123, MPI_COMM_WORLD, &iop, &status);
        // MPI_Get_count(&status, MPI_BYTE, &count);
        // char data[count];
        // MPI_Recv(data, count, MPI_BYTE, 0, 123, MPI_COMM_WORLD, &status);

        std::size_t count{};
        sts.logger.debug("Broadcasting distribution buffer size");
        MPI_Bcast(&count, 1, MPI_UNSIGNED_LONG, cs::mpi_root, MPI_COMM_WORLD);
        sts.logger.debug("Got distribution buffer size: {}", count);
        char data[count];
        sts.logger.debug("Broadcasting distribution");
        MPI_Bcast(&data, count, MPI_BYTE, cs::mpi_root, MPI_COMM_WORLD);
        sts.logger.debug("Got distribution");

        yas::intrusive_buffer asd = yas::intrusive_buffer(data, count);

        yas::mem_istream is(asd);
        yas::binary_iarchive<yas::mem_istream> ia(is);

        std::map<int, std::map<std::string, std::pair<int, int>>> distribution;
        std::map<int, std::map<fs::path, std::pair<int, int>>> _distribution;

        auto io = YAS_OBJECT("distribution", distribution);
        ia.serialize(io);
        sts.logger.debug("Parsed distribution");

        sts.logger.trace("Completed distribution:");
        for (const auto &[worker, stors] : distribution)
        {
            sts.logger./* data */ trace("Worker: {}", worker);
            for (const auto &[stor, bounds] : stors)
            {
                _distribution[worker].insert(std::make_pair(fs::path(stor), bounds));
                sts.logger.trace("Storage: '{}' [{}, {}]", stor, bounds.first, bounds.second);
            }
        }

        return _distribution;
    }

    RETURN_CODES entry(utils::MC &sts)
    {
        if (setup(sts) != RETURN_CODES::OK)
            return RETURN_CODES::ERROR;

        sts.logger.info("Initialized");

        sts.logger.debug("Getting distribution");
        std::map<int, std::map<fs::path, std::pair<int, int>>> _distribution = dns(sts);

        int Natoms{};
        MPI_Bcast(&Natoms, 1, MPI_INT, 0, MPI_COMM_WORLD);

        fs::path ss = sts.args.outfile / (sts.args.outfile.filename().string() + "." + std::to_string(sts.rank));

        logic::run(sts, ss, _distribution.at(sts.rank), Natoms);

        return RETURN_CODES::OK;
    }
}
#endif // !__MDN_NONROOT__