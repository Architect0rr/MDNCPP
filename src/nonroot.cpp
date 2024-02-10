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
#include "MDN.hpp"

// namespace mdn::nonroot{
    namespace fs = std::filesystem;
    RETURN_CODES MDN_nonroot::sanity()
    {
        int rnd;
        MPI_Bcast(&rnd, 1, MPI_INT, cs::mpi_root, wcomm);
        // sts.logger.info(std::string("Got number: " + std::to_string(rnd)));
        rnd += rank;
        // sts.logger.info(std::string("Sending back: " + std::to_string(rnd)));
        int *responces = NULL;
        MPI_Gather(&rnd, 1, MPI_INT, responces, 0, MPI_INT, cs::mpi_root, wcomm);
        SANITY_BARRIER(wcomm);

        return RETURN_CODES::OK;
    }

    RETURN_CODES MDN_nonroot::setup()
    {
        TRY
            sanity();
        CATCH_NOLOGGER("Error while checking sanity")

        SETUP_BARRIER(wcomm);

        fs::path log_folder = cwd / "logs";

        TRY
            setup_logger(log_folder, rank);
        CATCH("Error while setting up logger")

        return RETURN_CODES::OK;
    }

    uint64_t MDN_nonroot::dns(std::map<fs::path, std::pair<int, int>> *_distrib)
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
        logger.debug("Broadcasting distribution buffer size");
        MPI_Bcast(&count, 1, MPI_UNSIGNED_LONG, cs::mpi_root, wcomm);
        logger.debug("Got distribution buffer size: {}", count);
        char data[count];
        logger.debug("Broadcasting distribution");
        MPI_Bcast(&data, count, MPI_BYTE, cs::mpi_root, wcomm);
        logger.debug("Got distribution");

        yas::intrusive_buffer asd = yas::intrusive_buffer(data, count);

        yas::mem_istream is(asd);
        yas::binary_iarchive<yas::mem_istream> ia(is);

        std::map<int, std::map<std::string, std::pair<int, int>>> distribution;
        std::map<int, std::map<fs::path, std::pair<int, int>>> _distribution;

        auto io = YAS_OBJECT("distribution", distribution);
        ia.serialize(io);
        logger.debug("Parsed distribution");

        logger.trace("Completed distribution:");
        for (const auto &[worker, stors] : distribution)
        {
            logger./* data */ trace("Worker: {}", worker);
            for (const auto &[stor, bounds] : stors)
            {
                _distribution[worker].insert(std::make_pair(fs::path(stor), bounds));
                logger.trace("Storage: '{}' [{}, {}]", stor, bounds.first, bounds.second);
            }
        }

        uint64_t Natoms{};
        MPI_Bcast(&Natoms, 1, MPI_UINT64_T, 0, wcomm);

        *_distrib = _distribution.at(rank);
        return Natoms;
    }

    RETURN_CODES MDN_nonroot::entry()
    {
        TRY
            setup();
        CATCH_NOLOGGER("Error wahile setting up")
        logger.info("Initialized");


        logger.debug("Getting distribution");
        std::map<fs::path, std::pair<int, int>> *distribution = {};
        uint64_t Natoms{};
        TRY
            Natoms = dns(distribution);
        CATCH("Error while getting distribution")

        fs::path ss = args.outfile / (args.outfile.filename().string() + "." + std::to_string(rank));

        logger.info("Starting calculations...");
        TRY
            run(ss, *distribution, Natoms);
        CATCH("Error while doing caculations")

        return RETURN_CODES::OK;
    }
// } // namespace
#endif // !__MDN_NONROOT__