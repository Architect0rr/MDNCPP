#ifndef __MDN_NONROOT__
#define __MDN_NONROOT__

// #include <yas/serialize.hpp>
// #include <yas/std_types.hpp>
#include "../lib/yas/include/yas/mem_streams.hpp"
#include "../lib/yas/include/yas/binary_iarchive.hpp"
#include "../lib/yas/include/yas/binary_oarchive.hpp"

#include "mpi.h"

#include "utils.cpp"
#include "constants.cpp"
#include "barriers.cpp"
#include "logic.cpp"
#include "MDN.hpp"

namespace mdn{

    namespace fs = std::filesystem;
    RETURN_CODES MDN_nonroot::sanity(){
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

    RETURN_CODES MDN_nonroot::setup(){
        TRY
            sanity();
        CATCH_NOLOGGER("Error while checking sanity")

        SETUP_BARRIER(wcomm);

        TRY
            setup_logger(rank);
        CATCH("Error while setting up logger")

        return RETURN_CODES::OK;
    }

    uint64_t MDN_nonroot::dns(std::map<fs::path, std::pair<int, int>> &_distrib){
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
        std::unique_ptr<char> data(new char[count]);
        // char data[count];
        logger.debug("Broadcasting distribution");
        MPI_Bcast(data.get(), count, MPI_BYTE, cs::mpi_root, wcomm);
        logger.debug("Got distribution");

        yas::intrusive_buffer asd = yas::intrusive_buffer(data.get(), count);

        yas::mem_istream is(asd);
        yas::binary_iarchive<yas::mem_istream> ia(is);

        std::map<int, std::map<std::string, std::pair<int, int>>> distribution;
        std::map<int, std::map<fs::path, std::pair<int, int>>> _distribution;

        auto io = YAS_OBJECT("distribution", distribution);
        ia.serialize(io);
        logger.debug("Parsed distribution");

        logger.debug("Not dumping parsed distribution");

        ds2p(_distribution, distribution);

        uint64_t Natoms{};
        logger.debug("Bcasting Natoms");
        MPI_Bcast(&Natoms, 1, MPI_UINT64_T, cs::mpi_root, wcomm);
        logger.debug("Got Natoms");

        for (const auto &[k,v]: distribution.at(rank)){
            _distrib.emplace(k ,v);
        }
        logger.debug("Distribution complete");

        return Natoms;
    }

    RETURN_CODES MDN_nonroot::entry(){
        TRY
            setup();
        CATCH_NOLOGGER("Error wahile setting up")
        logger.info("Initialized");


        logger.debug("Getting distribution");
        std::map<fs::path, std::pair<int, int>> distribution;
        uint64_t summ = 0;
        uint64_t Natoms{};
        TRY
            Natoms = dns(distribution);
        CATCH("Error while getting distribution")
        for (const auto &[k,v]: distribution){
            logger.info("{}: {}({}) - {}", k.string(), v.first, v.second, v.first + v.second);
            summ += v.second;
        }
        logger.info("Total: {}", summ);

        fs::path ss = args.outfile.parent_path() / (args.outfile.filename().string() + "." + std::to_string(rank));
        logger.info("Output data will be written to {}", ss.string());

        uint64_t max_cluster_size = 0;
        logger.info("Starting calculations...");
        TRY
            run(ss, distribution, Natoms, max_cluster_size);
        CATCH("Error while doing caculations")
        logger.info("Calculations ended, sending result info");

        int *responces = NULL;
        MPI_Gather(&max_cluster_size, 1, MPI_UINT64_T, responces, 0, MPI_UINT64_T, cs::mpi_root, wcomm);

        const char * storpath = ss.string().c_str();
        MPI_Send(ss.string().c_str(), ss.string().length(), MPI_CHAR, cs::mpi_root, MPI_TAGS::STOR, wcomm);

        logger.info("Exiting entry point. NO RETURN");
        return RETURN_CODES::OK;
    }

} // namespace mdn

#endif // !__MDN_NONROOT__
