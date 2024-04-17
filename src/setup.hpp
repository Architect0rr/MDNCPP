#ifndef __MDN_SETUP__
#define __MDN_SETUP__

#include "nlohmann/json.hpp"

#include "mpi.h"

#include "utils.hpp"
#include "barriers.hpp"
#include "MDN.hpp"

namespace mdn {

    using json = nlohmann::json;
    namespace fs = std::filesystem;

    void MDN_root::sanity(){
        int rnd = std::rand();
        std::cout << "  ├─Bcasting: " << rnd << std::endl;
        MPI_Bcast(&rnd, 1, MPI_INT, cs::mpi_root, wcomm);
        int *responces = new int[size];
        std::cout << "  ├─Gathering..." << std::endl;
        MPI_Gather(&rnd, 1, MPI_INT, responces, 1, MPI_INT, cs::mpi_root, wcomm);

        if (responces[0] != rnd){
            std::cerr << "  └─Root sanity doesn't passed" << std::endl;
            throw std::runtime_error("Root sanity doesn't passed");
        }

        // logger.info(std::string("Random int: " + std::to_string(rnd)));
        // logger.info(std::string("Received arrray:"));
        std::cout << "  ├─Received array: ";
        for (int i = 0; i < size; ++i)
        {
            std::cout << responces[i] << " ";
            if (responces[i] - i != rnd)
            {
                std::cerr << "  └─Root sanity doesn't passed" << std::endl;
                throw std::runtime_error("Root sanity doesn't passed");
            }
        }
        std::cout << "  ├─Root sanity OK" << std::endl;
        std::cout << "  ├─Releasing sanity barrier" << std::endl;
        SANITY_BARRIER(wcomm);
        std::cout << "  └─Sanity barrier released" << std::endl;

        delete[] responces;
    }

    std::map<int, std::string> MDN_root::gather_storages(){
        std::map<int, std::string> storages;
        std::set<int> ncw;
        for (int i = 1; i < size; ++i) ncw.emplace(i);
        MPI_Status status;
        int flag = 0;
        int count{};
        while (!ncw.empty()){
            for (const int i : ncw){
                MPI_Iprobe(i, MPI_TAGS::STOR, wcomm, &flag, &status);
                if (flag){
                    MPI_Get_count(&status, MPI_CHAR, &count);
                    char* buf = new char[count];
                    MPI_Recv(buf, count, MPI_CHAR, i, MPI_TAGS::STOR, wcomm, &status);
                    storages.emplace(i, std::string(buf, count));
                    delete[] buf;
                    flag = 0;
                    status = MPI_Status();
                    ncw.erase(i);
                }
            }
        }
        return storages;
    }

    void MDN_root::setup(){
        std::cout << "├─Sanity:" << std::endl;
        TRY
            sanity();
        CATCH_NOLOGGER("Error while checking sanity")

        std::error_code ec;
        fs::path log_folder = cwd / folders::log;
        if (!fs::exists(log_folder)){
            if (!fs::create_directories(log_folder, ec)){
                std::cerr << "├─Cannot create non-existent directories: " << log_folder.string() << std::endl;
                std::cerr << "├─Error code: " << ec.value() << std::endl;
                std::cerr << "└─Message: " << ec.message() << std::endl;
                throw std::runtime_error("Cannot create non-existent directory: " + log_folder.string());
            }else{
                std::cout << "├─Created logging directory: " << log_folder << std::endl;
            }
        }else{
            std::cout << "├─Logging folder exists" << std::endl;
        }

        std::cout << "├─Releasing setup barrier" << std::endl;
        SETUP_BARRIER(wcomm);
        std::cout << "├─Setup barrier released" << std::endl;

        TRY
            setup_logger(rank);
        CATCH("Error while setting up logger")

        std::cout << "└─Root setup OK" << std::endl;
    }

    void MDN::sanity(){
        int rnd;
        MPI_Bcast(&rnd, 1, MPI_INT, cs::mpi_root, wcomm);
        // sts.logger.info(std::string("Got number: " + std::to_string(rnd)));
        rnd += rank;
        // sts.logger.info(std::string("Sending back: " + std::to_string(rnd)));
        int *responces = NULL;
        MPI_Gather(&rnd, 1, MPI_INT, responces, 0, MPI_INT, cs::mpi_root, wcomm);
        SANITY_BARRIER(wcomm);
    }

    void MDN::setup(){
        TRY
            sanity();
        CATCH_NOLOGGER("Error while checking sanity")

        SETUP_BARRIER(wcomm);

        TRY
            setup_logger(rank);
        CATCH("Error while setting up logger")
    }

} // namespace mdn

#endif // !__MDN_SETUP__