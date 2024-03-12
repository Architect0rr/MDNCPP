// STL libraries
#include <iostream>
#include <filesystem>

// Third-party librariesq
#include "argparse/argparse.hpp"
#include "spdlog/spdlog.h"

// Convenient libraries
#include "mpi.h"

// This project libraries
#include "utils.cpp"
#include "MDN.hpp"
#include "utils.hpp"
#include "root.cpp"
#include "nonroot.cpp"

namespace mdn{

    RETURN_CODES MDN::start(const int argc, const char ***argv){
        TRY
            parse_args(argc, argv);
        CATCH_NOLOGGER("Error while parsing args")
        TRY
        entry();
        CATCH_NOLOGGER("Error in entry point")
        logger.info("Worker shutdown");
        spdlog::shutdown();
        return RETURN_CODES::OK;
    }

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (size >= 2)
    {
        try
        {
            if (rank == mdn::cs::mpi_root)
            {
                mdn::MDN_root mmd(rank, size, MPI_COMM_WORLD, "root");
                mmd.start(argc, const_cast<const char***>(&argv));
            }
            else
            {
                mdn::MDN_nonroot mmd(rank, size, MPI_COMM_WORLD, "nonroot");
                mmd.start(argc, const_cast<const char***>(&argv));
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << "Some low-level std::exception. Aborting..." << std::endl;
            std::cerr << e.what() << std::endl;
            MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
        }
        catch (...)
        {
            std::cerr << "Some low-level unknown exception. Aborting..." << std::endl;
            MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
        }
    }
    else
    {
        std::cerr << "It is required at least 2 processes to be runned, shutting down" << std::endl;
        MPI_Abort(MPI_COMM_WORLD, mdn::RETURN_CODES::ERROR);
    }

    MPI_Finalize();
    return 0;
}
