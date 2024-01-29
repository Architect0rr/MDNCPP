#ifndef __MDN_ADIOS2WRITER__
#define __MDN_ADIOS2WRITER__
#include "adios2.h"

#include "mpi.h"

#include <filesystem>

namespace mdn
{
    namespace fs = std::filesystem;
    struct adser
    {
        adios2::ADIOS adios = adios2::ADIOS(MPI_COMM_WORLD);
        adios2::IO io = adios.DeclareIO("PReader");
        adios2::Engine reader;
        bool opened = false;
        adser()
        {
        }

        bool open(const fs::path &storage)
        {
            if (opened)
            {
                return false;
            }
            reader = io.Open(storage, adios2::Mode::Read);
            opened = true;
            return true;
        }
    };
}

#endif // !__MDN_ADIOS2WRITER__