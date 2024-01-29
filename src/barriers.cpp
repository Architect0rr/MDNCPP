#ifndef __MDN_BARRIERS__
#define __MDN_BARRIERS__

#include "mpi.h"

void SANITY_BARRIER(MPI_Comm comm)
{
    MPI_Barrier(comm);
}

void SETUP_BARRIER(MPI_Comm comm)
{
    MPI_Barrier(comm);
}

void DISTRIBUTION_BARRIER(MPI_Comm comm)
{
    MPI_Barrier(comm);
}

#endif // !__MDN_BARRIERS__