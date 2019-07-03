#ifndef NOMPI
#include <mpi.h>
#endif

double mpi_wtime_()
{
#ifndef NOMPI
    return MPI_Wtime();
#endif
#ifdef NOMPI
    return 0;
#endif
}
