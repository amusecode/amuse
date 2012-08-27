#include "copyright.h"
/*============================================================================*/
/*! \file baton.c
 *  \brief Stage independent processes by passing a baton.
 *
 * PURPOSE: Stage independent processes by passing a baton.  The
 * functions here assume that no MPI communication will proceed
 * between baton_start() and baton_stop() and that the arguments to
 * the functions will be identical.  This is easy enough to
 * accomplish, but can be error prone.  Use these functions with
 * caution!
 *
 * These functions were written and added to athena by T. A. Gardiner
 * in March, 2008.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 * - baton_start()
 * - baton_stop()							      */
/*============================================================================*/
#include "defs.h"
#include "prototypes.h"


#ifdef MPI_PARALLEL

/* ========================================================================== */

/*! \fn void baton_start(const int Nb, const int tag)
 *  \brief This function starts off a "Baton Passing" scheme for regulating
 * the maximum number of processes active at a given moment. 
 *
 * NOTES on baton_start() and baton_stop():
 *
 * The values for Nb and tag MUST be the same in the call to
 * baton_start() and in baton_stop().
 *
 * - Nb  -> Number of Batons = max number of ranks in action at a given time.
 * - tag -> MPI message tag to use for this baton passing.
 *
 * This function starts off a "Baton Passing" scheme for regulating
 * the maximum number of processes active at a given moment.  It works
 * by handing out Nb batons to the first (0 -> Nb-1) ranks.  Then
 * ranks (Nb -> 2Nb-1) wait for a signal from the ranks (0 -> Nb-1)
 * which is sent by the baton_stop() function and so the process
 * continues until all of the MPI ranks have completed. */
void baton_start(const int Nb, const int tag){

  MPI_Status stat;
  int status, isig, fr_id, my_id;

  status = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  if(status != MPI_SUCCESS)
    ath_error("[baton_start]: MPI_Comm_rank error = %d\n",status);

  fr_id = my_id - Nb;

  if(Nb > 0 && fr_id >= 0){
    /* Wait for another process to signal before beginning. */
    status = MPI_Recv(&isig, 0, MPI_INT, fr_id, tag, MPI_COMM_WORLD, &stat);

    if(status != MPI_SUCCESS)
      ath_error("[baton_start]: MPI_Recv error = %d\n",status);
  }

  return;
}


/* ========================================================================== */


/*! \fn void baton_stop(const int Nb, const int tag)
 *  \brief This function marks the end of a section of code which is regulated
 * by "Baton Passing".  It signals the next process to start, i.e. it
 * hands off the Baton. */
void baton_stop(const int Nb, const int tag){

  int status, isig, to_id, my_id, nproc;

  status = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
  if(status != MPI_SUCCESS)
    ath_error("[baton_stop]: MPI_Comm_rank error = %d\n",status);

  status = MPI_Comm_size(MPI_COMM_WORLD, &nproc);
  if(status != MPI_SUCCESS)
    ath_error("[baton_stop]: MPI_Comm_size error = %d\n",status);

  to_id = my_id + Nb;

  if(Nb > 0 && to_id < nproc){
    /* Signal the next process to begin. */
    status = MPI_Send(&isig, 0, MPI_INT, to_id, tag, MPI_COMM_WORLD);

    if(status != MPI_SUCCESS)
      ath_error("[baton_stop]: MPI_Send error = %d\n",status);
  }

  return;
}


/* ========================================================================== */


#endif /* MPI_PARALLEL */
