/*
 * -----------------------------------------------------------------
 * $Revision: 1.5 $
 * $Date: 2008/04/18 19:42:36 $
 * ----------------------------------------------------------------- 
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the header file for the CVODE dense linear solver, CVDENSE.
 * -----------------------------------------------------------------
 */

#ifndef _CVDENSE_H
#define _CVDENSE_H

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include "cvode_direct.h"
#include "sundials_dense.h"

/*
 * -----------------------------------------------------------------
 * Function: CVDense
 * -----------------------------------------------------------------
 * A call to the CVDense function links the main integrator with
 * the CVDENSE linear solver.
 *
 * cvode_mem is the pointer to the integrator memory returned by
 *           CVodeCreate.
 *
 * N is the size of the ODE system.
 *
 * The return value of CVDense is one of:
 *    CVDLS_SUCCESS   if successful
 *    CVDLS_MEM_NULL  if the cvode memory was NULL
 *    CVDLS_MEM_FAIL  if there was a memory allocation failure
 *    CVDLS_ILL_INPUT if a required vector operation is missing
 * -----------------------------------------------------------------
 */

SUNDIALS_EXPORT int CVDense(void *cvode_mem, int N);

#ifdef __cplusplus
}
#endif

#endif
