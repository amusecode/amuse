/*
 *    stardyn_util.h: transformation function from stellar evolution
 *		     package to dynamical package.
 *
 *.....................................................................
 *    version 1:  Aug 1996   Simon F. Portegies Zwart
 *    version 2:
 *...................................................................
 *     This file includes:
 *  1) Transformation routines from stellar evolution to n-body dynamcs.
 *
 *....................................................................
 */

#ifndef     _STARDYN_UTIL
#   define  _STARDYN_UTIL

#include  "starbase.h"
#include  "dyn.h"

#include "stdinc.h"
#include "starlab_constants.h"

#include "single_star.h"

//-----------------------------------------------------------------------------
//  stardyn_util  --  
//-----------------------------------------------------------------------------

// Copied from dyn.h...

void  compute_radii_percentiles(dyn *, bool);
void  compute_radii_quartiles(dyn *, bool);

#endif          // _STARDYN_UTIL





