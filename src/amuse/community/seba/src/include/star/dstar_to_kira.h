/*
 *    star_to_kira.h: transformation function from stellar evolution
 *		      package to dynamical package.
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

#ifndef     _DSTAR_TO_KIRA
#   define  _DSTAR_TO_KIRA

//#include  "starbase.h"
#include  "hdyn.h"

//#include "stdinc.h"
//#include "constants.h"

//#include "double_star.h"
#include "dstar_to_dyn.h"

#define MASS_UPDATE_LIMIT 0.0001
#define SMA_UPDATE_LIMIT  0.0001

/*-----------------------------------------------------------------------------
 *  dstar_to_hdyn  --  
 *-----------------------------------------------------------------------------
 */

bool create_or_delete_binary(hdyn* bi,
			     bool* update_dynamics,
			     int full_dump = 0);
bool binary_is_merged(dyn* bi);
bool binary_evolution(hdyn* b, int full_dump = 0);
void update_kepler_from_binary_evolution(hdyn* the_binary_node, 
					 real dm_fast);

#endif          // _DSTAR_TO_KIRA




