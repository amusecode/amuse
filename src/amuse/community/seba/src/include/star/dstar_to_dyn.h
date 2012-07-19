/*
 *    star_to_dyn.h: transformation function from stellar evolution
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

#ifndef     _DSTAR_TO_DYN
#   define  _DSTAR_TO_DYN

#include  "starbase.h"
#include  "hdyn.h"

#include "stdinc.h"
#include "starlab_constants.h"

#include "double_star.h"
#include "sstar_to_dyn.h"

/*-----------------------------------------------------------------------------
 *  dstar_to_dyn  --  
 *-----------------------------------------------------------------------------
 */

void update_dynamical_part_of_binary(dyn*);
void update_evolution_part_of_binary(dyn*);
bool check_binary_consistency(dyn*);
void rearrange_the_binary_dynamics_tree(dyn*);
void rearrange_binaries_and_update_tree(dyn*);
//bool evolve_the_stellar_system(dyn*, real); 

double_state make_state(dyn* b);

void dstar_stats(dyn*, bool mass_function = true, 
		 vec center = 0, 
		 bool verbose = true);

void print_binary_dstars(dyn*);

void adddouble(hdyn * b, real time=0, binary_type=Detached, 
               bool random_initialization=false, 
               real a_min=1, real a_max=1.e+6,
               real e_min=0, real e_max=1) ;

void update_dyn_from_binary_evolution(dyn* the_binary_node,
				      dyn* d1=NULL, real dm1_fast=0,
				                    real dm2_fast=0);

#endif          // _DSTAR_TO_DYN




