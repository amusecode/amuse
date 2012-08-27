#ifndef PARTICLE_H
#define PARTICLE_H 
#include "../copyright.h"
/*============================================================================*/
/*! \file particle.h
 *  \brief Global variables for all functions in in the src/particles
 *   directory.								      */
/*============================================================================*/

#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/*----------------------------------------------------------------------------*/
#ifdef PARTICLES

/*--------------------------- grid limit quantities -----------------------==-*/
/* left and right limit of grid indices */
int ilp,iup, jlp,jup, klp,kup;
/* left and right limit of grid boundary */
Real x1lpar, x1upar, x2lpar, x2upar, x3lpar, x3upar;

/*----------------- Quantities for Stopping time calculation -----------------*/
/*! \var Real *tstop0
 *  \brief Array of particle stopping time (for tstop=const) for each particle 
 *  type */
Real *tstop0;
/*! \var Real *grrhoa
 *  \brief an array of particle solid density times particle size in 
 *  normalized unit */
Real *grrhoa;
/*! \var Real alamcoeff
 *  \brief coefficient for the calculation of a/lambda */
Real alamcoeff;


/*! \var int ncell
 *  \brief number of neighbouring cells involved in 1D interpolation */
int ncell;

#ifdef SHEARING_BOX
/*! \var Real vshear
 *  \brief Shear velocity */
Real vshear;
#endif

#endif /* PARTICLES */

#endif /* PARTICLE_H */
