/*
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~
									     */

/** \brief @file c_stdinc.h  C version of the standard include file */

/*  version 1:  Sep 1995   Steve McMillan
 *  version 2:
 *
 *  This file includes:
 *  1) new naming conventions to add to or replace existing names in C
 *  2) a string manipulation macro
 *  3) mathematical constants
 *  4) functions  abs()  min(,)  max(,)
 *  5) macros to cast angular arguments in standard form
 */

#ifndef  STARLAB_STDINC_H
#  define  STARLAB_STDINC_H

#include <stdio.h>
#include  <stdlib.h>
#include  <math.h>
#include  <string.h>

/*=============================================================================
**  Starlab version specification :
**=============================================================================
*/

#include <version.h>


/*=============================================================================
**  New naming conventions to add to or replace existing names in C :
**=============================================================================
*/

/*-----------------------------------------------------------------------------
 *  real  --  a more general name for the standard floating-point data type
 *-----------------------------------------------------------------------------
 */

typedef double  real;

/*-----------------------------------------------------------------------------
 *  bool  --  another name for int, to indicate use in logical operations
 *-----------------------------------------------------------------------------
 */

typedef int bool;

/* Convenient definitions: */

#define  false  0
#define  FALSE  0
#define  true   1
#define  TRUE   1

/*-----------------------------------------------------------------------------
 *  local  --  a more descriptive name for variables or functions which
 *             are invisible outside the file in which they are defined.
 *-----------------------------------------------------------------------------
 */

#define  local      static


/*=============================================================================
**  A  string manipulation macro :
**=============================================================================
*/
/*-----------------------------------------------------------------------------
 *  streq  --  a macro which returns 1 if two strings are equal, 0 otherwise
 *-----------------------------------------------------------------------------
 */

#define  streq(x,y)  (strcmp((x), (y)) == 0)


/*=============================================================================
**  Mathematical constants : 
**=============================================================================
*/

/*-----------------------------------------------------------------------------
 *  pi, etc.  --  mathematical constants, as well as "infinity"
 *-----------------------------------------------------------------------------
 */

#ifndef PI
#  define   PI         3.14159265358979323846
#endif
#define   TWO_PI     (2 * (PI))
#define   HALF_PI    (0.5 * (PI))
#define   ONE_THIRD  0.33333333333333333333
#define   ONE_SIXTH  0.16666666666666666667

#define VERY_LARGE_NUMBER 1e300


/*=============================================================================
**  Macros to cast angular arguments in standard form :
**=============================================================================
*/

/*-----------------------------------------------------------------------------
 *  pos_angle  --  recasts an angular variable into the range [0, TWO_PI)
 *  sym_angle  --  recasts an angular variable into the range [-PI, PI)
 *                   (recasting: transforming modulo 2 pi)
 *         example:
 *                 to map an angular variable 'phi' into the smallest positive
 *                 value, use
 *
 *                     phi = pos_angle(phi);
 *
 *                 to map an angular variable 'phi' into the smallest value,
 *                 positive or negative, use
 *
 *                     phi = sym_angle(phi);
 *
 *-----------------------------------------------------------------------------
 */

#define  pos_angle(phi)    ((phi) - TWO_PI * floor((phi)/TWO_PI ))
#define  sym_angle(phi)    ((phi) - TWO_PI * floor(((phi)+PI)/TWO_PI ))

#if defined USE_XREAL
#  include "xreal.h"
#else
   typedef real xreal;
#endif

#endif
