
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/// @file util_math.h  Extra math functions.
//
//  version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
//  version 2:
//
//  This file includes:
//  1) ....

#ifndef  STARLAB_UTIL_MATH_H
#  define  STARLAB_UTIL_MATH_H

#include  "stdinc.h"

bool twiddles(real a, real b, real eps = 1.e-12);

real adjust_number_to_power(real newstep, real max_step_size);
 
real asinh(real);  // inverse of sinh()
int sign(real);
real acosh(real);  // inverse of cosh()

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/util_math.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~

