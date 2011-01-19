/*************************************************************************
file:         Main_loop.h
author:       Jan-Pieter Paardekooper
mail:         jppaarde@strw.leidenuniv.nl
version:      0.1
last change:  04.04.2008
---------------------------------------------------------------------
description:
This file contains the header file for the main routine
**************************************************************************/
/*
 * Date: Name
 * Put additional comments here
 *
 * 03.04.08 Jan-Pieter Paardekooper
 * Put comments suitable for doxygen in the relevant places
 *
*/

/***** TO DO *****
 *****************/

#ifndef MAINLOOP_H
#define MAINLOOP_H

#include "Common.h"
#include "SimpleX.h"
#include "mpi.h"

using namespace std;

int main_loop(int argc, char** argv,SimpleX* SimpleXGrid);


#endif
