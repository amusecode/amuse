
//this header contains the rates and cross sections needed by SimpleX.cpp

#ifndef RATES_H
#define RATES_H

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <float.h>

using namespace std;

// Recombination coefficient HII for case B from Hui & Gnedin 1997
double recombCoeffHIICaseB(  double& tempGas ) {  
  return 2.753e-14 * pow( 2.0 * 1.57807e5/tempGas, 1.5 ) / pow( ( 1.0 + pow( 0.729927007 * 1.57807e5/tempGas, 0.407 )), 2.242 );
}

// Recombination cooling coefficient HII for case B from Hui & Gnedin 1997
double recombCoolingHIICaseB(  double& tempGas ) {  
  return 3.435e-30 * tempGas * pow( 2.0 * 1.57807e5/tempGas, 1.970 ) / pow( ( 1.0 + pow( 0.88888888889 * 1.57807e5/tempGas, 0.376 )), 3.720 );// Hui & Gnedin 1997  
}



#endif
