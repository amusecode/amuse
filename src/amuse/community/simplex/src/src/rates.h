
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

/************************** Source spectra ***************************************/

//! Planck curve divided by frequency
double PlanckDivNu( double x, double T );

//! Planck curve
double Planck( double x, double T );

//! Integral of Planck curve
double PlanckInt( double _upperBound, double _T);


/************************** Cross sections ***************************************/

// The Planck curve and cross section have been scaled for stable numerics
// constants in Planck curve do not matter as they are divided out by normalisation
// Multiply any cross section with 1e-18 to obtain physical values


//! cross-section for H from Verner et al. 1996
double cross_HI( const double& _nu, const double& nu0HI );

//! Planck curve times cross section and divided by nu
double fHI( double x, double T );

//! Integral of Planck curve times cross section
double crossIntHI( double _upperBound, double _T );


/************************** Ionisations ***************************************/


//! collisional ionisation coefficient of hydrogen from Theuns et al. 1998
double coll_ion_coeff_HI( const double& _T );



/************************** Recombinations ************************************/

//! Recombination coefficient HII for case B from Hui & Gnedin 1997
double recomb_coeff_HII_caseB( const double& tempGas );

//! Recombination coefficient HII for case A from Hui & Gnedin 1997
double recomb_coeff_HII_caseA( const double& tempGas );


/************************** Heating *******************************************/


// The Planck curve and cross section have been scaled for stable numerics
// ants in Planck curve do not matter as they are divided out by normalisation 
// except for one factor of h
// Multiply any cross section with 1e-18 to obtain physical values

//! Planck curve times H cross section times h*(nu-nu_i) and divided by h*nu
//! Verner et al. 1996 cross section is used
double enerfHI( double x, double T );



/************************** Cooling ******************************************/


//! Recombination cooling coefficient HII for case B from Hui & Gnedin 1997
double recomb_cooling_coeff_HII_caseB( const double& tempGas );

//! Recombination cooling coefficient HII for case A from Hui & Gnedin 1997
double recomb_cooling_coeff_HII_caseA( const double& tempGas );

//! Collisional ionisation cooling coefficient of hydrogen from Theuns et al. 1998
double coll_ion_cooling_coeff_HI( const double& _T );

//! Collisional excitation cooling coefficient of hydrogen from Cen 1992
double coll_excit_cooling_coeff_HI( const double& _T );

//! Free-free cooling coefficient with gaunt factor from Theuns et al. 1998
double ff_cooling_coeff( const double& tempGas );


#endif
