/*************************************************************************
file:         Common.h
author:       Jan-Pieter Paardekooper
mail:         jppaarde@strw.leidenuniv.nl
version:      0.1
last change:  03.04.2008
---------------------------------------------------------------------
description:
This file contains useful information needed globally for the simulation,
like important constants, standard functions etc. It also contains the
createOptions array containing the simulation parameters
**************************************************************************/
/*
 * Date: Name
 * Put additional comments here
 *
*/

/***** To Do *******
 *
 *
 *******************/


#ifndef COMMON_H
#define COMMON_H

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

/// Maximum number of characters in a character array.
const int MAX_LINE = 256;		

const double ONE_HALF   = 0.50;
const double ONE_THIRD  = 0.333333333333333333333333333;
const double ONE_FOURTH = 0.25;

//! An enum for several initialisation options
enum    choices	{
	POISSON, 	    //!< Generate Poisson distribution 
	AUTOMATIC, 	    //!< Generate distribution automatically
        READ,    //!< Read generated distribution from file (structured grid)
	LONGCHAR,	    //!< Use long characteristic method  
	SHORTCHAR, 	    //!< Use short characteristic method
	IONIZED, 	    //!< Use ionisation routine
	REGULAR, 	    //!< Use normal radiation routine
	SMC,    //!< SMC dust model
	LMC,    //!< LMC dust model
	MW,     //!< Milky Way dust model
	NO_DUST, //!< No dust
};

const double alpha_A=4.18e-13;          //!< Case A recombination coefficient 
const double alpha_B=2.59e-13;          //!< Case B recombination coefficient 
const double alpha_1=1.59e-13;          //!< Directly to the ground state recombination coefficient 
const double recombCoefficient=0.7;     //!< Recombination coefficient power law coefficient 
const double A0=6.3e-18;                //!< v0 Photoionisation cross section 
const double h=6.6262e-34;              //!< Planck's constant 
const double k=1.380658e-16;              //!< Boltzmann's constant cgs
const double v0=3.2e15;                 //!< Hydrogen ionisation frequency 
const double tempIonGas=1e4;            //!< Assumed temperature for ionised gas 
const double parsecToCm=3.08568025e18;  //!< Number of cm in one parsec 
const double secondsPerMyr=3.1556926e13;//!< Number of seconds in 1 Myr 
const double grPerHParticle=1.673534e-24;//!< Number of grams per hydrogen particle 
 
const double planckConstant=6.6261e-27; //!< Planck's constant [erg s] 
const double speed_of_light=2.9979e10;    //!< Speed of light[cm/s] 
const double thomsonCross=6.6525e-25;  //!< Thomson cross section [cm^2] 
const double megaParsecToCm=3.086e24; //!< Number of centimeters in one Megaparsec [cm] 
const double yrInSec=3.1558e7;        //!< Number of seconds in one year [s] 
const double massSunInGr=1.989e33;    //!< Mass of the Sun [g] 
const double massProton=1.6726e-24;   //!< Mass of a proton [g] (proton mass) 
const double H_0=1e7/megaParsecToCm;  //!< Hubble's constant [s^{-1}] (100 km/s/Mpc) 

const double nu0HI=3.28798e15;          //!< Hydrogen ionisation frequency 
const double crossFactor = 1.0e-18;// Common factor which transforms the cross sections into physical magnitude
const double sigma_dust_SMC = 1.0e-22; //SMC in cm^2
const double sigma_dust_LMC = 3.0e-22; //LMC in cm^2
const double sigma_dust_MW  = 5.0e-22; //MW in cm^2


const int BBevalNum=10;									//!< Number of temperatures at which the blackbody integrals have been evaluated 
 
const double BBevalTemps[10]=						//!< The BBevalNum temperatues at which the blackbody integrals have been evaluated 
{ 
  1e4, 2e4, 3e4, 4e4, 5e4, 6e4, 7e4, 8e4, 9e4, 1e5 
};

// Chaels waardes
const double BBevalAbove[10]=        //!< The fraction of the Planck spectrum above the Lyman limit at a certain temperature
{
  0.000103496,                // T=1e4K
  0.0422316,                  // T=2e4K
  0.213447,                   // T=3e4K
  0.413349,                   // T=4e4K
  0.57335,                    // T=5e4K
  0.687975,                   // T=6e4K
  0.768041,                   // T=7e4K
  0.824229,                   // T=8e4K
  0.864264,                   // T=9e4K
  0.893322                    // T=1e5K
};

const double BBevalEffAbs[10]=        //!< The effective (weighted by a normalized Planck spectrum) absorption coefficient at a certain temperature
{
  0.831343,                        // T=1e4K
  0.697735,                        // T=2e4K
  0.590817,                        // T=3e4K
  0.505007,                        // T=4e4K
  0.435861,                        // T=5e4K
  0.379716,                        // T=6e4K
  0.333684,                        // T=7e4K
  0.295553,                        // T=8e4K
  0.263646,                        // T=9e4K
  0.236691                         // T=1e5K
};

double rec_rad_escape(const double& tau);



//! Function determining the inproduct of two vectors in arbitrary dimension \f$d\ge 1\f$
/**
  \param v1 Vector number 1
  \param v2 Vector number 2
  \param dimension Dimension of vector number 1 and 2
  \return \f$d\f$-dimensional inproduct of vectors 1 and 2
*/
double inproduct(double v1[], double v2[], int dimension);
float inproduct(float v1[], float v2[], int dimension);


//! Compute the determinant of a 4x4 array (used to check if point is within tetrahedron)
/**
  \param a[] 1D array containing all 16 entries of the matrix
*/
double det_4by4 (double a[]);


//! Routine to compute the circumcenter of a tetrahedron.
/**
  \param xTet The coordinates of the four vertices of a tetrahedron
  \param xCC The x-coordinate of the circumcenter
  \param yCC The y-coordinate of the circumcenter
  \param zCC The z-coordinate of the circumcenter
*/
void CalcCircumCenter(double xTet[4][3], double& xCC, double& yCC, double& zCC);

void quickSortPerm(vector<float>& inVec, vector<int>& permVec, const int& left, const int& right);
void quickSortPerm(vector<unsigned long long int>& inVec, vector<unsigned long long int>& permVec, const int& left, const int& right);



void polint(vector<double> xa, double *ya, const double x, double &y, double &dy);
double trapzd(double func(const double, const double), const double a, const double b, const double T, const int n);
double qromb(double func(const double, const double), double a, double b, double T);

//hydrogen cross section times Planck curve
double fHI( double x, double T );
//dust cross section times Planck curve for SMC
double f_dust_SMC( double x, double T );
//dust cross section times PLanck curve for LMC
double f_dust_LMC( double x, double T );
//dust cross section times PLanck curve for Milky Way
double f_dust_MW( double x, double T ); 
//Planck curve divided by frequency
double PlanckDivNu( double x, double T );
//Planck curve
double Planck( double x, double T );
// Fitting function for dust cross sections
double fitting_function( const double& x, const double& a, const double& b, const double& p, const double& q );
// Dust extinction
double ksi_dust( double x, short int dust_model );

int site_compare_key(const void *a, const void *b);

unsigned int pow( const unsigned int& a, const int& b);
int pow( const int& a, const int& b);

#endif
