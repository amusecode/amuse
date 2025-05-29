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
#include <limits>

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
  RESTART, //!< Run was restarted from previous simulation
	SMC,    //!< SMC dust model
	LMC,    //!< LMC dust model
	MW,     //!< Milky Way dust model
	NO_DUST, //!< No dust
	NO_WEIGHTS,
	LOG_WEIGHTS,
	ENERGY_WEIGHTS,
	IONISATION_WEIGHTS,
};

const double alpha_A=4.18e-13;          //!< Case A recombination coefficient 
const double alpha_B=2.59e-13;          //!< Case B recombination coefficient 
const double alpha_1=1.59e-13;          //!< Directly to the ground state recombination coefficient 
const double recombCoefficient=0.7;     //!< Recombination coefficient power law coefficient 
const double A0=6.3e-18;                //!< v0 Photoionisation cross section 
const double h=6.6262e-34;              //!< Planck's constant 
const double k=1.380658e-16;            //!< Boltzmann's constant cgs
const double k_B=1.380658e-16;          //!< Boltzmann's constant cgs
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

const double m_H = 1.67372e-24;       //!< Mass of hydrogen atom [g]

const double nu0HI=3.28798e15;          //!< Hydrogen ionisation frequency 
const double crossFactor = 1.0e-18;// Common factor which transforms the cross sections into physical magnitude
const double sigma_dust_SMC = 1.0e-22; //SMC in cm^2
const double sigma_dust_LMC = 3.0e-22; //LMC in cm^2
const double sigma_dust_MW  = 5.0e-22; //MW in cm^2

//! Tells you whether the interaction is with a H atoms or not.
//! If yes, the n_e factor must not be used.
//! List of cooling curves H, He, C, C-H, N, O, O-H, Ne, Si, Si-H, Fe, Fe-H
const bool withH[12] = {0,0,0,1,0,0,1,0,0,1,0,1};


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
void quickSortPerm(vector<unsigned long int>& inVec, vector<unsigned long int>& permVec, const int& left, const int& right);

void polint(vector<double> xa, double *ya, const double x, double &y, double &dy);
double trapzd(double func(const double, const double), const double a, const double b, const double T, const int n);
double qromb(double func(const double, const double), double a, double b, double T);
double zbrent(double func( double,  double),  double x2,  double tol,  double Temp,  double perBin);

inline float SIGN(const float &a, const double &b);

double rec_rad_escape(const double& tau);


int site_compare_key(const void *a, const void *b);

unsigned int pow( const unsigned int& a, const int& b);
int pow( const int& a, const int& b);



#endif
