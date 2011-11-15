//Implementation of rates


#include "rates.h"
#include "Common.h"

using namespace std;



/************************** Source spectra ***************************************/

// Planck curve divided by nu
double PlanckDivNu( double x, double T ){
  return pow( x, 2.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );
}

//Planck curve
double Planck( double x, double T ){
  return pow( x, 3.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );
}

// Integral of Planck curve
double PlanckInt( double _upperBound, double _T){
  return qromb( Planck, 1.0, _upperBound, _T);
}

/************************** Cross sections ***************************************/

// cross-section for H from Verner et al. 1996
double cross_HI( const double& _nu, const double& nu0HI ){

  return (_nu < nu0HI) ? 0.0 : 1.16032e43 * pow( 9.62231e-15 * _nu - 1.0, 2.) * 
    pow( _nu, -4.0185) * pow( 1.0 + 1.7107e-8 * pow( _nu, 0.5), -2.963 ) ;

}

// hydrogen cross section from Verner et al. 1996 times black body curve
double fHI( double x, double T ){

  //  return ( x < 1.0) ? 0.0 : 6.3 * pow( x, -1.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );

  return ( x < 1.0) ? 0.0 : 
    1.16032e61 * 
    pow( 31.6379628 * x - 1.0, 2.) * 
    pow( nu0HI * x, -4.0185) * 
    pow( 1.0 + 1.7107e-8 * pow( nu0HI * x , 0.5), -2.963 ) * 
    pow( x, 2.0 ) / 
    ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );

}

// Integral of Planck curve times cross section
double crossIntHI( double _upperBound, double _T ){
  return qromb(fHI, 1.0, _upperBound, _T);
}

/************************** Ionisations ***************************************/

//collisional ionisation coefficient from Theuns et al. 1998
double coll_ion_coeff_HI( const double& _T ){

   return 1.17e-10 * sqrt( _T ) * exp( -1.578091e5 / _T ) * pow(1.0 + sqrt(_T/1e5), -1.0);

}


/************************** Recombinations ************************************/

// Recombination coefficient HII for case B from Hui & Gnedin 1997
double recomb_coeff_HII_caseB( const double& tempGas ) {  

  // If no temperature is given, 3K is assumed
  return (tempGas>0.0) ? 2.753e-14 * pow( 2.0 * 1.57807e5/tempGas, 1.5 ) / pow( ( 1.0 + pow( 0.729927007 * 1.57807e5/tempGas, 0.407 )), 2.242 ) : 5.97797e-11;

}

// Recombination coefficient HII for case A from Hui & Gnedin 1997
double recomb_coeff_HII_caseA( const double& tempGas ) {  

  return 1.269e-13 * pow( 2.0 * 1.57807e5/tempGas, 1.503 ) / pow( ( 1.0 + pow( 3.831417625 * 1.57807e5/tempGas, 0.47 )), 1.923 );

}

/************************** Heating *******************************************/

//Planck curve times H cross section times h*(nu-nu_0) and divided by h*nu
//Verner et al. 1996 cross section is used
double enerfHI( double x, double T ){

  return ( x < 1.0) ? 0.0 :
    ( x - 1.0) *
    1.16032e61 * 
    pow( 31.6379628 * x - 1.0, 2.0) * 
    pow( nu0HI * x, -4.0185) * 
    pow( 1.0 + 1.7107e-8 * pow( nu0HI * x , 0.5), -2.963 ) * 
    pow( x, 2.0 ) / 
    ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 ); 

}


/************************** Cooling ******************************************/

// Recombination cooling coefficient HII for case B from Hui & Gnedin 1997
double recomb_cooling_coeff_HII_caseB( const double& tempGas ) {  

  return 3.435e-30 * tempGas * pow( 2.0 * 1.57807e5/tempGas, 1.970 ) / 
    pow( ( 1.0 + pow( 0.88888888889 * 1.57807e5/tempGas, 0.376 )), 3.720 );

}

// Recombination cooling coefficient HII for case A from Hui & Gnedin 1997
double recomb_cooling_coeff_HII_caseA( const double& tempGas ) {  

  return 1.778e-29 * tempGas * pow( 2.0 * 1.57807e5/tempGas, 1.965 ) / 
    pow( ( 1.0 + pow( 3.696857671 * 1.57807e5/tempGas, 0.502 )), 2.697 );

}

// Collisional ionisation cooling coefficient from Theuns et al. 1998
double coll_ion_cooling_coeff_HI( const double& _T ){
  
  return 2.54e-21 * sqrt( _T ) * exp( -1.578091e5 / _T ) * pow(1.0 + sqrt(_T/1e5), -1.0);

}

// Collisional excitation cooling from Cen 1992
double coll_excit_cooling_coeff_HI( const double& _T ){
  
  return 7.5e-19 * pow(1.0 + sqrt(_T/1.e5), -1.0) * exp( -1.18348e5 / _T);

}

//free-free cooling coefficient with gaunt factor from Theuns et al. 1998
double ff_cooling_coeff( const double& tempGas ){  

  //gaunt factor  
  double g_ff = 1.1 + 0.34 * exp( -0.33333333333 * pow(5.5-log(tempGas)/log(10), 2.0) );

  return 1.42e-27 * g_ff * sqrt( tempGas );// This must be multiplied with (n_HII + n_HeII + 4 n_HeIII) * n_e

}
