/*************************************************************************
file:         Common.cpp
author:       Jan-Pieter Paardekooper
mail:         jppaarde@strw.leidenuniv.nl
version:      0.1
last change:  03.04.2008
---------------------------------------------------------------------
description:
This file contains the read in routine and often used functions like error
messages and random number generators
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


#include "Common.h"
#include "Structs.h"

using namespace std;

double rec_rad_escape(const double& tau){

  //  double tmp1 = 1.0 - 2.*pow(tau, 2);
  //double tmp2 = 1.0 + 2.*tau;
  //double tmp3 = exp(-2.*tau);
  //double tmp4 = 3./(8.*tau);

  //double f_esc = tmp4*( tmp3*tmp2 - tmp1 );
 
  //double f_esc = (tau > 0.0 ) ? 3.*( exp(-tau) * (1.0 + tau) + 2.*pow(0.5*tau, 2) - 1.0 )/(8.*pow(0.5*tau,3)) : 1.0;

  //for very small tau solution is unstable
  double f_esc = (tau > 1.e-3) ? 3.*( exp(-tau) * (1.0 + tau) + 2.*pow(0.5*tau, 2.) - 1.0 )/(8.*pow(0.5*tau,3.)) : 1.0;


  // if(f_esc < 0.0){
  //   f_esc = 1.0;
  // }

  return f_esc;

}

double inproduct(double v1[], double v2[], int dimension) {
						// Omgeschreven tot generiek d dimensionaal inproduct...
  double v1temp[dimension];
  double v2temp[dimension];

  double length1 = 0.0, length2=0.0, result=0.0;
  
  for(int i=0; i<dimension; i++) {
    length1+=pow(v1[i],2.);
    length2+=pow(v2[i],2.);
  }
  double one_over_length1 = 1.0/( sqrt(length1) );
  double one_over_length2 = 1.0/( sqrt(length2) );
  
  for(int i=0; i<dimension; i++) {
    v1temp[i]=v1[i]*one_over_length1;
    v2temp[i]=v2[i]*one_over_length2;
    result+=v1temp[i]*v2temp[i];
  }

  return result;

}

// overloaded for floats
float inproduct(float v1[], float v2[], int dimension) {
						// Omgeschreven tot generiek d dimensionaal inproduct...
  float v1temp[dimension]; 
  float v2temp[dimension];

  float length1 = 0.0, length2=0.0, result=0.0;
  
  for(int i=0; i<dimension; i++) {
    length1+=pow(v1[i],2.);
    length2+=pow(v2[i],2.);
  }
  
  float one_over_length1 = 1.0/( sqrt(length1) );
  float one_over_length2 = 1.0/( sqrt(length2) );
  
  for(int i=0; i<dimension; i++) {
    v1temp[i]=v1[i]*one_over_length1;
    v2temp[i]=v2[i]*one_over_length2;
    result+=v1temp[i]*v2temp[i];
  }
  
  return result;

}

//calculate the determinant of a 4x4 matrix
double det_4by4 (double a[]) {

  double det;

  det =
      a[0+0*4] * (
          a[1+1*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        + a[1+3*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] ) )
    - a[0+1*4] * (
          a[1+0*4] * ( a[2+2*4] * a[3+3*4] - a[2+3*4] * a[3+2*4] )
        - a[1+2*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] ) )
    + a[0+2*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+3*4] - a[2+3*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+3*4] - a[2+3*4] * a[3+0*4] )
        + a[1+3*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) )
    - a[0+3*4] * (
          a[1+0*4] * ( a[2+1*4] * a[3+2*4] - a[2+2*4] * a[3+1*4] )
        - a[1+1*4] * ( a[2+0*4] * a[3+2*4] - a[2+2*4] * a[3+0*4] )
        + a[1+2*4] * ( a[2+0*4] * a[3+1*4] - a[2+1*4] * a[3+0*4] ) );

  return det;

}

//Calculate the coordinates of the circumcentre of a simplex
void CalcCircumCenter(double xTet[4][3], double& xCC, double& yCC, double& zCC) {

  double det0[16]={xTet[0][0], xTet[0][1], xTet[0][2], 1,
                   xTet[1][0], xTet[1][1], xTet[1][2], 1,
                   xTet[2][0], xTet[2][1], xTet[2][2], 1,
                   xTet[3][0], xTet[3][1], xTet[3][2], 1};

  double detx[16]={pow(xTet[0][0],2.)+pow(xTet[0][1],2.)+pow(xTet[0][2],2.), xTet[0][1], xTet[0][2], 1,
                   pow(xTet[1][0],2.)+pow(xTet[1][1],2.)+pow(xTet[1][2],2.), xTet[1][1], xTet[1][2], 1,
                   pow(xTet[2][0],2.)+pow(xTet[2][1],2.)+pow(xTet[2][2],2.), xTet[2][1], xTet[2][2], 1,
                   pow(xTet[3][0],2.)+pow(xTet[3][1],2.)+pow(xTet[3][2],2.), xTet[3][1], xTet[3][2], 1};

  double dety[16]={pow(xTet[0][0],2.)+pow(xTet[0][1],2.)+pow(xTet[0][2],2.), xTet[0][0], xTet[0][2], 1,
                   pow(xTet[1][0],2.)+pow(xTet[1][1],2.)+pow(xTet[1][2],2.), xTet[1][0], xTet[1][2], 1,
                   pow(xTet[2][0],2.)+pow(xTet[2][1],2.)+pow(xTet[2][2],2.), xTet[2][0], xTet[2][2], 1,
                   pow(xTet[3][0],2.)+pow(xTet[3][1],2.)+pow(xTet[3][2],2.), xTet[3][0], xTet[3][2], 1};

  double detz[16]={pow(xTet[0][0],2.)+pow(xTet[0][1],2.)+pow(xTet[0][2],2.), xTet[0][0], xTet[0][1], 1,
                   pow(xTet[1][0],2.)+pow(xTet[1][1],2.)+pow(xTet[1][2],2.), xTet[1][0], xTet[1][1], 1,
                   pow(xTet[2][0],2.)+pow(xTet[2][1],2.)+pow(xTet[2][2],2.), xTet[2][0], xTet[2][1], 1,
                   pow(xTet[3][0],2.)+pow(xTet[3][1],2.)+pow(xTet[3][2],2.), xTet[3][0], xTet[3][1], 1};


  double a=0.5/det_4by4(det0);
  double Dx=det_4by4(detx);
  double Dy=-det_4by4(dety);
  double Dz=det_4by4(detz);

  xCC=Dx*a;
  yCC=Dy*a;
  zCC=Dz*a;

}

/*******************************************************************/
void quickSortPerm(vector<float>& inVec, vector<int>& permVec, const int& left, const int& right) {

  int i = left, j = right;
  int intBuff;
  float buff;
  float pivot = inVec[(left+right)/2];
  while(i <= j){
    while (inVec[i]>pivot) i++;
    while (inVec[j]<pivot) j--;
    if (i<=j){
      // swap vector elements
      buff = inVec[i];
      inVec[i] = inVec[j];
      inVec[j] = buff;

      // swap elements of other vector in same manner
      intBuff = permVec[i];
      permVec[i] = permVec[j];
      permVec[j] = intBuff;
      i++; j--;
    }
  }

  if(left < j) quickSortPerm(inVec, permVec, left, j);
  if(right > i) quickSortPerm(inVec, permVec, i, right);
}

/*******************************************************************/
void quickSortPerm(vector<unsigned long long int>& inVec, vector<unsigned long long int>& permVec, const int& left, const int& right) {

  int i = left, j = right;
  int intBuff;
  unsigned long long int buff;
  unsigned long long int pivot = inVec[(left+right)/2];
  while(i <= j){
    while (inVec[i]>pivot) i++;
    while (inVec[j]<pivot) j--;
    if (i<=j){
      // swap vector elements
      buff = inVec[i];
      inVec[i] = inVec[j];
      inVec[j] = buff;

      // swap elements of other vector in same manner
      intBuff = permVec[i];
      permVec[i] = permVec[j];
      permVec[j] = intBuff;
      i++; j--;
    }
  }

  if(left < j) quickSortPerm(inVec, permVec, left, j);
  if(right > i) quickSortPerm(inVec, permVec, i, right);
}

/*******************************************************************/
void quickSort(vector<float>& inVec, vector<int>& permVec, const int& left, const int& right) {

  int i = left, j = right;
  int intBuff;
  float buff;
  float pivot = inVec[(left+right)/2];
  while(i <= j){
    while (inVec[i]>pivot) i++;
    while (inVec[j]<pivot) j--;
    if (i<=j){
      // swap vector elements
      buff = inVec[i];
      inVec[i] = inVec[j];
      inVec[j] = buff;

      // swap elements of other vector in same manner
      intBuff = permVec[i];
      permVec[i] = permVec[j];
      permVec[j] = intBuff;
      i++; j--;
    }
  }

  if(left < j) quickSortPerm(inVec, permVec, left, j);
  if(right > i) quickSortPerm(inVec, permVec, i, right);
}

/***************************************************************/
/*               Romberg integration routines                  */
/***************************************************************/
void polint(vector<double> xa, double *ya, const double x, double &y, double &dy)
{
  int i,m,ns=0;
  double den,dif,dift,ho,hp,w;

  int n=xa.size();
  double *c, *d;
  c = new double[n];
  d = new double[n];
  dif=fabs(x-xa[0]);
  for (i=0;i<n;i++) {
    if ((dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=0;i<n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ((den=ho-hp) == 0.0){
	cerr << "Error in routine polint" << endl;
	exit(1);
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    y += (dy=(2*(ns+1) < (n-m) ? c[ns+1] : d[ns--]));
  }
  delete [] c;
  delete [] d;
}

double trapzd(double func(const double, const double), const double a, const double b, const double T, const int n)
{
  double x,tnm,sum,del;
  static double s;
  int it,j;

  if (n == 1) {
    return (s=0.5*(b-a)*(func(a, T)+func(b, T)));
  } else {
    for (it=1,j=1;j<n-1;j++) it <<= 1;
    tnm=it;
    del=(b-a)/tnm;
    x=a+0.5*del;
    for (sum=0.0,j=0;j<it;j++,x+=del) sum += func(x, T);
    s=0.5*(s+(b-a)*sum/tnm);
    return s;
  }
}

double qromb(double func(const double, const double), double a, double b, double T){

  const int JMAX=20, JMAXP=JMAX+1, K=5;
  const double EPS=1.0e-10;
  double ss, dss;
  double s[JMAX], h[JMAXP], s_t[K];
  vector<double> h_t(K, 0.0);
  int i,j;

  h[0]=1.0;
  for (j=1;j<=JMAX;j++) {
    s[j-1]=trapzd(func,a,b,T,j);
    if (j >= K) {
      for (i=0;i<K;i++) {
	h_t[i]=h[j-K+i];
	s_t[i]=s[j-K+i];
      }
      polint(h_t, s_t, 0.0, ss, dss);
      if (fabs(dss) <= EPS*fabs(ss)) return ss;
    }
    h[j]=0.25*h[j-1];
  }
  cerr << "Too many steps in routine qromb" << endl;
  h_t.clear();
  return 0.0;
}

/*******************************************************************/
// Planck curve times cross sections and divided by nu
// The Planck curve and cross section have been scaled for stable numerics
// Constants in Planck curve do not matter as they are divided out by normalisation
// Multiply any cross section with 1e-18 to obtain physical values
// x = nu/nu_0

//hydrogen
double fHI( double x, double T ){
  return 6.3 * pow( x, -1.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );
}

//dust
double f_dust_SMC( double x, double T ){

  double ksi = ksi_dust( x, SMC );
  double planck = pow( x, 3.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );

  return ksi*planck;

}
double f_dust_LMC( double x, double T ){

  double ksi = ksi_dust( x, LMC );
  double planck = pow( x, 3.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );

  return ksi*planck;

}
double f_dust_MW( double x, double T ){

  double ksi = ksi_dust( x, MW );
  double planck = pow( x, 3.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );

  return ksi*planck;

}
// Planck curve divided by nu
double PlanckDivNu( double x, double T ){
  return pow( x, 2.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );
}

//Planck curve
double Planck( double x, double T ){
  return pow( x, 3.0 ) / ( exp( planckConstant * x * nu0HI  / ( k * T ) ) - 1.0 );
}

//fitting function for dust extinction from Pei 1992 and Gnedin et al. 2008
double fitting_function( const double& x, const double& a, const double& b, const double& p, const double& q ){

  double denom = pow( x, p ) + pow( x, -1.0*q ) + b;

  return a/denom;

}

//determine dust extinction for scaled frequency x
double ksi_dust( double x, short int dust_model ){

  //change to real frequency
  double nu = x*nu0HI;

  double lambda = speed_of_light/nu; //wavelength in cm
  lambda *= 1e4; //in micrometer

  //values from gnedin et al. 2008 and Pei et al. 1992
  vector <float> a_arr, b_arr, p_arr, q_arr, l_arr;

  //short int dust_model = 1;
  short int n_param = 0;
  if(dust_model == SMC){
    n_param = 7;
    a_arr.resize(n_param);b_arr.resize(n_param);p_arr.resize(n_param);q_arr.resize(n_param);l_arr.resize(n_param);
    //SMC
    a_arr[0]=185.;a_arr[1]=27.;a_arr[2]=0.005;a_arr[3]=0.01;a_arr[4]=0.012;a_arr[5]=0.03;a_arr[6]=10.;
    b_arr[0]=90.;b_arr[1]=15.5;b_arr[2]=-1.95;b_arr[3]=-1.95;b_arr[4]=-1.8;b_arr[5]=0.;b_arr[6]=1.9;
    p_arr[0]=2.; p_arr[1]=4.; p_arr[2]=2.; p_arr[3]=2.; p_arr[4]=2.; p_arr[5]=2.; p_arr[6]=4.;
    q_arr[0]=2.; q_arr[1]=4.; q_arr[2]=2.; q_arr[3]=2.; q_arr[4]=2.; q_arr[5]=2.; q_arr[6]=15.;
    l_arr[0]=0.042; l_arr[1]=0.08; l_arr[2]=0.22; l_arr[3]=9.7; l_arr[4]=18.; l_arr[5]=25.; l_arr[6]=0.067;
  }else if(dust_model == LMC){
    n_param = 7;
    a_arr.resize(n_param);b_arr.resize(n_param);p_arr.resize(n_param);q_arr.resize(n_param);l_arr.resize(n_param);
    //LMC
    a_arr[0]=90.;a_arr[1]=19.;a_arr[2]=0.0023;a_arr[3]=0.005;a_arr[4]=0.006;a_arr[5]=0.02;a_arr[6]=10.;
    b_arr[0]=90.;b_arr[1]=21.;b_arr[2]=-1.95;b_arr[3]=-1.95;b_arr[4]=-1.8;b_arr[5]=0.;b_arr[6]=1.9;
    p_arr[0]=2.; p_arr[1]=4.5; p_arr[2]=2.; p_arr[3]=2.; p_arr[4]=2.; p_arr[5]=2.; p_arr[6]=4.;
    q_arr[0]=2.; q_arr[1]=4.5; q_arr[2]=2.; q_arr[3]=2.; q_arr[4]=2.; q_arr[5]=2.; q_arr[6]=15.;
    l_arr[0]=0.046; l_arr[1]=0.08; l_arr[2]=0.22; l_arr[3]=9.7; l_arr[4]=18.; l_arr[5]=25.; l_arr[6]=0.067;
  }else if(dust_model == MW){
    n_param = 6;
    a_arr.resize(n_param);b_arr.resize(n_param);p_arr.resize(n_param);q_arr.resize(n_param);l_arr.resize(n_param);
    //MW
    a_arr[0]=165.;a_arr[1]=14.;a_arr[2]=0.045;a_arr[3]=0.002;a_arr[4]=0.002;a_arr[5]=0.012;
    b_arr[0]=90.;b_arr[1]=4.;b_arr[2]=-1.95;b_arr[3]=-1.95;b_arr[4]=-1.8;b_arr[5]=0.;
    p_arr[0]=2.; p_arr[1]=6.5; p_arr[2]=2.; p_arr[3]=2.; p_arr[4]=2.; p_arr[5]=2.; 
    q_arr[0]=2.; q_arr[1]=6.5; q_arr[2]=2.; q_arr[3]=2.; q_arr[4]=2.; q_arr[5]=2.; 
    l_arr[0]=0.047; l_arr[1]=0.08; l_arr[2]=0.22; l_arr[3]=9.7; l_arr[4]=18.; l_arr[5]=25.;
  }else{
    cerr << " Error, no dust model, exiting!" << endl;
    exit(-1);
  }

  double dust_extinction = 0.0;
  for(short int i=0; i<n_param; i++){
    double l = lambda/l_arr[i];
    dust_extinction += fitting_function(l, a_arr[i], b_arr[i], p_arr[i], q_arr[i] );
  }

  return dust_extinction;

}

/**************************************************************************************************/
/*! This routine is a comparison kernel used in a sort routine to group
 *  vertices that are exported to the same processor.
 */
int site_compare_key(const void *a, const void *b)
{
  if( ( ( Send_Site *) a)->get_process() < ( ( ( Send_Site *) b)->get_process() ) )
    return -1;

  if( ( ( Send_Site *) a)->get_process() > ( ( ( Send_Site *) b)->get_process() ) )
    return +1;

  return 0;
}

unsigned int pow(const unsigned int& a, const int& b){

  double base = (double) a;
  double power = (double) b;

  double result = pow(base,power);

  unsigned int res = (unsigned int) floor(result);

  return res;

}

int pow(const int& a, const int& b){

  double base = (double) a;
  double power = (double) b;

  double result = pow(base,power);

  int res = (int) floor(result);

  return res;

}
