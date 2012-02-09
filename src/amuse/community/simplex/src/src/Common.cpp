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
//sorts ascending!!!!
void quickSortPerm(vector<unsigned long int>& inVec, vector<unsigned long int>& permVec, const int& left, const int& right) {

  int i = left, j = right;
  int intBuff;
  unsigned long int buff;
  unsigned long int pivot = inVec[(left+right)/2];
  while(i <= j){
    while (inVec[i]<pivot) i++;
    while (inVec[j]>pivot) j--;
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
//sorts descending!!!!!
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

//find the zero of function func
double zbrent(double func( double,  double),  double x2,  double tol,  double Temp,  double perBin){

   int ITMAX=100;
   double EPS=numeric_limits<double>::epsilon();
  int iter;
  double a=1.0,b=x2,c=x2,d=0.0,e=0.0,min1,min2;
  double fa=func(a, Temp)-perBin, fb=func(b, Temp)-perBin,fc,p,q,r,s,tol1,xm;

  if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)){
    cerr << "Root must be bracketed in zbrent" << endl;
    exit(1);
  }
  fc=fb;
  for (iter=0;iter<ITMAX;iter++) {
    if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
      c=a;
      fc=fa;
      e=d=b-a;
    }
    if (fabs(fc) < fabs(fb)) {
      a=b;
      b=c;
      c=a;
      fa=fb;
      fb=fc;
      fc=fa;
    }
    tol1=2.0*EPS*fabs(b)+0.5*tol;
    xm=0.5*(c-b);
    if (fabs(xm) <= tol1 || fb == 0.0) return b;
    if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
      s=fb/fa;
      if (a == c) {
	p=2.0*xm*s;
	q=1.0-s;
      } else {
	q=fa/fc;
	r=fb/fc;
	p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	q=(q-1.0)*(r-1.0)*(s-1.0);
      }
      if (p > 0.0) q = -q;
      p=fabs(p);
      min1=3.0*xm*q-fabs(tol1*q);
      min2=fabs(e*q);
      if (2.0*p < (min1 < min2 ? min1 : min2)) {
	e=d;
	d=p/q;
      } else {
	d=xm;
	e=d;
      }
    } else {
      d=xm;
      e=d;
    }
    a=b;
    fa=fb;
    if (fabs(d) > tol1)
      b += d;
    else
      b += SIGN(tol1,xm);
    fb=func(b, Temp)-perBin;
  }
  cerr<< "Maximum number of iterations exceeded in zbrent" << endl;
  return 0.0;
}

inline float SIGN( const float &a, const double &b)
{
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}



double rec_rad_escape(const double& tau){

  //  double tmp1 = 1.0 - 2.*pow(tau, 2.);
  //double tmp2 = 1.0 + 2.*tau;
  //double tmp3 = exp(-2.*tau);
  //double tmp4 = 3./(8.*tau);

  //double f_esc = tmp4*( tmp3*tmp2 - tmp1 );
 
  //double f_esc = (tau > 0.0 ) ? 3.*( exp(-tau) * (1.0 + tau) + 2.*pow(0.5*tau, 2.) - 1.0 )/(8.*pow(0.5*tau,3.)) : 1.0;

  //for very small tau solution is unstable
  double f_esc = (tau > 1.e-3) ? 3.*( exp(-tau) * (1.0 + tau) + 2.*pow(0.5*tau, 2.) - 1.0 )/(8.*pow(0.5*tau,3.)) : 1.0;


  // if(f_esc < 0.0){
  //   f_esc = 1.0;
  // }

  return f_esc;

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
