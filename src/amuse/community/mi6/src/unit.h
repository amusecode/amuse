#ifndef UNIT_H
#define UNIT_H

#include<sys/time.h>

/////////// astrophisical constant (MKS) //////////////
const double Gsi = 6.67259E-11; // m^3/kg/sec^2
const double csi = 2.99792458E08;  // m/sec
const double Pcsi = 3.08567802E16;   //1pc
const double AUsi = 1.49597892E11;  // 1AU
const double Rsolsi = 6.9599E08;  //1R_sol
const double Msolsi = 1.989E30; //1M_sol
const double Yearsi = 3.1558149984E07;  //1year


const double n2=0.5;
const double n3=1.0/3.0;
const double n4=1.0/4.0;
const double n5=1.0/5.0;
const double n6=1.0/6.0;
const double n7=1.0/7.0;
const double n8=1.0/8.0;
const double n9=1.0/9.0;
const double n10=1.0/10.0;
const double n12=1.0/12.0;
const double n24=1.0/24.0;
const double n30=1.0/30.0;
const double n32=1.0/32.0;
const double n60=1.0/60.0;
const double n40=1.0/40.0;
const double n42=1.0/42.0;
const double n84=1.0/84.0;
const double n120=1.0/120.0;
const double n280=1.0/280.0;
const double n720=1.0/720.0;
const double n1024=1.0/1024.0;
const double n1680=1.0/1680.0;
const double n32768=1.0/32768.0;
const double n3_28=3.0/28.0;
const double n208_15=208.0/15.0;
const double n24_5=24.0/5.0;
const double n12_5=12.0/5.0;
const double n8_5=8.0/5.0;
const double n32_5=32.0/5.0;
const double n4_5=4.0/5.0;

const double LARGE_DOUBLE = 99999999999999999.9;
const int LARGE_INT = 999999999;
#ifndef NOMPI
#include <mpi.h>
#else
#ifdef WIN32
#include <windows.h>
#endif
#endif

//#include"/data/iwasawa/work/amuse/prerequisites/include/mpi.h"

inline double gettimeofday_sec(){
  
#ifndef NOMPI
    return MPI_Wtime();
#else
#ifdef WIN32
    FILETIME filetime;
    GetSystemTimeAsFileTime(&filetime);
    unsigned long long longtime = filetime.dwHighDateTime;
    longtime <<=32;
    longtime |= filetime.dwLowDateTime;
    longtime /=10;
    longtime -= 11644473600000000ULL;
    return longtime / 1000000.0;
#else
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec;
#endif
#endif
}


#endif
