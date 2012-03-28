#ifndef UTILIS_H
#define UTILIS_H

#include <cmath>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <string>
#include <cstdlib>
#include <vector>
#include "./types.h"
#include <mpi.h>
#include "sys/time.h"
#include <stdexcept>
#include <cuda_runtime_api.h>

using namespace std;

class BadConversion : public std::runtime_error {
   public:
      BadConversion(std::string const& s)
      : std::runtime_error(s)
      { }
}; 

inline void get_times(double *tm __attribute__((unused))){
#ifdef CHECK_TIMES
   int local_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);

   if(local_rank == 0){
      cudaDeviceSynchronize();
      struct timeval tv;

      gettimeofday(&tv, NULL);
      int sec = tv.tv_sec;
      int microsec = tv.tv_usec;

      *tm = sec + microsec * 0.000001;
   }

   MPI_Barrier(MPI_COMM_WORLD);
#endif
}


inline void set_times(double time __attribute__((unused)), double *T __attribute__((unused))){
#ifdef CHECK_TIMES
   int local_rank;
   MPI_Comm_rank(MPI_COMM_WORLD, &local_rank);

   if(local_rank == 0)
      *T = time;

   MPI_Barrier(MPI_COMM_WORLD);
#endif

}
inline unsigned int countlines(std::ifstream &fbuf){

    std::string line;
    unsigned int count = 0;

    while(getline(fbuf, line))
      count++;

    ///////rewind file///////
    fbuf.clear();
    fbuf.seekg(0, ios::beg);
    ////////////////////////

    return count;
}


inline bool isNumeric_int(const char *str, int *number){
   std::istringstream test( str );
   *number = -1;
   test >> *number;

   return (test.rdbuf()->in_avail() == 0);

}

inline char* to_char(std::string s){
   char *buf = new char[s.length()+1];
   strcpy(buf,s.c_str());
   return buf;
}

inline std::string to_string( char* val ){
   std::stringstream s;
   s << val;
   return s.str();
}
inline std::string to_string( int val ){
   std::stringstream s;
        s << val;
        return s.str();
}

inline double to_double (std::string val){
   double number;
   std::istringstream stm;
   stm.str(val);
   if(!(stm>>number))
      throw BadConversion("convertToDouble(\"" + val + "\")");
   return number;
}

inline unsigned int to_uint (std::string val){
   unsigned int number;
   std::istringstream stm;
   stm.str(val);
   if(!(stm>>number))
      throw BadConversion("convertToInt(\"" + val + "\")");
   return number;
}



class Utilis
{

 public :

    vector<double> CdM(double *x, double *y, double *z, double *w, int num);
  vector<double> CdD(double *x, double *y, double *z, double *w, int num);

 private :

  //variables
  int k;
  double cmx, cmy, cmz, total_mass;

};

#endif

