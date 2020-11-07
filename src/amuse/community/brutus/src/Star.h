#include <iostream>
using namespace std;
#include <fstream>
#include <cstdlib>

#ifndef __Star_h
class DEBUG
{
public:
  static ofstream DEB;
};
#endif

#include <vector>

#include "mpreal.h"
using namespace mpfr;

#ifndef __Star_h
#define __Star_h

#define use_additional_acc

class Star {
public:
  size_t id;
  mpreal m,r;
  vector<mpreal> x;
  vector<mpreal> v;
  vector<mpreal> a, a0;
  #ifdef use_additional_acc
  vector<mpreal> a_step;
  #endif // use_additional_acc
  Star();
  Star(mpreal m, vector<mpreal> x, vector<mpreal> v
    #ifdef use_additional_acc
    , vector<mpreal> a_step
    #endif // use_additional_acc
  );

  friend ostream & operator << (ostream &so, const Star &si) {
    so << si.m << " " << si.x[0] << " "<< si.x[1] << " "<< si.x[2] << " "
                      << si.v[0] << " "<< si.v[1] << " "<< si.v[2] << endl;
    return so;
  }
};

#endif


