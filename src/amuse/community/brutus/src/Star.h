#include <iostream>
using namespace std;

#include <vector>

#include "mpreal.h"
using namespace mpfr;

#ifndef __Star_h
#define __Star_h

class Star {
public:
  mpreal m;
  vector<mpreal> r;
  vector<mpreal> v;
  vector<mpreal> a, a0;

  Star();
  Star(mpreal m, vector<mpreal> r, vector<mpreal> v);

  friend ostream & operator << (ostream &so, const Star &si) {
    so << si.m << " " << si.r[0] << " "<< si.r[1] << " "<< si.r[2] << " "
                      << si.v[0] << " "<< si.v[1] << " "<< si.v[2] << endl;
    return so; 
  }
};

#endif


