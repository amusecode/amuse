#include <iostream>
#include <vector>

#include "mpreal.h"
using namespace mpfr;

#ifndef __Star_h
#define __Star_h

class Star {
public:
  mpreal m;
  std::vector<mpreal> r;
  std::vector<mpreal> v;
  std::vector<mpreal> a, a0;

  Star();
  Star(mpreal m, std::vector<mpreal> r, std::vector<mpreal> v);

  friend std::ostream & operator << (std::ostream &so, const Star &si) {
    so << si.m << " " << si.r[0] << " "<< si.r[1] << " "<< si.r[2] << " "
                      << si.v[0] << " "<< si.v[1] << " "<< si.v[2] << std::endl;
    return so; 
  }
};

#endif


