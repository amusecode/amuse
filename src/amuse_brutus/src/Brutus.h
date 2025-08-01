#include "Star.h"
#include "Cluster.h"
#include "Bulirsch_Stoer.h"

#ifndef __Brutus_h
#define __Brutus_h

class Brutus {
  mpreal t;
  int N;  
  std::vector<mpreal> data;

  mpreal tolerance;
  int numBits;

  mpreal eta, dt;

  Cluster cl;
  Bulirsch_Stoer bs;

  public:

  Brutus();
  Brutus(std::vector<mpreal> &data);
  Brutus(mpreal &t, std::vector<mpreal> &data, mpreal &tolerance);
  Brutus(mpreal &t, std::vector<mpreal> &data, mpreal &tolerance, int &numBits);

  void set_data(std::vector<mpreal> &data);
  void set_eta(mpreal &eta);
  void set_tolerance(mpreal &tolerance);
  void set_numBits(int &numBits);
  void set_t_begin(mpreal &t_begin);

  mpreal get_eta(mpreal tolerance);
  mpreal get_tolerance();
  int get_numBits();
  int get_numBits(mpreal tolerance);
  mpreal fit_slope(std::vector<mpreal> &x, std::vector<mpreal> &y);

  void setup();
  void evolve(mpreal t_end);
  
  mpreal get_t();
  std::vector<mpreal> get_data();
  std::vector<double> get_data_double();
  std::vector<std::string> get_data_string();
};

#endif


