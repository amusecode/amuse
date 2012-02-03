using namespace std;

#include <fstream>

#include <cstdlib>
#include <vector>

#include "Star.h"
#include "Force.h"
#include "Dynamics.h"

#ifndef __CLUSTER_H
#define __CLUSTER_H

class Cluster
{
  mpreal t;
  int N;

  Dynamics dynamics;

  public:
  Force force;
  vector<Star> star;

  // Constructors
  Cluster();
  Cluster( string file, Force fo );

  void add_star(int id, mpreal m, mpreal radius, mpreal x, mpreal y, mpreal z, mpreal vx, mpreal vy, mpreal vz);

  // Set
  void set_t(mpreal T);
  void set_N(int N);

  // Get
  mpreal get_t();
  int get_N();
  Force* get_pointer_to_force();
  Star* get_pointer_to_star();
  Star* get_pointer_to_star(int index);
  mpreal get_E();

  mpreal get_dt();

  // Calculate
  void calc_a();
  void calc_a_dt();
  void leapfrog(mpreal dt);

  // Printers
  void print();
  void print( ofstream &data );
};

#endif


