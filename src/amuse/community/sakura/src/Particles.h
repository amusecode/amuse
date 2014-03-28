#include <iostream>
using namespace std;
#include <fstream>
#include <vector>

#include "Particle.h"

#ifndef __PARTICLES_H
#define __PARTICLES_H

class Particles {

  double t;
  int N;
  vector<Particle> particle;

  public:

  Particles();
  Particles(vector<Particle> particle);
  Particles(double t, vector<Particle> particle);
  Particles(double t, int N, vector<Particle> particle);

  void set_t(double t);
  void set_N(int N);
  void set_particles(vector<Particle> particle);
  void set_particle(int index, Particle particle);

  double get_t();
  int get_N();
  vector<Particle> get_particles();
  Particle get_particle(int index);
  Particle* get_pointer_to_star(int index);

  void set_data(vector<double> &data);
  vector<double> get_data();
  
  void add_to_t(double dt);
  void add_particle(Particle particle);
  void remove_particle(int index);

  void print();
  void print(double t_cpu);
  void print(ofstream &odata);
  void print(double t_cpu, ofstream &odata);
};

#endif


