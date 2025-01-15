#include <iostream>
#include <fstream>
#include <vector>

#include "Particle.h"

#ifndef __PARTICLES_H
#define __PARTICLES_H

class Particles {

  double t;
  int N;
  std::vector<Particle> particle;

  public:

  Particles();
  Particles(std::vector<Particle> particle);
  Particles(double t, std::vector<Particle> particle);
  Particles(double t, int N, std::vector<Particle> particle);

  void set_t(double t);
  void set_N(int N);
  void set_particles(std::vector<Particle> particle);
  void set_particle(int index, Particle particle);

  double get_t();
  int get_N();
  std::vector<Particle> get_particles();
  Particle get_particle(int index);
  Particle* get_pointer_to_star(int index);

  void set_data(std::vector<double> &data);
  std::vector<double> get_data();
  
  void add_to_t(double dt);
  void add_particle(Particle particle);
  void remove_particle(int index);

  void print();
  void print(double t_cpu);
  void print(std::ofstream &odata);
  void print(double t_cpu, std::ofstream &odata);
};

#endif


