#include "Communicator.h"

#include <string>
#include <cmath>
#include <cstdlib>

#include "Particle.h"
#include "Particles.h"

#include "Interaction.h"
#include "Timestep.h"

#include "Twobody.h"

#ifndef __SAKURA_H
#define __SAKURA_H

class Sakura {

  Particles particles;
  Interaction interaction;  
  Timestep timestep;

  Twobody twobody;

  unsigned int Nstep;
  double sum_dt, sum_inv_dt;

  public:

  Sakura();
  Sakura(Particles particles);
  Sakura(Particles particles, double dt);

  void set_particles(Particles particles);
  void set_dt(double dt);
  void set_tolerance(double tolerance);
  
  Particles get_particles();
  Particles* get_pointer_to_particles();
  double get_dt();
  double get_tolerance();
  unsigned int get_Nstep();
  double get_sum_dt();
  double get_sum_inv_dt();
  Particle* get_pointer_to_star(int index);

  void evolve(double t, Communicator &communicator);
  void step(vector<Particle> &particle, double dt, Communicator &communicator);

  void initial_calculations();
  double get_t();
  void set_t(double t);
  void step(double dt, Communicator &communicator);

  vector<double> get_coordinates(vector<Particle> &particle);
  void update_particles(vector<Particle> &particle, vector<double> &coordinates);
  void update_particles(vector<double> &coordinates);
  vector<double> get_data();
};

#endif


