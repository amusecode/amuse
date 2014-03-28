#include "Communicator.h"

using namespace std;
#include <vector>
#include <cmath>

#include "Particle.h"
#include "Particles.h"

#ifndef __INTERACTION_H
#define __INTERACTION_H

class Interaction {

  public:

  Interaction();
  void calc_a(vector<Particle> &particle, Communicator &communicator);
};

#endif
