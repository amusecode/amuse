using namespace std;
#include <vector>
#include <cmath>

#ifndef __Diagnostics_h
#define __Diagnostics_h

class Diagnostics {
  public:

  // Conserved quantities
  double get_mass(vector<double> &data);

  vector<double> get_rcm(vector<double> &data);
  vector<double> get_vcm(vector<double> &data);
  vector<double> get_lcm(vector<double> &data);

  double get_kinetic_energy(vector<double> &data);
  double get_potential_energy(vector<double> &data);
  double get_energy(vector<double> &data);

  // System properties
  double get_virial_radius(vector<double> &data);
  double get_harmonic_radius(vector<double> &data);
  double get_velocity_disperion(vector<double> &data);
};

#endif

