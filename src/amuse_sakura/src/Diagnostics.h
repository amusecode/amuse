#include <vector>
#include <cmath>

#ifndef __Diagnostics_h
#define __Diagnostics_h

class Diagnostics {
  public:

  // Conserved quantities
  double get_mass(std::vector<double> &data);

  std::vector<double> get_rcm(std::vector<double> &data);
  std::vector<double> get_vcm(std::vector<double> &data);
  std::vector<double> get_lcm(std::vector<double> &data);

  double get_kinetic_energy(std::vector<double> &data);
  double get_potential_energy(std::vector<double> &data);
  double get_energy(std::vector<double> &data);

  // System properties
  double get_virial_radius(std::vector<double> &data);
  double get_harmonic_radius(std::vector<double> &data);
  double get_velocity_disperion(std::vector<double> &data);
};

#endif

