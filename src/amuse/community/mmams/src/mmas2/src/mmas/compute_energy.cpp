#include "mmas.h"
#include "include/units.h"

real mmas::compute_stellar_energy(usm &model) {
  real energy = 0;
  real m_prev = 0.0;
  for (int i = 0; i < model.get_num_shells(); i++) {
    mass_shell &shell = model.get_shell(i);
    if (shell.radius == 0) {
      m_prev  = shell.mass;
      continue;
    }
    real de_grav  = -shell.mass/shell.radius;
    real de_therm = 1.5*uK*shell.temperature/shell.mean_mu/uM_U + 
      uA_RAD*pow(shell.temperature, 4.0)/shell.density;
    de_therm *= 1.0/(uG*uMSUN/uRSUN);
    energy += (de_grav + de_therm) * (shell.mass - m_prev);
    m_prev  = shell.mass;
  }
  return energy;
}
