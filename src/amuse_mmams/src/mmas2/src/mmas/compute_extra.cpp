#include "mmas.h"
#include "eos/eos.h"

real mmas::get_lagrad(usm &model, real frac) {
  PRL(model.star_mass);
  for (int i = 1; i < model.get_num_shells(); i++)
    if (model.get_shell(i).mass > frac*model.star_mass)
      return model.get_shell(i-1).radius;
  return -1;
}

void mmas::compute_extra() {
#define am(x) (1.0+Amass[x]/2.0)/Amass[x]
  real Amass[] = {1, 4, 16, 14, 12, 20, 24, 28, 56};
  real mass, radius;

  mass = radius = 0;
  for (int i = 0; i < model_a->get_num_shells(); i++) {
    mass_shell &shell = model_a->get_shell(i);
    real mean_mu = 2 * shell.composition.H1 + 
      am(1) * shell.composition.He4 +
      am(2) * shell.composition.O16 +
      am(3) * shell.composition.N14 +
      am(4) * shell.composition.C12 + 
      am(5) * shell.composition.Ne20 + 
      am(6) * shell.composition.Mg24 +
      am(7) * shell.composition.Si28 + 
      am(8) * shell.composition.Fe56;
    mean_mu = 1.0/mean_mu;

    shell.mean_mu = mean_mu;
    shell.temperature = compute_temperature(shell.density, shell.pressure, shell.mean_mu);
    shell.entropy     = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
    real Pgas = shell.density/(shell.mean_mu*uM_U) * uK * shell.temperature;
    shell.beta = Pgas / (shell.pressure);

    mass   = max(mass, shell.mass);
    radius = max(mass, shell.radius);
  }
  model_a->star_mass   = mass;
  model_a->star_radius = radius;

  if (model_a->get_shell(0).composition.H1 < 0.1) {
    model_a->star_age = 1.0;
  } else {
    model_a->star_age = 0.0;
  }

  mass = radius = 0;
  for (int i = 0; i < model_b->get_num_shells(); i++) {
    mass_shell &shell = model_b->get_shell(i);
    real mean_mu = 2 * shell.composition.H1 + 
      am(1) * shell.composition.He4 +
      am(2) * shell.composition.O16 +
      am(3) * shell.composition.N14 +
      am(4) * shell.composition.C12 + 
      am(5) * shell.composition.Ne20 + 
      am(6) * shell.composition.Mg24 +
      am(7) * shell.composition.Si28 + 
      am(8) * shell.composition.Fe56;
    mean_mu = 1.0/mean_mu;

    shell.mean_mu = mean_mu;
    shell.temperature = compute_temperature(shell.density, shell.pressure, shell.mean_mu);
    shell.entropy     = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
    real Pgas = shell.density/(shell.mean_mu*uM_U) * uK * shell.temperature;
    shell.beta = Pgas / (shell.pressure);

    mass   = max(mass, shell.mass);
    radius = max(mass, shell.radius);
  }
  model_b->star_mass   = mass;
  model_b->star_radius = radius;
  if (model_b->get_shell(0).composition.H1 < 0.1) {
    model_b->star_age = 1.0;
  } else {
    model_b->star_age = 0.0;
  }


  usm *temp_model = model_a;
  if (model_a->star_mass < model_b->star_mass) {
    model_a = model_b;
    model_b = temp_model;
  }
  PRL(model_a->star_age);
  PRL(model_b->star_age);
  
}

real mmas::compute_orbital_energy() {
  return 0;
}
