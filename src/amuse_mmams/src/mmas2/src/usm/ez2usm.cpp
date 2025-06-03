#include "usm.h"

int const n_mesh = 199;
real radii[n_mesh], radii2r[n_mesh], radii2m[n_mesh];
real tempr[n_mesh], tempr2r[n_mesh], tempr2m[n_mesh];
real mass[n_mesh], mass2r[n_mesh], mass2m[n_mesh];
real lum[n_mesh], lum2r[n_mesh], lum2m[n_mesh];
real psi[n_mesh], psi2r[n_mesh], psi2m[n_mesh];
real X_H[n_mesh], X_H2r[n_mesh], X_H2m[n_mesh];
real X_He[n_mesh], X_He2r[n_mesh], X_He2m[n_mesh];
real X_C[n_mesh], X_C2r[n_mesh], X_C2m[n_mesh];
real X_N[n_mesh], X_N2r[n_mesh], X_N2m[n_mesh];
real X_Ne[n_mesh], X_Ne2r[n_mesh], X_Ne2m[n_mesh];
real X_O[n_mesh], X_O2r[n_mesh], X_O2m[n_mesh];
real X_Mg[n_mesh], X_Mg2r[n_mesh], X_Mg2m[n_mesh];
real pressure[n_mesh], pressure2r[n_mesh], pressure2m[n_mesh];
real density[n_mesh], density2r[n_mesh], density2m[n_mesh];
real opacity[n_mesh], opacity2r[n_mesh], opacity2m[n_mesh];
real e_thermal[n_mesh], e_thermal2r[n_mesh], e_thermal2m[n_mesh];
real entropy[n_mesh], entropy2r[n_mesh], entropy2m[n_mesh];
real mean_mu[n_mesh], mean_mu2r[n_mesh], mean_mu2m[n_mesh];
real mu_e[n_mesh], mu_e2r[n_mesh], mu_e2m[n_mesh];

int main(int argc, char *argv[])
{
  cerr << "Reading EZ log file from stdin" << endl;

  real rtemp;
  
  for (int i=0; i<8 ; i++) cin >> rtemp;
  for (int i=0; i<24; i++) cin >> rtemp;
  for (int i=0; i<27; i++) cin >> rtemp;

  for (int i=0; i<n_mesh; i++) cin >> radii[n_mesh-1-i];           // 1
  for (int i=0; i<n_mesh; i++) cin >> tempr[n_mesh-1-i];           // 2
  for (int i=0; i<n_mesh; i++) cin >> mass[n_mesh-1-i];           // 3
  for (int i=0; i<n_mesh; i++) cin >> lum[n_mesh-1-i];           // 4
  for (int i=0; i<n_mesh; i++) cin >> psi[n_mesh-1-i];           // 5
  for (int i=0; i<n_mesh; i++) cin >> X_H[n_mesh-1-i];           // 6
  for (int i=0; i<n_mesh; i++) cin >> X_He[n_mesh-1-i];           // 7
  for (int i=0; i<n_mesh; i++) cin >> X_C[n_mesh-1-i];           // 8
  for (int i=0; i<n_mesh; i++) cin >> X_N[n_mesh-1-i];           // 9
  for (int i=0; i<n_mesh; i++) cin >> X_O[n_mesh-1-i];           // 10
  for (int i=0; i<n_mesh; i++) cin >> X_Ne[n_mesh-1-i];           // 11
  for (int i=0; i<n_mesh; i++) cin >> X_Mg[n_mesh-1-i];           // 12
  for (int i=0; i<n_mesh; i++) cin >> pressure[n_mesh-1-i];           // 13
  for (int i=0; i<n_mesh; i++) cin >> density[n_mesh-1-i];           // 14
  for (int i=0; i<n_mesh; i++) cin >> opacity[n_mesh-1-i];           // 15
  for (int i=0; i<n_mesh; i++) cin >> e_thermal[n_mesh-1-i];           // 16
  for (int i=0; i<n_mesh; i++) cin >> entropy[n_mesh-1-i];           // 17
  for (int i=0; i<n_mesh; i++) cin >> mean_mu[n_mesh-1-i];           // 18
  for (int i=0; i<n_mesh; i++) cin >> mu_e[n_mesh-1-i];           // 19

  mass_shell zone;
  usm star;

  real max_mass = 0.0;
  real max_radius = 0.0;

  for (int i = 0; i < n_mesh; i++) {
    zone.id = i;

    zone.radius = radii[i];
    zone.mass = mass[i];

    max_mass = max(max_mass, mass[i]);
    max_radius = max(max_radius, radii[i]);

    zone.density = density[i];
    zone.pressure = pressure[i];
    zone.e_thermal = e_thermal[i];
    zone.temperature = tempr[i];
    zone.entropy = entropy[i];
    
    zone.mean_mu = mean_mu[i];
    
    zone.composition.H1 = X_H[i];
    zone.composition.He4 = X_He[i];
    zone.composition.O16 = X_O[i];
    zone.composition.N14 = X_N[i];
    zone.composition.C12 = X_C[i];
    zone.composition.Ne20 = X_Ne[i];
    zone.composition.Mg24 = X_Mg[i];
    
    star.add_shell(zone);

  }
  star.build_hashtable();

  star.star_mass = max_mass;
  star.star_radius = max_radius;

  cerr << "Dumping stellar model into stdout " << endl;
  star.write(stdout);
  
  cerr << "end-of-code" << endl;
  return 0;

}
