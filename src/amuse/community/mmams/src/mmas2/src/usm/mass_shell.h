#ifndef _MASS_SHELL_H_
#define _MASS_SHELL_H_

#include "chemical_composition.h"

class mass_shell {
protected:
public:
  int  id;

  real dm;
  real mass;
  real radius;

  real density;
  real pressure;
  real e_thermal;
  real entropy;
  real temperature;
  real beta;
  real buoyancy;

  real luminosity;
  real mean_mu;

  real angular_momentum;
  real opacity;

  real potential;
  
  chemical_composition composition;

  mass_shell() {
    id = 0;
    dm = mass = radius = 0;
    density = pressure = e_thermal = entropy = temperature = 0;
    beta = 0;
    luminosity = 0;
    mean_mu = 0;
    angular_momentum = 0;
    opacity = 0;
    potential = 0;
    buoyancy = 0.0;
  }
  ~mass_shell() {
  }

  void operator=(mass_shell &shell) {
    id = shell.id;

    dm = shell.dm;
    mass = shell.mass;
    radius = shell.radius;

    density = shell.density;
    pressure = shell.pressure;
    e_thermal = shell.e_thermal;
    entropy = shell.entropy;
    temperature = shell.temperature;
    beta = shell.beta;

    luminosity = shell.luminosity;
    mean_mu = shell.mean_mu;

    angular_momentum = shell.angular_momentum;
    opacity = shell.opacity;

    potential = shell.potential;
    buoyancy = shell.buoyancy;

    composition = shell.composition;
  }

  void write(FILE *fout) {
    fprintf(fout, "(mass_shell\n");

    writei(id, "id");

    writef(dm, "dm");
    writef(mass, "mass");
    writef(radius, "radius");

    writef(density, "density");
    writef(pressure, "pressure");
    writef(e_thermal, "e_thermal");
    writef(entropy, "entropy");
    writef(temperature, "temperature");
    writef(beta, "beta");
    

    writef(luminosity, "luminosity");
    writef(mean_mu, "molecular_weight");

    writef(angular_momentum, "angular_momentum");
    writef(opacity, "opacity");
    writef(potential, "potential");
    writef(buoyancy, "buoyancy");
    
    composition.write(fout);

    fprintf(fout, ")mass_shell\n");
  }

  void read(FILE *fin, int self_check) {
    char line[LINESIZE], variable[128], equal[10], value1[128], value2[128], value3[128];

    if (self_check == 1) {
      fgets(line, LINESIZE, fin);
      sscanf(line, "%s", variable);
      if (strcmp(variable, "(mass_shell") != 0) {
	cerr << " I tried to read, but it seems for me that it is not a mass_shell data " << endl;
	cerr << "       ..... I am sorry, but I have to give up. bye bye ... " << endl;
	PRL(sqrt(-1.0));
	exit(-1);
      }
    }

    while(1) {
      fgets(line, LINESIZE, fin);
      sscanf(line, "%s %s  %s %s %s", variable, equal,  value1, value2, value3);

      readi(id, "id");

      readf(dm, "dm");
      readf(mass, "mass");
      readf(radius, "radius");

      readf(density, "density");
      readf(pressure, "pressure");
      readf(e_thermal, "e_thermal");
      readf(entropy, "entropy");
      readf(temperature, "temperature");
      readf(beta, "beta");
    
      readf(luminosity, "luminocity");
      readf(mean_mu, "mean_mu");
      readf(mean_mu, "molecular_weight");
      
      readf(angular_momentum, "angular_momentum");
      readf(opacity, "opacity");
      readf(potential, "potential");
      readf(buoyancy, "buoyancy");

      check("(chemical_composition") composition.read(fin, 0);

      check(")mass_shell") break;
    }
  }

};

#endif // _MASS_SHELL_H_
