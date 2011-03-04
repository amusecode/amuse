// Equation of State

#ifndef __EOS__
#define __EOS__

#include "include/stdinc.h"
#include "include/units.h"


double calc_temp(double, double);
double find_beta(double, double);

real compute_beta(real pressure, real entropy, real mean_mu);
real compute_entropy(real density, real temperature, real mean_mu);
real compute_temperature(real density, real pressure, real mean_mu);
real compute_density(real entropy, real pressure, real mean_mu);
real compute_pressure(real density, real ethermal, real mean_mu);
real compute_energy(real density, real temperature, real mean_mu);

// class equation_of_state {
// protected:

//   /* macrosopic thermodynamic parameters */

//   real _gamma_gas;
//   real _density;
//   real _e_thermal;
//   real _pressure;
//   real _temperature;
//   real _mean_mu;
//   real _entr_var;

//   /* microscopic parameters, like chemical composition, etc ... */

//   int _n_elements;
//   real *_elements;

// public:

//   real& gamma_gas() {return _gamma_gas;}
//   real& density() {return _density;}
//   real& e_thermal() {return _e_thermal;}
//   real& pressure() {return _pressure;}
//   real& temperature() {return _temperature;}
//   real& mean_mu() {return _mean_mu;}
//   real& entr_var() {return _entr_var;}
//   int&  n_elements() {return _n_elements;}
//   real& elements(int i) {return _elements[i];}
  
  
//   equation_of_state() {
//     _gamma_gas = 5.0/3.0;
//     _n_elements = 0;
//     _elements = NULL;
//     _density = _e_thermal = _pressure = _temperature = 0;
//     _mean_mu = _entr_var = 0;
//   };
//   virtual ~equation_of_state() {
//     if (_elements != NULL) delete[] _elements;
//   };

//   virtual eos_type get_eos_type() = 0;
  
//   virtual void compute_state() = 0;
//   virtual real get_sound_speed() = 0;
//   virtual real get_entropy() = 0;
//   virtual real get_A() = 0;
//   virtual real get_beta() = 0;
//   virtual real compute_e_thermal() = 0;
//   virtual real compute_temperature(real,real,real) = 0;
//   virtual real compute_e_thermal(real, real, real) = 0;
//   virtual real compute_mean_mu() = 0;

//   void set_elements(int n, real el[]) {
//     if (_elements != NULL) delete[] _elements;
//     _n_elements = n;
//     if (_n_elements > 0) {
//       _elements = new real[_n_elements];
//       for (int i = 0; i < _n_elements; i++) _elements[i] = el[i];
//     }
//   }
  
//   void operator=(equation_of_state &eos) {
//     _gamma_gas   = eos._gamma_gas;
//     _density     = eos._density;
//     _e_thermal   = eos._e_thermal;
//     _entr_var    = eos._entr_var;
//     _pressure    = eos._pressure;
//     _temperature = eos._temperature;
//     _mean_mu     = eos._mean_mu;

//     if (_elements != NULL) {
//       delete[] _elements;
//       _elements = NULL;
//       _n_elements = 0;
//     }
//     if (eos._n_elements > 0) {
//       _n_elements = eos._n_elements;
//       _elements = new real[_n_elements];
//       for (int i = 0; i < _n_elements; i++) _elements[i] = eos._elements[i];
//     }
//   }

//   void write(FILE *fout) {
//     fprintf(fout, "(equation_of_state\n");
    
//     writef(_gamma_gas, "gamma_gas");
//     writef(_density, "density");
//     writef(_e_thermal, "e_thermal");
//     writef(_entr_var, "entr_var");
//     writef(_pressure, "pressure");
//     writef(_temperature, "temperature");
//     writef(_mean_mu, "mean_mu");

//     writei(_n_elements, "n_elements");
//     fprintf(fout, "   elements = ");
//     for (int i = 0; i < _n_elements; i++) {
//       double val = _elements[i];
//       fprintf(fout, "%.17le  ", val);
//     }
//     fprintf(fout, "\n");

//     fprintf(fout, ")equation_of_state\n");
//   }

//   void read(FILE *fin, int self_check) {
//     char line[LINESIZE], variable[128], equal[10], value1[128], value2[128], value3[128];
    
//     if (self_check == 1) {
//       fgets(line, LINESIZE, fin);
//       sscanf(line, "%s", variable);
//       if (strcmp(variable, "(equation_of_state") != 0) {
// 	cerr << " I tried to read, but it seems for me that it is not a <equation_of_state> data " << endl;
// 	cerr << "       ..... I am sorry, but I have to give up. bye bye ... " << endl;
// 	PRL(sqrt(-1.0));
// 	exit(-1);
//       }
//     }

//     while(1) {
//       fgets(line, LINESIZE, fin);
//       sscanf(line, "%s %s  %s %s %s", variable, equal,  value1, value2, value3);

//       readf(_gamma_gas, "gamma_gas");
//       readf(_density, "density");
//       readf(_e_thermal, "e_thermal");
//       readf(_entr_var, "entr_var");
//       readf(_pressure, "pressure");
//       readf(_temperature, "temperature");
//       readf(_mean_mu, "mean_mu");

//       if (strcmp(variable, "n_elements") == 0) {
// 	_n_elements = atoi(value1);
// 	if (_elements != NULL) delete[] _elements;
// 	fscanf(fin, "%s %s", value1, value1);
// 	if (_n_elements > 0) _elements = new real[_n_elements];
// 	else _elements = NULL;
// 	for (int i = 0; i < _n_elements; i++) {
// 	  double val = 0;
// 	  fscanf(fin, "%lf" , &val);
// 	  _elements[i] = val;
// 	}
// 	fgets(line, LINESIZE, fin);
//       }
      
//       check(")equation_of_state") break;
//     }
//   }
// };


// /* --------------------------- */
// /* adiabatic equation os state */
// /* --------------------------- */

// class adiabatic : public equation_of_state {
// protected:
//   /* inhertied from the base class */

// public:
//   adiabatic() : equation_of_state() {
//     _gamma_gas = 5.0/3.0;
//   };
//   ~adiabatic() {}
  
//   void compute_state() {
//     _pressure = _e_thermal * _density * (_gamma_gas - 1);
//     _temperature = _pressure/_density * _mean_mu;
//   }
//   real get_sound_speed() {
//     real sound_speed = sqrt( _gamma_gas * _pressure/_density);
//     return sound_speed;
//   }
//   real get_entropy() { return 0;}
//   real get_A() { return 0;}
//   real get_beta() {return 1.0;}
//   real compute_e_thermal() {return _e_thermal;}

//   eos_type get_eos_type() { return ADIABATIC;}
//   real compute_temperature(real density, real pressure, real mean_mu) {
//     return -1.0;
//   };
//   real compute_e_thermal(real, real, real) {
//     return -1.0;
//   }
//   real compute_mean_mu() {
//     return -1;
//   };

// };

// /* ---------------------------- */
// /* mixture of gas and radiation */
// /* ---------------------------- */

// class gas_and_radiation : public equation_of_state {
// protected:
//   /* inhertied from the base class */

// public:
//   gas_and_radiation() : equation_of_state() {
//     _gamma_gas = 5.0/3.0;
//   };
//   ~gas_and_radiation() {}
  
//   void compute_state() {
// #ifdef _ENERGY_EQUATION_
//     real rho = _density*rhounit;

//     real alpha = 1.5*k_boltz/(_mean_mu*m_u);
//     real betax = a_rad/rho;

//     real q = alpha/betax;
//     real r = -(_e_thermal*eunit) / betax;

//     real Tcrit = 1.0e6;
//     q = 3.0/2.0 * rho*k_boltz/(a_rad * (_mean_mu*m_proton) * pow(Tcrit, 3.0));
//     r = rho * (_e_thermal*eunit)/(a_rad * pow(Tcrit, 4.0));

//     real r4 = pow(r, 1.0/4.0);
//     real eps = q * r4/r;
//     if (eps < 1.0e-6)
//       _temperature = r4 * (1 - eps/4.0);
//     else
//       _temperature = calc_temp(q, -r);
    
//     _temperature *= Tcrit;
    
//     /* pressure computations */
    
//     real p_gas = rho/(_mean_mu*m_proton) * k_boltz * _temperature;
//     real p_rad = a_rad/3.0 * pow(_temperature,4.0);
//     real p_tot = p_rad + p_gas;
//     _pressure = p_tot/punit;
// #else
//     compute_e_thermal();
// #endif // _ENERGY_EQUATION_
//   }

//   real get_entropy() {
//     real rho = _density*rhounit;
//     real mu = _mean_mu * m_proton;
//     real gas_entropy = 
//       3.0/2.0 * k_boltz/mu * log(3.0/2.0 * k_boltz*_temperature/mu * pow(rho, -2.0/3.0));
//     real radiation_entropy = 4.0/3.0 * a_rad/rho * pow(_temperature, 3.0);

//     real entropy = gas_entropy + radiation_entropy;
//     return entropy;
//   }
//   real get_A() {
//     real rho = _density * rhounit;
//     real Pgas = rho/(_mean_mu*m_proton) * k_boltz * _temperature;
//     real Prad = a_rad/3.0 * pow(_temperature, 4.0);
//     real Ptot = Pgas + Prad;
//     real beta = Pgas/Ptot;
//     real A = log(Pgas) - 5.0/3.0*log(rho) + 8.0/3.0/beta - log(aunit);
//     return A;
//   }

//   real get_beta() {
//     real p_gas = _density*rhounit/(_mean_mu*m_proton) * k_boltz * _temperature/punit;
//     real beta = p_gas/_pressure;
//     return beta;
//   }

//   real get_sound_speed() {
//     real beta = get_beta();
//     real gamma1 = pow( 4.0 - 3.0 * beta, 2.0 ) * (_gamma_gas - 1.0);
//     gamma1 = beta + gamma1/( beta + 12.0 * (1.0 - beta) * (_gamma_gas - 1.0) );
//     real sound_speed = sqrt(_gamma_gas * _pressure/_density);
//     return sound_speed;
//   }

//   real compute_e_thermal() {
//     double rho = _density * rhounit;
//     double mmu = _mean_mu * m_proton;
    
//     double delta = 3.0*pow(k_boltz,4.0)/a_rad * 1.0/(pow(mmu*aunit,3.0)*mmu);
//     double rho_tilde = log(rho/delta) + 3.0 * double(_entr_var);

//     double tol = 1.0e-6;
//     double beta = find_beta(tol, rho_tilde);

//     if (abs(1.0 - beta) > tol) {
//       double temp = 3.0*rho*k_boltz/(a_rad * mmu) * (1.0-beta)/beta;
//       _temperature = pow(temp, 1.0/3.0);
      
//       double p_gas = rho/(_mean_mu*m_proton) * k_boltz * _temperature;
//       double p_rad = a_rad/3.0 * pow(_temperature,4.0);
//       _pressure = (p_rad + p_gas) / punit;
//     } else {
//       double p_tot = aunit * exp(double(_entr_var)) * pow(rho, 5.0/3.0) * exp(-8.0/3.0);
//       _temperature = p_tot * (_mean_mu*m_proton)/(rho*k_boltz);
//       _pressure = p_tot / punit;
//     }

//     _e_thermal = (3.0/2.0 * k_boltz * _temperature/(_mean_mu*m_proton) + 
// 		 a_rad * pow(_temperature, 4.0)/rho) / eunit;
//     return _e_thermal;
//   }

//   eos_type get_eos_type() {return GAS_AND_RADIATION;}

//   real compute_temperature(real density, real pressure, real mean_mu) {
//     double t0 = 1.0e6;
//     double r = - 3*pressure/(a_rad*pow(t0,4));
//     double q = + (3*density/(mean_mu*m_proton*t0))*
//       k_boltz/(a_rad*t0*t0);
//     double temp;
    
//     if (pressure < 0) {
//       real energy = -pressure;
//       r = - density*energy/a_rad * 1.0/pow(t0, 4.0);
//       q = (k_boltz/a_rad)/(5.0/3.0 - 1) * density/(mean_mu*m_proton) * 1.0/pow(t0, 3.0);
//     }
    
//     calc_temp(q, r, temp);
//     temp *= t0;
// //     PRC(density); PRC(pressure); PRL(mean_mu);
// //     PRC(r); PRC(q); PRL(temp);

//     return temp;
//   };

//   real compute_e_thermal(real density, real pressure, real mean_mu) {
// //     real temp = compute_temperature(density, pressure, mean_mu);
// //     real e_th = 3.0/2.0 * k_boltz * temp/(mean_mu*m_proton) + 
// //       a_rad * pow(temp, 4.0)/density;
//     real e_th = 1;
//     return e_th;
//   }

//   real compute_mean_mu() {
//     return -1.0;
//   }

// };

#endif // __EOS__
