#include "mmas.h"
#include "usm/usm.h"


#define _MACOSX_
#ifdef _MACOSX_
// extern "C" void u_fpu_setup(void);
#else
#include <fenv.h>
#endif // _MACOSX_


void is_file(char *fn) {
  ifstream fin(fn, ifstream::in);
  fin.close();
  if (fin.fail() != 0) {
    cerr << "File \"" << fn << "\" does not seem to exist ! " << endl;
    exit(-1);
  };
}

int main(int argc, char *argv[]) {

#ifdef _MACOSX_
  // trap FPE on Mac OS X
// u_fpu_setup();
#else
  // trap FPE on Linux PC
  feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

  cerr << endl << "   exec line:   ";
  for (int i = 0; i < argc; i++) cerr << argv[i] << " ";
  cerr << endl << endl;
  
  usm model_a, model_b;
  real r_p   = 0.0;
  real v_inf = 0.0;
  FILE *fn_model_a = NULL;
  FILE *fn_model_b = NULL;
  int dump_mixed = 0;
  
  extern int  pgetopt(int argc, char ** argv,  char * optstr);
  extern char *poptarg;
  extern char *poparr[];
  char* param_string = "hm::xf:n:";
  int ic;
  int n_shells = int(1e4);
  real ff = 1;
  while ((ic = pgetopt(argc, argv, param_string)) != -1) {
    switch(ic) {
    case 'h':
      cerr << " Make Me a [Massive] Star \n";
      cerr << endl;
      cerr << "Usage: setup_two_stars OPTIONS > OUTPUT\n";
      cerr << endl;
      cerr << "Options: \n";
      cerr << "            -m    pass two .usm models  \n";
      cerr << "            -x    generate mixed product [no] \n";
      cerr << "            -f    param for shock heating [1] \n";
      cerr << "            -n    number of shells in final product [1e5] \n";
//       cerr << "            -p    r_p [0] \n";
//       cerr << "            -v    v_infty [0] \n";
	cerr << endl;
      return 0;
      break;
    case 'm':
      is_file(poparr[0]);
      fn_model_a = fopen(poparr[0], "r");
      is_file(poparr[1]);
      fn_model_b = fopen(poparr[1], "r");
      break;
    case 'x':
      dump_mixed = 1;
      break;
    case 'f':
      ff = atof(poptarg);
      break;
    case 'n':
      n_shells = atoi(poptarg);
      break;
    }
  }
  
  if (fn_model_a == NULL || fn_model_b == NULL) {
    cerr << "no .usm models have been provided ... quit \n";
    exit(-1);
  }
  
  model_a.read(fn_model_a, 1);
  model_b.read(fn_model_b, 1);

  mmas merger(model_a, model_b, r_p, v_inf);
  merger.merge_stars_consistently(n_shells, 1);

  
//   real rh1 = merger.get_lagrad(model_a, 0.5);
//   real rh2 = merger.get_lagrad(model_b, 0.5);
//   PRC(rh1); PRL(rh2);
  
//   real rf1 = merger.get_lagrad(model_a, 0.85);
//   real rf2 = merger.get_lagrad(model_b, 0.85);
//   PRC(rf1); PRL(rf2);

//   real rho_1 = model_a.star_mass/pow(model_a.star_radius, 3.0)/4.0/PI*3*5.89994;
//   real rho1c = model_a.get_shell(0).density;
//   real rho_2 = model_b.star_mass/pow(model_b.star_radius, 3.0)/4.0/PI*3*5.89994;
//   real rho2c = model_b.get_shell(0).density;
//   PRC(log(rho_1/rho1c)/log(10.0)); PRL(log(rho_2/rho2c)/log(10));
  


//    merger.mixing();
//    merger.shock_heating_3();
//    merger.shock_heating_4();

//   PRL(merger.mass_loss());
//   merger.shock_heating(ff);
//   merger.merge_stars(1.0, n_shells); 

  real energy_a = merger.compute_stellar_energy(merger.get_model_a());
  real energy_b = merger.compute_stellar_energy(merger.get_model_b());
  real energy_p = merger.compute_stellar_energy(merger.get_product());
  PRL(energy_a);
  PRL(energy_b);
  PRC(energy_a+energy_b); PRL(energy_p);
  
  real de = (energy_p - (energy_a + energy_b))/(energy_a + energy_b);
  de *= 1.0/(merger.mass_loss()/100.0);
  PRL(de);
  
  cerr << "Smoothing \n";
//   merger.smooth_product();
  merger.mixing_product(200);
  if (!dump_mixed) {
    cerr << "Dumping unmixed product in stdout \n";
    merger.get_product().write(stdout);
  } else {
    cerr << "Dumping mixed product in stdout \n";
    merger.get_mixed_product().write(stdout);
  }
  
  cerr << " --- \n";
  cerr << "It seems that I've successfully merged two stars. Enjoy the product ;). " << endl;
  cerr << " --- \n";
  return 0;
}
