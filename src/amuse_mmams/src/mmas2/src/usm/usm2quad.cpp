#include "usm.h"
#include "eos/eos.h"
#include <iomanip>

real rate_PP(real density, real temperature, real X_H) {
  real t6 = temperature / 1.0e6;
  real e_pp = 2.38e6 * density * pow(X_H, 2.0) * pow(t6, -2.0/3.0) * exp(-33.80/pow(t6,1.0/3.0));
  return e_pp;
}

real rate_CNO(real density, real temperature, real X_H, real X_CNO) {
  real t6 = temperature / 1.0e6;
  real e_CNO = 8.67e27 * X_CNO*X_H*density * pow(t6, -2.0/3.0) * exp(-152.28/pow(t6, 1.0/3.0));
  return e_CNO;
}

int main(int argc, char *argv[])
{
  cerr << "Reading .usm format from stdin" << endl;

  usm model;
  model.read(stdin, 1);

  /*
    write down star in columned format 
    suitable for using for plotting packages
  */
  
  int n = model.get_num_shells();

  real luminosity = 0.0;
  cout << n << endl;

  for (int ix = 0; ix < n; ix++) {
    int i = n-1 - ix;
    mass_shell &shell = model.get_shell(i);
    
    real dm;
    if (i < n-1) {
      dm = model.get_shell(i+1).mass - shell.mass;
    } else {
      dm = shell.mass - model.get_shell(i-1).mass;
    }

#define am(x) (1.0+Amass[x]/2.0)/Amass[x]
    real Amass[] = {1, 4, 16, 14, 12, 20, 24, 28, 56};
    real mean_mu = 2 * shell.composition.H1 + 
      am(1) * shell.composition.He4 +
      am(2) * shell.composition.O16 +
      am(3) * shell.composition.N14 +
      am(4) * shell.composition.C12 + 
      am(5) * shell.composition.Ne20 + 
      am(6) * shell.composition.Mg24;
    mean_mu = 1.0/mean_mu;
    
    shell.mean_mu = mean_mu;
    
    real temperature = compute_temperature(shell.density, shell.pressure, shell.mean_mu);
    shell.temperature = temperature;
    
    real e_thermal = 3.0/2.0 * uK * temperature/(mean_mu*uM_U) + uA_RAD * pow(temperature, 4.0)/shell.density;
    shell.e_thermal = e_thermal;

    /* compute beta parameter */
//     real Pgas = shell.density/(shell.mean_mu*uM_U) * uK * shell.temperature;
//     real beta = Pgas / (shell.pressure);
    chemical_composition che = shell.composition;
    
    real gas_entropy = 
      3.0/2.0 * uK/(shell.mean_mu*uM_P) *
      log(3.0/2.0 * uK*shell.temperature/(shell.mean_mu*uM_P) * pow(shell.density, -2.0/3.0));
    real radiation_entropy = 4.0/3.0 * uA_RAD/shell.density * pow(shell.temperature, 3.0);
//     real buoyancy =  log(Pgas) - 5.0/3.0*log(shell.density) + 8.0/3.0/beta - 8.0/3.0;

//     real e_pp  = rate_PP(shell.density, shell.temperature, che.H1);
//     real e_cno = rate_CNO(shell.density, shell.temperature, che.H1, che.O16+che.N14+che.C12);

    real cH = pow(1.0/shell.mean_mu, 5.0/2.0)/che.H1;
    cH = pow(cH, che.H1/shell.mean_mu);

    real cHe = pow(4.0/shell.mean_mu, 5.0/2.0)/(1 - che.H1);
    cHe = pow(cHe, (1-che.H1)*4/shell.mean_mu);

//     real cc = cH * cHe;
    
    real mmu = shell.mean_mu * uM_U;

    radiation_entropy = 4.0/3.0 * mmu * 
      uA_RAD * pow(shell.temperature, 3.0)/(shell.density * uK);

    gas_entropy = pow(mmu*uK*shell.temperature/uH/uH, 3.0/2.0) * mmu/shell.density;
//     PRC(mmu); PRC(shell.temperature); PRC(uH); PRC(mmu*u_k*shell.temperature/u_h); PRL(gas_entropy);
    gas_entropy = log(gas_entropy);


    real Mg24 = 1.0 - (che.H1 + che.He4 + che.O16 + che.N14 + che.C12 + che.Ne20);

    cout  << scientific << setprecision(15)
	  << i << " "                        // 1
	  << shell.mass*1.9891e33   << " "   // 2
	  << shell.radius*6.9599e10 << " "   // 3
	  << shell.temperature << " "        // 4
	  << shell.density << " "            // 5
	  << luminosity    << " "            // 6
	  << che.H1 << " "               // 7
	  << che.He4 << " "             // 8
	  << che.O16 << " "             // 9
	  << che.N14 << " "             // 10
	  << che.C12 << " "             // 11
	  << che.Ne20 << " "            // 12
	  << Mg24 << " "            // 13
// 	  << che.Si28 << " "            // 14
// 	  << che.Fe56 << " "            // 15
	 << endl;

  }

  cerr << "end-of-code" << endl;
  return 0;

}

