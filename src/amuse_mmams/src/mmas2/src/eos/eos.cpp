#include "eos.h"

inline real beta_func(real beta, real y) {
  real zeta = 5.0*log(1.0 - beta) - 8.0*log(beta) + 32.0/beta;
  zeta = y - zeta;
  return zeta;
}

real compute_beta(real pressure, real entropy, real mean_mu) {
  real eps = 1.0e-7;
  int max_iter = -int(log(eps)/log(2.0)) + 1;
  
  real beta_min = eps;
  real beta_max = 1.0 - eps;

  real delta = 3*pow(uK, 4.0)/uA_RAD;
  real y = 3.0*log(pressure) - 5*log(delta) + 12*log(entropy)
    + 20*log(mean_mu*uM_U);

  real beta = 0;
  for (int i = 0; i < max_iter; i++) {
    beta = 0.5 * (beta_min + beta_max);
    real zeta = beta_func(beta, y);
    
    if (zeta < 0)   beta_min = beta;
    else            beta_max = beta;
  }
  beta = 0.5 * (beta_min + beta_max);
  return beta;
}

double calc_temp(double q, double r) {
  double k,b,piece1,piece2;
  double y1,y2,yy,aa,b2,c2,kh;
  double x3;

  /*
    Subroutine to solve 4th order equations to determine the temperature x3
    for an equation of state with both ideal gas and radiation pressure.
    Written by Scott Fleming 10/04/02 and James Lombardi 2002-2003
    
    The fourth order equation comes from u_gas+ u_rad = u, with
    u_gas proportional to T and u_rad proportional to T^4
    
    In general, we can transform a 4th order equation to x^4+px^2+qx+r=0
    (see pages 57-58 of Stillwell's "Mathematics and its history" text)
    but we fortunately don't even have an x^2 term (that is, p=0).
    
    Follow Stillwell, we can transform this into a cubic equation:
    First solve for y by using the fact that B^2-4AC=0
    equation is then:  y^3=ry+q^2/8
    using the solution of cubic equations found in Stillwell page 55:
  */

  /* Modified E.Gaburov 4-Jan'06 */

  k=0.125*pow(q,2);
  kh=0.5*k;
  if(pow(kh,2)-pow(r/3.,3)<=0){
    PRL(k); PRL(r); PRL(pow(kh,2)-pow(r/3.,3));
    cerr<<"bad input: imaginary results?\n";
    exit(-1);
  }

  piece1=kh+pow(pow(kh,2)-pow(r/3.,3),0.5);
  piece2=kh-pow(pow(kh,2)-pow(r/3.,3),0.5);

  y1=pow(piece1,1./3.);
      
  // C++ doesn't like cube roots of neg. #'s
  y2=-pow(fabs(piece2),1./3.);

  yy=y1+y2;

  aa=2.*yy;
  b=-q;
  //c=-r+pow(y,2);

  b2=pow(aa,0.5);
  c2=0.5*b/b2+yy;

  x3=0.5*(-b2+pow(pow(b2,2)-4.*c2,0.5));

  if(piece1<0) {
    cerr << "piece 1 lt 0"<< endl;
    PRL(k); PRL(r); PRL(piece1); PRL(piece2);
  }
  if(piece1==-piece2){
    cerr<<"piece 1 eq -piece 2 (rad pressure dominates)" << endl;
    PRL(k); PRL(r); PRL(piece1); PRL(piece2);
    PRL(b2); PRL(c2); PRL(x3); PRL(pow(-r, 0.25));
    x3=pow(-r-q*pow(-r-q*pow(-r,0.25),0.25),0.25);
    PRL(x3);
  }
  if(piece2>=0) x3=-(r+pow(r/q,4))/q;
  return x3;
}

real compute_entropy(real density, real temperature, real mean_mu) {
  real Pgas = density/(mean_mu*uM_U) * uK * temperature;
  real Prad = uA_RAD/3.0 * pow(temperature, 4.0);
  real Ptot = Pgas + Prad;
  real beta = Pgas/Ptot;
  real A = log(Pgas) - 5.0/3.0*log(density) + 8.0/3.0/beta;
  A = exp(A);
  return A;
}

real compute_density(real pressure, real entropy, real mean_mu) {
  real beta    = compute_beta(pressure, entropy, mean_mu);
  real density = (beta*pressure)/entropy * exp(8.0/(3.0*beta));
  density = pow(density, 3.0/5.0);
  return density;
}

real compute_temperature(real density, real pressure, real mean_mu) {
  real t0 = 1.0e6;
  real r = - 3*pressure/uA_RAD * 1.0/pow(t0, 4.0); 
  real q = + 3*density * uK/(uA_RAD*mean_mu*uM_U) * 1.0/pow(t0, 3.0); 
  
  real temperature = t0 * calc_temp(q, r);
  return temperature;
}

real compute_pressure(real density, real energy, real mean_mu) {
  
  mean_mu *= uM_U;

  real t0 = 1.0e6;
  real r = - density*energy/uA_RAD * 1.0/pow(t0, 4.0);
  real q = (uK/uA_RAD)/(5.0/3.0 - 1) * density/mean_mu * 1.0/pow(t0, 3.0);
  
  real temperature = t0 * calc_temp(q, r);

  real p_gas = density/mean_mu * uK * temperature;
  real p_rad = uA_RAD/3.0 * pow(temperature, 4.0);
//   real beta = p_gas/(p_gas + p_rad);
  real pressure = (p_gas + p_rad);
  return pressure;
  
  
}

real compute_energy(real density, real temperature, real mean_mu) {
  real ugas = 1.5*uK/uM_U * temperature/mean_mu;
  real urad = uA_RAD*pow(temperature, 4.0)/density;
  return (ugas + urad);
}

#ifdef _TOOLBOX_

int main(int argc, char *argv[]) {
  gas_and_radiation eos;

  real rho, e_th;
  eos.density() =  rho = 5.0/rhounit;
  eos.e_thermal() = e_th = 4e+15/eunit;
  eos.mean_mu() = 1.0;

  eos.compute_state();
  real A = eos.get_A();
  PRL(A);
  PRC(eos.density()); PRC(eos.temperature()); PRC(eos.pressure()); PRL(eos.e_thermal());
  A += 0.0;
  eos.compute_e_thermal(A);
  PRC(eos.density()); PRC(eos.temperature()); PRC(eos.pressure()); PRL(eos.e_thermal());
  eos.compute_state();
  PRC(eos.density()); PRC(eos.temperature()); PRC(eos.pressure()); PRL(eos.e_thermal());
  PRL(eos.get_A());
  
  return 0;
}
#endif // _TOOLBOX_

