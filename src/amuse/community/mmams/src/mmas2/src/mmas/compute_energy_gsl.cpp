#include "mmas.h"
#include "include/units.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>

gsl_interp       *interp_radius;
gsl_interp_accel *acc_radius;
gsl_interp       *interp_temp, *interp_m_mu, *interp_dens;
gsl_interp_accel *acc_temp,    *acc_m_mu,    *acc_dens;

int n_shells;
double *arr_mass, *arr_radius;
double *arr_temp, *arr_dens, *arr_m_mu;

double integrand(double mass, void *params) {
  double radius = gsl_interp_eval(interp_radius, arr_mass, arr_radius, mass, acc_radius);
  double temp   = gsl_interp_eval(interp_temp, arr_mass, arr_temp, mass, acc_temp);
  double dens   = gsl_interp_eval(interp_dens, arr_mass, arr_dens, mass, acc_dens);
  double m_mu   = gsl_interp_eval(interp_m_mu, arr_mass, arr_m_mu, mass, acc_m_mu);

  real f = -mass/radius;
  real de_therm = 1.5*uK*temp/m_mu/uM_U + 
    uA_RAD*pow(temp, 4.0)/dens;
  de_therm *= 1.0/(uG*uMSUN/uRSUN);

  return (f + de_therm);
}

real mmas::compute_stellar_energy(usm &model) {
  n_shells = model.get_num_shells();

  arr_mass   = new double[n_shells];
  arr_radius = new double[n_shells];
  arr_temp   = new double[n_shells];
  arr_dens   = new double[n_shells];
  arr_m_mu   = new double[n_shells];
  for (int i = 0; i < n_shells; i++) {
    mass_shell &shell = model.get_shell(i);
    if (i > 0 && shell.mass <= arr_mass[i-1]){
        n_shells = i;
        break;
    }
    arr_mass[i]   = shell.mass;
    arr_radius[i] = shell.radius;
    arr_temp[i]   = shell.temperature;
    arr_dens[i]   = shell.density;
    arr_m_mu[i]   = shell.mean_mu;
//     fprintf(stderr, "m= %lg: r= %lg, t= %lg, rho= %lg, mu= %lg\n",
// 	    shell.mass, shell.radius, shell.temperature, shell.density, shell.mean_mu);
  }

//   fprintf(stderr,"radius\n");
  acc_radius    = gsl_interp_accel_alloc();
  interp_radius = gsl_interp_alloc(gsl_interp_linear, n_shells);
  int status = gsl_interp_init(interp_radius, arr_mass, arr_radius, n_shells);
  if (status != GSL_SUCCESS) return -1e99;
  
//   fprintf(stderr,"temp\n");
  acc_temp    = gsl_interp_accel_alloc();
  interp_temp = gsl_interp_alloc(gsl_interp_linear, n_shells);
  gsl_interp_init(interp_temp, arr_mass, arr_temp, n_shells);

//   fprintf(stderr,"dens\n");
  acc_dens    = gsl_interp_accel_alloc();
  interp_dens = gsl_interp_alloc(gsl_interp_linear, n_shells);
  gsl_interp_init(interp_dens, arr_mass, arr_dens, n_shells);
 
//   fprintf(stderr,"mu\n");
  acc_m_mu    = gsl_interp_accel_alloc();
  interp_m_mu = gsl_interp_alloc(gsl_interp_linear, n_shells);
  gsl_interp_init(interp_m_mu, arr_mass, arr_m_mu, n_shells);


  int n_limit = 2*n_shells;

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(n_limit);
  double result, error;
  
  gsl_function F;
  F.function = &integrand;
  
  double eps_rel = 1.0e-2;
  double eps_abs = 0.0;
  double m_0 = 0.0;
  double m_1 = model.get_shell(n_shells-1).mass;

  gsl_integration_qag(&F, 
		      m_0, m_1,
		      eps_abs, eps_rel,
		      n_limit, 1,
		      w, &result, &error);
  
  gsl_integration_workspace_free(w);
		       
  gsl_interp_accel_free(acc_radius);
  gsl_interp_free(interp_radius);
  gsl_interp_accel_free(acc_temp);
  gsl_interp_free(interp_temp);
  gsl_interp_accel_free(acc_dens);
  gsl_interp_free(interp_dens);
  gsl_interp_accel_free(acc_m_mu);
  gsl_interp_free(interp_m_mu);
  delete arr_mass;
  delete arr_radius;
  delete arr_temp;
  delete arr_dens;
  delete arr_m_mu;
  
  
  return result;
}
