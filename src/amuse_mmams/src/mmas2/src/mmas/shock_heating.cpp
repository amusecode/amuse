#include "mmas.h"
#include "eos/eos.h"
#include <gsl/gsl_interp.h>
#include <gsl/gsl_errno.h>

// int mmas::shock_heating(usm &model, real strength, real da_const) {
//   for (int i = 0; i < model.get_num_shells(); i++) {
//     mass_shell &shell = model.get_shell(i);
//     shell.entropy  = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
//     real de  = pow(10, strength) * pow(shell.pressure, da_const);
//     shell.entropy += pow(10, strength) * pow(shell.pressure, da_const);
//   }
// }

extern bool error_occurred;

int mmas::shock_heating_4(usm& model, real a, real b, real c, real d) {
  real density_max = 0, pressure_max = 0, entropy_min = 1e300;
  for (int i = 0; i < model.get_num_shells(); i++) {
    mass_shell &shell = model.get_shell(i);
    shell.entropy  = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
    density_max = max(density_max, shell.density);
    pressure_max = max(pressure_max, shell.pressure);
    entropy_min =  min(entropy_min, shell.entropy);
  }

  for (int i = 0; i < model.get_num_shells(); i++) {
    mass_shell &shell = model.get_shell(i);
    shell.entropy  = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
//     real x = log(shell.pressure/pressure_max)/log(10.0);
//     real de = 10**(a+b*b*x) + 10**(c + d*d*x))/log(10.0);;
    real de = 1;
    if (de > 100) de = 100;
    shell.entropy *= 1 + pow(10.0, de);
  }
  return 0;
}


int mmas::shock_heating_4() {
  /* TAMS */
//   shock_heating_3(*model_a, +0.62, +0.50, +0.34);   /* 80 */
//   shock_heating_3(*model_b, -0.49, -0.58, +0.00);   /* 8 */

  shock_heating_3(*model_a, +0.67, +0.34, -0.29);   /* 40 */
  shock_heating_3(*model_b, -0.70, -0.46, -0.16);   /* 8 */
  
  return 0;
}
int mmas::shock_heating_3(usm& model, real a, real b, real c) {
  real density_max = 0, pressure_max = 0, entropy_min = 1e300;
  for (int i = 0; i < model.get_num_shells(); i++) {
    mass_shell &shell = model.get_shell(i);
    shell.entropy  = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
    density_max = max(density_max, shell.density);
    pressure_max = max(pressure_max, shell.pressure);
    entropy_min =  min(entropy_min, shell.entropy);
  }

  for (int i = 0; i < model.get_num_shells(); i++) {
    mass_shell &shell = model.get_shell(i);
    shell.entropy  = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
    real x = log(shell.pressure/pressure_max)/log(10.0);
    real de = a+b*x + c*c * x*x;
    if (de > 100) de = 100;
    shell.entropy *= 1 + pow(10.0, de);
  }
  return 0;
}


int mmas::shock_heating_3() {
  /* TAMS */
//   shock_heating_3(*model_a, +0.62, +0.50, +0.34);   /* 80 */
//   shock_heating_3(*model_b, -0.49, -0.58, +0.00);   /* 8 */

  shock_heating_3(*model_a, +0.67, +0.34, -0.29);   /* 40 */
  shock_heating_3(*model_b, -0.70, -0.46, -0.16);   /* 8 */
  
  return 0;
}

int mmas::shock_heating(usm &model, real a, real b) {
  real density_max = 0, pressure_max = 0, entropy_min = 1e300;
  for (int i = 0; i < model.get_num_shells(); i++) {
    mass_shell &shell = model.get_shell(i);
    shell.entropy  = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
    density_max = max(density_max, shell.density);
    pressure_max = max(pressure_max, shell.pressure);
    entropy_min =  min(entropy_min, shell.entropy);
  }

  for (int i = 0; i < model.get_num_shells(); i++) {
    mass_shell &shell = model.get_shell(i);
    shell.entropy  = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
    real de = a + b*log(shell.pressure/pressure_max)/log(10.0);
    if (de > 100) de = 100;
    shell.entropy += entropy_min * pow(10.0, de);
//     shell.entropy *= 1.0 +  pow(10.0, de);
  }

  return 0;
}

int mmas::shock_heating(real ff) {

  struct fit_params{
    real a, b, c, d;
  };
  
  int status, n = 4;
  double tams_dm[]         = {8.3, 8.9, 5.0, 1.9};
  double tams_q[]          = {0.8, 0.4, 0.2, 0.1};
  fit_params tams_params[] = { {-0.12, -0.63, -1.025, -0.77},          // 10+8
			       {0.25,  -0.37, -1.32,  -0.95},          // 20+8
			       {0.35,  -0.22, -0.84,  -0.77},          // 40+8
			       {0.22,  -0.17, -0.59,  -0.80} };        // 80+8


  double hams_dm[]         = {6.8, 4.7, 2.1, 0.8};
  double hams_q[]          = {0.8, 0.4, 0.2, 0.1};
  fit_params hams_params[] = { {-0.48, -0.79, -0.89, -0.85},          // 10+8
			       {-0.22, -0.66, -0.88, -0.79},          // 20+8
			       {-0.09, -0.44, -0.77, -0.85},          // 40+8
			       {-0.13, -0.28, -0.80, -1.03} };        // 80+8

  double *dm_cur, *q_cur;
  fit_params *fit_cur;
  if (model_a->star_age == 0) {
    dm_cur  = hams_dm;
    q_cur   = hams_q;
    fit_cur = hams_params;
  } else {
    dm_cur  = tams_dm;
    q_cur   = tams_q;
    fit_cur = tams_params;
  }

  real q = model_b->star_mass/model_a->star_mass;
  if (q > 1) {
    cerr << "Something gone wrong ... q = " << q << " > 1 \n";
    exit(-1);
  }
  q = min(max(q, 0.05), 1.0);
//   PRL(q);

  double *xq = new double[n];
  double *xa = new double[n];
  double *xb = new double[n];
  double *xc = new double[n];
  double *xd = new double[n];
  for (int j = 0; j < n; j++) {
    int i = n-1 - j;
    xq[j] = q_cur[i];
    xa[j] = fit_cur[i].a;
    xb[j] = fit_cur[i].b;
    xc[j] = fit_cur[i].c;
    xd[j] = fit_cur[i].d;
  }

  /* get a */

  gsl_interp *interpolation = gsl_interp_alloc (gsl_interp_linear,n);
  status = gsl_interp_init(interpolation, xq, xa, n);
  if (status != GSL_SUCCESS) {
      cerr << "Error in gsl_interp_init for parameter a, in mmas::shock_heating(real ff) ";
      PRL(status);
      return -1;
  }
  gsl_interp_accel * accelerator =  gsl_interp_accel_alloc();
  real a = gsl_interp_eval(interpolation, xq, xa, q, accelerator);
  if (error_occurred) {
      cerr << "Error in gsl_interp_eval for parameter a, in mmas::shock_heating(real ff)" << endl;
      error_occurred = false;
      return -1;
  }
  gsl_interp_accel_free(accelerator);
  gsl_interp_free(interpolation);
//   PRL(a);

  /* get b */

  interpolation = gsl_interp_alloc (gsl_interp_linear,n);
  status = gsl_interp_init(interpolation, xq, xb, n);
  if (status != GSL_SUCCESS) {
      cerr << "Error in gsl_interp_init for parameter b, in mmas::shock_heating(real ff) ";
      PRL(status);
      return -1;
  }
  accelerator =  gsl_interp_accel_alloc();
  real b = gsl_interp_eval(interpolation, xq, xb, q, accelerator);
  if (error_occurred) {
      cerr << "Error in gsl_interp_eval for parameter b, in mmas::shock_heating(real ff)" << endl;
      error_occurred = false;
      return -1;
  }
  gsl_interp_accel_free(accelerator);
  gsl_interp_free(interpolation);
//   PRL(b);
 
  /* get c */

  interpolation = gsl_interp_alloc (gsl_interp_linear,n);
  status = gsl_interp_init(interpolation, xq, xc, n);
  if (status != GSL_SUCCESS) {
      cerr << "Error in gsl_interp_init for parameter c, in mmas::shock_heating(real ff) ";
      PRL(status);
      return -1;
  }
  accelerator =  gsl_interp_accel_alloc();
  real c = gsl_interp_eval(interpolation, xq, xc, q, accelerator);
  if (error_occurred) {
      cerr << "Error in gsl_interp_eval for parameter c, in mmas::shock_heating(real ff)" << endl;
      error_occurred = false;
      return -1;
  }
  gsl_interp_accel_free(accelerator);
  gsl_interp_free(interpolation);
//   PRL(c);

  /* get d */

  interpolation = gsl_interp_alloc (gsl_interp_linear,n);
  status = gsl_interp_init(interpolation, xq, xd, n);
  if (status != GSL_SUCCESS) {
      cerr << "Error in gsl_interp_init for parameter d, in mmas::shock_heating(real ff) ";
      PRL(status);
      return -1;
  }
  accelerator =  gsl_interp_accel_alloc();
  real d = gsl_interp_eval(interpolation, xq, xd, q, accelerator);
  if (error_occurred) {
      cerr << "Error in gsl_interp_eval for parameter d, in mmas::shock_heating(real ff)" << endl;
      error_occurred = false;
      return -1;
  }
  gsl_interp_accel_free(accelerator);
  gsl_interp_free(interpolation);
//   PRL(d);
 
  a += log(ff)/log(10.0);
  c += log(ff)/log(10.0);
  shock_heating(*model_a, a, b);  
  shock_heating(*model_b, c, d);  

  return 0;
}
