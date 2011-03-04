#include "mmas.h"
#include "eos/eos.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_integration.h>

inline real gij(real mi, real mj, real M, real alpha) {
  real c = alpha/(M*M);
  real gij1 = exp(-c*square(mi-mj));
  real gij2 = exp(-c*square(mi+mj));
  real gij3 = exp(-c*square(mi+mj - 2*M));
  return gij1 + gij2 + gij3;
}

void mmas::mixing(usm &model) {
  int n = model.get_num_shells();
  
  real entr_min = 1e100, entr_max = 0;
  vector<real> dm;
  real m_prev = 0;
  for (int i = 0; i < n; i++) {
    mass_shell &shell  = model.get_shell(i);
    entr_min = min(entr_min, shell.entropy);
    entr_max = max(entr_max, shell.entropy);
    dm.push_back(shell.mass - m_prev);
    m_prev = shell.mass;
  }

  real M = model.star_mass;

  real alpha = 10 * square(log(entr_max/entr_min)/log(10));
  PRL(alpha);
  vector<real> m_mu;
  for (int k = 0; k < n; k++) {
    mass_shell &shell_k = model.get_shell(k);
    real X_H = 0.0;
    for (int i = 0; i < n; i++) {
      mass_shell &shell_i = model.get_shell(i);
      real fac = 0;
      for (int j = 0; j < n; j++) {
	mass_shell &shell_j = model.get_shell(j);
	fac += gij(shell_i.mass, shell_j.mass, M, alpha)*dm[j];
      }
      real xH = (4.0/shell_i.mean_mu - 3.0)/5.0;
      X_H += xH * gij(shell_i.mass, shell_k.mass, M, alpha) * dm[i] / fac;
    }
    real mmu = 4.0/(5*X_H + 3);
    m_mu.push_back(mmu);
  }

  for (int i = 0; i < n; i++) {
    mass_shell &shell = model.get_shell(i);
//     PRC(shell.mass); PRC(shell.mean_mu); PRL(m_mu[i]);
    shell.mean_mu = m_mu[i];
  }

}

void mmas::mixing() {

  cerr << "Mixing star A\n";
  mixing(*model_a);

  cerr << "Mixing star B\n";
  mixing(*model_b);

  cerr << "done mixing \n";
  
  
}


void mmas::mixing_product(int n) {
  delete mixed_product;
  mixed_product = new usm;

  real dm_bin = product->star_mass / n;
  
  real m_prev = 0, m_bin = 0;
  real Mtot = 0, Utot = 0, Vtot = 0, r_mean = 0;
  int  n_in_bin = 0;
  int  j = 0;
  real H1tot = 0, He4tot = 0, O16tot = 0, N14tot = 0, C12tot = 0, Ne20tot = 0, Mg24tot = 0, Si28tot = 0, Fe56tot = 0;
  
  
  for (int i = 0; i < product->get_num_shells(); i++) {
    mass_shell &shell_i = product->get_shell(i);
    real dm = shell_i.mass - m_prev;
    m_prev  = shell_i.mass;

    r_mean += shell_i.radius;
    n_in_bin += 1;
 
    m_bin += dm;
    Mtot  += dm;
    Vtot  += dm/shell_i.density;
    H1tot   += dm*shell_i.composition.H1;
    He4tot  += dm*shell_i.composition.He4;
    O16tot  += dm*shell_i.composition.O16;
    N14tot  += dm*shell_i.composition.N14;
    C12tot  += dm*shell_i.composition.C12;
    Ne20tot += dm*shell_i.composition.Ne20;
    Mg24tot += dm*shell_i.composition.Mg24;
    Si28tot += dm*shell_i.composition.Si28;
    Fe56tot += dm*shell_i.composition.Fe56;
    Utot  += compute_energy(shell_i.density, shell_i.temperature, shell_i.mean_mu) * dm;

    if (m_bin > dm_bin) {
//       PRC(j); PRC(n); PRC(Mtot); PRC(m_bin); PRC(dm_bin); PRL(n_in_bin);

      mass_shell shell_j; j++;
//       mass_shell &shell_j = mixed_product->get_shell(j++);
      shell_j.radius      = r_mean/n_in_bin;
      shell_j.mass        = Mtot;
      shell_j.density     = m_bin/Vtot;
      shell_j.composition.H1   = H1tot/m_bin;
      shell_j.composition.He4  = He4tot/m_bin;
      shell_j.composition.O16  = O16tot/m_bin;
      shell_j.composition.N14  = N14tot/m_bin;
      shell_j.composition.C12  = C12tot/m_bin;
      shell_j.composition.Ne20 = Ne20tot/m_bin;
      shell_j.composition.Mg24 = Mg24tot/m_bin;
      shell_j.composition.Si28 = Si28tot/m_bin;
      shell_j.composition.Fe56 = Fe56tot/m_bin;

#define am(x) (1.0+Amass[x]/2.0)/Amass[x]
      real Amass[] = {1, 4, 16, 14, 12, 20, 24, 28, 56};
      shell_j.mean_mu = 2 * shell_j.composition.H1 + 
	am(1) * shell_j.composition.He4 +
	am(2) * shell_j.composition.O16 +
	am(3) * shell_j.composition.N14 +
	am(4) * shell_j.composition.C12 + 
	am(5) * shell_j.composition.Ne20 + 
	am(6) * shell_j.composition.Mg24 +
	am(7) * shell_j.composition.Si28 + 
	am(8) * shell_j.composition.Fe56;
      shell_j.mean_mu = 1.0/shell_j.mean_mu;

      shell_j.e_thermal   = Utot/m_bin;
      shell_j.pressure    = compute_pressure(shell_j.density, shell_j.e_thermal, shell_j.mean_mu);
      shell_j.temperature = compute_temperature(shell_j.density, shell_j.pressure, shell_j.mean_mu);
      shell_j.entropy     = compute_entropy(shell_j.density, shell_j.temperature, shell_j.mean_mu);
      mixed_product->add_shell(shell_j);

      m_bin -= dm_bin;
      m_bin = Utot = Vtot = r_mean = 0;
      H1tot = He4tot = O16tot = N14tot = C12tot = Ne20tot = Mg24tot = Si28tot = Fe56tot = 0;
      n_in_bin = 0;
    }
  }

  mixed_product->build_hashtable();
}

/* =================== */


inline real smoothing_kernel(real x, real h) {
  real q = fabs(x)/h;
  real wk = 0.0;
  if (q < 1.0) {
    wk = 1 - 3.0/2 * q*q + 3.0/4 * q*q*q;
  } else if (q < 2) {
    q = 2.0 - q;
    wk = 1.0/4 * q*q*q;
  }
  wk *= 2.0/3.0;
//   wk = 0.25; if (q > 2) wk = 0;
  return wk/h;
}

struct smoothing_params {
  gsl_interp       *int_x, *int_y;
  gsl_interp_accel *acc_x, *acc_y;
  double           *arr_x, *arr_y, *smoothed_y;
  double            h, x;
  int               n;
};

double smoothing_integrand(double x, void *params) {
  struct smoothing_params *p = (struct smoothing_params* )params;
  double y;
  if      (x < 0)                   y = p->arr_y[0];
  else if (x >= p->arr_x[p->n-1])   y = p->arr_y[p->n-1];
  else                              y = gsl_interp_eval(p->int_y, p->arr_x, p->arr_y, x, p->acc_y);
  
  real   wk = smoothing_kernel(p->x - x, p->h);
  return y * wk;
}

void smoothing_integrate(smoothing_params &params, int n) {
  int n_int = 2*n;

  params.acc_y = gsl_interp_accel_alloc();
  params.int_y = gsl_interp_alloc(gsl_interp_linear, n);
  params.n     = n;
  gsl_interp_init(params.int_y, params.arr_x, params.arr_y, n);

  /* setup integration workspace */

  gsl_integration_workspace *w = gsl_integration_workspace_alloc(n_int);
  gsl_function F; 
  F.function = &smoothing_integrand;
  F.params   = &params;
  
  double eps_abs = 0.0, eps_rel = 0.01;
  double result, error;
  gsl_set_error_handler_off();
 
  for (int i = 0; i < n; i++) {
    params.x = params.arr_x[i];
    params.h = 0.1; //0.1*params.x + 1.0e-4;

    gsl_integration_qag(&F, 
			params.x - 2*params.h, params.x + 2*params.h,
			eps_abs, eps_rel,
			n_int, 1,
			w, &result, &error);
    params.smoothed_y[i] = result;
  }
  
  gsl_integration_workspace_free(w);
  gsl_interp_accel_free(params.acc_y);
  gsl_interp_free(params.int_y);
}

void mmas::smooth_product() {
  int n_shells = product->get_num_shells();   /* number of shells in the product */

  smoothing_params params;

  params.arr_x      = new double[n_shells];
  params.arr_y      = new double[n_shells];
  params.smoothed_y = new double[n_shells];

  /* composition */
  cerr << "Smoothing composition\n";
  for (int i = 0; i < n_shells; i++) {
    mass_shell &shell = product->get_shell(i);
    params.arr_x[i] = shell.radius;
    params.arr_x[i] = shell.mass;
    params.arr_y[i] = (4.0/shell.mean_mu - 3.0)/5.0;
  }
  smoothing_integrate(params, n_shells);
  for (int i = 0; i < n_shells; i++) {
    mass_shell &shell = product->get_shell(i);
    shell.mean_mu = 4.0/(5.0*params.smoothed_y[i] + 3);
  }

  /* thermal energy */
  cerr << "Smoothing thermal energy\n";
  for (int i = 0; i < n_shells; i++) {
    mass_shell &shell = product->get_shell(i);
    params.arr_y[i] = compute_energy(shell.density, shell.temperature, shell.mean_mu);
  }
  smoothing_integrate(params, n_shells);
  for (int i = 0; i < n_shells; i++) {
    mass_shell &shell = product->get_shell(i);
    shell.e_thermal =  params.smoothed_y[i];
//     real x = params.arr_x[i];
//     real y = params.arr_y[i];
//     real ys = params.smoothed_y[i];
//     PRC(x); PRC(y); PRL(ys);
  }

  /* density */
//   cerr << "Smoothing density\n";
//   for (int i = 0; i < n_shells; i++) {
//     mass_shell &shell = product->get_shell(i);
//     params.arr_y[i] = shell.density;
//   }
//   smoothing_integrate(params, n_shells);
//   for (int i = 0; i < n_shells; i++) {
//     mass_shell &shell = product->get_shell(i);
//     shell.density = params.smoothed_y[i];
//   }
  
  for (int i = 0; i < n_shells; i++) {
    mass_shell &shell = product->get_shell(i);
    shell.pressure    = compute_pressure(shell.density, shell.e_thermal, shell.mean_mu);
    shell.temperature = compute_temperature(shell.density, shell.pressure, shell.mean_mu);
    shell.entropy     = compute_entropy(shell.density, shell.temperature, shell.mean_mu);
  }

  delete[] params.arr_x;
  delete[] params.arr_y;
  delete[] params.smoothed_y;

}
