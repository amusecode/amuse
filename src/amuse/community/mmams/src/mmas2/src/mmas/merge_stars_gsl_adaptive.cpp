#include "mmas.h"
#include "eos/eos.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include "gsl/gslSpline.h"
#include "gsl/gslInterp.h"

#define ERROR_TRY0 1.0e-6
#define NTRY 3000
enum {NEXCEED=10, FAILED, NEGATIVE_PRESSURE};
#define ERROR_SCALE_INC 3.9
#define ERROR_SCALE_DEC 2.0

#define NUMBER_OF_CHEMICAL_SPECIES     9

#define randx (1.0)
//#define randx (1.0+1.0e-4*drand48())

/* This routine builds the merger product of the two parent stars */

/* 
 #
 #  (a)  dP/dm = - 1.0/(4*pi) * Gm/r^4
 #  (b)  dr/dm = 1.0/(4*pi*r^2*rho)
 #
 #  Eggleton'71 ->  mu = m^(2/3)
 #                   x = r^2
 #
 #  (c)  dP/dmu = - 3*G/(8*pi) * (mu/x)^2
 #  (d)  dx/dmu =   3/(4*pi) * (mu/x)^(1/2)/rho
 # 
 #  solving (c) & (d)
 #  with the following boundary conditions
 #  P(m = Mstar) = g/chi ~ 0
 #  r(m = 0)     = 0
 #
 #  (c) E.Gaburov, 2007, Sterrenkundig Instituut, Universiteit van Amsterdam
 #
 #
 #  Limitations: No angular momentum is taken into account. Only non-rotating products
 #
 # TODO LIST:      
 #   Rotating stars
 #
*/

double p_unit, rho_unit;

struct hse_func_params {
  gslInterp *entr_A, *entr_B;
  gslInterp *mmu_A,  *mmu_B;
  vector<double> *Morig_temp, *eggleton_mu;
  double mcur_A, mcur_B;
  double m_max_A, m_max_B;
  int negative_pressure;
  int which_star;
};

int hse_func(double mu, const double y[], double dydmu[], void *params) {

  hse_func_params *p = (hse_func_params*)params;

  double pressure = y[0];
  double x        = y[1]; 

  if (pressure < 0) {
    p->negative_pressure = 1;
    return GSL_FAILURE;
  }

  double dm = pow(mu, 1.5);
  if (p->eggleton_mu->size() > 0) 
    dm -= pow(p->eggleton_mu->back(), 1.5);
  
  double mc_at = p->mcur_A + dm; 
  mc_at = min(mc_at, p->m_max_A);
  double rho_A = compute_density(pressure*p_unit, 
				 p->entr_A->eval(mc_at), 
				 p->mmu_A->eval(mc_at))/rho_unit;
  
  double mc_bt = p->mcur_B + dm;
  mc_bt = min(mc_bt, p->m_max_B);
  double rho_B = compute_density(pressure*p_unit, 
				 p->entr_B->eval(mc_bt), 
				 p->mmu_B->eval(mc_bt))/rho_unit;

  if (p->mcur_A == p->m_max_A) {
    rho_A = -1;
  } else if (p->mcur_B == p->m_max_B) {
    rho_B = -1;
  }

  double density;
  
  if (rho_A > rho_B) {
    density = rho_A;
    p->Morig_temp->push_back(mc_at);     /* track location of the shell */
    p->which_star = 1;
  } else { 
    density = rho_B;
    p->Morig_temp->push_back(-mc_bt);    /* track location of the shell */
    p->which_star = 2;
  }
  
  double mu_x = pow(4.0*PI/3 * density, 2.0/3.0);
  if (x > 0) mu_x = mu/x;
  
  // dP/dmu
  dydmu[0] = -3.0/(8*PI) * mu_x*mu_x;

  // dx/dmu
  dydmu[1] = +3.0/(4*PI) * sqrt(mu_x)/density;

  return GSL_SUCCESS;
}

int hse_jac(double mu, const double y[], double *dfdy, double dfdt[], void *params) {
  return GSL_SUCCESS;
}

int solve_HSE(double p_centre, double m_product,
	      gslInterp &mmu_A,  gslInterp &mmu_B,
	      gslInterp &entr_A, gslInterp &entr_B,
	      double m_min_A, double m_min_B,
	      double m_max_A, double m_max_B,
	      vector<double> &Morig, 
	      vector<double> &mass, vector<double> &radius, vector<double> &pressure,
	      double relative_error, int n_shells) {

  int n_dim = 2;

  double eps_abs = 1.0e-4;
  double eps_rel = 0;
  eps_abs = 0; 
  eps_rel = relative_error;

  const gsl_odeiv_step_type *T = gsl_odeiv_step_rkf45;
  gsl_odeiv_step            *s = gsl_odeiv_step_alloc(T, n_dim);
  gsl_odeiv_control         *c = gsl_odeiv_control_y_new(eps_abs, eps_rel);
  gsl_odeiv_evolve          *e = gsl_odeiv_evolve_alloc(n_dim);

  hse_func_params params;

  vector<double> eggleton_mu;			
  eggleton_mu.push_back(0);
  params.eggleton_mu = &eggleton_mu;

  params.mcur_A  = m_min_A;
  params.m_max_A = m_max_A;
  params.mmu_A   = &mmu_A;
  params.entr_A  = &entr_A;
  
  params.mcur_B  = m_min_B;
  params.m_max_B = m_max_B;
  params.mmu_B   = &mmu_B;
  params.entr_B  = &entr_B;
  
  gsl_odeiv_system sys = {hse_func, hse_jac, n_dim, &params};

  double m_c = 0.0;
  double m_1 = pow(m_product, 2.0/3.0);
  double dm  = 1.0e-6;
  double y[2] = {p_centre/p_unit, 0.0};   // {pressure, x}
  
  Morig.clear();
  mass.clear();
  pressure.clear();
  radius.clear();
  
  params.negative_pressure = 0;
  while (m_c < m_1) {
    vector<double> Morig_temp;
    params.Morig_temp = &Morig_temp;
    
    params.which_star = 0;
    int status  = gsl_odeiv_evolve_apply(e, c, s,
					 &sys,
					 &m_c, m_1,
					 &dm, y);
    double dm = pow(m_c, 1.5);
    if  (eggleton_mu.size() > 0) {
      dm -= pow(eggleton_mu.back(), 1.5);
      eggleton_mu.push_back(m_c);
    } else {
      eggleton_mu.push_back(m_c);
    }


    if (params.which_star == 1) {
      params.mcur_A += dm;
      params.mcur_A = min(params.mcur_A, m_max_A);
    } else {
      params.mcur_B += dm;
      params.mcur_B = min(params.mcur_B, m_max_B);
    }    
    real mpr = pow(m_c, 1.5);
    real ma = params.mcur_A;
    real mb = params.mcur_B;
//     PRC(dm); PRC(mpr); PRC(ma+mb); PRC(ma); PRC(m_max_A); PRC(mb); PRL(m_max_B); 
  
    if (params.negative_pressure == 1)  {
      pressure.push_back(-p_unit);
      return NEGATIVE_PRESSURE;
    }
    

    /* if more than desired number of shells
       was generated, reduce relative accuracy of the solver
    */
    if (mass.size() + 1 > n_shells*3) {
      return NEXCEED;
    }

    if (status != GSL_SUCCESS)
      return FAILED;
    
    mass.push_back(pow(m_c, 1.5));
    pressure.push_back(y[0]*p_unit);
    radius.push_back(sqrt(y[1]));
    Morig.push_back(Morig_temp.back());

  }
  return 0;
}

struct merge_stars_params {
  mmas *mmas_ptr;
  vector<double> Morig, mass, radius, pressure;
  double m_product, error_try, initial_error_try;
  double p_old;
  gslInterp *mmu_A, *entr_A;
  gslInterp *mmu_B, *entr_B;
  double m_min_A, m_min_B;
  double m_max_A, m_max_B;
  int n_shells;
  int status;
};

double merge_stars_eq(double p_centre, void *params) {
  merge_stars_params *p = (struct merge_stars_params*)params;
  p->error_try = p->initial_error_try;
  int do_loop = 1;
  int pc = 0;
  while (do_loop == 1) {
    if (pc > 12) {
       PRC(p_centre);
       PRC(pc); PRC(p->n_shells); PRC(p->error_try); PRL(p->mass.size());
    }
    p->status = solve_HSE(p_centre*p_unit, p->m_product,
			  *p->mmu_A,  *p->mmu_B,
			  *p->entr_A, *p->entr_B,
			  p->m_min_A, p->m_min_B,
			  p->m_max_A, p->m_max_B,
			  p->Morig,
			  p->mass, p->radius, p->pressure, p->error_try,
			  p->n_shells);
    pc++;
    if (p->status == NEXCEED)
        p->error_try *= ERROR_SCALE_INC;
    else 
        if (p->mass.size() < p->n_shells)
            p->error_try *= 1.0/ERROR_SCALE_DEC;
        else
            do_loop = 0;
  }
  double dp = p->pressure.back()/p_centre/p_unit;
//  fprintf(stderr,  "merge_stars_eq: dp= %g \n", dp);
  return dp;
}

int mmas::merge_stars(double f_lost, int n_desired_shells) {
  p_unit = uG * pow(uMSUN,2.0)/pow(uRSUN,4.0);
  rho_unit = uMSUN/pow(uRSUN, 3.0);

  vector<double> Mass, entr, m_mu, id_sA;
  sort_model(*model_a, Mass, entr, m_mu, id_sA);
  for (size_t i = 0; i < entr.size(); i++) {
    entr[i] *=  randx;
    m_mu[i] *=  randx;
  }
  gslInterp mmu_A (Mass, m_mu);
  gslInterp entr_A(Mass, entr);
  double m_min_A = max(Mass[1], Mass[2]*1.0e-4);
  double m_max_A = (1.0 - 1.0e-10)*Mass.back();

  vector<double> id_sB;
  sort_model(*model_b, Mass, entr, m_mu, id_sB); 
  for (size_t i = 0; i < entr.size(); i++) {
    entr[i] *=  randx;
    m_mu[i] *=  randx;
  }
  gslInterp mmu_B (Mass, m_mu);
  gslInterp entr_B(Mass, entr);
  double m_min_B = max(Mass[1], Mass[2]*1.0e-4);
  double m_max_B = (1.0 - 1.0e-10)*Mass.back();

  /* ------------------- */
  /*  merging two stars */
  /* ------------------- */

  merge_stars_params p;
  p.mmas_ptr = this;

  p.m_min_A = m_min_A;
  p.m_max_A = m_max_A;
  p.mmu_A   = &mmu_A;
  p.entr_A  = &entr_A;

  p.m_min_B = m_min_B;
  p.m_max_B = m_max_B;
  p.mmu_B   = &mmu_B;
  p.entr_B  = &entr_B;

  if (n_desired_shells > 10000)
    p.n_shells = 10000;
  else
    p.n_shells = n_desired_shells;

  
  gsl_function F;
  F.function = &merge_stars_eq;
  F.params   = &p;
  
//   const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
  const gsl_root_fsolver_type *T = gsl_root_fsolver_falsepos;
//   const gsl_root_fsolver_type *T = gsl_root_fsolver_bisection;
  gsl_root_fsolver *s            = gsl_root_fsolver_alloc(T);

  /* bracketing central pressure */

  p.m_product = p.m_max_A + p.m_max_B;
  f_lost = (1 - f_lost*mass_loss()/100.0);
  p.m_product *= f_lost;           // mass loss

  cerr << "Bracketing central pressure: \n";
  double p_0 = (model_a->get_shell(0).pressure + model_b->get_shell(0).pressure)/p_unit;

  p.initial_error_try = ERROR_TRY0;
  merge_stars_eq(p_0, &p);
  p.initial_error_try = p.error_try;
  
  double sgn = p.pressure.back();

  int iterb_max = 100;
  int iterb = 0;
  int factor = 2;
  int stop = 0;
//  double p_last_last = p_0;
  double p_last = p_0;
  while(!stop) {
    if (iterb++ > iterb_max) {
      cerr << "Failed to bracket root. Quit\n";
      exit(-1);
    }
//    p_last_last = p_last;
    p_last = p_0;
    if (sgn < 0) p_0 *= factor;
    else         p_0 *= 1.0/factor;

    merge_stars_eq(p_0, &p);
//    cerr << "try p_centre= " << p_0 << "  ";
//    cerr << ", p_last_shell= " << p.pressure.back()/p_unit << "  "; // << endl;
//    PRC(p.Morig.size()); PRL(p.mass.size());

    if (sgn*p.pressure.back() < 0) stop = 1;
    
//     PRC(p.error_try); PRC(p.mass.size());
//     if (p.mass.size() < NTRY)   p.error_try *= 1.0/5;
//     if (p.mass.size() > 5*NTRY) p.error_try *= 10;
//     PRL(p.error_try);
  }

  const double p_min = min(p_last, p_0);
  const double p_max = max(p_last, p_0);
  cerr << "p is in [" << p_min << ", " << p_max << "] \n";

  p.n_shells = n_desired_shells;
    
  fprintf(stderr, "\n----------------------------------------------------\n");
  fprintf(stderr, "Computing p_centre using %s method:\n", 
	  gsl_root_fsolver_name (s));

  gsl_root_fsolver_set (s, &F, p_min, p_max);

  fprintf (stderr, "%5s [%9s, %9s] %9s %9s\n",
	   "iter", "lower", "upper", "root", "err(est)");
  
  int status;
  int iter = 0, iter_max = 1000;
  double rel_err = 1.0e-4;
  do {
    status      = gsl_root_fsolver_iterate (s);
    double r    = gsl_root_fsolver_root    (s);
    double x_lo = gsl_root_fsolver_x_lower (s);
    double x_hi = gsl_root_fsolver_x_upper (s);

//     if ((fabs(x_hi - x_lo)/r < rel_err*10) && (p.n_shells != n_desired_shells)) {
//       cerr << "Now use n_desired_shells ... \n";
//       p.n_shells = n_desired_shells;
//     }

    status = gsl_root_test_interval (x_lo, x_hi,
				     0, rel_err);
    
    if (status == GSL_SUCCESS)
      fprintf (stderr, "Converged in %u iterations:\n", iter+1);
    
//     if (iter%100 == 0 || status == GSL_SUCCESS)
      fprintf (stderr, "%5d [%.7f, %.7f] %.7f %.7f\n",
	       iter, x_lo, x_hi,
	       r, x_hi - x_lo);
    iter++;
  }
  while (status == GSL_CONTINUE && iter < iter_max);
  fprintf(stderr, "----------------------------------------------------\n");

  double p_centre = gsl_root_fsolver_root(s);
  gsl_root_fsolver_free (s);

  /* --- */

  /* warning, this routine is unaware of possible
     shell-sorting perfomed above.
     Currently, no sorting is performed, and
     sort_model is a stub routine. Nevertheless,
     in future one may change this, and therefore please
     do not forget to modify the subsequent routine too in order
     to make it aware of the sorting.
  */
  
  vector<double> mass, chemicals[NUMBER_OF_CHEMICAL_SPECIES];

  int n = model_a->get_num_shells();
  for (int i = 0; i < n; i++) {
    mass_shell &shell = model_a->get_shell(i);
    chemical_composition &che = shell.composition;
    mass.push_back(shell.mass);
    chemicals[0].push_back(che.H1*randx);
    chemicals[1].push_back(che.He4*randx);
    chemicals[2].push_back(che.O16*randx);
    chemicals[3].push_back(che.N14*randx);
    chemicals[4].push_back(che.C12*randx);
    chemicals[5].push_back(che.Ne20*randx);
    chemicals[6].push_back(che.Mg24*randx);
    chemicals[7].push_back(che.Si28*randx);
    chemicals[8].push_back(che.Fe56*randx);
//     PRC(shell.mass); PRL(chemicals[0].back());
  }
  gslInterp *cheA[NUMBER_OF_CHEMICAL_SPECIES];
  for (int i = 0; i < NUMBER_OF_CHEMICAL_SPECIES; i++) {
    cheA[i] = new gslInterp(mass, chemicals[i]);
    chemicals[i].clear();
  }
  mass.clear();

  n = model_b->get_num_shells();
  for (int i = 0; i < n; i++) {
    mass_shell &shell = model_b->get_shell(i);
    chemical_composition &che = shell.composition;
    mass.push_back(shell.mass);
    chemicals[0].push_back(che.H1*randx);
    chemicals[1].push_back(che.He4*randx);
    chemicals[2].push_back(che.O16*randx);
    chemicals[3].push_back(che.N14*randx);
    chemicals[4].push_back(che.C12*randx);
    chemicals[5].push_back(che.Ne20*randx);
    chemicals[6].push_back(che.Mg24*randx);
    chemicals[7].push_back(che.Si28*randx);
    chemicals[8].push_back(che.Fe56*randx);
  }
  gslInterp *cheB[NUMBER_OF_CHEMICAL_SPECIES];
  for (int i = 0; i < NUMBER_OF_CHEMICAL_SPECIES; i++) {
    cheB[i] = new gslInterp(mass, chemicals[i]);
    chemicals[i].clear();
  }
  mass.clear();

  delete product;
  product = new usm;
  
  mass_shell shell;
  int mass_size = p.mass.size();
  for (int i = 0; i < mass_size-1; i++) {
    product->add_shell(shell);
  }
  product->build_hashtable();
  product->star_mass   = p.mass[p.mass.size()-1];
  product->star_radius = p.radius[p.mass.size()-1];

  PRL(mass_size);
  
  for (int i = 0; i < mass_size-1; i++) {
    mass_shell &shell = product->get_shell(i);
    chemical_composition &che = shell.composition;
    shell.radius      = p.radius[i];
    shell.mass        = p.mass[i];
    shell.pressure    = p.pressure[i];
    if (p.Morig[i] > 0) {
      real m = p.Morig[i];
      shell.entropy   = entr_A.eval(m);
      shell.mean_mu   = mmu_A.eval(m);
      che.H1   = cheA[0]->eval(m);
      che.He4  = cheA[1]->eval(m);
      che.O16  = cheA[2]->eval(m);
      che.N14  = cheA[3]->eval(m);
      che.C12  = cheA[4]->eval(m);
      che.Ne20 = cheA[5]->eval(m);
      che.Mg24 = cheA[6]->eval(m);
      che.Si28 = cheA[7]->eval(m);
      che.Fe56 = cheA[8]->eval(m);
    } else {
      real m = -p.Morig[i];
      shell.entropy   = entr_B.eval(m);
      shell.mean_mu   = mmu_B.eval(m);
      che.H1   = cheB[0]->eval(m);
      che.He4  = cheB[1]->eval(m);
      che.O16  = cheB[2]->eval(m);
      che.N14  = cheB[3]->eval(m);
      che.C12  = cheB[4]->eval(m);
      che.Ne20 = cheB[5]->eval(m);
      che.Mg24 = cheB[6]->eval(m);
      che.Si28 = cheB[7]->eval(m);
      che.Fe56 = cheB[8]->eval(m);
    }
    
    shell.density     = compute_density(shell.pressure, shell.entropy, shell.mean_mu);
    shell.temperature = compute_temperature(shell.density, shell.pressure, shell.mean_mu);
    if (i%100 == 101) {
      PRC(shell.mass);
      PRC(shell.pressure);
      PRC(p.Morig[i]);
      PRL(che.H1);
    }
  }
  cerr << " Done merging stars \n";

  for (int i = 0; i < NUMBER_OF_CHEMICAL_SPECIES; i++) {
    delete cheA[i];
    delete cheB[i];
  }
  return 0;
}
