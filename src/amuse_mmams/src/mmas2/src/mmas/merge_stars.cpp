#include "mmas.h"
#include "eos/eos.h"

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
 #  The second order predictor-corrector method is invoked.
 #
 #  Limitations: No angular momentum is taken into account. Only non-rotating products
 #
 # TODO LIST:      
 #   Add rotation
 #
*/

real p_unit, rho_unit;
real max_product_mass;

int solve_HSE(real p_centre,
	      vector<real> &dm_a,   vector<real> &dm_b,   vector<real> &dm_p,   vector<int> &id_a,
	      vector<real> &entr_a, vector<real> &entr_b, vector<real> &entr_p, vector<int> &id_b,
	      vector<real> &mu_a,   vector<real> &mu_b,   vector<real> &mu_p,   vector<int> &id_p,
	      vector<real> &mass,   vector<real> &radius, vector<real> &pressure) {
  int n_a = dm_a.size();
  int n_b = dm_b.size();
  int n_p = n_a + n_b;

  dm_p.resize(n_p);
  entr_p.resize(n_p);
  mu_p.resize(n_p);
  id_p.resize(n_p);

  int nc_a = 0;
  int nc_b = 0;
  int nc_p = 0;
  
  real m_product = 0;

  vector<real> mu, x, pres;
  mu.push_back(0.0);
  x.push_back(0.0);
  pres.push_back(p_centre/p_unit);

  while ((nc_a + nc_b) < n_p-2) {
//      PRC(nc_a); PRC(nc_b); PRC(nc_p); PRL(nc_a+nc_b);

    /* predictor step */
    
    real density_a = compute_density(pres[nc_p]*p_unit, entr_a[nc_a], mu_a[nc_a]);
    real density_b = compute_density(pres[nc_p]*p_unit, entr_b[nc_b], mu_b[nc_b]);
//      PRC(pow(mu[nc_p],1.5)); PRC(entr_a[nc_a]); PRL(entr_b[nc_b]);
//      PRC(pres[nc_p]); PRC(density_a); PRL(density_b);

    real density_p = -1;
    real entrh, mmuh;
    if (nc_a == n_a - 1) density_a = -1;
    if (nc_b == n_b - 1) density_b = -1;
    if (density_a > density_b) {
      density_p = density_a;
      dm_p  [nc_p] = dm_a  [nc_a];
      entr_p[nc_p] = entr_a[nc_a];
      mu_p  [nc_p] = mu_a  [nc_a];
      id_p  [nc_p] = id_a  [nc_a];
      entrh = entr_a[nc_a];
      mmuh  = mu_a  [nc_a];
      if (nc_a < n_a - 1) {
	entrh = (entrh + entr_a[nc_a+1])/2.0;
	mmuh  = (mmuh   + mu_a  [nc_a+1])/2.0;
      }
      nc_a += 1;
    } else { //if (density_b > density_a) {
      density_p = density_b;
      dm_p  [nc_p] = dm_b  [nc_b];
      entr_p[nc_p] = entr_b[nc_b];
      mu_p  [nc_p] = mu_b  [nc_b];
      id_p  [nc_p] = id_b  [nc_b];
      entrh = entr_b[nc_b];
      mmuh  = mu_b  [nc_b];
      if (nc_b < n_b - 1) {
	entrh = (entrh + entr_b[nc_b+1])/2.0;
	mmuh  = (mmuh  + mu_b  [nc_b+1])/2.0;
      }

      nc_b += 1;
//     } else {
//       cerr << "cannot happen ... quit mmas::merge_stars() \n";
//       exit(-1);
    }
    
    density_p *= 1.0/rho_unit;

    real dmu = pow(m_product + dm_p[nc_p], 2.0/3.0) - pow(m_product, 2.0/3.0);
    m_product += dm_p[nc_p];
    real dmuh = dmu * 0.5;
    
    real muh = mu[nc_p] + dmuh;
    real mu_x = pow(4.0*PI/3.0 * density_p, 2.0/3.0);
    if (x[nc_p] > 0) mu_x = mu[nc_p]/x[nc_p];

    real dp_dmu = -3.0/8.0/PI * (mu_x*mu_x);
    real dx_dmu = +3.0/4.0/PI * sqrt(mu_x)/density_p;

    real ph     = pres[nc_p] + dp_dmu * dmuh;
    real xh     = x   [nc_p] + dx_dmu * dmuh;
    if (ph < 0) {
      pres.push_back(ph);
      break;
    }

    /* corrector step */
    
    density_p = compute_density(ph*p_unit, entrh, mmuh)/rho_unit;

    mu_x = muh/xh;
    dp_dmu = -3.0/8.0/PI * (mu_x*mu_x);
    dx_dmu = +3.0/4.0/PI * sqrt(mu_x)/density_p;
    
    real p1     = pres[nc_p] + dp_dmu * dmu;
    real x1     = x   [nc_p] + dx_dmu * dmu;
    real mu1    = mu  [nc_p] + dmu;

    mu.push_back  (mu1);
    pres.push_back(p1);
    x.push_back   (x1);
    if (p1 < 0) {
      break;
    }

    if (pow(mu1, 1.5) > max_product_mass) {
      break;
    }

//     PRC(dmu); PRL(dp_dmu);
    nc_p += 1;
  }

  mass.resize(mu.size());
  for (int i = 0; i < mu.size(); i++) mass[i] = pow(mu[i], 1.5);

  pressure.resize(pres.size());
  for (int i = 0; i < pres.size(); i++) pressure[i] = pres[i]*p_unit;

  radius.resize(x.size());
  for (int i = 0; i < x.size(); i++) radius[i] = sqrt(x[i]);

  return pres.size();
}

int mmas::merge_stars(real f_loss) {
  p_unit = uG * pow(uMSUN,2.0)/pow(uRSUN,4.0);
  rho_unit = uMSUN/pow(uRSUN, 3.0);


  max_product_mass = f_loss *  (model_a->star_mass + model_b->star_mass);

  vector<real> dm_a, entr_a, mu_a;
  vector<real> dm_b, entr_b, mu_b;
  vector<real> dm_p, entr_p, mu_p;
  vector<int> id_a, id_b, id_p;

  vector<real> mass, radius, pressure;

  real m_prev = 0;
  for (int i = 0; i < model_a->get_num_shells(); i++) {
    mass_shell &shell = model_a->get_shell(i);
    real dm;
    if (i == 0) dm = shell.mass;
    else        dm = shell.mass - m_prev;
    m_prev = shell.mass;
    
    dm_a.push_back(dm);
    mu_a.push_back(shell.mean_mu);
    entr_a.push_back(shell.entropy);
    id_a.push_back(-(i+1));
  }

  m_prev = 0;
  for (int i = 0; i < model_b->get_num_shells(); i++) {
    mass_shell &shell = model_b->get_shell(i);

    real dm;
    if (i == 0) dm = shell.mass;
    else        dm = shell.mass - m_prev;
    m_prev = shell.mass;

    dm_b.push_back(dm);
    mu_b.push_back(shell.mean_mu);
    entr_b.push_back(shell.entropy);
    id_b.push_back(i+1);
  }

  cerr << "Bracketing central pressure: \n";
  real p_0 = 1;
  int sz = solve_HSE(p_0*p_unit,
		     dm_a,   dm_b,   dm_p,   id_a,
		     entr_a, entr_b, entr_p, id_b,
		     mu_a,   mu_b,   mu_p,   id_p,
		     mass,   radius, pressure);
  while(pressure[pressure.size()-1] < 0) {
    p_0 *= 2;
    cerr << "try p_centre= " << p_0 << endl;
    sz = solve_HSE(p_0*p_unit,
		   dm_a,   dm_b,   dm_p,   id_a,
		   entr_a, entr_b, entr_p, id_b,
		   mu_a,   mu_b,   mu_p,   id_p,
		   mass,   radius, pressure);
  }

  real p_1 = p_0;
  p_0 = p_1/2.0;
  cerr << "Central pressure is in the following range: [" << p_0 << ", " << p_1 << "]\n";

  cerr << "Iterating until converged: \n";
  int iter = 0;
  bool converged = false;
  real p_init = 0;
  while (!converged) {
    real p_init_old = p_init;
    p_init = (p_0 + p_1)/2.0;
    sz = solve_HSE(p_init*p_unit,
		   dm_a,   dm_b,   dm_p,   id_a,
		   entr_a, entr_b, entr_p, id_b,
		   mu_a,   mu_b,   mu_p,   id_p,
		   mass,   radius, pressure); 
    cerr << "iteration= " << iter << ": try p_centre= " << p_init
	 << ", p_last_shell= " << pressure[pressure.size()-1]/p_unit << endl;
    if (fabs(pressure[pressure.size() - 1])/p_unit < 1.0e-7) {
      cerr << "  ********** \n";
      cerr << "Congratulations! converged @ p_centre= " << p_init << endl;
      cerr << "  ********** \n";
      break;
    }
    if (pressure[pressure.size() - 1] < 0) p_0 = p_init;
    else                                   p_1 = p_init;
    iter++;
  }

  delete product;
  product = new usm;

  mass_shell shell;
  for (int i = 0; i < mass.size()-1; i++) {
    product->add_shell(shell);
  }
  product->build_hashtable();
  product->star_mass   = mass[mass.size()-1];
  product->star_radius = radius[mass.size()-1];

  for (int i = 0; i < mass.size()-1; i++) {
    mass_shell &shell = product->get_shell(i);
    shell.radius      = radius[i];
    shell.mass        = mass[i];
    shell.density     = compute_density(pressure[i], entr_p[i], mu_p[i]);
    shell.pressure    = pressure[i];
    shell.entropy     = entr_p[i];
    shell.mean_mu     = mu_p[i];
    shell.temperature = compute_temperature(shell.density, shell.pressure, shell.mean_mu);
  }

  return 0;
}
