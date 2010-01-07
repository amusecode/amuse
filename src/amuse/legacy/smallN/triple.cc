#include "amuse_interface.h"

local void create_binary(real m1, real m2, hdyn *b1, hdyn *b2,
			 kepler &k_in, hdyn *b12)
{
  b12->set_oldest_daughter(b1);
  b1->set_parent(b12);
  b2->set_parent(b12);
  b1->set_younger_sister(b2);
  b2->set_elder_sister(b1);

  real mtot = m1+m2;
  b12->set_mass(mtot);

  b1->set_mass(m1);
  b1->set_pos(-m2*k_in.get_rel_pos()/mtot);
  b1->set_vel(-m2*k_in.get_rel_vel()/mtot);
  b2->set_mass(m2);
  b2->set_pos(m1*k_in.get_rel_pos()/mtot);
  b2->set_vel(m1*k_in.get_rel_vel()/mtot);
}

local hdyn *create_triple(real sigma, real e_in, real e_out, real inclination,
			  real m1, real m2, real m3)
{
  // Without loss of generality choose the outer period to be 1.

  real p_out = 1, p_in = 1/sigma;
  vec l(1,0,0), t(0,1,0), n(0,0,1);

  // Create the inner orbit, wlog in the x-y plane with semi-major
  // axis along the x-axis.

  kepler k_in;
  k_in.set_time(0);
  k_in.set_total_mass(m1+m2);
  real a_in = semi_from_period(m1, m2, p_in);
  k_in.set_semi_major_axis(a_in);
  k_in.set_eccentricity(e_in);
  k_in.set_mean_anomaly(randinter(-M_PI, M_PI));
  k_in.set_orientation(l, t, n);
  k_in.initialize_from_shape_and_phase();

  // Create the outer orbit, with orientation to be specified.

  kepler k_out;
  k_out.set_time(0);
  k_out.set_total_mass(m1+m2+m3);
  real a_out = semi_from_period(m1+m2, m3, p_out);
  k_out.set_semi_major_axis(a_out);
  k_out.set_eccentricity(e_out);
  k_out.set_mean_anomaly(randinter(-M_PI, M_PI));
  l = random_unit_vec();
  t = l^random_unit_vec();
  t /= abs(t);
  n = l^t;
  k_out.set_orientation(l, t, n);
  k_out.initialize_from_shape_and_phase();

  // Combine into a triple system.

  hdyn *b = new hdyn(), *b12 = new hdyn();
  hdyn *b1 = new hdyn(), *b2 = new hdyn(), *b3 = new hdyn();

  b1->set_index(1);
  b2->set_index(2);
  b3->set_index(3);

  create_binary(m1, m2, b1, b2, k_in, b12);
  create_binary(m1+m2, m3, b12, b3, k_out, b);

  my_sys_stats(b);
  b->flatten_node();
  return b;
}

local bool my_integrate(hdyn *b, real t_end, real dt_struc,
			real eps2 = 0,
			real break_r2 = VERY_LARGE_NUMBER,
			real dt_log = VERY_LARGE_NUMBER,
			real dt_energy = VERY_LARGE_NUMBER,
			real dt_snap = VERY_LARGE_NUMBER)
{
  // SmallN_evolve expects a flat tree.  Do that now.

  b->flatten_node();

  // Note: every time we return from smallN_evolve, we change the data
  // in the simulation, even though the function as written should be
  // exactly restartable.  The reason appars to be tha difference in
  // precision between data that remain in the processor (80 bits?)
  // and thoae that get restored from memory (64 bits).  As a
  // workaround, we always call kira_smallN for a specified length of
  // time (1 time unit here, but this should be adaptive); so long as
  // dt_struc is a multiple of this, we should have reproducible
  // results...

  real dt_dyn = 1;
  bool over = false;
  while (b->get_time() < t_end) {
    real t_struc = b->get_time() + dt_struc;
    while (b->get_time() < t_struc) {
      // int status =
      smallN_evolve(b, b->get_time() + dt_dyn, break_r2,
		    false, dt_log, dt_energy, dt_snap);
    }
    if (over = check_structure(b, eps2)) break;
  }
  return over;
}

int main()
{
  // Set up the interaction.

  real sigma = 10, e_in = 0, e_out = 0.1, inclination = 1,
       m1 = 1, m2 = 1, m3 = 1;
  hdyn *b = create_triple(sigma, e_in, e_out, inclination, m1, m2, m3);

  kira_options ko;
  ko.perturber_criterion = 2;
  b->set_kira_options(&ko);

  b->set_system_time(0);
  for_all_nodes(hdyn, b, bi)
    bi->set_time(0);

  // The time unit is the outer orbital period.

  real t_end = 1000;
  real t_struc = 50;
  if (my_integrate(b, t_end, t_struc)) {
    dyn *d = get_tree(b);
    my_sys_stats(d);
  } else
    cerr << "interaction not over!" << endl;
}
