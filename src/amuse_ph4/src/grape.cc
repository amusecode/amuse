
#include "grape.h"

// Local lookalikes for Sapporo interface functions.

extern "C" {

  int g6_open(int id)
  {
      return g6_open_(&id);
  }

  int g6_close(int id)
  {
      return g6_close_(&id);
  }

  int g6_npipes()
  {
      return g6_npipes_();
  }

  int g6_set_tunit(double t)
  {
      return g6_set_tunit_(&t);
  }

  int g6_set_xunit(double x)
  {
      return g6_set_xunit_(&x);
  }

  int g6_set_ti(int id, double ti)
  {
      return g6_set_ti_(&id, &ti);
  }

  int g6_set_j_particle(int cluster_id,
			int address,
			int index,
			double tj, double dtj,
			double mass,
			double k18[3], double j6[3],
			double a2[3], double v[3], double x[3])
  {
      return g6_set_j_particle_(&cluster_id,
				&address,
				&index,
				&tj, &dtj,
				&mass, k18, j6,
				a2, v, x);
  }

  void g6calc_firsthalf(int cluster_id,
			int nj, int ni,
			int index[], 
			double xi[][3], double vi[][3],
			double aold[][3], double j6old[][3],
			double phiold[3], 
			double eps2, double h2[])
  {
      return g6calc_firsthalf_(&cluster_id,
			       &nj, &ni,
			       index,
			       xi, vi,
			       aold, j6old,
			       phiold,
			       &eps2, h2);
  }

  int g6calc_lasthalf(int cluster_id,
		      int nj, int ni,
		      int index[],
		      double xi[][3], double vi[][3],
		      double eps2, double h2[],
		      double acc[][3], double jerk[][3], double pot[])
  {
      return g6calc_lasthalf_(&cluster_id,
			      &nj, &ni,
			      index,
			      xi, vi,
			      &eps2, h2,
			      acc, jerk, pot);
  }

  int g6calc_lasthalf2(int cluster_id,
		       int nj, int ni,
		       int index[],
		       double xi[][3], double vi[][3],
		       double eps2, double h2[],
		       double acc[][3], double jerk[][3], double pot[],
		       int inn[])
  {
      return g6calc_lasthalf2_(&cluster_id,
			       &nj, &ni,
			       index,
			       xi, vi,
			       &eps2, h2,
			       acc, jerk, pot,
			       inn);
  }

  int g6_initialize_jp_buffer(int cluster_id, int buf_size)
  {
      return g6_initialize_jp_buffer_(&cluster_id, &buf_size);
  }

  int g6_flush_jp_buffer(int cluster_id)
  {
      return g6_flush_jp_buffer_(&cluster_id);
  }

  int g6_reset(int cluster_id)
  {
      return g6_reset_(&cluster_id);
  }

  int g6_reset_fofpga(int cluster_id)
  {
      return g6_reset_fofpga_(&cluster_id);
  }

  int g6_read_neighbour_list(int cluster_id)
  {
      return g6_read_neighbour_list_(&cluster_id);
  }

  int g6_get_neighbour_list(int cluster_id,
			    int ipipe,
			    int maxlength,
			    int *n_neighbours,
			    int neighbour_list[])
  {
      return g6_get_neighbour_list_(&cluster_id,
				    &ipipe,
				    &maxlength,
				    n_neighbours,
				    neighbour_list);
  }

  void g6_set_neighbour_list_sort_mode(int mode) {}

  int g6_set_overflow_flag_test_mode(int aflag, int jflag, int pflag)
  {
      return 0;
  }
}
