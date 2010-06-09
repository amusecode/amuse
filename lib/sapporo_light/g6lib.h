#ifndef _G6LIB_
#define _G6LIB_

extern "C" {
  int g6_open_(int *id);
  int g6_close_(int *id);
  int g6_npipes_();
  int g6_set_tunit_(double*);
  int g6_set_xunit_(double*);
  int g6_set_ti_(int *id, double *ti);
  int g6_set_j_particle_(int *cluster_id,
			 int *address,
			 int *index,
			 double *tj, double *dtj,
			 double *mass,
			 double k18[3], double j6[3],
			 double a2[3], double v[3], double x[3]);
  void g6calc_firsthalf_(int *cluster_id,
			 int *nj, int *ni,
			 int index[], 
			 double xi[][3], double vi[][3],
			 double aold[][3], double j6old[][3],
			 double phiold[3], 
			 double *eps2, double h2[]);
  int g6calc_lasthalf_(int *cluster_id,
		       int *nj, int *ni,
		       int index[], 
		       double xi[][3], double vi[][3],
		       double *eps2, double h2[],
		       double acc[][3], double jerk[][3], double pot[]);
  int g6calc_lasthalf2_(int *cluster_id,
			int *nj, int *ni,
			int index[], 
			double xi[][3], double vi[][3],
			double *eps2, double h2[],
			double acc[][3], double jerk[][3], double pot[],
			int *inn);
  
  int g6_initialize_jp_buffer_(int* cluster_id, int* buf_size);
  int g6_flush_jp_buffer_(int* cluster_id);
  int g6_reset_(int* cluster_id);
  int g6_reset_fofpga_(int* cluster_id);

  int g6_read_neighbour_list_(int* cluster_id);

  int g6_get_neighbour_list_(int *cluster_id,
			     int *ipipe,
			     int *maxlength,
			     int *n_neighbours,
			     int neighbour_list[]);
}

#endif
