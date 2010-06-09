#include "sapporo.h"

sapporo grav;

extern "C" {
  int g6_open_(int *id) {  return grav.open(*id); }
  int g6_close_(int *id) {  std::cerr<<"close!"<<std::endl; return grav.close(*id); }
  int g6_npipes_() { return grav.get_n_pipes(); }
  int g6_set_tunit_(double*) {return 0;}
  int g6_set_xunit_(double*) {return 0;}
  int g6_set_ti_(int *id, double *ti) {return  grav.set_ti(*id, *ti); }
  int g6_set_j_particle_(int *cluster_id,
			 int *address,
			 int *index,
			 double *tj, double *dtj,
			 double *mass,
			 double k18[3], double j6[3],
			 double a2[3], double v[3], double x[3]) {
    return grav.set_j_particle(*cluster_id, *address, *index, *tj, *dtj,
			       *mass, k18, j6, a2, v, x);
  }
  void g6calc_firsthalf_(int *cluster_id,
			 int *nj, int *ni,
			 int index[], 
			 double xi[][3], double vi[][3],
			 double aold[][3], double j6old[][3],
			 double phiold[3], 
			 double *eps2, double h2[]) {
    grav.calc_firsthalf(*cluster_id, *nj, *ni,
			index, xi ,vi, aold, j6old, phiold,
			*eps2, h2);
  }
  int g6calc_lasthalf_(int *cluster_id,
		       int *nj, int *ni,
		       int index[], 
		       double xi[][3], double vi[][3],
		       double *eps2, double h2[],
		       double acc[][3], double jerk[][3], double pot[]) {
    return grav.calc_lasthalf(*cluster_id, *nj, *ni,
			      index, xi, vi, *eps2, h2, acc, jerk, pot);
  }
  int g6calc_lasthalf2_(int *cluster_id,
			int *nj, int *ni,
			int index[], 
			double xi[][3], double vi[][3],
			double *eps2, double h2[],
			double acc[][3], double jerk[][3], double pot[],
			int *inn) {
    return grav.calc_lasthalf2(*cluster_id, *nj, *ni,
			       index, xi, vi, *eps2, h2, acc, jerk, pot, inn);
  }
  
  int g6_initialize_jp_buffer_(int* cluster_id, int* buf_size) {return 0;}
  int g6_flush_jp_buffer_(int* cluster_id) {return 0;}
  int g6_reset_(int* cluster_id) {return 0;}
  int g6_reset_fofpga_(int* cluster_id) {return 0;}

#ifdef NGB
  int g6_read_neighbour_list_(int* cluster_id) {return grav.read_ngb_list(*cluster_id);}
#else
  int g6_read_neighbour_list_(int* cluster_id) {return 0;}
#endif

  int g6_get_neighbour_list_(int *cluster_id,
			     int *ipipe,
			     int *maxlength,
			     int *n_neighbours,
			     int neighbour_list[]) {
#ifdef NGB
    return grav.get_ngb_list(*cluster_id,
			     *ipipe,
			     *maxlength,
			     *n_neighbours,
			     neighbour_list);
#else
    return 0;
#endif
  }
			     
			     

}

