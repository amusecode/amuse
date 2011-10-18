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

  int g6_close(int clusterid){ return g6_close_(&clusterid);}
  int g6_open(int clusterid){  return g6_open_(&clusterid);}
  int g6_npipes(){ return g6_npipes_();}
  void g6_set_ti(int cluster_id, double ti){  g6_set_ti_(&cluster_id, &ti);}
  int 	g6_reset(int cluster_id){  return g6_reset_(&cluster_id);}
  int 	g6_reset_fofpga(int cluster_id){  return g6_reset_fofpga_(&cluster_id);}
  int 	g6_set_j_particle(int clusterid, int address, int index, double tj, double dtj, double mass, double k18[3], double j6[3], double a2[3], double v[3], double x[3]){   return g6_set_j_particle_(&clusterid, &address,&index,&tj,&dtj,&mass,k18,j6,a2,v,x);}
  void 	g6calc_firsthalf(int cluster_id, int nj, int ni, int index[], double xi[][3], double vi[][3], double aold[][3], double j6old[][3], double phiold[], double eps2, double h2[]){   g6calc_firsthalf_ (&cluster_id, &nj, &ni, index, xi, vi, aold, j6old, phiold, &eps2,h2);}
  int 	g6calc_lasthalf(int cluster_id, int nj, int ni, int index[], double xi[][3], double vi[][3], double eps2, double h2[], double acc[][3], double jerk[][3], double pot[]){ return g6calc_lasthalf_ (&cluster_id, &nj, &ni, index, xi, vi, &eps2, h2, acc, jerk, pot); }
  void  g6_set_tunit(double newtunit){g6_set_tunit_(&newtunit);}
  void  g6_set_xunit(double newxunit){g6_set_xunit_(&newxunit);}




int g6_set_j_particle_multisend_mxfast_(int * clusterid,
                                 int * nclusters,
                                 int *address,
                                 int *index,
                                 double *mass,
                                 double x[3] ) //position
{
  double v[3] = {0,0,0};
  double tj = 0; 
  double dtj = 0;
  double k18[3] = {0,0,0};
  double j6[3] = {0,0,0};
  double a2[3] = {0,0,0};
   
  //Call the usual set_j_particle function with 0's for the v value
  return g6_set_j_particle_(clusterid, address, index, &tj, &dtj, mass, k18,j6,a2,v,x);

}

int g6_reinitialize(int clusterid){}
int g6_print_chip_status(int clusterid){}
int g6_flush_jp_buffer_and_multisend(int clusterid, int something){}
int g6_initialize_jp_buffer(int clusterid, int list_max){}
int g6_set_overflow_flag_test_mode_(int *clusterid){}

}

#endif
