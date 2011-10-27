#include "sapporohostclass.h"

sapporo grav;



//4th order specific for GRAPE6 compatability



extern "C" {

  const char *kernelFile = "CUDA/kernels4th.ptx";	
  double *dsmin_i;        //Distance of nearest neighbour
  
  double acc_i[3];      //To store the multiplied acc
  double jrk_i[3];      //To store the multiplied jrk
  


  int g6_open_special(int ndev, int *list)
  {  
    //Now we have the number of devices and the (possible) devices
    //if how_many is negative we let the Driver decide which
    //devices to use. Otherwise they should be specified in the config file
    
    //Open the GPUs
    int res     = grav.open(kernelFile, list, ndev, FOURTH);
    
    //Read properties and allocate memory buffers
    int n_pipes = grav.get_n_pipes();
    dsmin_i     = (double*)malloc(sizeof(double)*n_pipes);      
    
    return res;    
  }  
  
  int g6_open_(int *id)
  {
    //Check for a config file if its there use it  
    int *devList = NULL;
    int how_many = 0;
    FILE *fd;
    if ((fd = fopen("sapporo2.config", "r"))) {
      char line[256];
      fprintf(stderr, "sapporo2::open - config file is found\n");
      fgets(line, 256, fd);
      sscanf(line, "%d", &how_many);

      //Read the devices we want to use
      if(how_many > 0)
      {
        devList = new int[how_many];
        for (int i = 0; i < how_many; i++) {
            fgets(line, 256, fd);
            sscanf(line, "%d", &devList[i]);
        }
      }
    } else {
      fprintf(stderr," sapporo2::open - no config file is found \n");
      how_many = 0;
    }
    
    int res = g6_open_special(how_many, devList);
    
    delete[] devList;
    
    return res;
  }


  int g6_close_(int *id) {  return grav.close(); }
  int g6_npipes_() { return grav.get_n_pipes(); }
  int g6_set_tunit_(double*) {return 0;}
  int g6_set_xunit_(double*) {return 0;}
  int g6_set_ti_(int *id, double *ti) {return  grav.set_time(*ti); }
  int g6_set_j_particle_(int *cluster_id,
			 int *address,
			 int *index,
			 double *tj, double *dtj,
			 double *mass,
			 double k18[3], double j6[3],
			 double a2[3], double v[3], double x[3]) {
    //Grape6 is 4th order so set snp and crk to null and set individual softening
    //value to 0
    
    //multiply the acc and jrk by a factor 2 and 6, for GRAPE compatability    
    acc_i[0] = a2[0]*2; acc_i[1] = a2[1]*2; acc_i[2] = a2[2]*2;
    jrk_i[0] = j6[0]*6; jrk_i[1] = j6[1]*6; jrk_i[2] = j6[2]*6;
        
    return grav.set_j_particle(*address, *index, *tj, *dtj,
			       *mass, k18, jrk_i, acc_i, v, x, NULL, NULL, 0); 
  }
  void g6calc_firsthalf_(int *cluster_id,
			 int *nj, int *ni,
			 int index[], 
			 double xi[][3], double vi[][3],
			 double aold[][3], double j6old[][3],
			 double phiold[3], 
			 double *eps2, double h2[]) {    
    grav.startGravCalc(*nj, *ni,
			index, xi ,vi, aold, j6old, phiold,
			*eps2, h2, NULL);
  }
  int g6calc_lasthalf_(int *cluster_id,
		       int *nj, int *ni,
		       int index[], 
		       double xi[][3], double vi[][3],
		       double *eps2, double h2[],
		       double acc[][3], double jerk[][3], double pot[]) {
    
    return grav.getGravResults(*nj, *ni,
                               index, xi, vi, *eps2, h2,
                               acc, jerk, NULL, NULL, pot, NULL, NULL, false);
    
//     return grav.calc_lasthalf(*cluster_id, *nj, *ni,
// 			      index, xi, vi, *eps2, h2, acc, jerk, pot);
  }
  int g6calc_lasthalf2_(int *cluster_id,
			int *nj, int *ni,
			int index[], 
			double xi[][3], double vi[][3],
			double *eps2, double h2[],
			double acc[][3], double jerk[][3], double pot[],
			int *inn) {
    
    return grav.getGravResults(*nj, *ni,
                               index, xi, vi, *eps2, h2,
                               acc, jerk, NULL, NULL, pot, inn, dsmin_i, true);    
    
//     return grav.calc_lasthalf2(*cluster_id, *nj, *ni,
// 			       index, xi, vi, *eps2, h2, acc, jerk, pot, inn);
  }
  
    
  void get_j_part_data(int addr, int nj,
                         double *pos, 
                         double *vel, 
                         double *acc, 
                         double *jrk, 
                         double *ppos, 
                         double *pvel,
                         double &mass,
                         double &eps2,
                         int    &id) {    

    
    //Retrieve predicted properties for particle at address addr
    double idt;
    
    grav.retrieve_j_particle_state(addr,     mass, 
                                       idt,  eps2,
                                       pos, vel,
                                       acc, jrk, ppos, pvel, NULL);    
     id = (int) idt;
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
			     
			     

} //end extern C

