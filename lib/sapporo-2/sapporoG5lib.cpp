#include "sapporohostclass.h"

sapporo grav;



//GRAPE5 compatable interface

/*
  Implemented functions
void g5_set_range(double xmin, double xmax, double mmin)  
void g5_open()
void g5_close();
int g5_get_number_of_pipelines() 
void g5_set_n(int nj)
void g5_set_eps_to_all(double eps)

void g5_set_xmj(int addr, int nj, double xj[][3], double jmass[])
void g5_calculate_force_on_x(double xi[][3], double a[][3], double *phi, int ni)
void g5_set_range(double xmin, double xmax, double mmin) //No use on GPU

//Sapporo specific
int g5_open_special(int ndev, int *list)

*/



extern "C" {

//   const char *kernelFile = "CUDA/kernelsG5DS.ptx";
  const char *kernelFile = "CUDA/kernelsG5SP.ptx";        
  
  int nj_to_use = -1;
  double eps_to_use = -1;
  
  double *h2_i;

  int g5_open_special(int ndev, int *list)
  {  
    //Now we have the number of devices and the (possible) devices
    //if how_many is negative we let the Driver decide which
    //devices to use. Otherwise they should be specified in the config file
    
    //Open the GPUs
    int res     = grav.open(kernelFile, list, ndev, GRAPE5);
    
    //Read properties and allocate memory buffers
    int n_pipes = grav.get_n_pipes();
    
    h2_i        = (double*)malloc(sizeof(double)*n_pipes);      
    return res;    
  }  
  
  void g5_open()
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

    g5_open_special(how_many, devList);
    
    delete[] devList;  
  }


  void g5_close() {  grav.close(); }
  int g5_get_number_of_pipelines() { return grav.get_n_pipes(); }
  
  void g5_set_n(int nj)
  {
    nj_to_use = nj;
  }
  
  void g5_set_eps2_to_all(double eps)
  {
    eps_to_use = eps;
  }
  
  void g5_set_eps_to_all(double eps)
  {
    eps_to_use = eps*eps;
  }
  

  //Set the j-particles
  void g5_set_xmj(int addr, int nj, double xj[][3], double jmass[])
  {
    //Address is used to set ranged of particles
    for(int i = addr; i < addr+nj; i++)
    {
       grav.set_j_particle(i, i, 0, 0, jmass[i], NULL, NULL, NULL, NULL,
                                 xj[i], NULL, NULL, 0);
    }
  }  
  
  void g5_set_jp(int addr, int nj, double jmass[], double xj[][3])
  {
    g5_set_xmj(addr, nj, xj, jmass);
  }
  
  //Force calculation functions, a bit different than in g5 library
  //g5_set_xi will copy particles AND start the runs
  //These 3 are used if the user handles the number of pipe limit him/herself
  void g5_set_xi(int ni, double xi[][3])
  {
    int n_pipes = grav.get_n_pipes();
    double eps2 = eps_to_use;
    
    if(eps2 < 0)
    {
      fprintf(stderr, "SAPPORO ERROR: eps not set using 'g5_set_eps_to_all'exit now \n");
      exit(-1);
    }
    
    if(ni > n_pipes)
    {
      fprintf(stderr, "SAPPORO ERROR: ni > n_pipes exit now \n");     
      exit(0);
    }
    //Start the gravity computation
    grav.startGravCalc(nj_to_use, ni, NULL, xi, NULL, NULL, NULL, NULL, eps2, h2_i, NULL);    
  }
  void g5_run(){} //Already started when we set xi
  
  void g5_get_force(int ni, double a[][3], double *phi)
  {
    //Retrieve the results
    grav.getGravResults(nj_to_use, ni,
                          NULL,  NULL, NULL, 0, NULL,             
                          a, NULL, NULL, NULL, 
                          phi, NULL, NULL, false);  
  }
  
  
  //Set the i particles, calculate force and return the results
  void g5_calculate_force_on_x(double xi[][3], double a[][3], double *phi, int ni)
  { 
    //ni is number of i-particles, can be equal or larger than NJ! 
    //Have to walk over the ni particles in steps equal to number of pipes
    int n_pipes = grav.get_n_pipes();
    
    double eps2 = eps_to_use;
    
    if(eps2 < 0)
    {
      fprintf(stderr, "SAPPORO ERROR: eps not set using 'g5_set_eps_to_all'exit now \n");
      exit(-1);
    }
    
    if(nj_to_use <= 0)
    {
      fprintf(stderr, "SAPPORO ERROR: nj not set using 'g5_set_n'exit now \n");
      exit(-1);
    }
      
    //We calc the grav on i particles in blocks    
    for(int j=0; j < ni; j+= n_pipes)
    {              
      int tempni = min(n_pipes, ni-j);
     
      //Start the gravity computation
      grav.startGravCalc(nj_to_use, tempni, NULL, &xi[j], NULL, NULL, NULL, NULL, eps2, h2_i, NULL);
  
      //Retrieve the results
      grav.getGravResults(nj_to_use, tempni,
                          NULL,  NULL, NULL, 0, NULL,             
                          &a[j], NULL, NULL, NULL, 
                          &phi[j], NULL, NULL, false);  
    }//End ni loop 
    

  }//end calculate force
  
  
  
  void g5_set_range(double xmin, double xmax, double mmin)
  {
    //Does nothing in the GPU library
  }
  
  //Some non relevant functions fo the GPU
  int g5_get_jmemsize()
  {
    return INT_MAX;
  }
  
  int g5_get_number_of_real_pipelines()
  {
    return 1;
  }
  
  int g5_get_pcibus_freq()
  {
    return 1;
  }

} //end extern C

