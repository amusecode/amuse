#include "sapporohostclass.h"

sapporo grav;

extern "C" {
  
   const char *kernelFile = "CUDA/kernels6thDP.ptx";   
   
   double *h2_i;

  /*
  Initialize the GPU library

  SHOULD BE CALLED FIRST! AND ONLY ONCE!

  */

  void initialize_special(int ndev, int *list)
  { 
      
      static int initCount = 0;
      
      if(initCount > 0)
        return;
        
        
      //Open the sapporo library using the specified kernel file, gpu devices and
      //integration order
      grav.open(kernelFile, list, ndev, SIXTH);

      int n_pipes = grav.get_n_pipes();
      
      //Allocate arrays 
       h2_i        = (double*)malloc(sizeof(double)*n_pipes);      
       
       for(int i=0; i < n_pipes; i++)
       {
         h2_i[i] = 0;
       }
  }

  void initialize()
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
      }
      
      initialize_special(how_many, devList);
      
      delete[] devList;   
  }



  /*
  add = address
  pos, vel, acc, jrk, snp, crk, 
  mass, 
  time = current particle time
  id = unique particle id
  eps2 = softening of j-particle


  */  
  void set_j_particle(int add, double pos[3], double vel[3], double acc[3],
              double jrk[3], double snp[3], double crk[3], double mass, double time, int id, double eps2)
  {
    //Not used variables
    double k18[3];
    double dt = 0;
    
    grav.set_j_particle(add, id, time, dt, mass, k18, jrk, acc, vel, pos,
                        snp, crk, eps2);
  }


  /*
  Set time of the prediction

  time = time to which particles are predicted
  nj = amount of particles that are predicted

  */
  void predict_all(double time, int nj)
  {
    //Set the next time
    grav.set_time(time); 
    
    //Force predictor to run, it will only be run if required
    //but this way we ensure that data is there for pick_up_predictor
    grav.forcePrediction(nj);    
  }

  /*

  Do not execute prediction, but only copy the particles
  into the predicted buffers.

  */
  void no_predict_all(double time, int nj)
  {
    //Set the next time
      grav.set_no_time(); 
      
      //Force predictor to run, it will only be run if required
      //but this way we ensure that data is there for pick_up_predictor
      grav.forcePrediction(nj);
  }

  /*

  Return the predicted values for a particle at an address

  addr = address of the particle
  
  id   = the particle id
  mass = the mass of the particle
  eps2 = the softening value of the particle
  pos = buffer to store predicted position
  vel = buffer to store predicted velocity
  acc = buffer to store predicted acceleration

  */
  void pick_up_predictor_2(int addr, int &id, double &mass, double &eps2, double pos[3],  double vel[3],  double acc[3])
  {
    double idTemp;
    //Retrieve predicted properties for particle at address addr
    grav.retrieve_predicted_j_particle(addr, mass, idTemp, eps2, pos, vel, acc);
    id = (int) idTemp;
  }


  /*
  Calculate the gravity on the i-particles

  //Input
  ni = number of particles to be integrated
  nj = number of sources
  pos, vel, acc, mass, eps2

  //Output
  acc, jrk, snp, potential (phi)
  nnb = nearest neighbour ID
  nnb_r2 = distance to the nearest neighbour. (Squared distance + softening)
  nnb_r2 =   double r2 = EPS2 + dx*dx + dy*dy + dz*dz;
  */
    
  void calc_force_on_predictors(int ni, int nj, int ids[], double pos[][3], double vel[][3],double acc[][3], 
                                double mass[], double eps2[], 
                                double accNew[][3], double jrkNew[][3], 
                                double snpNew[][3], double crkNew[][3],
                                double phi[], int nnb[], double nnb_r2[])
  {
    //ni is number of i-particles, can be equal or larger than NJ!
    //pred[] is an array containing the i-particles
    //Have to walk over the ni particles in steps equal to number of pipes
    int n_pipes = grav.get_n_pipes(); 
    
    //We calc the grav on i particles in blocks    
    for(int j=0; j < ni; j+= n_pipes)
    { 
      int    tempni = min(n_pipes, ni-j);
      double eps2_t   = eps2[0];
      
      //TODO check eps2, and make sure we start a kernel without neighbours, since h2 is not 
      //set by the calling function, but we want the nearest neighbour info, so we should use
      //some kind of neighbour info
      grav.startGravCalc(nj, tempni, &ids[j], &pos[j], &vel[j], &acc[j], NULL, NULL, eps2_t, h2_i, &eps2[j]);
      
      
      //Call results, reuse the previous used arrays as much as possible
      //id_i, pos_i, vel_i, 0, h2_i ->  are not used in the function we call
      //Store the results first in our temporary continous arrays and then 
      //loop over them to put them in the correct arrays
      grav.getGravResults(nj, tempni,
                          &ids[j],  &pos[j], &vel[j], 0, &h2_i[0],             
                          &accNew[j], &jrkNew[j], &snpNew[j], &crkNew[j], 
                          &phi[j], &nnb[j], &nnb_r2[j], true);  
                                
      
    } //for j  to n_pipes
  } //end grav calc             
                              
} //end extern C                              

