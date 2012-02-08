#ifndef _6THLIB_
#define _6THLIB_


extern "C" {
  /*
  Initialize the GPU library

  SHOULD BE CALLED FIRST! AND ONLY ONCE!

  */

  void initialize_special(int ndev, int *list);

  void initialize();

  /*
  add = address
  pos, vel, acc, jrk, snp, crk, 
  mass, 
  time = current particle time
  id = unique particle id
  eps2 = softening of j-particle
  */  
  void set_j_particle(int add, double pos[3], double vel[3], double acc[3],
              double jrk[3], double snp[3], double crk[3], double mass, double time, int id, double eps2);


  /*
  Set time of the prediction

  time = time to which particles are predicted
  nj = amount of particles that are predicted

  */
  void predict_all(double time, int nj);  

  /*
  Do not execute prediction, but only copy the particles
  into the predicted buffers.
  */
  void no_predict_all(double time, int nj);

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
  void pick_up_predictor_2(int addr, int &id, double &mass, double &eps2, double pos[3],  double vel[3],  double acc[3]);


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
                                double phi[], int nnb[], double nnb_r2[]);                              
} //end extern C                              

#endif
