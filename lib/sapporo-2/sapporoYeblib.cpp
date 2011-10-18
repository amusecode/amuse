#include "sapporohostclass.h"

#include "vec.h" 

sapporo grav;
sapporo grav2;



//6th order specific for Michiko's code
typedef double DB3[3];


extern "C" {
  

  
struct predictor{
                double pos [3]; //  6
                double vel [3]; // 12
                double acc [3]; // 18
                double mass;    // 20
                double id  ;    // 22
                //double pad [2];
                double eps2;
                    predictor(){
                            assert(sizeof(*this) == 12*8);
                    }
            };  

struct Force{
    //static const int nword = 9;
    vec acc;
    vec jrk;
    vec snp;
    double phi;
    int nnb_id;
    float  nnb_r2;
    Force() : acc(0.0), jrk(0.0), snp(0.0), phi(0.0), nnb_id(-1), nnb_r2(HUGE){}    void operator += (Force &rhs){
        acc += rhs.acc;
        jrk += rhs.jrk;
        snp += rhs.snp;
        phi += rhs.phi;
    }
    void clear(){
        acc[0] = acc[1] = acc[2] = 0.0;
        jrk[0] = jrk[1] = jrk[2] = 0.0;
        snp[0] = snp[1] = snp[2] = 0.0;
        phi = 0.0;
        nnb_id = -1;
        nnb_r2 = HUGE;

    }
    
    
  //// for BH particles ////
    void add_BH_force(double eps2, double mBH, vec pos_BH, vec vel_BH, vec acc_BH,
                      vec pos_i, vec vel_i, vec acc_i, int iBH){

            vec pos0 = pos_BH - pos_i;
            vec vel0 = vel_BH - vel_i;
            vec acc0 = acc_BH - acc_i;
            double dist2 = pos0*pos0;
            double r2 = dist2 + eps2;
            double rinv2 = 1.0 / r2;
            double rinv = sqrt(rinv2);
            double rinv3 = rinv*rinv2;

            double alpha = rinv2*(vel0*pos0);
            double beta = rinv2*(vel0*vel0 + pos0*acc0) + alpha*alpha;

            vec acc1 = mBH * rinv3 * pos0;
            vec jerk1 = mBH*rinv3*vel0 - 3.0*alpha*acc1;
            double phi1 = -mBH * rinv;

            snp += mBH*rinv3*acc0 - 6.0*alpha*jerk1 - 3.0*beta*acc1;
            acc += acc1;
            jrk += jerk1;
            phi += phi1;

            if(dist2<nnb_r2){
                nnb_id = iBH;
                nnb_r2 = dist2;
                //cout << index << endl;
            }
    }
    
};  

namespace yebisu{
struct predictor{
                double pos [3]; //  6
                double vel [3]; // 12
                double acc [3]; // 18
                double mass;    // 20
                double id  ;    // 22
                //double pad [2];
                double eps2;
                    predictor(){
                            assert(sizeof(*this) == 12*8);
                    }
            };  

struct Force{
    //static const int nword = 9;
    vec acc;
    vec jrk;
    vec snp;
    double phi;
    int nnb_id;
    float  nnb_r2;
    Force() : acc(0.0), jrk(0.0), snp(0.0), phi(0.0), nnb_id(-1), nnb_r2(HUGE){}    void operator += (Force &rhs){
        acc += rhs.acc;
        jrk += rhs.jrk;
        snp += rhs.snp;
        phi += rhs.phi;
    }
    void clear(){
        acc[0] = acc[1] = acc[2] = 0.0;
        jrk[0] = jrk[1] = jrk[2] = 0.0;
        snp[0] = snp[1] = snp[2] = 0.0;
        phi = 0.0;
        nnb_id = -1;
        nnb_r2 = HUGE;

    }
    
    
  //// for BH particles ////
    void add_BH_force(double eps2, double mBH, vec pos_BH, vec vel_BH, vec acc_BH,
                      vec pos_i, vec vel_i, vec acc_i, int iBH){

            vec pos0 = pos_BH - pos_i;
            vec vel0 = vel_BH - vel_i;
            vec acc0 = acc_BH - acc_i;
            double dist2 = pos0*pos0;
            double r2 = dist2 + eps2;
            double rinv2 = 1.0 / r2;
            double rinv = sqrt(rinv2);
            double rinv3 = rinv*rinv2;

            double alpha = rinv2*(vel0*pos0);
            double beta = rinv2*(vel0*vel0 + pos0*acc0) + alpha*alpha;

            vec acc1 = mBH * rinv3 * pos0;
            vec jerk1 = mBH*rinv3*vel0 - 3.0*alpha*acc1;
            double phi1 = -mBH * rinv;

            snp += mBH*rinv3*acc0 - 6.0*alpha*jerk1 - 3.0*beta*acc1;
            acc += acc1;
            jrk += jerk1;
            phi += phi1;

            if(dist2<nnb_r2){
                nnb_id = iBH;
                nnb_r2 = dist2;
                //cout << index << endl;
            }
    }
    
};  
}


  const char *kernelFile = "CUDA/kernels6thDP.ptx";       

  DB3    *pos_i;
  DB3    *vel_i;
  DB3    *acc_i;
  DB3    *jrk_i;
  DB3    *snp_i;
  DB3    *crk_i;
  double *mass_i;
  int *id_i;
  int *nnb_i;             //Nearest neighbour
  double *dsmin_i;        //Distance of nearest neighbour
  double *eps2_i;
  double *h2_i;
  double *pot_i;
  

  void initialize_special(int ndev, int *list)
  { 
    
    static int initCount = 0;
    
    if(initCount > 0)
      return;
   
    initCount++;
    if(initCount == 2)
      grav = grav2;
//     if(initCalled)
//       return;
//     initCalled = true;
      
      
    //Open the sapporo library using the specified kernel file, gpu devices and
    //integration order
    grav.open(kernelFile, list, ndev, SIXTH);

    int n_pipes = grav.get_n_pipes();
    //Allocate arrays 
    pos_i       = (DB3*)malloc(sizeof(DB3)*n_pipes);    
    vel_i       = (DB3*)malloc(sizeof(DB3)*n_pipes);    
    acc_i       = (DB3*)malloc(sizeof(DB3)*n_pipes);    
    jrk_i       = (DB3*)malloc(sizeof(DB3)*n_pipes);    
    snp_i       = (DB3*)malloc(sizeof(DB3)*n_pipes);    
    crk_i       = (DB3*)malloc(sizeof(DB3)*n_pipes);    
    mass_i      = (double*)malloc(sizeof(double)*n_pipes);    
    id_i        = (int*)malloc(sizeof(int)*n_pipes);            
    nnb_i       = (int*)malloc(sizeof(int)*n_pipes);  
    eps2_i      = (double*)malloc(sizeof(double)*n_pipes);    
    h2_i        = (double*)malloc(sizeof(double)*n_pipes);      
    dsmin_i     = (double*)malloc(sizeof(double)*n_pipes);      
    pot_i       = (double*)malloc(sizeof(double)*n_pipes);   
  }
  
  
  void initialize(int n)
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
  
  
//(int ndev, int *list)

  void set_jp_2(int add, double pos[3], double vel[3], double acc[3],
            double jrk[3], double snp[3], double crk[3], double mass, double time, int id,
            double eps2)
  {
    //Not used variables
    double k18[3];
    double dt = 0;
    
    grav.set_j_particle(add, id, time, dt, mass, k18, jrk, acc, vel, pos,
                        snp, crk, eps2);
  }
  
  void predict_all(double time, int nj)
  {
    //Set the next time
    grav.set_time(time); 
    
    //Force predictor to run, it will only be run if required
    //but this way we ensure that data is there for pick_up_predictor
    grav.forcePrediction(nj);
  }
  
  void no_predict_all(double time, int nj)
  {
    //Set the next time
    grav.set_no_time(); 
    
    //Force predictor to run, it will only be run if required
    //but this way we ensure that data is there for pick_up_predictor
    grav.forcePrediction(nj);
  }  
  
  
  void pick_up_predictor_2(int addr, yebisu::predictor &pred)
  {
    //Retrieve predicted properties for particle at address addr
    grav.retrieve_predicted_j_particle(addr,     pred.mass, 
                                       pred.id,  pred.eps2,
                                       pred.pos, pred.vel,
                                       pred.acc);
  }
  
  //Get the force on i-particles, using one common eps2 value
  void calc_force_on_predictors(int ni, predictor pred[], Force force[], int nj)
  {
    //ni is number of i-particles, can be equal or larger than NJ!
    //pred[] is an array containing the i-particles
    //Have to walk over the ni particles in steps equal to number of pipes
    int n_pipes = grav.get_n_pipes();
    
    
    //We calc the grav on i particles in blocks    
    for(int j=0; j < ni; j+= n_pipes)
    {  
      for(int i = 0; i < n_pipes; i++)  //Copy data in arrays
      {
        if(j+i < ni)
        {
          pos_i[i][0] = pred[j+i].pos[0]; 
          pos_i[i][1] = pred[j+i].pos[1];
          pos_i[i][2] = pred[j+i].pos[2];
          vel_i[i][0] = pred[j+i].vel[0]; 
          vel_i[i][1] = pred[j+i].vel[1];
          vel_i[i][2] = pred[j+i].vel[2];       
          acc_i[i][0] = pred[j+i].acc[0]; 
          acc_i[i][1] = pred[j+i].acc[1];
          acc_i[i][2] = pred[j+i].acc[2]; 
          mass_i[i]   = pred[j+i].mass;
          id_i[i]     = (int)pred[j+i].id;
          eps2_i[i]   = pred[j+i].eps2;
          h2_i[i]     = 0;
          
        } //end if j+i < ni
      }//end for i < n_pipes
      
      int tempni = min(n_pipes, ni-j);

      double eps2 = pred[0].eps2;
      
      //TODO check eps2, and make sure we start a kernel without neighbours, since h2 is not 
      //set by the calling function, but we want the nearest neighbour info, so we should use
      //some kind of neighbour info
      grav.startGravCalc(nj, tempni, id_i, pos_i, vel_i, acc_i, NULL, NULL, eps2, h2_i, eps2_i);
      
      //Call results, reuse the previous used arrays as much as possible
      //id_i, pos_i, vel_i, 0, h2_i ->  are not used in the function we call
      //Store the results first in our temporary continous arrays and then 
      //loop over them to put them in the correct arrays
      
      grav.getGravResults(nj, tempni,
                          id_i,  pos_i, vel_i, 0, h2_i,             
                          acc_i, jrk_i, snp_i, crk_i, 
                          pot_i, nnb_i, dsmin_i, true);  
                          
      //Put results in the force array of the calling program
      for(int i = 0; i < n_pipes; i++)  //Copy data in arrays
      {
        if(j+i < ni)
        {
          force[j+i].acc[0] = acc_i[i][0];
          force[j+i].acc[1] = acc_i[i][1];
          force[j+i].acc[2] = acc_i[i][2];
          force[j+i].jrk[0] = jrk_i[i][0];
          force[j+i].jrk[1] = jrk_i[i][1];
          force[j+i].jrk[2] = jrk_i[i][2];
          force[j+i].snp[0] = snp_i[i][0];
          force[j+i].snp[1] = snp_i[i][1];
          force[j+i].snp[2] = snp_i[i][2];          
          force[j+i].phi    = -pot_i[i];   
          force[j+i].nnb_id = nnb_i[i];   
          force[j+i].nnb_r2 = dsmin_i[i];   
          
/*          if(i < 5){
              fprintf(stderr, "in yeblib acc: %lf %lf %lf \tphi: %lf\tjrk:%lf %lf %lf\tsnp:%lf %lf %lf\tr2: %lf\n",                                                                             
                  force[i].acc[0], force[i].acc[1], force[i].acc[2], force[i].phi,                                                                                                
                  force[i].jrk[0], force[i].jrk[1], force[i].jrk[2],                                                                                                              
                  force[i].snp[0], force[i].snp[1], force[i].snp[2],                                                                                                              
                  force[i].nnb_r2 );
          }
                  */                                                   

          
          
        } //end if j+i < ni              
      }//end for i < n_pipes
    }//End ni loop     
  }// end calc_force_on_predictors
  
  
 void calc_force_on_predictors_epsj(int ni, predictor pred[], Force force[], int nj)
 {
   //Does exactly the same as calc_force_on_predictors, except that now we use individual softening
   //namely the softening value of the J-particles
   cerr << "TODO: Not implemented yet, have to check how I do this in the compute kernels, before "; 
   cerr << "I finalise this function! \n";
   
   exit(0);   
 }
  
}