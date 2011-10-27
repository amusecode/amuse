/****************************
 *                          
 *     libSAPPORO  v2.0     
 *  a multiGPU GRAPE6 and beyond library !
 *        library           
 *                         
 * (c) 2010                 
 *
 * TODO License
 *
 ********************************/

#ifndef _SAPPORO_H_
#define _SAPPORO_H_

#include <builtin_types.h>
#include <iostream>
#include <vector>
#include <map>

#include <math.h>
using namespace std;

#include "sapdevclass.h"


#define INT_AS_FLOAT(x) (*((float*)&(x)))
#define FLOAT_AS_INT(x) (*((int*)&(x)))

//All DS stuff will be taken care on _on the device_
typedef float2 DS;  // double single;

struct DS4 {
  DS x, y, z, w;
  DS4 operator=(DS4 val) {
    x = val.x;
    y = val.y;
    z = val.z;
    w = val.w;
    return val;
  }
};
struct DS2 {
  DS x, y;
  DS2 operator=(DS2 val) {
    x = val.x;
    y = val.y;
    return val;
  }
};

inline DS to_DS(double a) {
  DS b;
  b.x = (float)a;
  b.y = (float)(a - b.x);
  return b;
}

inline double to_double(DS a) {
  double b;
  b = (double)((double)a.x + (double)a.y);
  return b;
}



struct memPointerJstruct
{
  int     *address;
  double2 *t_j;
  double4 *pos_j;
  double4 *vel_j;
  double4 *acc_j;       //Acceleration
  double4 *jrk_j;       //Jerk
  double4 *snp_j;       //Snap
  double4 *crk_j;       //Crack
  int     *id_j;  
  int     count;
  int     toCopy;
} typedef memPointerJstruct;


double get_time();

  
class sapporo {
protected:
   
  int n_pipes;
  int nCUDAdevices;

  int nj_modified;     //Modified number of particles that has to be send to the device  
  int nj_max;          //Total number of allocated memory particles
  int nj_total;        //Total number of j particles in the system (same as nj max?)
//   

  bool nj_updated;      //Indicates if there are updated j particles
  bool predJOnHost;     //Indicate if the predicted J-particles are on the host or not

  double EPS2;
  
  map<int, int4> mappingFromIndexToDevIndex;
  //int4.x = device
  //int4.y = arraylocation
  //int4.z = device address
 
  vector<int>     address_j;
  vector<int>     id_j;
  vector<double2> t_j;
  vector<double4> pos_j;
  vector<double4> vel_j;
  vector<double4> acc_j;
  vector<double4> jrk_j;
  vector<double4> snp_j;
  vector<double4> crk_j;

  
  double t_i;
  vector<int>     id_i; 
  vector<double4> pos_i;
  vector<double4> vel_i;
  vector<double4> accin_i;
  vector<double4> acc_i;
  vector<double4> jrk_i;
  vector<double4> snp_i;
  vector<double4> crk_i;
  vector<double>   ds_i;
  vector<int>     ngb_list_i;
  
  vector<memPointerJstruct> jMemAddresses;
  
  
  bool ngb_list_copied;
  bool predict;
  
  void cleanUpDevice();
  void free_cuda_memory(int);
  void allocate_cuda_memory(int);
  void send_j_particles_to_device(int);
  void send_i_particles_to_device(int, int);
  void fetch_data_from_device(int, int );
  void retrieve_i_particle_results(int ni);
  int  fetch_ngb_list_from_device();
  void initialize_firstsend();
  void increase_jMemory(); 
  
  void copyJInDev(int nj);
  void predictJParticles(int nj)  ;

  double evaluate_gravity(int, int);
  
  
  bool isFirstSend;             //Used to check if this is the first time we sent particles to the
                                //device, so we have to allocate memory
  int integrationOrder;         //Order of the ingtegrator we use, should be set during open call, default is fourth order
                        

  
//TODO ask EG why we increase the i particles 
public:
  sapporo() {
//     n_pipes = NTHREADS;
    n_pipes = NPIPES;
    pos_i.resize(n_pipes);
    vel_i.resize(n_pipes);
    t_i = 0.0;

    ngb_list_copied = false;
    
    //Clear the J particles
    address_j.clear();    
    t_j.clear();
    pos_j.clear();
    vel_j.clear();
    acc_j.clear();
    jrk_j.clear();
    id_j.clear();
    
    predict     = false;
    isFirstSend = true;
    nj_updated  = false;
    
    
    integrationOrder = FOURTH;
    

  };
  ~sapporo() {
     cleanUpDevice();
  };
  
  
  //Device communication functions
  void send_j_particles_to_device();
  void send_i_particles_to_device(int i);
    
  //Library interface functions
  int open(std::string kernelFile, int *devices, int nprocs, int order);
  int close();
  int get_n_pipes();
  int set_time(double ti);
  int set_no_time();

  int set_j_particle(int address,
                     int index,
                     double tj, double dtj,
                     double mass,
                     double k18[3], double j6[3],
                     double a2[3], double v[3], double x[3],
                     double snp[3], double crk[3], double eps);
                     
  void startGravCalc(int nj, int ni,
                     int index[], 
                     double xi[][3], double vi[][3],
                     double aold[][3], double j6old[][3],
                     double phiold[3], 
                     double eps2, double h2[],
                     double eps2_i[]);
  /*int calc_lasthalf(int cluster_id,
                    int nj, int ni,
                    int index[], 
                    double xi[][3], double vi[][3],
                    double eps2, double h2[],
                    double acc[][3], double jerk[][3], double pot[]);
  int calc_lasthalf2(int cluster_id,
                     int nj, int ni,
                     int index[], 
                     double xi[][3], double vi[][3],
                     double eps2, double h2[],
                     double acc[][3], double jerk[][3], double pot[],
                     int nnbindex[]);*/
                     
  int getGravResults(int nj, int ni,
                     int index[], 
                     double xi[][3],      double vi[][3],
                     double eps2,         double h2[],
                     double acc[][3],     double jerk[][3], 
                     double snp[][3],     double crk[][3],
                     double pot[],        int nnbindex[],
                     double dsmin_i[],    bool ngb);
                     
  void forcePrediction(int nj);            
  void retrieve_predicted_j_particle(int addr,       double &mass, 
                                     double &id,     double &eps2,
                                     double pos[3],  double vel[3],
                                     double acc[3]);
                                     
  void retrieve_j_particle_state(int addr,       double &mass, 
                                 double &id,     double &eps2,
                                 double pos[3],  double vel[3],
                                 double acc[3],  double jrk[3], double ppos[3],
                                 double pvel[3], double pacc[3]);
  
  int fetch_ngb_list_from_device(int);
  int read_ngb_list(int);
  int get_ngb_list(int cluster_id,
                   int ipipe,
                   int maxlength,
                   int &nblen,
                   int nbl[]);
};

#endif 
