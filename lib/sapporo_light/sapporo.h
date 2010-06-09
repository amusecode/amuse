/****************************
 *                          
 *     libSAPPORO  v1.0     
 *  a multiGPU GRAPE6-like  
 *        library           
 *                         
 * (c) 2008                 
 *     Evghenii Gaburov     
 *     Sterrenkunding Instituut "Anton Pannekoek"
 *     Section Computational Science
 *     Universiteit van Amsterdam
 *
 ********************************/

#ifndef _SAPPORO_H_
#define _SAPPORO_H_

#include "sapporo_defs.h"
#include <builtin_types.h>
#include <iostream>
#include <vector>
#include <math.h>
using namespace std;

// #include <cutil.h>
#include <cuda_runtime.h>


#define INT_AS_FLOAT(x) (*((float*)&(x)))
#define FLOAT_AS_INT(x) (*((int*)&(x)))

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



struct dev_struct {
  int ni, nj;

  /*********** j-particles **********/

  int    *address_j;

  DS2    *t_j;
  DS4    *Ppos_j;     // predicted position
  float4 *Pvel_j;     // predicted velocity

  DS4    *pos_j;
  float4 *vel_j;
  float4 *acc_j;
  float4 *jrk_j;

  /*********** i-particles **********/

  DS4    *pos_i;
  float4 *vel_i;
  float4 *acc_i;
  float4 *jrk_i;
  int    *ngb_list_i;
  float  *ds_i;
};

#include "sapporo_multi.h"

double get_time();

extern "C"
{
  int  get_device_count();
  cudaError_t host_evaluate_gravity(sapporo_multi_struct);
}

class sapporo {
protected:
  int n_pipes, nj_max;
  int nCUDAdevices;
  int device_id;

  bool predict;
  double EPS2;

  int            nj_modified;
  vector<int>    address_j;
  vector<DS2>    t_j;
  vector<DS4>    pos_j;
  vector<float4> vel_j;
  vector<float4> acc_j;
  vector<float4> jrk_j;

  
  DS t_i;
  vector<DS4>    pos_i;
  vector<float4> vel_i;
  vector<float4> acc_i;
  vector<float4> jrk_i;
  vector<float>  ds_i;
  vector<int>    ngb_list_i;
  
  dev_struct device;
  
  bool ngb_list_copied;
  

  void free_cuda_memory(int);
  void allocate_cuda_memory(int);
  void send_j_particles_to_device(int);
  void send_i_particles_to_device(int, int);
  void fetch_data_from_device(int, int );

  double evaluate_gravity(int, int);

public:
  sapporo() {
    n_pipes = NTHREADS;
    pos_i.resize(n_pipes);
    vel_i.resize(n_pipes);
    acc_i.resize(n_pipes*MAXCUDADEVICES);
    jrk_i.resize(n_pipes*MAXCUDADEVICES);
    ds_i.resize(n_pipes*MAXCUDADEVICES);
    t_i = (DS){0,0};

    ngb_list_i.resize(n_pipes*NGB_PP*MAXCUDADEVICES);
    ngb_list_copied = false;
    
    address_j.clear();
    
    t_j.clear();
    pos_j.clear();
    vel_j.clear();
    acc_j.clear();
    jrk_j.clear();

    predict = false;
    nj_modified = 0;
    
    device.address_j = NULL;

    device.t_j   = NULL;
    device.pos_j = NULL;
    device.acc_j = NULL;
    device.vel_j = NULL;
    device.jrk_j = NULL;

    device.pos_i = NULL;
    device.acc_i = NULL;
    device.vel_i = NULL;
    device.jrk_i = NULL;
    
  };
  ~sapporo() {};
  
  int open(int cluster_id);
  int close(int cluster_id);
  int get_n_pipes();
  int set_ti(int cluster_id, double ti);

  int set_j_particle(int cluster_id,
		     int address,
		     int index,
		     double tj, double dtj,
		     double mass,
		     double k18[3], double j6[3],
		     double a2[3], double v[3], double x[3]);
  void calc_firsthalf(int cluster_id,
		      int nj, int ni,
		      int index[], 
		      double xi[][3], double vi[][3],
		      double aold[][3], double j6old[][3],
		      double phiold[3], 
		      double eps2, double h2[]);
  int calc_lasthalf(int cluster_id,
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
		     int nnbindex[]);

  
  int fetch_ngb_list_from_device(int);
  int  read_ngb_list(int);
  int get_ngb_list(int cluster_id,
		   int ipipe,
		   int maxlength,
		   int &nblen,
		   int nbl[]);
};

#endif 
