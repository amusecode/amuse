#include "sapporo.h"
#include <sys/time.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>

sapporo_multi_struct sapporo_multi_data[MAXCUDADEVICES];


double get_time() {
  struct timeval Tvalue;
  struct timezone dummy;
  
  gettimeofday(&Tvalue,&dummy);
  return ((double) Tvalue.tv_sec +
          1.e-6*((double) Tvalue.tv_usec));
}

int sapporo::open(int cluster_id) {
  fprintf(stderr, "sapporo::open --- ver 1.5 --- \n");
  nCUDAdevices = get_device_count();
  if (nCUDAdevices > 0) {
    fprintf(stderr, "sapporo::open - found %d CUDA device(s)\n",
            nCUDAdevices);
  } else {
    fprintf(stderr, "sapporo::open - FATAL! No CUDA device found\n");
    exit(-1);
  }
  
  nj_max = 131072;
  if(cluster_id >= nCUDAdevices) {
    fprintf(stderr, "sapporo::open - FATAL! Cluster id is too high, no cuda device with given id\n");
    return -1;
  }
  
  device_id = cluster_id;
  cudaSetDevice(device_id);
  allocate_cuda_memory(device_id);
  
  return 0;
}

int sapporo::close(int cluster_id) {
  printf("sapporo::close --- ver 1.5 --- \n");
  free_cuda_memory(device_id);
  if(address_j.size() > 0) {
    mapping_from_address_j_to_index_in_update_array.clear();
    address_j.clear();
    t_j.clear();
    pos_j.clear();
    vel_j.clear();
    acc_j.clear();
    jrk_j.clear();
  }
  return 0;
}

int sapporo::set_ti(int cluster_id, double ti) {
  t_i = to_DS(ti);
  predict = true;
  return 0;
}

int sapporo::get_n_pipes() {
  return n_pipes;
}

int sapporo::set_j_particle(int cluster_id,
                           int address,
                           int id,
                           double tj, double dtj,
                           double mass,
                           double k18[3], double j6[3],
                           double a2[3], double v[3], double x[3]) {

  if (address > nj_max) {
    fprintf(stderr, "Fatal! address= %d > nj_max= %d. I am giving up.\n",
            address, nj_max);
    exit(-1);
  }

    DS  Dmass = (DS){mass, INT_AS_FLOAT(id)};
    map<int,int>::iterator iterator = mapping_from_address_j_to_index_in_update_array.find(address);
    map<int,int>::iterator end = mapping_from_address_j_to_index_in_update_array.end();
    
    if(iterator != end)
    {
        int index = (*iterator).first;
        //printf("found index: %d for address: %d\n", index, address);
        t_j[index] = (DS2){to_DS(tj), to_DS(dtj)};
        pos_j[index] = ( (DS4){to_DS(x[0]), to_DS(x[1]), to_DS(x[2]), Dmass} );
        vel_j[index] = ( (float4){v[0],    v[1],    v[2],    0.0} );
        acc_j[index] = ( (float4){a2[0]*2, a2[1]*2, a2[2]*2, 0.0} );
        jrk_j[index] = ( (float4){j6[0]*6, j6[1]*6, j6[2]*6, 0.0} );
    }
    else
    {
        mapping_from_address_j_to_index_in_update_array[address] = address_j.size();
        //printf("set index: %d for address: %d\n",  address_j.size(), address);
        address_j.push_back(address);

        t_j.push_back( (DS2){to_DS(tj), to_DS(dtj)} );

        pos_j.push_back( (DS4){to_DS(x[0]), to_DS(x[1]), to_DS(x[2]), Dmass} );
        vel_j.push_back( (float4){v[0],    v[1],    v[2],    0.0} );
        acc_j.push_back( (float4){a2[0]*2, a2[1]*2, a2[2]*2, 0.0} );
        jrk_j.push_back( (float4){j6[0]*6, j6[1]*6, j6[2]*6, 0.0} );
        nj_modified = address_j.size();
    }
  return 0;
};

void sapporo::calc_firsthalf(int cluster_id,
                            int nj, int ni,
                            int id[], 
                            double xi[][3], double vi[][3],
                            double aold[][3], double j6old[][3],
                            double phiold[3], 
                            double eps2, double h2[]) {
                                
    
    
  
  for (int i = 0; i < ni; i++) {
    DS Dmass = (DS){h2[i], INT_AS_FLOAT(id[i])};
    pos_i[i] = (DS4) { to_DS(xi[i][0]),
                       to_DS(xi[i][1]),
                       to_DS(xi[i][2]),
                       Dmass };
    vel_i[i] = (float4){ vi[i][0], vi[i][1], vi[i][2], eps2};
    EPS2 = eps2;

  }

  if (address_j.size() > 0) {
    send_j_particles_to_device(device_id);
  }

  send_i_particles_to_device(device_id, ni);
  if(address_j.size() > 0) {
    mapping_from_address_j_to_index_in_update_array.clear();
    address_j.clear();
    t_j.clear();
    pos_j.clear();
    vel_j.clear();
    acc_j.clear();
    jrk_j.clear();
  }
  evaluate_gravity(ni, nj);

}

int sapporo::calc_lasthalf(int cluster_id,
                          int nj, int ni,
                          int index[], 
                          double xi[][3], double vi[][3],
                          double eps2, double h2[],
                          double acc[][3], double jerk[][3], double pot[]) {

  for (int i = 0; i < ni; i++) {
    pot[i] = 0;
    acc[i][0] =  acc[i][1] =  acc[i][2] = 0;
    jerk[i][0] = jerk[i][1] = jerk[i][2] = 0;
  }
  fetch_data_from_device(device_id, ni);
  
  for (int i = 0; i < ni; i++) {
      float4 acci = acc_i[i];
      float4 jrki = jrk_i[i];
      pot[i]    += acci.w;
      
      acc[i][0] += acci.x;
      acc[i][1] += acci.y;
      acc[i][2] += acci.z;
      
      jerk[i][0] += jrki.x;
      jerk[i][1] += jrki.y;
      jerk[i][2] += jrki.z;
  }
  return 0;
};


int sapporo::calc_lasthalf2(int cluster_id,
                           int nj, int ni,
                           int index[], 
                           double xi[][3], double vi[][3],
                           double eps2, double h2[],
                           double acc[][3], double jerk[][3], double pot[],
                           int nnbindex[]) {
  

    float ds_min[NTHREADS];
    for (int i = 0; i < ni; i++) {
        pot[i] = 0;
        acc[i][0]  = acc[i][1]  = acc[i][2]  = 0;
        jerk[i][0] = jerk[i][1] = jerk[i][2] = 0;
        nnbindex[i] = 0;
        ds_min[i] = 1.0e10;
    }

    double t1 = get_time();
    fetch_data_from_device(device_id, ni);

    for (int i = 0; i < ni; i++) {

        float4 acci = acc_i[i];
        float4 jrki = jrk_i[i];
        float  ds   = ds_i[i];

        //       fprintf(stdout, "device= %d, ni= %d pot = %g\n",
        //               dev, i, acci.x);

        pot[i]    += acci.w;

        acc[i][0] += acci.x;
        acc[i][1] += acci.y;
        acc[i][2] += acci.z;

        jerk[i][0] += jrki.x;
        jerk[i][1] += jrki.y;
        jerk[i][2] += jrki.z;

        if (ds < ds_min[i]) {
        int nnb = (int)(jrki.w);
        nnbindex[i] = nnb; 
        ds_min[i] = ds;
    }
  }
  return 0;
};

int sapporo::read_ngb_list(int cluster_id) {
    bool overflow = false;
    int ni = fetch_ngb_list_from_device(device_id);

    for (int j = 0; j < ni; j++) {
        if (ngb_list_i[device_id*NGB_PP*n_pipes + j*NGB_PP] >= NGB_PP) {
            overflow = true;
        }
    }
    

    return overflow;
}

int sapporo::get_ngb_list(int cluster_id,
                         int ipipe,
                         int maxlength,
                         int &nblen,
                         int nbl[]) {

    if (ipipe >= device.ni) {
        fprintf(stderr, "Fatal! ipipe= %d >= dev.ni= %d. I give up.\n",
                ipipe, device.ni);
        exit(-1);
    }
    
    bool overflow = false;
    nblen = 0;
    int offset = device_id*NGB_PP*n_pipes + NGB_PP*ipipe;
    int len = ngb_list_i[offset];
    memcpy(nbl+nblen, &ngb_list_i[offset+1], sizeof(int)*min(len, maxlength - len));
    nblen += len;
    if (nblen >= maxlength) {
        overflow = true;
    }

    sort(nbl, nbl + min(nblen, maxlength));
    return overflow;
}
