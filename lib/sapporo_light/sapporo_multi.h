#ifndef _SAPPORO_MULTI_
#define _SAPPORO_MULTI_

struct sapporo_multi_struct {
  int device;
  int offset;

  float EPS2;

  bool ngb;
  int nj, ni;
  int nj_total;
  int  nj_max;
  int     nj_modified;
  bool    predict;
  float   t_i_x, t_i_y;
  
  int     *address_j;
  
  DS2     *t_j;
  DS4     *Ppos_j;
  float4  *Pvel_j;
  
  DS4     *pos_j;
  float4  *vel_j;
  float4  *acc_j;
  float4  *jrk_j;
  
  DS4     *pos_i;
  float4  *vel_i;
  float4  *acc_i;
  float4  *jrk_i;
  
  float   *ds_i;
  
  int *ngb_list;
};


#endif
