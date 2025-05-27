#ifndef SCHEDULE_H
#define SCHEDULE_H

#include"sort.h"

inline void select_ip(double *value, int *index, int &Nip, const int &Ntot){
  double ref_value = value[index[0]];
  int i=0;
  while( ref_value == value[index[i]] ){
    i++;
    if(i >= Ntot) break;
  }
  Nip = i;
}

inline void sort_select_iparticle(int &Nip, const int &Ntot, int index[], double next_time[]){
  if(Nip >=1 ){
    Qsort_index(next_time, index, 0, Nip-1);
    select_ip(next_time, index, Nip, Ntot);
  }
}
#endif //SCHEDULE_H
