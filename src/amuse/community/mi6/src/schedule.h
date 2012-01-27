#ifndef SCHEDULE_H
#define SCHEDULE_H

//#include"sort/sort.h"
#include"sort.h"

inline void select_ip(double value[], 
		      int address[], 
		      int &Nip, 
		      const int &Ntot){
  double ref_value = value[address[0]];
  int i=0;
  while( ref_value == value[address[i]] ){
    i++;
    if(i >= Ntot) break;
  }
  Nip = i;
}

inline void sort_select_iparticle(Particle prt[], 
				  int &Nip, 
				  const int &Ntot, 
				  int address[], 
				  double next_time[]){
  Qsort_index(next_time, address, 0, Nip-1);
  select_ip(next_time, address, Nip, Ntot);
}

#endif //SCHEDULE_H
