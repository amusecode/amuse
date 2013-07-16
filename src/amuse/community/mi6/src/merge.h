#ifndef MERGE_H
#define MERGE_H

#include"external_field.h"
#include"Particle.h"

int merge_check(Particle prt[],
		int address[],
		const int &Nip,
		const int &Nmerge);


int get_merge_candidates(int index, Particle **candidate1, Particle **candidate2);

void merge_prt();

void merge(Particle *prt0,
	   Particle *prt1);

void destroy(Particle *prt0, 
	     Particle *prt1);


void accrete(Particle *prt0);

void get_merged_prt(Particle *(prt_merged[]), int &Nmerge_loop);

void get_accreted_prt(Particle *(prt_accreted[]), int &Naccrete_loop);

void set_Ndead(const int &_Nmerge,
	       const int &_Naccrete);

#endif // MERGE_H
