#ifndef ENERGY_H
#define ENERGY_H
#include"Particle.h"

inline void calc_energy(Particle prt[],
			const int &Ntot,
			double &E,
			double &Ek,
			double &Ep,
			const int &mode=0){
  if(mode == 0){
    Ek=0.0; Ep=0.0; E=0.0;
  }
  for(int i=0; i<Ntot; i++){
    Ek += 0.5*prt[i].mass*prt[i].vel*prt[i].vel;
    Ep += 0.5*prt[i].mass*(prt[i].phi - prt[i].phi_ex);
    Ep += prt[i].mass*prt[i].phi_ex;
  }
  E = Ek+Ep;
}

inline void calc_energy(Particle prt[],
			int address[],
			const int &Ntot,
			double &E,
			double &Ek,
			double &Ep,
			const int &mode=0){
  if(mode == 0){
    Ek=0.0; Ep=0.0; E=0.0;
  }
  for(int i=0; i<Ntot; i++){
    int id = address[i];
    Ek += 0.5*prt[id].mass*prt[id].vel*prt[id].vel;
    Ep += 0.5*prt[id].mass*(prt[id].phi - prt[id].phi_ex);
    Ep += prt[id].mass*prt[id].phi_ex;
  }
  E = Ek+Ep;
}

/*
void energy0(Particle *prt, 
	     const int &Ntot, 
	     double &E, 
	     double &Ek, 
	     double &Ep, 
	     const double &eps2=1e-6){
  Ek=0.0; Ep=0.0; E=0.0;
  for(int i=0; i<Ntot; i++){
    prt[i].phi=0.0;
    for(int j=0; j<Ntot; j++){
      if(i != j){
	Vector3 rij=prt[i].pos-prt[j].pos;
	double r2=rij*rij+eps2;
	double R = 1.0/sqrt(r2);
	prt[i].phi -= prt[j].mass*R;
      }
    }
  }
  for(int i=0; i<Ntot; i++){
    Ek += 0.5*prt[i].mass*prt[i].vel*prt[i].vel;
    Ep += 0.5*prt[i].mass*prt[i].phi;
    E=Ek+Ep;
  }
}
*/
#endif //ENERGY_H
