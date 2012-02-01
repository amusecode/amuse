#ifndef EVOLVE_H
#define EVOLVE_H

#include"Particle.h"
#include"merge.h"
#include"schedule.h"
#include"mpi_interface.h"
#include"external_field.h"

#ifdef SAP
#include "6thorder.h"
#else
#include "6thorder_dummy.h"
#endif

extern int EX_FLAG;

void set_NSTEP(const int &_NSTEP);

int get_NSTEP();

void set_NBH_NJP(const int &_NBH, 
		 const int &_NJP);

void set_NJP(const int &_NJP);

void set_eta(const double &eta_s,
	     const double &eta_fs,
	     const double &eta_smbh,
	     const double &eta_imbh);

void setj_to_sapporo(Particle prt[],
		     int address[],
		     const int &Nip);

void set_eps2(const double &eps2_fs_fs,
	      const double &eps2_bh_fs,
	      const double &eps2_bh_bh);



void calc_force_using_sapporo(Particle prt[],
			      int address[],
			      const int &Nip,
			      const double &Tsys,
			      const int &mode);

void evolve_initialize(Particle prt[],
		       int address[],
		       const int &Ntot,
		       const int &NBH,
		       const int &Njp,
		       const double &Tsys);

/*
int evolve_onestep(Particle prt[],
		   int address[],
		   int &Nip,
		   const int &Ntot,
		   const int &NBH,
		   double &Tsys,
		   const double &Tmerge,
		   const double &maxdt,
		   const int &first_loop,
		   double &Egr);
*/

int evolve_onestep(Particle prt[],
		   int address[],
		   int &Nip,
		   const int &Ntot,
		   const int &NBH,
		   double &Tsys,
		   const double &Tmerge,
		   const double &maxdt,
		   const int &first_loop,
		   double &Egr,
		   const int &itr);

void sort_time_all(Particle prt[], 
		   int address[], 
		   int Ntot);

#endif //EVOLVE_H
