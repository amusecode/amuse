#ifndef EXTERNAL_FIELD_H
#define EXTERNAL_FIELD_H

#include"Vector3.h"
#include"force.h"
#include"Particle.h"
#include"mpi_interface.h"

/////////////////////////////
///// forces from SMBH //////
/////////////////////////////

void calc_force_from_point_mass(Particle &prti,
				const int &mode);

void calc_force_from_point_mass_to_array(Particle prti[],
					 int address[],
					 const int &Nip,
					 const int &mode);

void accrete_mass(const double &mass);

double calc_rtide_cu(const double &mass,
		     const double &radius);



void get_SMBH(double &mass, 
	      Vector3 &pos,
	      Vector3 &vel);

void set_SMBH(const double &mass, 
	      const Vector3 &pos,
	      const Vector3 &vel);


void copy_SMBH_OLD_TO_NEW();

void copy_SMBH_NEW_TO_OLD();

void set_speed_of_light(const double &_c);

#endif //EXTERNAL_FIELD_H
