#ifndef EXTERNAL_FIELD_H
#define EXTERNAL_FIELD_H

#include"mpi_interface.h"
#include"Vector3.h"
#include"force.h"
#include"Particle.h"

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

int set_speed_of_light(double value);
int get_speed_of_light(double *value);

#ifdef __cplusplus
extern "C" {
#endif

int get_calculate_postnewtonian(int *value);
int set_calculate_postnewtonian(int value);
int get_calculate_postnewtonian_only_first_order(int *value);
int set_calculate_postnewtonian_only_first_order(int value);

#ifdef __cplusplus
}
#endif
#endif //EXTERNAL_FIELD_H
