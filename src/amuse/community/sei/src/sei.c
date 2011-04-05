/********************************************
 *
 * Implementation of the SEI integrator
 * Hanno Rein & Scott Tremaine 2011
 *
 * The functions correspond to the operators
 * H0, Phi and SEI as described in the paper.
 * Each of them evolves the particles p for 
 * one timestep dt. 
 *
 * Compile with:
 * gcc -Wall -o sei -lm sei.c
 *
 ********************************************/ 
#include <math.h>
#include <stdio.h>
#include "sei.h"
#define OMEGA 1.

// Cache sin() tan() values.
double lastdt=0;
double sindt, tandt;

// This function evolves a particle under 
// Hamiltonian H0 exactly up to machine precission.
void operator_H0(double dt, struct particle* p){
	if (lastdt!=dt){
		// Only calculate sin() and tan() if timestep changed
		sindt = sin(OMEGA*(-dt));
		tandt = tan(OMEGA*(-dt/2.));
		lastdt = dt;
	}
		
	// Integrate vertical motion
	const double zx = p->z * OMEGA;
	const double zy = p->vz;
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double zt1 =  zx - tandt*zy;			
	const double zyt =  sindt*zt1 + zy;
	const double zxt =  zt1 - tandt*zyt;	
	p->z  = zxt/OMEGA;
	p->vz = zyt;

	// Integrate motion in xy directions
	const double aO = 2.*p->vy + 4.*p->x*OMEGA;	// Center of epicyclic motion
	const double bO = p->y*OMEGA - 2.*p->vx;	

	const double ys = (p->y*OMEGA-bO)/2.; 		// Epicycle vector
	const double xs = (p->x*OMEGA-aO); 
	
	// Rotation implemeted as 3 shear operators
	// to avoid round-off errors
	const double xst1 =  xs - tandt*ys;			
	const double yst  =  sindt*xst1 + ys;
	const double xst  =  xst1 - tandt*yst;	

	p->x  = (xst+aO)    /OMEGA;			
	p->y  = (yst*2.+bO) /OMEGA - 3./2.*aO*dt;	
	p->vx = yst;
	p->vy = -xst*2. -3./2.*aO;
}

// This function evolves a particle under the operator
// Phi exactly up to machine precission.
void operator_phi(double dt, struct particle* p){
	// The force used here is for test cases 2 and 3 
	// in Rein & Tremaine 2011. 
	const double r = sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
	const double forcex = -1./(r*r*r)*p->x;
	const double forcey = -1./(r*r*r)*p->y;
	const double forcez = -1./(r*r*r)*p->z;

	p->vx += forcex * dt;
	p->vy += forcey * dt;
	p->vz += forcez * dt;
}

// This function evolves a particle under the operator
// H_Kin exactly up to machine precission (drift).
void operator_HKin(double dt, struct particle* p){
	p->x += p->vx * dt;
	p->y += p->vy * dt;
	p->z += p->vz * dt;
}

// This function is the SEI integrator.
// It is symplectic, second order accurate and time-reversible.
void operator_sei(double dt, struct particle* p){
	operator_H0(dt/2.,p);
	operator_phi(dt,p);
	operator_H0(dt/2.,p);
}

// This function corresponds to the modified 
// Kick operator for leap-frog. It is not governed by a Hamiltonian.
void operator_kick(double dt, struct particle* p){
	const double r = sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
	const double forcex = -1./(r*r*r)*p->x;
	const double forcey = -1./(r*r*r)*p->y;
	const double forcez = -1./(r*r*r)*p->z;

	const double vxn = p->vx;
	p->vx += (3.*OMEGA*OMEGA*p->x + 2.*OMEGA*p->vy + forcex)* dt;
	p->vy += (-2.*OMEGA*vxn + forcey) * dt;
	p->vz += (-OMEGA*OMEGA*p->z + forcez) * dt;
}

// This function is the classical leap-frog integrator.
// It is NOT symplectic, second order accurate and time-reversible.
void operator_lf(double dt, struct particle* p){
	operator_HKin(dt/2.,p);
	operator_kick(dt,p);
	operator_HKin(dt/2.,p);
}

// This function is the Quinn et al. integrator.
void operator_quinn(double dt, struct particle* p){
	//Kick 1
	double r = sqrt(p->x*p->x + p->y*p->y + p->z*p->z);
	double forcex = -1./(r*r*r)*p->x;
	double forcey = -1./(r*r*r)*p->y;
	double forcez = -1./(r*r*r)*p->z;
	const double vxn14 	= p->vx - 0.5*dt*(OMEGA*OMEGA*p->x-forcex);
	const double Pyn	= p->vy + 2.*OMEGA*p->x + 0.5*dt*forcey; // This variable needs to be saved for Kick 2.
	p->vx	= vxn14 + dt*OMEGA*Pyn;
	p->vy 	= Pyn - OMEGA*p->x - OMEGA*(p->x + dt*p->vx);
	p->vz  += 0.5*dt*(-OMEGA*OMEGA*p->z + forcez);

	// Drift
	operator_HKin(dt,p);

	//Kick 2
	r = sqrt(p->x*p->x + p->y*p->y + p->z*p->z); // This doesn't have to be calculated twice per timestep.
	forcex = -1./(r*r*r)*p->x;
	forcey = -1./(r*r*r)*p->y;
	forcez = -1./(r*r*r)*p->z;
	const double vxn34	= p->vx + dt*OMEGA*Pyn;
	p->vx 	= vxn34 -0.5*dt*(OMEGA*OMEGA*p->x - forcex);
	p->vy 	= Pyn - 2.*OMEGA*p->x + 0.5*dt*forcey;
	p->vz  += 0.5*dt*(-OMEGA*OMEGA-p->z + forcez);
}

// Setup a test case (perturbed epicyclic motion)
void setup(struct particle* p){
	const double impactparameter = 8.0*0.69336127435063; // 8 Hill radii
	p->x  = impactparameter;
	p->y  = M_PI*300./2.*OMEGA*impactparameter;
	p->z  = 0.; 
	p->vx = 0;
	p->vy = -3./2.*OMEGA*impactparameter;
	p->vz = 0.; 
}

// Run a test case and output the positions
/*
int main(){
	double time = 0;
	const double dt = 1e-1*2.*M_PI;
	struct particle p_sei;
	struct particle p_lf;
	struct particle p_quinn;
	setup(&p_sei);
	setup(&p_lf);
	setup(&p_quinn);
	while (time<200.*M_PI){
		// Output all positions every timestep.
		printf("%e\t",time);
		printf("%e\t%e\t%e\t",p_sei.x,p_sei.y,p_sei.z);
		printf("%e\t%e\t%e\t",p_lf.x,p_lf.y,p_lf.z);
		printf("%e\t%e\t%e\t",p_quinn.x,p_quinn.y,p_quinn.z);
		printf("\n");
		// Integrate
		operator_sei(dt,&p_sei);
		operator_lf(dt,&p_lf);
		operator_quinn(dt,&p_quinn);
		time+=dt;
	}
	return 1;
}

*/
