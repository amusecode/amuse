#include <iostream>

#include "src/common.hpp"
#include "src/scf.hpp"
#include "src/integrate.hpp"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/inner_product.h>

using namespace std;
using namespace etics;

Integrator IntegratorObj;

// GLOBAL VARIABLES
Real ConstantStep = 0.001953125;
Real T, Step, dT1, dT2, Tcrit;
int N;
extern Real mass;
extern int k3gs, k3bs, k4gs, k4bs;

/*extern*/ Particle *hostP;
thrust::host_vector<Particle> PPP;
/*extern*/ thrust::device_vector<Particle> PPPPP;

// /*extern*/   thrust::device_vector<vec3>   F0xxxxxx;
// /*extern*/   thrust::device_vector<Real>   PotPotPot; // ugly name
// /*extern*/   thrust::device_vector<vec3>   F1xxxxxx;
/*extern*/ vec3 *F1_ptr;



// extern Particle *P_h;
// extern thrust::device_vector<Particle> P;
// 
// extern thrust::device_vector<vec3>   F0;
// extern thrust::device_vector<Real>   Potential;
// extern thrust::device_vector<vec3>   F1;


void CommitParticles();
// void InitSCF(int N);
// void ForceSCF(int N, Real *Potential, Particle *PPPPP, vec3 *F);
void DriftStep();
void KickStep();
void CommitForces();
int InitilizeIntegratorMemory();
#define PTR(x) (thrust::raw_pointer_cast((x).data()))


int initialize_code() {
#warning initscf should be here!!! just the problem is that N is requires, so fix it
    return 0;
}

int recommit_parameters() {
    return 0;
}

int commit_parameters() {
    return 0;
}

int new_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius) {
    Particle p;
    p.m = mass;
    p.pos = vec3(x, y, z);
    p.vel = vec3(vx, vy, vz);
    PPP.push_back(p);
    *id = N;
    N++;
    return 0;
}

int commit_particles() {
//     cerr << "calling commit_particles" << endl;
    cerr << "we did commit_particles()" << endl;
    etics::scf::Init(N, 180, 64, 2605, 384);
#warning hardcoded launch configuration
    IntegratorObj = Integrator(&PPP[0], N);
    return 0;
}

struct CenterMassFunctor {
    Real ConstantMass;
    __host__ __device__ CenterMassFunctor(Real _ConstantMass) : ConstantMass(_ConstantMass) {}
    __host__ __device__ vec3 operator() (const Particle &p) const {return p.pos;}
};

struct ShiftFunctor {
    vec3 Shift;
    __host__ __device__ ShiftFunctor(vec3 _Shift) : Shift(_Shift) {}
    __host__ __device__ Particle operator() (Particle &p) const {
        p.pos += Shift;
        p.CalculateR2();
        return p;
    }
};

int counttt = 0;
bool FirstStep = true;
int evolve_model(double t) {
//     PPPPP = PPP;
    cerr << "call evolve_model t_end = " << t << " dt = " << t - T  << "****************" << endl;

//     vec3 CenterMass = thrust::transform_reduce(PPPPP.begin(), PPPPP.end(), CenterMassFunctor(mass), vec3(0,0,0), thrust::plus<vec3>());
//     CenterMass = CenterMass * (1.0/N); //ugly should divide by the total mass
//     cerr << "CENTER OF MASS " << CenterMass.x << endl;
//     
// //     thrust::transform(PPPPP.begin(), PPPPP.end(), PPPPP.begin(), ShiftFunctor(-CenterMass));
//     
//     vec3 CenterMass2 = thrust::transform_reduce(PPPPP.begin(), PPPPP.end(), CenterMassFunctor(mass), vec3(0,0,0), thrust::plus<vec3>());
//     CenterMass2 = CenterMass2 * (1.0/N); //ugly should divide by the total mass
//     cerr << "CENTER OF MASS after correction " << CenterMass2.x << endl;
// 

    Step = ConstantStep;
    while (T <= t) {
        // Take the drift step.
        IntegratorObj.DriftStep(Step);

        // Calculate the forces in the new positions.
//         ForceSCF(N, PTR(PotPotPot), PTR(PPPPP), PTR(F1xxxxxx));
        IntegratorObj.CalculateGravity();

        // Finish by taking the kick step.
        // The kick functor also "commits" the predicted forces into the "acc" member.
        IntegratorObj.KickStep(Step);

        // N particles were implicitly propagated in this iteration.

        // Advance global time.
        T += Step;
    }
// 
//     vec3 CenterMass3 = thrust::transform_reduce(PPPPP.begin(), PPPPP.end(), CenterMassFunctor(mass), vec3(0,0,0), thrust::plus<vec3>());
//     CenterMass3 = CenterMass3 * (1.0/N); //ugly should divide by the total mass
//     cerr << "CENTER OF MASS after evolve " << CenterMass3.x << endl;
//     
//     cerr << "done evolve; transform" << endl;
// //     thrust::transform(PPPPP.begin(), PPPPP.end(), PPPPP.begin(), ShiftFunctor(+CenterMass)); // antishift
//     
//     vec3 CenterMass4 = thrust::transform_reduce(PPPPP.begin(), PPPPP.end(), CenterMassFunctor(mass), vec3(0,0,0), thrust::plus<vec3>());
//     CenterMass4 = CenterMass4 * (1.0/N); //ugly should divide by the total mass
//     cerr << "CENTER OF MASS after antishift " << CenterMass4.x << endl;
// 
//     cerr << "done transform; download to RAM" << endl;
    IntegratorObj.CopyParticlesToHost(&PPP[0]);
//     
//     cerr << "done download; return" << endl;
    return 0;
}

int set_begin_time(double time_begin) {
//     cerr << "called set_begin_time(" << time_begin << endl;
    return 0;
}


int get_begin_time(double *time_begin) {
    *time_begin = 0;
    return 0;
}

int get_mass(int index_of_the_particle, double *mass) {
    *mass = PPP[index_of_the_particle].m;
    return 0;
}

int get_time(double *time) {
    *time = T;
    return 0;
}

int set_mass(int index_of_the_particle, double mass) {
//     cerr << "calling set_mass" << endl;
    PPP[index_of_the_particle].m = mass;
    return 0;
}

int get_index_of_first_particle(int *index_of_the_particle) {
//     cerr << "calling get_index_of_first_particle" << endl;
    *index_of_the_particle = 0;
    return 0;
}

int get_total_radius(double *radius) {
    return -2;
}

int get_potential_at_point(double soft, double x, double y, double z, double *phi) {
    return -2;
}

int get_total_mass(double *mass) {
    return -2;
}

int set_eps2(double epsilon_squared) {
    return -1;
}

int get_eps2(double *epsilon_squared) {
    *epsilon_squared = 0;
    return -1;
}

int get_number_of_particles(int *number_of_particles) {
//     cerr << "calling get_number_of_particles" << endl;
    *number_of_particles = PPP.size();
    return 0;
}

int get_index_of_next_particle(int index_of_the_particle, int *index_of_the_next_particle) {
    *index_of_the_next_particle = index_of_the_particle + 1;
    return 0;
}

int delete_particle(int index_of_the_particle) {
    return -2;
}

int get_potential(int index_of_the_particle, double *potential) {
    return -2;
}

int synchronize_model() {
//     cerr << "calling synchronize_model" << endl;
    return 0;
}

int set_state(int index_of_the_particle, double mass, double radius, double x, double y, double z, double vx, double vy, double vz) {
    cerr << "calling set_state" << endl;
//     cerr << "calling set_state" << endl;
    PPP[index_of_the_particle].pos = vec3(x, y, z);
    PPP[index_of_the_particle].vel = vec3(vx, vy, vz);
    return 0;
}

int get_state(int index_of_the_particle, double *mass, double *radius, double *x, double *y, double *z, double *vx, double *vy, double *vz) {
//     cerr << "calling get_state" << endl;
    Particle p = PPP[index_of_the_particle];
    *mass = index_of_the_particle;
    *x = p.pos.x;
    *y = p.pos.y;
    *z = p.pos.z;
    *vx = p.vel.x;
    *vy = p.vel.y;
    *vz = p.vel.z;
    return 0;
}

int get_time_step(double *time_step) {
//     cerr << "calling get_time_step" << endl;
    *time_step = ConstantStep;
    return 0;
}

int set_time_step(double time_step) {
    cerr << "calling set_time_step" << endl;
    ConstantStep = time_step;
    return 0;
}

int get_launch_config(int **launch_config) {
    return 0;
}

int set_launch_config(int *launch_config) {
//     k3gs = launch_config[0];
//     k3bs = launch_config[1];
//     k4gs = launch_config[2];
//     k4bs = launch_config[3];
    return -2;
}

int recommit_particles() {
//     cerr << "calling recommit_particles" << endl;
#warning put something here
    cerr << "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhh" << endl;
    PPPPP = PPP;
    return -2;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, double az) {
    return -2;
}

int get_center_of_mass_position(double *x, double *y, double *z) {
//     vec3 CenterMass = thrust::transform_reduce(PPPPP.begin(), PPPPP.end(), CenterMassFunctor(mass), vec3(0,0,0), thrust::plus<vec3>());
//     CenterMass = CenterMass * (1.0/N); //ugly should divide by the total mass
//     *x = CenterMass.x;
//     *y = CenterMass.y;
//     *z = CenterMass.z;
    return 0;
}

int get_center_of_mass_velocity(double *vx, double *vy, double *vz) {
    return -2;
}

int get_radius(int index_of_the_particle, double *radius) {
    *radius = 0;
    return 0;
}

int set_radius(int index_of_the_particle, double radius) {
    // should store the radius somewhere but completely ignored by code
//     cerr << "calling set_radius" << endl;
    return 0;
}

int cleanup_code() {
    IntegratorObj.~Integrator();
    cerr << "bye" << endl;
    return 0;
}

int get_gravity_at_point(double soft, double x, double y, double z, double *forcex, double *forcey, double *forcez) {
   return -2;
}

int get_velocity(int index_of_the_particle, double *vx, double *vy, double *vz) {
    *vx = PPP[index_of_the_particle].vel.x;
    *vy = PPP[index_of_the_particle].vel.y;
    *vz = PPP[index_of_the_particle].vel.z;
    return 0;
}

int get_position(int index_of_the_particle, double *x, double *y, double *z) {
    *x = PPP[index_of_the_particle].pos.x;
    *y = PPP[index_of_the_particle].pos.y;
    *z = PPP[index_of_the_particle].pos.z;
    return 0;
}

bool already_printed = false;

int set_position(int index_of_the_particle, double x, double y, double z) {
    if (already_printed == false) {
        cerr << "calling set_position" << endl;
        cerr << "---------index_of_the_particle=" << index_of_the_particle << endl;
        cerr << "--------- x" << PPP[index_of_the_particle].pos.x << "--->" << x << endl;
        cerr << "--------- y" << PPP[index_of_the_particle].pos.y << "--->" << y << endl;
        cerr << "--------- z" << PPP[index_of_the_particle].pos.z << "--->" << z << endl;
        already_printed = true;
    }
    PPP[index_of_the_particle].pos = vec3(x, y, z);
    counttt++;
    return 0;
}

int get_acceleration(int index_of_the_particle, double *ax, double *ay, double *az) {
    return -2;
}

int set_velocity(int index_of_the_particle, double vx, double vy, double vz) {
//     cerr << "calling set_velocity" << endl;
    PPP[index_of_the_particle].vel = vec3(vx, vy, vz);
    return 0;
}


int get_kinetic_energy(double *kinetic_energy) {
    *kinetic_energy = IntegratorObj.KineticEnergy();
    return 0;
}

int get_potential_energy(double *potential_energy) {
    *potential_energy = IntegratorObj.PotentialEnergy();
    return 0;
}

int update_force_potential_arrays(double tttt) {
#warning time shouldnt be a parameter to this one
//     ForceSCF(N, PTR(PotPotPot), PTR(PPPPP), PTR(F0xxxxxx));
    return 0;
}