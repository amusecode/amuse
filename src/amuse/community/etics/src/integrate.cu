#include <thrust/device_vector.h>
#include <thrust/inner_product.h>
#include "common.hpp"
#include "scf.hpp"
#include "integrate.hpp"

#define PTR(x) (thrust::raw_pointer_cast((x).data()))

namespace etics
{

// The 'drift' step is performed using the 'acc' member.
struct DriftFunctor {
    Real Step, Step2;
    __host__ __device__ DriftFunctor(Real _Step) : Step(_Step), Step2(_Step*_Step) {}
    __host__ __device__ Particle operator() (Particle &p) const {
        p.pos += p.vel*Step + p.acc*0.5*Step2;
        p.CalculateR2(); // needed for both MEX and SCF, but not generally needed in leapfrog
        return p;
    }
};

// The 'kick' step is performed using the 'acc' member and also the force F,
// calculated at the new (predicted) position.
struct KickFunctor {
    Real Step;
    __host__ __device__ KickFunctor(Real _Step) : Step(_Step) {}
    __host__ __device__ Particle operator() (Particle& p, const vec3& F) const {
        p.vel += (p.acc + F)*0.5*Step;
        p.acc = F;
        return p;
    }
};

struct KineticEnergyFunctor {
    __host__ __device__ Real operator() (const Particle &p) const {return 0.5*p.m*p.vel.abs2();}
};

struct PotentialEnergyFunctor {
    __host__ __device__ Real operator() (const Particle &p, const Real &Potential) const {return p.m*Potential;}
};

Integrator::Integrator() {
    N = 0;
    Time = 0;
}

Integrator::Integrator(Particle *P_h, int _N) {
    N = _N;
    Time = 0;
    P = thrust::device_vector<Particle>(P_h, P_h+N);
    Potential = thrust::device_vector<Real>(N);
    Force = thrust::device_vector<vec3>(N);
    Method = &etics::scf::CalculateGravity;
    CalculateGravity();
    KickStep(0); // Just to "commit" the forces to the particle list.
}

Integrator::~Integrator() {
    N = 0;
    // Weird way to free Thrust memory.
    P.clear();           P.shrink_to_fit();
    Potential.clear();   Potential.shrink_to_fit();
    Force.clear();       Force.shrink_to_fit();
}

void Integrator::CalculateGravity() {
    (*Method)(PTR(P), N, PTR(Potential), PTR(Force));
}

void Integrator::DriftStep(Real Step) {
    thrust::transform(P.begin(), P.end(), P.begin(), DriftFunctor(Step));
}

void Integrator::KickStep(Real Step) {
    thrust::transform(P.begin(), P.end(), Force.begin(), P.begin(), KickFunctor(Step));
}

Real Integrator::GetTime() {
    return Time;
}

int Integrator::GetN() {
    return N;
}

Real Integrator::KineticEnergy() {
    return thrust::transform_reduce(
      P.begin(), P.end(),
      KineticEnergyFunctor(),
      (Real)0, // It must be clear to the function that this zero is a Real.
      thrust::plus<Real>()
    );
}

Real Integrator::PotentialEnergy() {
    return 0.5*thrust::inner_product(
      P.begin(), P.end(),
      Potential.begin(),
      (Real)0,
      thrust::plus<Real>(),
      PotentialEnergyFunctor()
    );
}

void Integrator::CopyParticlesToHost(Particle *P_h) {
    thrust::copy(P.begin(), P.end(), P_h);
}

void Integrator::CopyParticlesToHost(Particle **P_h, int *_N) {
    Particle *LocalList = new Particle[N];
    thrust::copy(P.begin(), P.end(), LocalList);
    *P_h = LocalList;
    *_N = N;
}

} // namespace etics
