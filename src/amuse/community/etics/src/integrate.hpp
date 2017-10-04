#pragma once
#include <thrust/device_vector.h>

namespace etics {
    class Integrator {
      public:
        Integrator();
        Integrator(Particle *P_h, int _N);
        ~Integrator();
        void CalculateGravity();
        void DriftStep(Real Step);
        void KickStep(Real Step);
        Real GetTime();
        int GetN();
        Real KineticEnergy();
        Real PotentialEnergy();
        void CopyParticlesToHost(Particle *P_h);
        void CopyParticlesToHost(Particle **P_h, int *_N);
      private:
        int N;
        double Time;
        thrust::device_vector<Particle> P;
        thrust::device_vector<Real> Potential;
        thrust::device_vector<vec3> Force;
        void (*Method)(Particle*, int, Real*, vec3*);
    };
}
