#pragma once
#include <string>
#include <time.h>
// A_ON_SHARED_MEMORY moves the A structure (after it has been calculated) from constant to shared memory. This might be faster, or not.
// #define A_ON_SHARED_MEMORY
#define ETICS_HDF5

#ifndef __CUDACC__
    #define __host__
    #define __device__
#endif

#if defined ETICS_SINGLE_PRECISION && defined ETICS_DOUBLE_PRECISION
    #error Contradictory precision flags!
#endif

#ifndef ETICS_SINGLE_PRECISION
    #define ETICS_DOUBLE_PRECISION
    #define Real double
    #define MPI_ETICS_REAL MPI_DOUBLE
#else
    #define Real float
    #define MPI_ETICS_REAL MPI_FLOAT
#endif

#if defined MEX && defined SCF
    #error Contradictory method flags!
#endif

#if defined MEX && !defined LMAX
    #error LMAX not defined.
#endif

#if defined SCF && ((!defined NMAX) || (!defined LMAX))
    #error Both NMAX and LMAX should be defined!
#endif

#if defined MEX && defined NMAX
    #warning NMAX is defined, but will be ignored since the method is MEX!
#endif

struct vec3 {
    Real x, y, z;

    __host__ __device__ vec3() : x(0), y(0), z(0) {}
    __host__ __device__ vec3(Real _x, Real _y, Real _z) : x(_x), y(_y), z(_z) {}
    __host__ __device__ Real abs2() const {return x*x + y*y + z*z;}
    __host__ __device__ vec3 operator+ (const vec3& V)   const {return vec3(this->x + V.x, this->y + V.y, this->z + V.z);}
    __host__ __device__ vec3 operator- (const vec3& V)   const {return vec3(this->x - V.x, this->y - V.y, this->z - V.z);}
    __host__ __device__ vec3 operator* (const Real& C) const {return vec3(this->x*C, this->y*C, this->z*C);}
    __host__ __device__ vec3& operator+= (const vec3& V) {this->x += V.x; this->y += V.y; this->z += V.z; return *this;}
    __host__ __device__ vec3& operator-= (const vec3& V) {this->x -= V.x; this->y -= V.y; this->z -= V.z; return *this;}
    __host__ __device__ vec3 operator- () const {return vec3(-this->x, -this->y, -this->z);}
    __host__ __device__ vec3 operator+ () const {return vec3(this->x, this->y, this->z);}
};


class Particle {
    public:
        int ID;
        unsigned char Status;
        Real m;
        vec3 pos, vel, acc;
        Real R2;

        __host__ __device__ Particle() {}
        __host__ __device__ void CalculateR2() {R2 = pos.abs2();}
        __host__ __device__ bool operator< (const Particle& p) const {return (this->R2 < p.R2);}
};

#define ETICS_MY_CLOCK CLOCK_MONOTONIC
class HostTimer {
    timespec StartTime, StopTime;
  public:
    void Start() {
        clock_gettime(ETICS_MY_CLOCK, &StartTime);
    }
    void Stop() {
        clock_gettime(ETICS_MY_CLOCK, &StopTime);
    }
    double Difference() {
        timespec TimeDiff;
        if ((StopTime.tv_nsec - StartTime.tv_nsec) < 0) {
            TimeDiff.tv_sec  = StopTime.tv_sec - StartTime.tv_sec - 1;
            TimeDiff.tv_nsec = 1000000000 + StopTime.tv_nsec - StartTime.tv_nsec;
        } else {
            TimeDiff.tv_sec  = StopTime.tv_sec - StartTime.tv_sec;
            TimeDiff.tv_nsec = StopTime.tv_nsec - StartTime.tv_nsec;
        }
        return (double)TimeDiff.tv_sec + ((double)TimeDiff.tv_nsec)*1.0e-9;
    }
};

#ifdef __CUDACC__
class DeviceTimer {
    cudaEvent_t StartTime, StopTime;
  public:
    DeviceTimer() {
        cudaEventCreate(&StartTime);
        cudaEventCreate(&StopTime);
    }
    ~DeviceTimer() {
        cudaEventDestroy(StartTime);
        cudaEventDestroy(StopTime);
    }
    void Start() {
        cudaEventRecord(StartTime);
    }
    void Stop() {
        cudaEventRecord(StopTime);
        cudaEventSynchronize(StopTime);
    }
    double Difference(){
        float TimeDiff;
        cudaEventElapsedTime(&TimeDiff, StartTime, StopTime);
        return (double)TimeDiff * 1.0e-3;
    }
};
#endif
