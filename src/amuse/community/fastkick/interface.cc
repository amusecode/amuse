#ifndef NOMPI
#include <mpi.h>
#endif

#include <iostream>
#include <string>
#include <cmath>
#include <vector>
#include "worker_code.h"

#ifdef GPU
#include "cuda_fastkick.h"
#endif


using namespace std;

#ifndef NOMPI
int mpi_rank = 0;
int mpi_size = 1;
#endif

int n_local, n_total;
vector<double> m, x, y, z;

#ifdef GPU
#endif


// Control parameters:

double eps2 = 0;


int initialize_code() {
#ifndef NOMPI
    int error;
    error = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if(error) {
        cerr << "MPI_Comm_rank returned: " << error << endl;
        return -1;
    }
    error = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    if(error) {
        cerr << "MPI_Comm_size returned: " << error << endl;
        return -1;
    }
#endif
    n_total = 0;
    n_local = 0;
#ifdef GPU
    return cuda_initialize_code();
#else
    return 0;
#endif
}

int commit_parameters() {
    return 0;
}
int recommit_parameters() {
    return 0;
}
int commit_particles() {
#ifdef GPU
    return cuda_commit_particles(m, x, y, z);
#else
    return 0;
#endif
}
int recommit_particles() {
#ifdef GPU
    return cuda_recommit_particles(m, x, y, z);
#else
    return 0;
#endif
}

int get_eps2(double *_epsilon_squared){
    *_epsilon_squared = eps2;
    return 0;
}
int set_eps2(double _epsilon_squared){
    eps2 = _epsilon_squared;
    return 0;
}

int cleanup_code() {
    m.clear();
    x.clear();
    y.clear();
    z.clear();
#ifdef GPU
    return cuda_cleanup_code();
#else
    return 0;
#endif
}

int new_particle(int *id, double *mass_in, double *x_in, double *y_in, double *z_in, int length) {
    int i, i_in, next_id;
    int n_total_old, n_local_old;
    n_total_old = n_total;
    n_local_old = n_local;
    n_total += length;
#ifdef NOMPI
    n_local = n_total;
#else
    n_local += length / mpi_size + (length % mpi_size > mpi_rank);
#endif
    m.resize(n_local);
    x.resize(n_local);
    y.resize(n_local);
    z.resize(n_local);
#ifdef NOMPI
    for (i=n_local_old, i_in=0; i<n_local; i++, i_in++) {
#else
    for (i=n_local_old, i_in=mpi_rank; i<n_local; i++, i_in+=mpi_size) {
#endif
        m[i] = mass_in[i_in];
        x[i] = x_in[i_in];
        y[i] = y_in[i_in];
        z[i] = z_in[i_in];
    }
    for (i=0, next_id=n_total_old; i<length; i++, next_id++) {
        id[i] = next_id;
    }
    return 0;
}

int delete_particle(int *id, int length) {
    return -1;
}

int get_potential_at_point(double *eps_in, double *x_in, double *y_in, double *z_in, 
        double *phi, int length){
    double dx, dy, dz, r;
    
#ifdef GPU
    int result = cuda_get_potential_at_point(eps2, eps_in, x_in, y_in, z_in, phi, length);
    // communicate errors for MPI
#else
    for (int j = 0; j < length; j++) {
        phi[j] = 0;
        for (int i = 0; i < n_local; i++) {
            dx = x[i] - x_in[j];
            dy = y[i] - y_in[j];
            dz = z[i] - z_in[j];
            r = sqrt(dx*dx + dy*dy + dz*dz + eps2 + eps_in[j]*eps_in[j]);
            phi[j] -= m[i]/r;
        }
    }
#endif
    
#ifndef NOMPI
    if (mpi_rank) {
        MPI_Reduce(phi, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, phi, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

#ifdef GPU
    return result;
#else
    return 0;
#endif
}

int get_gravity_at_point(double *eps_in, double *x_in, double *y_in, double *z_in, 
        double *ax, double *ay, double *az, int length){
    double dx, dy, dz, r2, tmp;
    
#ifdef GPU
    int result = cuda_get_gravity_at_point(eps2, eps_in, x_in, y_in, z_in, ax, ay, az, length);
    // communicate errors for MPI
#else
    for (int j = 0; j < length; j++) {
        ax[j] = 0;
        ay[j] = 0;
        az[j] = 0;
        for (int i = 0; i < n_local; i++) {
            dx = x[i] - x_in[j];
            dy = y[i] - y_in[j];
            dz = z[i] - z_in[j];
            r2 = (dx*dx + dy*dy + dz*dz + eps2 + eps_in[j]*eps_in[j]);
            tmp = m[i] / (r2 * sqrt(r2));
            ax[j] += tmp * dx;
            ay[j] += tmp * dy;
            az[j] += tmp * dz;
        }
    }
#endif
    
#ifndef NOMPI
    if (mpi_rank) {
        MPI_Reduce(ax, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(ay, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(az, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, ax, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, ay, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, az, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

#ifdef GPU
    return result;
#else
    return 0;
#endif
}

