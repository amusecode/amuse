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

int n_local, n_total, offset;
vector<double> m, x, y, z; // Only contains particles assigned to this process
vector<double> m_all, x_all, y_all, z_all;


// Control parameters:

double eps2 = 0;


int handle_result(int result){
    int local_result = result;
#ifndef NOMPI
    if (mpi_rank) {
        MPI_Reduce(&local_result, NULL, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, &local_result, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif
    return local_result;
}

int initialize_code() {
    cerr << "initialize_code" << endl;
#ifndef NOMPI
    int error;
    error = handle_result(MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank));
    if(error) {
        cerr << "MPI_Comm_rank returned: " << error << endl;
        return -1;
    }
    error = handle_result(MPI_Comm_size(MPI_COMM_WORLD, &mpi_size));
    if(error) {
        cerr << "MPI_Comm_size returned: " << error << endl;
        return -1;
    }
#endif
    n_total = 0;
    n_local = 0;
#ifdef GPU
    return handle_result(cuda_initialize_code());
#else
    return 0;
#endif
}

int commit_parameters() {
    cerr << "commit_parameters" << endl;
    return 0;
}
int recommit_parameters() {
    return 0;
}
int commit_particles() {
    cerr << "commit_particles" << endl;
    n_total = m_all.size();
#ifdef NOMPI
    n_local = n_total;
    m = m_all;
    x = x_all;
    y = y_all;
    z = z_all;
#else
    n_local = n_total / mpi_size + (n_total % mpi_size > mpi_rank);
    offset = (n_total / mpi_size) * mpi_rank + min(n_total % mpi_size, mpi_rank);
    m.assign(m_all.begin()+offset, m_all.begin()+offset+n_local);
    x.assign(x_all.begin()+offset, x_all.begin()+offset+n_local);
    y.assign(y_all.begin()+offset, y_all.begin()+offset+n_local);
    z.assign(z_all.begin()+offset, z_all.begin()+offset+n_local);
#endif
#ifdef GPU
    int result = handle_result(cuda_commit_particles(m, x, y, z));
#else
    int result = 0;
#endif
    cerr << "commit_particles done" << endl;
    return result;
}

int recommit_particles() {
    cerr << "recommit_particles" << endl;
#ifdef GPU
    int result = handle_result(cuda_cleanup_code());
    if (result < 0) {
        cerr << "recommit_particles: cuda_cleanup_code returned: " << result << endl;
        return result;
    }
#endif
    return commit_particles();
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
    m_all.clear();
    x_all.clear();
    y_all.clear();
    z_all.clear();
#ifdef GPU
    return handle_result(cuda_cleanup_code());
#else
    return 0;
#endif
}

int new_particle(int *id, double *mass_in, double *x_in, double *y_in, double *z_in, int length) {
    int i, next_id = m_all.size();
    for (i=0; i<length; i++) {
        id[i] = next_id++;
    }
    m_all.insert(m_all.end(), mass_in, mass_in+length);
    x_all.insert(x_all.end(), x_in, x_in+length);
    y_all.insert(y_all.end(), y_in, y_in+length);
    z_all.insert(z_all.end(), z_in, z_in+length);
    return 0;
}

int delete_particle(int *id, int length) {
    return -1;
}
int get_mass(int *id, double * out, int length) {
    int i;
    for (i=0; i<length; i++) {
	    if (id[i] < 0 || id[i] > m_all.size()) {
		    return -1;
	    }
    	    out[i] = m_all[id[i]];
    }
    return 0;
}
int set_mass(int *id, double * in, int length) {
    int i;
    for (i=0; i<length; i++) {
	    if (id[i] < 0 || id[i] > m_all.size()) {
		    return -1;
	    }
            m_all[id[i]] = in[i];
    }
    return 0;
}

int local_get_potential_at_point(double *eps_in, double *x_in, double *y_in, double *z_in, 
        double *phi, int length){
#ifdef GPU
    int result = handle_result(cuda_get_potential_at_point(eps2, eps_in, x_in, y_in, z_in, phi, length));
#else
    int result = 0;
    double dx, dy, dz, r;
    double dr2, eps2_total;
    for (int j = 0; j < length; j++) {
        eps2_total = eps2 + eps_in[j]*eps_in[j];
        phi[j] = 0;
        for (int i = 0; i < n_local; i++) {
            dx = x[i] - x_in[j];
            dy = y[i] - y_in[j];
            dz = z[i] - z_in[j];
            dr2 = dx*dx + dy*dy + dz*dz;
            if (dr2 > 0 && m[i] > 0) {
                r = sqrt(dr2 + eps2_total);
                phi[j] -= m[i]/r;
            }
        }
    }
#endif
    return result;
}

int get_potential_at_point(double *eps_in, double *x_in, double *y_in, double *z_in, 
        double *phi, int length){
    int result = local_get_potential_at_point(eps_in, x_in, y_in, z_in, phi, length);
#ifndef NOMPI
    if (mpi_rank) {
        MPI_Reduce(phi, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, phi, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif
    return result;
}

int get_gravity_at_point(double *eps_in, double *x_in, double *y_in, double *z_in, 
        double *ax, double *ay, double *az, int length){
    double dx, dy, dz, r2, tmp;
    
#ifdef GPU
    cerr << "get_gravity_at_point: start cuda_get_gravity_at_point with " << length << " points and " << n_local << " particles." << endl;
    int result = handle_result(cuda_get_gravity_at_point(eps2, eps_in, x_in, y_in, z_in, ax, ay, az, length));
#else
    double dr2, eps2_total;
    for (int j = 0; j < length; j++) {
        eps2_total = eps2 + eps_in[j]*eps_in[j];
        ax[j] = 0;
        ay[j] = 0;
        az[j] = 0;
        for (int i = 0; i < n_local; i++) {
            dx = x[i] - x_in[j];
            dy = y[i] - y_in[j];
            dz = z[i] - z_in[j];
            dr2 = dx*dx + dy*dy + dz*dz;
            if (dr2 > 0) {
                r2 = dr2 + eps2_total;
                tmp = m[i] / (r2 * sqrt(r2));
                ax[j] += tmp * dx;
                ay[j] += tmp * dy;
                az[j] += tmp * dz;
            }
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

int get_potential_energy(double *potential_energy){
    cerr << "get_potential_energy" << endl;
    int i;
    
    // First calculate local potential energy:
#ifdef GPU
    int result = handle_result(cuda_get_potential_energy(eps2, potential_energy));
    if (result < 0) return -1;
#else
    double dx, dy, dz, dr2, r;
    *potential_energy = 0;
    for (int j = 1; j < n_local; j++) {
        for (i = 0; i < j; i++) {
            dx = x[i] - x[j];
            dy = y[i] - y[j];
            dz = z[i] - z[j];
            dr2 = dx*dx + dy*dy + dz*dz;
            r = sqrt(dr2 + eps2);
            *potential_energy -= m[i]*m[j]/r;
        }
    }
#endif

// Now calculate contributions to the potential energy from particles on other processes
#ifndef NOMPI
    int n_others = n_total - n_local;
    if (n_others > 0){
        double *eps_in = new double[n_others];
        double *x_in = new double[n_others];
        double *y_in = new double[n_others];
        double *z_in = new double[n_others];
        double *mass_others = new double[n_others];
        double *phi = new double[n_others];
        for (i = 0; i < offset; i++) {
            eps_in[i] = 0;
            x_in[i] = x_all[i];
            y_in[i] = y_all[i];
            z_in[i] = z_all[i];
            mass_others[i] = m_all[i];
        }
        int j = offset;
        for (i = offset+n_local; i < n_total; i++, j++) {
            eps_in[j] = 0;
            x_in[j] = x_all[i];
            y_in[j] = y_all[i];
            z_in[j] = z_all[i];
            mass_others[j] = m_all[i];
        }
        int result = local_get_potential_at_point(eps_in, x_in, y_in, z_in, phi, n_others);
        if (result < 0) return -1;
        double potential_energy_others = 0;
        for (i = 0; i < n_others; i++) {
            potential_energy_others += mass_others[i]*phi[i];
        }
        *potential_energy += 0.5 * potential_energy_others;
        delete[] eps_in;
        delete[] x_in;
        delete[] y_in;
        delete[] z_in;
        delete[] mass_others;
        delete[] phi;
    }

    if (mpi_rank) {
        MPI_Reduce(potential_energy, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, potential_energy, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif
    cerr << "get_potential_energy done" << endl;
    return 0;
}

