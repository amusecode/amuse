// Standalone Stromgren sphere test.

#include <cassert>
#include <iostream>
#include <vector>
#include <cmath>

#include <mpi.h>
#include "aton_cpp.h"


double time_step(double dx, double c_light) {
  const double cfl = 0.9;
  return cfl*dx/c_light/3;
}

double calc_dx(double L) {
  int ncellx, ncelly, ncellz, nbound;
  aton::get_grid_size(&ncellx, &ncelly, &ncellz, &nbound);
  assert(ncellx == ncelly && ncelly == ncellz);
  return L / ncellx;
}

// Coordinates are in the range [0,1].
int cell_index(double x, double y, double z) {
  int ncellx, ncelly, ncellz, nbound;
  aton::get_grid_size(&ncellx, &ncelly, &ncellz, &nbound);
  return aton::cell_index4(int(x*ncellx), int(y*ncelly), int(z*ncellz), 0);
}

void output_slice(aton::State state) {
  int ncellx, ncelly, ncellz, nbound;
  aton::get_grid_size(&ncellx, &ncelly, &ncellz, &nbound);

  for (int i = 0; i < ncellx; i++) {
    for (int j = 0; j < ncelly; j++) {
      int index = aton::cell_index4(i, j, ncellz/2, 0);
      std::cout << i << " " << j << " " << state.xHII[index] << std::endl;
    }
    std::cout << std::endl;
  }
}

double find_half_point(const std::vector<double> xHI) {
  const double threshold = 0.5;
  
  for (int i = 0; i < (int)xHI.size(); i++) {
    if (xHI[i] >= threshold) {
      if (i > 0) {
        return i-1 + (threshold - xHI[i-1]) / (xHI[i] - xHI[i-1]);
        break;
      } else {
        return 0.0;
      }
    }
  }
}


void average_spherical(aton::State state,
                       std::vector<double>* xHI,
                       std::vector<double>* xHII,
                       std::vector<double>* T) {
  int ncellx, ncelly, ncellz, nbound;
  aton::get_grid_size(&ncellx, &ncelly, &ncellz, &nbound);
  int N = ncellx;
  assert(ncellx == N && ncelly == N && ncellz == N);

  std::vector<double> one(N);
  xHI->clear();
  xHI->resize(N);
  xHII->clear();
  xHII->resize(N);
  T->clear();
  T->resize(N);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      for (int k = 0; k < N; k++) {
      	int index = aton::cell_index4(i, j, k, 0);
        double Dx = double(i) / N - 0.5;
        double Dy = double(j) / N - 0.5;
        double Dz = double(k) / N - 0.5;
        double r = std::sqrt(Dx*Dx + Dy*Dy + Dz*Dz);
        if (r > 0.5) {
          continue;
        }
        int bin = int(N*r);
        
        one[bin] += 1.0;
        (*xHI)[bin] += 1.0 - state.xHII[index];
        (*xHII)[bin] += state.xHII[index];
 	(*T)[bin] += state.T[index];
      }
    }
  }
  
  for (int i = 0; i < N; i++) {
    if (one[i] <= 0.0) {
      continue;
    }
    (*xHII)[i] /= one[i];
    (*xHI)[i] /= one[i];
    (*T)[i] /= one[i];
  }
}


void output_spherical(aton::State state) {
  std::vector<double> xHI, xHII, T;
  average_spherical(state, &xHI, &xHII, &T);  
  int N = (int)xHI.size();
  for (int i = 0; i < N; i++) {
    double r = 2 * i / double(N);
    std::cout << r << " " << xHI[i] << " " << xHII[i] << " " << T[i] << std::endl;
  }
}

double ionized_radius(aton::State state) {
  std::vector<double> xHI, xHII, T;
  average_spherical(state, &xHI, &xHII, &T);  
  int N = (int)xHI.size();
  return find_half_point(xHI) * 2 / N;
}

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);

  aton::init_cuda_device(false);
  aton::gpu_allocate();

  aton::State state;
  aton::cpu_allocate(&state);
  
  const double Myr = 3.156e13; // 1 Myr [s]
  const double kpc = 3.086e19; // 1 kpc [m] 
  
  double c_light = 2.99792458e8;
  double fudgecool = 0.1;
  double nH = 1e3;  // [m^-3]
  double T = 1e4;  // [K]
  double xHII = 1.2e-3;
  double E = 1e-10;
  double source = 5e48;  // [photons / s]
  
  std::cerr << "Enter nH, c: ";
  double c_light_factor;
  std::cin >> nH >> c_light_factor;
  c_light *= c_light_factor;
  std::cout << "# nH = " << nH << " [m^-3]" << std::endl;
  std::cout << "# c = " << c_light_factor << std::endl;
 
  const double alphaB = 2.59e-19; // [m^3 s^-1]
  double t_rec = 1.0 / (alphaB * nH); // [s]
  std::cout << "# t_rec = " << t_rec/Myr << " [Myr]" << std::endl;
  const double G = 6.67e-11; // [m^3 kg^-1 s^-2]
  const double mH = 1.67e-27; // [kg]
  double t_hubble = sqrt(3.0 / (8*M_PI*G*mH*nH)); // [s]
  std::cout << "# t_hubble = " << t_hubble/Myr << " [Myr]" << std::endl;  
  double r_s = pow(3*source / (4*M_PI*alphaB*nH*nH), 1.0/3.0); // [m]
  std::cout << "# r_s = " << r_s/kpc << " [kpc]" << std::endl;
  
  double L = 1.2 * 2 * r_s; // [m]
  std::cout << "# L = " << L/kpc << " [kpc]" << std::endl;
  double dx = calc_dx(L);
  double final_t = 5 * std::min(t_rec, t_hubble); // [s]
  std::cout << "# final_t = " << final_t/Myr << " [Myr]" << std::endl;

  double dt = time_step(dx, c_light);
  
  state.Init(E, nH, T, xHII, 0.0);
  int index = cell_index(0.5, 0.5, 0.5);
  state.photon_source[index] = source / dx / dx / dx;
  
  aton::validate(state, c_light);

  aton::cpu_to_gpu(state);

  int num_iterations = int(final_t / dt);
  std::cerr << "# num_iterations: " << num_iterations << std::endl;

  int n = aton::get_boundary_buffer_size();
  double *boundary_values = new double[6*n];

  for (int i = 0; i < num_iterations; i++) {
    if (i % 1000 == 0) {
      //std::cerr << i << std::endl;
      aton::gpu_to_cpu(state);
      double r = ionized_radius(state);
      double t = i*dt / Myr; // [Myr]
      std::cout << t << " " << r << std::endl;
    }

    aton::gpu_transport(c_light, dx, dt);
    aton::gpu_add_sources(dx, dt);
    aton::gpu_cooling(c_light, dx, dt, 1.0, 0.0, fudgecool);
   
    // Zero-gradient boundary conditions.
    aton::gpu_to_cpu_boundary_values(boundary_values);
    aton::cpu_to_gpu_boundary_values(boundary_values);
  }

  aton::gpu_to_cpu(state);

  //output_slice(state);
  //double r = output_spherical(state);
  //std::cout << r << std::endl;

  MPI_Finalize();

  return 0;
}

