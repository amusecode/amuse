#include "Diagnostics.h"

// Conserved quantities

double Diagnostics::get_mass(vector<double> &data) {
  int N = data.size()/7;
  double M = 0;
  for(int i=0; i<N; i++) {
    M += data[i*7];
  }
  return M;
}
vector<double> Diagnostics::get_rcm(vector<double> &data) {
  double M = get_mass(data);
  vector<double> rcm(3,0);
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      rcm[j] += data[i*7]*data[i*7+(j+1)];
    }
  }
  for(int i=0; i<3; i++) rcm[i] /= M;
  return rcm;
}
vector<double> Diagnostics::get_vcm(vector<double> &data) {
  double M = get_mass(data);
  vector<double> vcm(3,0);
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      vcm[j] += data[i*7]*data[i*7+(j+4)];
    }
  }
  for(int i=0; i<3; i++) vcm[i] /= M;
  return vcm;
}
vector<double> Diagnostics::get_lcm(vector<double> &data) {
  double M = get_mass(data);
  vector<double> lcm(3,0);
  int N = data.size()/7;
  for(int i=0; i<N; i++) {
    lcm[0] += data[i*7]*(data[i*7+2]*data[i*7+6]-data[i*7+3]*data[i*7+5]);
    lcm[1] += data[i*7]*(data[i*7+3]*data[i*7+4]-data[i*7+1]*data[i*7+6]);
    lcm[2] += data[i*7]*(data[i*7+1]*data[i*7+5]-data[i*7+2]*data[i*7+4]);
  }
  for(int i=0; i<3; i++) lcm[i] /= M;
  return lcm;
}
double Diagnostics::get_kinetic_energy(vector<double> &data) {
  int N = data.size()/7;
  double ek = 0;
  for(int i=0; i<N; i++) {
    double m  = data[i*7];
    double vx = data[i*7+4];
    double vy = data[i*7+5];
    double vz = data[i*7+6];
    double v2 = vx*vx + vy*vy + vz*vz;
    ek += 0.5*m*v2;
  }
  return ek;
}
double Diagnostics::get_potential_energy(vector<double> &data) {
  int N = data.size()/7;
  double ep = 0;
  for(int i=0; i<N-1; i++) {
    double mi = data[i*7];
    double xi = data[i*7+1];
    double yi = data[i*7+2];
    double zi = data[i*7+3];
    for(int j=i+1; j<N; j++) {
      double mj = data[j*7];
      double xj = data[j*7+1];
      double yj = data[j*7+2];
      double zj = data[j*7+3];

      double dx = xj - xi;
      double dy = yj - yi;
      double dz = zj - zi;
      double dr2 = dx*dx + dy*dy + dz*dz;
      ep -= mi*mj/sqrt(dr2);
    }
  }
  return ep;
}
double Diagnostics::get_energy(vector<double> &data) {
  double ek = get_kinetic_energy(data);
  double ep = get_potential_energy(data);
  return ek+ep;
}

// System properties

double Diagnostics::get_virial_radius(vector<double> &data) {
  double M = get_mass(data);
  double ep = get_potential_energy(data);
  double rv = -1*M*M / (2*ep);
  return rv;
}
double Diagnostics::get_harmonic_radius(vector<double> &data) {
  int N = data.size()/7;

  double ep = 0;
  for(int i=0; i<N-1; i++) {
    double xi = data[i*7+1];
    double yi = data[i*7+2];
    double zi = data[i*7+3];
    for(int j=i+1; j<N; j++) {
      double xj = data[j*7+1];
      double yj = data[j*7+2];
      double zj = data[j*7+3];

      double dx = xj - xi;
      double dy = yj - yi;
      double dz = zj - zi;
      double dr2 = dx*dx + dy*dy + dz*dz;
      ep -= 1.0/sqrt(dr2);
    }
  }
  ep /= (N*N);

  double M = get_mass(data);
  double rh = -1*M*M / (2*ep);
  return rh;
}
double Diagnostics::get_velocity_disperion(vector<double> &data) {
  double M = get_mass(data);
  double ek = get_kinetic_energy(data);
  double sigma = sqrt(2*ek / M);
  return sigma;  
}





