#include <iostream>
//~#include <string>
//~#include <cmath>
#include <vector>

using namespace std;

int cuda_initialize_code();
int cuda_commit_particles(vector<double>& m, vector<double>& x, vector<double>& y, vector<double>& z);
int cuda_recommit_particles(vector<double>& m, vector<double>& x, vector<double>& y, vector<double>& z);
int cuda_cleanup_code();
int cuda_get_potential_at_point(double eps2, double *eps, double *x, double *y, double *z, double *phi, int N);
int cuda_get_gravity_at_point(double eps2, double *eps, double *x, double *y, double *z, double *ax, double *ay, double *az, int N);
