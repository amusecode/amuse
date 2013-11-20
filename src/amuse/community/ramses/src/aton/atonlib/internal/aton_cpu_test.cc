#include <cstdio>
#include <cstdlib>
using namespace std;

void ComputeTempSingleCell(double *cuxion,
			   double *cudensity,
			   double *cutemperature,
			   double *cuegy_new,
			   double *cuflx_new_x,
			   double *cuflx_new_y,
			   double *cuflx_new_z,
			   double fudgecool,
			   double c,
			   double dt,
			   double aexp,
			   double hubblet,
			   int *num_iterations);

int main() {
  double nH = 1.008625e+03;
  double T = 1.088057e+03;
  double x = 2.279228e-02;
  double N = 8.190702e+02;
  double Fx = 7.021692e+08;
  double Fy = -9.628663e-04;
  double Fz = -1.285492e+09;
  int num_iterations;

  const double fudgecool = 0.1;
  const double c = 0.01 * 2.99792458e8;
  const double dt = 35805249191221.930;
  const double aexp = 8.68396887149750774E-002;
  const double hubblet = 0.0;

  ComputeTempSingleCell(&x, &nH, &T, &N, &Fx, &Fy, &Fz,
			fudgecool, c, dt, aexp, hubblet,
			&num_iterations);

  return 0;
}
