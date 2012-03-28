#include <utilis.h>
#include <cstdlib>
#include <iomanip>

using namespace std;

vector<double> Utilis::CdM(double *x, double *y, double *z, double *w, int num){

  vector<double> CMV;

  cmx = cmy = cmz = 0.0;
  total_mass = 0.0;

  for(k = 0; k < num; k++)
    total_mass += w[k];

  total_mass = 1.0/total_mass;

  for(k = 0; k < num; k++){
    cmx += x[k] * w[k];
    cmy += y[k] * w[k];
    cmz += z[k] * w[k];
  }

  cmx *= total_mass;
  cmy *= total_mass;
  cmz *= total_mass;

  CMV.push_back(cmx);
  CMV.push_back(cmy);
  CMV.push_back(cmz);

  for(k = 0; k < num; k++){
    x[k] -= cmx;
    y[k] -= cmy;
    z[k] -= cmz;
  }

  return CMV;

}

