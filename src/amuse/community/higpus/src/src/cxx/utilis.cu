#include <utilis.h>
#include <cstdlib>
#include <iomanip>
#include <functions.h>

using namespace std;

vector<double> Utilis::CdM(double4 *a, int num){

  vector<double> CMV;

  cmx = cmy = cmz = 0.0;
  total_mass = 0.0;

  for(k = 0; k < num; k++)
    total_mass += a[k].w;

  total_mass = 1.0/total_mass;

  for(k = 0; k < num; k++){
    cmx += a[k].x * a[k].w;
    cmy += a[k].y * a[k].w;
    cmz += a[k].z * a[k].w;
  }

  cmx *= total_mass;
  cmy *= total_mass;
  cmz *= total_mass;

  CMV.push_back(cmx);
  CMV.push_back(cmy);
  CMV.push_back(cmz);

  for(k = 0; k < num; k++){
    a[k].x -= cmx;
    a[k].y -= cmy;
    a[k].z -= cmz;
  }

  return CMV;

}

