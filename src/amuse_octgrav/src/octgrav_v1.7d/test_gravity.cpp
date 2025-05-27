#include "octgrav.h"
#include <math.h>

int main(int argc, char* argv[]) {

  cerr << "n= " << argc << ": " << argv[0] << "  " << argv[1] << endl;
  int n_bodies = 32768;
  if (argc > 1) {
    n_bodies = atoi(argv[1]);
    cerr << "n_bodies= " << n_bodies << endl;
  } 

  float eps   = 0.0;
  double theta = 0.5;

  vector<float4> bodies_pos(n_bodies);
  vector<float4> bodies_grav(n_bodies);
  vector<float4> bodies_grav_ex(n_bodies);


  for (int i = 0; i < n_bodies; i++) {

    double tm = pow(drand48(),2./3.);
    double r  = sqrt(tm);
    double P  =-1./sqrt(tm+1);
    double fr = P*P*P;
    double rh = 0.75/M_PI*pow(-P,5u);
    double phi           = 2*M_PI*drand48();
    double cth = 2*drand48()-1;
    double R   = r*sqrt(1.-cth*cth);
    double sph           = sin(phi);
    double cph           = cos(phi);

    bodies_pos[i].z = r * cth;
    bodies_pos[i].x = R * sph;
    bodies_pos[i].y = R * cph;
    bodies_pos[i].w = 1.0/n_bodies;
    


  }  

  if (argc > 2)
    theta = atof(argv[2]);
  fprintf(stderr, "theta= %g\n", theta);
  
 
  for (int kk = 0; kk < 1; kk++) {
    octgrav system;
    
    system.set_softening(eps);
    system.set_opening_angle(theta);
    double t00 = get_time();
    system.evaluate_gravity(bodies_pos, bodies_grav);
    double dt = get_time() - t00;
    
    cerr << "It took " << dt << " seconds\n";
  }

//   system.set_opening_angle(0.00001);
//   system.evaluate_gravity(bodies_pos, bodies_grav_ex);

  double fcom[3] = {0,0,0};
  double tcom[3] = {0,0,0}; 
  double tot_mass = 0;

  double fcomt[3] = {0,0,0};
  double tcomt[3] = {0,0,0}; 
  double4 tot_force  = {0,0,0,0};
  double4 tot_torque = {0,0,0,0};

  double tot_pot = 0.0;
  for (int i = 0; i < n_bodies; i++) {
    fprintf(stdout, "%g %g\n",
// #ifdef EXACT
//   	    (float)acc[i].w,
// 	    (float)sqrt(acc[i].x*acc[i].x +
// 			acc[i].y*acc[i].y +
// 			acc[i].z*acc[i].z));
// #else
 	    (float)bodies_grav[i].w,
 	    (float)sqrt(bodies_grav[i].x*bodies_grav[i].x + 
 			bodies_grav[i].y*bodies_grav[i].y + 
 			bodies_grav[i].z*bodies_grav[i].z));
// #endif
    float4 pos  = bodies_pos[i];
    float4 acc = bodies_grav[i];
    
    tot_pot += pos.w * acc.w;
      
    tot_force.x += pos.w * acc.x;
    tot_force.y += pos.w * acc.y;
    tot_force.z += pos.w * acc.z;
    
    tot_torque.x += pos.w * (acc.y*pos.z - acc.z*pos.y);
    tot_torque.y += pos.w * (acc.z*pos.x - acc.x*pos.z);
    tot_torque.z += pos.w * (acc.x*pos.y - acc.y*pos.x);
  }
  fprintf(stderr, "tot_force=  [ %lg %lg %lg ]\n",
	  tot_force.x, tot_force.y, tot_force.z);
  fprintf(stderr, "tot_torque= [ %lg %lg %lg ]\n",
	  tot_torque.x, tot_torque.y, tot_torque.z);
  fprintf(stderr, "tot_pot= %lg\n", tot_pot);

  fprintf(stderr, "end-of-program\n");
  
  fprintf(stderr, "end-of-program\n");
  return 0;
}
