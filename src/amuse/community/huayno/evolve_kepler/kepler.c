/* 
Kepler integrator, two-body problem

Marcell Marosvolgyi 2010

code based on:
Fundamentals of Celestial Mechanics, J.M.A. Danby 2nd Edition
*/

#include <stdio.h>
#include <math.h>
//#include <gmp.h>
//#include <mpfr.h>

#define PR 10

DOUBLE x, y, z;
DOUBLE vx, vy, vz;
//DOUBLE time = 0.0;
DOUBLE mu;

int factorial (int k) {
  int r = 1;
  int i;

  for (i=1; i<=k; i++) {
    r *=i;
  }
  return r;
}

DOUBLE norm(DOUBLE x) {
  if (x<0) return -x;
  return x;
}

DOUBLE sign(DOUBLE x) {
  if (norm(x)<1e-12) return 0.0;
  if (x<0) return -1.0;
  else return 1.0;
}

void stumpff(DOUBLE s, DOUBLE *c0, DOUBLE *c1, DOUBLE *c2, DOUBLE *c3) {
  DOUBLE sqrt_s;

  if (s>0.0) {
    sqrt_s = pow(s, 0.5);
    *c0 = cos(sqrt_s);
    *c1 = sin(sqrt_s)/sqrt_s;
  }
  else if (s<0.0){
    sqrt_s = pow(-s, 0.5);
    *c0 = cosh(sqrt_s);
    *c1 = sinh(sqrt_s)/sqrt_s;
  }
  else printf("Error in stumpff s = 0\n");
  *c2 = 1.0/s * (1.0 - *c0);
  *c3 = 1.0/s * (1.0 - *c1);
}

DOUBLE initial_guess_for_s (DOUBLE dt,
			    DOUBLE r0, DOUBLE u, DOUBLE alpha) {
  DOUBLE A, En, Ec, Es, E, X, Y, Dm, sigma;
  DOUBLE s;

  if (norm(dt/r0) <=.2) {
    s = dt/r0 - (dt*dt*u)/(2.0*r0*r0*r0);
  }
  else if (alpha>0) {
    //elliptic motion initial guess
    A = mu/alpha;
    En = pow(mu/A/A/A, 0.5);
    Ec = 1.0 - r0/A;
    Es = u/En/A/A;
    E = pow(Ec*Ec + Es*Es, 0.5);
    dt = dt - floor(En*dt/2.0/3.14159265358979) * (2.0*3.14159265358979)/En;
    Y = En*dt-Es;
    sigma = sign(Es*cos(Y) + Ec*sin(Y));
    X = Y + sigma * 0.85*E;
    s = X/pow(alpha, 0.5);
  }
  else {
    //hyperbolic motion
    A = mu/alpha;
    En = pow(-mu/A/A/A, 0.5);
    Ec = 1.0 - r0/A;
    Es = u/pow(-A*mu, 0.5);
    E = pow(Ec*Ec - Es*Es, 0.5);
    Dm = En*dt;
    if (Dm<0) s = -log((-2.0 * Dm + 1.8 * E)/(Ec - Es))/pow(-alpha, 0.5);
    else s = log((2.0 * Dm + 1.8 * E)/(Ec + Es))/pow(-alpha, 0.5);
  }
  return s;
}

int kepler_solve (DOUBLE dt,
		  DOUBLE *F, DOUBLE *G, DOUBLE *Fdot, DOUBLE *Gdot) {
  DOUBLE f, fp;
  DOUBLE ds, s;
  DOUBLE r0, v0s, u, alpha;
  DOUBLE c0, c1, c2, c3;
  int i=0;

  r0 = pow(x*x +  y*y + z*z, 0.5);
  v0s = vx * vx + vy*vy + vz * vz;
  u = x*vx + y*vy + z*vz;
  alpha = 2.0*mu/r0 - v0s;

  s = initial_guess_for_s (dt, r0, u, alpha);
  ds = 1.0;
  while (norm(ds) > 1.0e-12) {
    stumpff(s*s*alpha, &c0, &c1, &c2, &c3);
    c1 *= s; c2 *= s*s; c3 *= s*s*s;
    f   = r0 * c1 + u * c2 + mu * c3 - dt;
    fp  = r0 * c0 + u * c1 + mu * c2;
    ds = f/fp;
    s -= ds;
    if (i++>50) {
      printf("Convergence error in Newton method\n");
      return -1; 
    }
  }
  *F = 1.0 - (mu/r0) * c2;
  *G = dt - mu * c3;
  *Fdot = - (mu/fp/r0) * c1;
  *Gdot = 1.0 - (mu/fp) * c2;
  return 0;
}

int kepler_evolve (DOUBLE dt) {
  DOUBLE x_new, y_new, z_new;
  DOUBLE vx_new, vy_new, vz_new;
  DOUBLE F, G, Fdot, Gdot;
  //DOUBLE dt;

  //dt = time_new - time;

  if (kepler_solve(dt, &F, &G, &Fdot, &Gdot) == 0) {  
    x_new = x * F + vx * G;
    y_new = y * F + vy * G;
    z_new = z * F + vz * G;
    vx_new = x * Fdot + vx * Gdot;
    vy_new = y * Fdot + vy * Gdot;
    vz_new = z * Fdot + vz * Gdot;
    x = x_new; y = y_new; z = z_new; vx = vx_new; vy = vy_new; vz = vz_new;
    //time = time_new;
    return 0;
  }
  else return -1;
}

int set_mu(DOUBLE mu_) {
  mu = mu_;
  return 0;
}

int kepler_set_position(DOUBLE r[3]) {
  x = r[0]; y = r[1]; z = r[2];
  return 0;
}

int kepler_get_position(DOUBLE r[3]) {
  r[0] = x; r[1] = y; r[2] = z;
  return 0;
}

int kepler_set_velocity(DOUBLE v[3]) {
  vx = v[0]; vy = v[1]; vz = v[2];
  return 0;
}

int kepler_get_velocity(DOUBLE v[3]) {
  v[0] = vx; v[1] = vy; v[2] = vz;
  return 0;
}
