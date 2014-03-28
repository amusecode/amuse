#include "Twobody.h"

#define SIGN(a) ((a) < 0 ? -1 : 1)
#define alpha_factor 0.0

Twobody::Twobody() {
  // tolerance = sqrt(pow(2.0, -53)); //   1e-6;
  tolerance = 1e-12;
}
Twobody::Twobody(double tolerance) {
  this->tolerance = tolerance;
}

void Twobody::set_tolerance(double tolerance) {
  this->tolerance = tolerance;
}
double Twobody::get_tolerance() {
  return tolerance;
}
////////////////////////////////////////////////////////////////////////////
// Laguerre's method Solver
////////////////////////////////////////////////////////////////////////////
bool Twobody::solve_by_leapfrog(double mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double dt) {
  double dr2 = x*x + y*y + z*z;
  double dr = sqrt(dr2);
  double dr3 = dr2*dr;
  double ax = -mu/dr3*x;
  double ay = -mu/dr3*y;
  double az = -mu/dr3*z;

  x += vx*dt + ax*dt*dt/2;
  y += vy*dt + ay*dt*dt/2;
  z += vz*dt + az*dt*dt/2;

  vx += ax*dt/2;
  vy += ay*dt/2;
  vz += az*dt/2;

  dr2 = x*x + y*y + z*z;
  dr = sqrt(dr2);
  dr3 = dr2*dr;

  ax = -mu/dr3*x;
  ay = -mu/dr3*y;
  az = -mu/dr3*z;  

  vx += ax*dt/2;
  vy += ay*dt/2;
  vz += az*dt/2;  
}

bool Twobody::solve(double mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double dt) {
  double Cm, Cr, Cv;
  normalize(mu, x, y, z, vx, vy, vz, dt, Cm, Cr, Cv);

  double mymu = mu;
  double myx = x;
  double myy = y;
  double myz = z;
  double myvx = vx;
  double myvy = vy;
  double myvz = vz;
  double mydt = dt;

  bool converged = try_solve(mymu, myx, myy, myz, myvx, myvy, myvz, mydt);

  if(!converged) {
    int counter = 2;
    while(!converged) {
      double dt_trial = dt/counter;
      mymu = mu;
      myx = x;
      myy = y;
      myz = z;
      myvx = vx;
      myvy = vy;
      myvz = vz;
      for(int i=0; i<counter; i++) {
        converged = try_solve(mymu, myx, myy, myz, myvx, myvy, myvz, dt_trial);
	if(!converged) break;
      }
      counter += 2;
    }
  }

  mu = mymu;
  x = myx;
  y = myy;
  z = myz;
  vx = myvx;
  vy = myvy;
  vz = myvz;

  //denormalize(mu, x, y, z, vx, vy, vz, Cm, Cr, Cv);

  return converged;
}
bool Twobody::try_solve(double mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double dt) {
  double M1 = mu;
  double sqmu = sqrt(mu);
  double r0 = sqrt(x*x+y*y+z*z); // current radius
  double v2 = (vx*vx+vy*vy+vz*vz);  // current velocity
  double r0dotv0 = (x*vx + y*vy + z*vz);
  double alpha = (2.0/r0 - v2/M1);  // inverse of semi-major eqn 2.134 MD
// here alpha=1/a and can be negative
  double v0r = r0dotv0/r0;

  double x_p;

  //bool converged = calc_universal_anomaly_newton(sqmu, r0, v0r, alpha, dt, x_p); // solve universal kepler eqn
  bool converged = calc_universal_anomaly_laguerre(sqmu, r0, v0r, alpha, dt, x_p); // solve universal kepler eqn
  //bool converged = calc_universal_anomaly_halley(sqmu, r0, v0r, alpha, dt, x_p); // solve universal kepler eqn
  //bool converged = calc_universal_anomaly_chebyshev(sqmu, r0, v0r, alpha, dt, x_p); // solve universal kepler eqn

  if(converged) {
    double smu = sqrt(M1);
    double foo = 1.0 - r0*alpha;
    double sig0 = r0dotv0/smu;
    double x2,x3,alx2,Cp,Sp,r;
    x2 = x_p*x_p;
    x3 = x2*x_p;
    alx2 = alpha*x2;
    Cp = calc_C(alx2);
    Sp = calc_S(alx2);
    r = sig0*x_p*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42  PC
  // f,g functions equation 2.38a  PC
    double f_p= 1.0 - (x2/r0)*Cp;
    double g_p= dt - (x3/smu)*Sp;
  // dfdt,dgdt function equation 2.38b PC
    double dfdt = x_p*smu/(r*r0)*(alx2*Sp - 1.0);
    double dgdt = 1.0 - (x2/r)*Cp;
    double xn, yn, zn, vxn, vyn, vzn;
    if (r0 > 0.0){ // error catch if a particle is at Sun
      xn = x*f_p + g_p*vx; // eqn 2.65 M+D
      yn = y*f_p + g_p*vy; 
      zn = z*f_p + g_p*vz;
      vxn = dfdt*x + dgdt*vx; //eqn 2.70 M+D
      vyn = dfdt*y + dgdt*vy;
      vzn = dfdt*z + dgdt*vz;
    }
    else {
      xn = x; 
      yn = y;
      zn = z;
      vxn = vx; 
      vyn = vy; 
      vzn = vz;
    }
    x = xn;
    y = yn;
    z = zn;
    vx = vxn;
    vy = vyn;
    vz = vzn;
  }

  return converged;
}

bool Twobody::calc_universal_anomaly_laguerre(double &smu, double &r0, double &v0r, double &alpha, double dt, double &x) {
  int N_LAG = 5.0;
  double M1 = smu*smu;
  double r0dotv0 = v0r*r0;
  double foo = 1.0 - r0*alpha;
  double sig0 = r0dotv0/smu;
  x = M1*dt*dt/r0; // initial guess could be improved 
  //double x = sqmu*fabs(alpha)*dt;
  bool converged = false;	
  int counter = 0;
  double u = 1.0;
  while(!converged) {
    double x2,x3,alx2,Cp,Sp,F,dF,ddF,z,u1,u2;
    x2 = x*x;
    x3 = x2*x;
    alx2 = alpha*x2;
    if(fabs(alx2) > 1e3) {
      return false;
    }
    Cp = calc_C(alx2);
    Sp = calc_S(alx2);
    F = sig0*x2*Cp + foo*x3*Sp + r0*x - smu*dt; // eqn 2.41 PC
    dF = sig0*x*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42 PC
    ddF = sig0*(1.0-alx2*Cp) + foo*x*(1.0 - alx2*Sp);
    z = fabs((N_LAG - 1.0)*((N_LAG - 1.0)*dF*dF - N_LAG*F*ddF));
    z = sqrt(z);
    u1 = N_LAG*F;
    u2 = dF + SIGN(dF)*z + tolerance;
    if(u2 == 0) {
      return false;
    }
    u = u1/u2; 
    x -= u;
    converged = true;
    if(fabs(u) > tolerance) converged = false;
    counter++;
    if(counter > 1e2) {
      return false;
    }
  } 
/*
  double x2,x3,alx2,Cp,Sp,F,dF,ddF,z,u1,u2;
  x2 = x*x;
  x3 = x2*x;
  alx2 = alpha*x2;
  if(fabs(alx2) > 1e3) {
    return false;
  }
  Cp = calc_C(alx2);
  Sp = calc_S(alx2);
  F = sig0*x2*Cp + foo*x3*Sp + r0*x - smu*dt; // eqn 2.41 PC
  dF = sig0*x*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42 PC
  ddF = sig0*(1.0-alx2*Cp) + foo*x*(1.0 - alx2*Sp);
  z = fabs((N_LAG - 1.0)*((N_LAG - 1.0)*dF*dF - N_LAG*F*ddF));
  z = sqrt(z);
  u1 = N_LAG*F;
  u2 = dF + SIGN(dF)*z + tolerance;
  if(u2 == 0) {
    return false;
  }
  u = u1/u2; 
  x -= u;
  converged = true;
  if(fabs(u) > tolerance) converged = false;
*/
  return converged;
}
bool Twobody::calc_universal_anomaly_newton(double &smu, double &r0, double &v0r, double &alpha, double dt, double &x) {
  double M1 = smu*smu;
  double r0dotv0 = v0r*r0;
  double foo = 1.0 - r0*alpha;
  double sig0 = r0dotv0/smu;
  x = M1*dt*dt/r0; // initial guess could be improved 
  //double x = sqmu*fabs(alpha)*dt;
  bool converged = false;	
  int counter = 0;
  double u = 1.0;
  while(!converged) {
    double x2,x3,alx2,Cp,Sp,F,dF,ddF,z,u1,u2;
    x2 = x*x;
    x3 = x2*x;
    alx2 = alpha*x2;
    Cp = calc_C(alx2);
    Sp = calc_S(alx2);

    F = sig0*x2*Cp + foo*x3*Sp + r0*x - smu*dt; // eqn 2.41 PC
    dF = sig0*x*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42 PC

    double my_factor;
    if(F*dF >= 0) my_factor = alpha_factor;
    else my_factor = -alpha_factor;

    u1 = F;
    u2 = dF + my_factor*F; 
    if(u2 == 0) {
      return false;
    }

    u = u1/u2 ;
    if( !is_valid_number(u) ) {
      return false;
    }

    x -= u;

    converged = true;

    if(fabs(u) > tolerance) converged = false;

    counter++;
    if(counter > 1e2) {
      return false;
    }
  } 
  return converged;
}
bool Twobody::calc_universal_anomaly_halley(double &smu, double &r0, double &v0r, double &alpha, double dt, double &x) {
  double M1 = smu*smu;
  double r0dotv0 = v0r*r0;
  double foo = 1.0 - r0*alpha;
  double sig0 = r0dotv0/smu;
  x = M1*dt*dt/r0; // initial guess could be improved 
  //double x = sqmu*fabs(alpha)*dt;
  bool converged = false;	
  int counter = 0;
  double u = 1.0;
  while(!converged) {
    double x2,x3,alx2,Cp,Sp,F,dF,ddF,z1,z2;
    x2 = x*x;
    x3 = x2*x;
    alx2 = alpha*x2;
    Cp = calc_C(alx2);
    Sp = calc_S(alx2);

    F = sig0*x2*Cp + foo*x3*Sp + r0*x - smu*dt; // eqn 2.41 PC
    dF = sig0*x*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42 PC
    ddF = sig0*(1.0-alx2*Cp) + foo*x*(1.0 - alx2*Sp);

    double myfactor1 = alpha_factor;
    double myfactor2 = -alpha_factor;
    
    double dF1 = dF + myfactor1*F;
    double dF2 = dF + myfactor2*F;

    double z2a = 2.0*dF1*dF1 - F*(ddF + 2.0*myfactor1*dF1);
    double z2b = 2.0*dF2*dF2 - F*(ddF + 2.0*myfactor2*dF2);

    if(fabs(z2a) > fabs(z2b)) {
      if(z2a == 0) return false;
      z1 = 2.0*F*dF1;
      u = z1/z2a;
    }
    else {
      if(z2b == 0) return false;
      z1 = 2.0*F*dF2;
      u = z1/z2b;
    }
    if( !is_valid_number(u) ) {
      return false;
    }

    x -= u;

    converged = true;

    if(fabs(u) > tolerance) converged = false;

    counter++;
    if(counter > 1e2) {
      return false;
    }
  } 
  return converged;
}
bool Twobody::calc_universal_anomaly_chebyshev(double &smu, double &r0, double &v0r, double &alpha, double dt, double &x) {
  double M1 = smu*smu;
  double r0dotv0 = v0r*r0;
  double foo = 1.0 - r0*alpha;
  double sig0 = r0dotv0/smu;
  x = M1*dt*dt/r0; // initial guess could be improved 
  //double x = sqmu*fabs(alpha)*dt;
  bool converged = false;	
  int counter = 0;
  double u = 1.0;
  while(!converged) {
    double x2,x3,alx2,Cp,Sp,F,dF,ddF,z;
    x2 = x*x;
    x3 = x2*x;
    alx2 = alpha*x2;
    Cp = calc_C(alx2);
    Sp = calc_S(alx2);

    F = sig0*x2*Cp + foo*x3*Sp + r0*x - smu*dt; // eqn 2.41 PC
    dF = sig0*x*(1.0 - alx2*Sp)  + foo*x2*Cp + r0; // eqn 2.42 PC
    ddF = sig0*(1.0-alx2*Cp) + foo*x*(1.0 - alx2*Sp);

    double z1 = dF + alpha_factor*F;
    double z2 = dF - alpha_factor*F;
    if(fabs(z1) > fabs(z2)) z = z1;
    else z = z2;

    double u1 = F / z;
    double u2 = 0.5*F*F*(ddF + 2.0*alpha_factor*dF) / pow(z, 3);
    u = u1 + u2;

    x -= u;

    converged = true;

    if(fabs(u) > tolerance) converged = false;

    else if(fabs(u) == 0) {
      return false;
    }

    else if( !is_valid_number(u) ) {
      return false;
    }

    counter++;
    if(counter > 1e2) {
      return false;
    }
  } 
  return converged;
}

double Twobody::calc_C(double &z) {
  if(fabs(z)<1e-4) {
    return 1.0/2.0*(1.0 - z/12.0*(1.0 - z/30.0*(1.0 - z/56.0)));
  }
  double u = sqrt(fabs(z));
  if(z>0.0) return (1.0- cos(u))/ z;
  else      return (cosh(u)-1.0)/-z;
}
double Twobody::calc_S(double &z) {
  if (fabs(z)<1e-4) {
    return 1.0/6.0*(1.0 - z/20.0*(1.0 - z/42.0*(1.0 - z/72.0)));
  }
  double u = sqrt(fabs(z));
  double u3 = u*u*u;
  if(z>0.0) return (u -  sin(u))/u3;
  else      return (sinh(u) - u)/u3;
}
////////////////////////////////////////////////////////////////////////////
// Normalizer
////////////////////////////////////////////////////////////////////////////
void Twobody::normalize(double &mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double &dt, double &Cm, double &Cr, double &Cv) {
  //Cm = 1.0/mu;
  Cr = 1.0/sqrt(x*x + y*y + z*z);
  double v2 = vx*vx + vy*vy + vz*vz;
  double e = v2/2 - mu*Cr;
  //Cv = sqrt(Cm/Cr);

  if(e < -tolerance) {
    double a_semi_major = -mu/(2.0*e);
    double P = 2.0*acos(-1.0)*a_semi_major*sqrt(a_semi_major/mu);
    int numP = dt/P;
    dt = dt - numP*P;
  }

  /*mu *= Cm;
  x *= Cr;
  y *= Cr;
  z *= Cr;
  vx *= Cv;
  vy *= Cv;
  vz *= Cv;
  dt *= Cr/Cv;*/
}
void Twobody::denormalize(double &mu, double &x, double &y, double &z, double &vx, double &vy, double &vz, double &Cm, double &Cr, double &Cv) {
  /*mu /= Cm;
  x /= Cr;
  y /= Cr;
  z /= Cr;
  vx /= Cv;
  vy /= Cv;
  vz /= Cv;*/
}
bool Twobody::is_valid_number(double &x) {
  bool valid = true;
  if(isinf(x)) {
    valid = false;
  }
  else if(isnan(x)) {
    valid = false;
  }
  return valid;
}




