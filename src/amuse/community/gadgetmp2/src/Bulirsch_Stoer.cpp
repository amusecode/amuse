#include "Bulirsch_Stoer.h"

Bulirsch_Stoer::Bulirsch_Stoer() {
  tolerance = "1e-6";
  n_max = 32;
  k_max = 128;
}
Bulirsch_Stoer::Bulirsch_Stoer(mpreal tolerance) {
  this->tolerance = tolerance;
  n_max = 32;
  k_max = 128;
}
Bulirsch_Stoer::Bulirsch_Stoer(mpreal tolerance, int n_max, int k_max) {
  this->tolerance = tolerance;
  this->n_max = n_max;
  this->k_max = k_max;
}

void Bulirsch_Stoer::set_tolerance(mpreal tolerance) {
  this->tolerance = tolerance;
}
void Bulirsch_Stoer::set_n_max(int n_max) {
  this->n_max = n_max;
}
void Bulirsch_Stoer::set_k_max(int k_max) {
  this->k_max = k_max;
}

mpreal Bulirsch_Stoer::get_tolerance() {
  return tolerance;
}
int Bulirsch_Stoer::get_n_max() {
  return n_max;
}
int Bulirsch_Stoer::get_k_max() {
  return k_max;
}

bool Bulirsch_Stoer::integrate(Cluster &cl, mpreal &dt) {
  Cluster c0 = cl;
  Cluster c  = cl;	
  mpreal timestep = dt;

  bool flag = step(c, timestep);

  if(flag == false) {
    int k = 2;
    while( flag == false && k <= k_max ) {
      timestep = dt/k;
      c = c0;
      flag = step(c, timestep);
      k += 2;
    }
  }  

  if(flag == true) {
    cl = c;
    dt = timestep;
  }

  return flag;
}
bool Bulirsch_Stoer::step(Cluster &cl, mpreal &dt) {
  bool flag;
  int n;
  vector<mpreal> h;
  vector<Cluster> c;
  Cluster cl_exp0 = cl;
  Cluster cl_exp  = cl;

  // n=1
  n=1;
  h.push_back( dt/n );
  c.push_back( cl );
  c[0].step( h[0] );
  cl_exp0 = c[0];

  // n=2
  n=2;
  h.push_back( dt/n );
  c.push_back( cl );
  for(int i=0; i<n; i++) c[1].step( h[1] );
  extrapol( cl_exp, h, c );
  
  flag = error_control(cl_exp0, cl_exp);

  if(flag == false) {
    while(flag == false && n <= n_max) {
      n += 2;
      h.push_back( dt/n );
      c.push_back( cl );
      for(int i=0; i<n; i++) c[n/2].step( h[n/2] );
      cl_exp0 = cl_exp;
      extrapol( cl_exp, h, c );   

      flag = error_control(cl_exp0, cl_exp);   
    }    
  }

  cl = cl_exp;

  return flag;
}
void Bulirsch_Stoer::extrapol(Cluster &cl_exp, vector<mpreal> &dt, vector<Cluster> &c) {
  int M = dt.size();
  int N = c[0].s.size();

  vector<mpreal> x_sample(M), y_sample(M), z_sample(M), vx_sample(M), vy_sample(M), vz_sample(M);
  for(int i=0; i<N; i++) {
    for(int j=0; j<M; j++) {
      x_sample[j]  = c[j].s[i].r[0];
      y_sample[j]  = c[j].s[i].r[1];
      z_sample[j]  = c[j].s[i].r[2];
      vx_sample[j] = c[j].s[i].v[0];
      vy_sample[j] = c[j].s[i].v[1];
      vz_sample[j] = c[j].s[i].v[2];
    }
    cl_exp.s[i].r[0]  = extrapolate(dt, x_sample, "0.0");
    cl_exp.s[i].r[1]  = extrapolate(dt, y_sample, "0.0");
    cl_exp.s[i].r[2]  = extrapolate(dt, z_sample, "0.0");
    cl_exp.s[i].v[0] = extrapolate(dt, vx_sample, "0.0");
    cl_exp.s[i].v[1] = extrapolate(dt, vy_sample, "0.0");
    cl_exp.s[i].v[2] = extrapolate(dt, vz_sample, "0.0");
  }
}
mpreal Bulirsch_Stoer::extrapolate(vector<mpreal> x, vector<mpreal> y, mpreal x0) {
  int N = x.size();
  if(N == 1) {
    return y[0];
  }
  else {
    for(int i=1; i<N; i++) {
      for(int j=0; j<N-i; j++) {
	y[j] = ( (x0-x[j+i])*y[j] + (x[j]-x0)*y[j+1] ) / ( x[j]-x[j+i] );
      }
    }
    return y[0];
  }  
}
bool Bulirsch_Stoer::error_control(Cluster &c1, Cluster &c2) {
  int N = c1.s.size();
  bool flag = true;
  for(int i=0; i<N; i++) {
    if( fabs(c2.s[i].r[0]-c1.s[i].r[0]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].r[1]-c1.s[i].r[1]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].r[2]-c1.s[i].r[2]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].v[0]-c1.s[i].v[0]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].v[1]-c1.s[i].v[1]) > tolerance ) {
      flag = false;
      break;
    }
    if( fabs(c2.s[i].v[2]-c1.s[i].v[2]) > tolerance ) {
      flag = false;
      break;
    }
  }

  return flag;
}


