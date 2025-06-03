#include "Bs_integrator.h"

/////////////////////////////////////////
// Constructors
/////////////////////////////////////////
Bs_integrator::Bs_integrator() {
  epsilon = "1e-6";
  n_max = 64;
  k_max = 64;
}
Bs_integrator::Bs_integrator(mpreal e) {
  epsilon = e;
  n_max = 64;
  k_max = 64;
}
Bs_integrator::Bs_integrator(mpreal e, int n, int k) {
  epsilon = e;
  n_max = n;
  k_max = k;
}
/////////////////////////////////////////
// Set
/////////////////////////////////////////
void Bs_integrator::set_epsilon(mpreal e) {
  epsilon = e;
}
void Bs_integrator::set_n_max(int n) {
  n_max = n;
}
void Bs_integrator::set_k_max(int k) {
  k_max = k;
}
/////////////////////////////////////////
// Get
/////////////////////////////////////////
mpreal Bs_integrator::get_epsilon() {
  return epsilon;
}
int Bs_integrator::get_n_max() {
  return n_max;
}
int Bs_integrator::get_k_max() {
  return k_max;
}
/////////////////////////////////////////
// Functions
/////////////////////////////////////////
void Bs_integrator::integrate(Cluster &cluster, mpreal &dt) {
  Cluster cl = cluster;
  mpreal timestep = dt;

  step(cl, timestep);

  if(flag == 0) {
    int k = 2;
    while( flag == 0 && k <= k_max ) {
      timestep = dt/k;
      cl = cluster;
      step(cl, timestep);
      k += 2;
    }
  }  

  if(flag == 0) {
    flag = 0;
  }
  else {
    flag = 1;

    cluster = cl;
    dt = timestep;
  }
}
void Bs_integrator::step(Cluster &cluster, mpreal dt) {
  flag = 1;

  int n;
  vector<mpreal> h;
  vector<Cluster> cl;
  Cluster cl_exp0 = cluster;
  Cluster cl_exp = cluster;

  // n=1
  n=1;
  h.push_back( dt/n );
  cl.push_back( cluster );
  cl[0].leapfrog( h[0] );
  cl_exp0 = cl[0];

  // n=2
  n=2;
  h.push_back( dt/n );
  cl.push_back( cluster );
  for(int i=0; i<n; i++) cl[1].leapfrog( h[1] );
  extrapol( cl_exp, h, cl );
  
  flag = error_control(cl_exp0, cl_exp);

  if(flag == 0) {
    while(flag == 0 && n <= n_max) {
      n += 2;
      h.push_back( dt/n );
      cl.push_back( cluster );
      for(int i=0; i<n; i++) cl[n/2].leapfrog( h[n/2] );
      cl_exp0 = cl_exp;
      extrapol( cl_exp, h, cl );   

      flag = error_control(cl_exp0, cl_exp);   
    }    
  }

  cluster = cl_exp;
}
void Bs_integrator::extrapol( Cluster &cl_exp, vector<mpreal> &h, vector<Cluster> &cl )
{
  int M = h.size();
  int N = cl[0].get_N();

  vector< Star* > st(M);
  for(int i=0; i<M; i++) st[i] = cl[i].get_pointer_to_star();  
  Star* st_exp = cl_exp.get_pointer_to_star(); 

  for(int i=0; i<N; i++)
  {
    vector<mpreal> x_sample(M), y_sample(M), z_sample(M), vx_sample(M), vy_sample(M), vz_sample(M);

    for(int j=0; j<M; j++)
    {
      x_sample[j] = st[j]->x;
      y_sample[j] = st[j]->y;
      z_sample[j] = st[j]->z;
      vx_sample[j] = st[j]->vx;
      vy_sample[j] = st[j]->vy;
      vz_sample[j] = st[j]->vz;
    }

    st_exp->x = extrapolate(h, x_sample, "0.0");
    st_exp->y = extrapolate(h, y_sample, "0.0");
    st_exp->z = extrapolate(h, z_sample, "0.0");
    st_exp->vx = extrapolate(h, vx_sample, "0.0");
    st_exp->vy = extrapolate(h, vy_sample, "0.0");
    st_exp->vz = extrapolate(h, vz_sample, "0.0");

    for(int j=0; j<M; j++) st[j]++;   
    st_exp++; 
  }
}
mpreal Bs_integrator::extrapolate(vector<mpreal> h, vector<mpreal> A, mpreal x)
{
  int N = h.size();	// # data points
  if(N == 1)		// 0th order fit
  {
    return A[0];
  }
  else
  {
    for(int i=1; i<N; i++)    
    {
      for(int j=0; j<N-i; j++)
      {
	A[j] = ( (x-h[j+i])*A[j] + (h[j]-x)*A[j+1] ) / ( h[j]-h[j+i] );
      }
    }
    return A[0];
  }  
}
int Bs_integrator::error_control(Cluster &cl0, Cluster &cl)
{
  int flag = 1;

  int N = cl.get_N();
  Star* st0 = cl0.get_pointer_to_star(); 
  Star* st = cl.get_pointer_to_star();  

  for(int i=0; i<N; i++)
  {
    if( abs(st0->x-st->x) > epsilon )
    {
      flag = 0;
      break;
    }
    if( abs(st0->y-st->y) > epsilon )
    {
      flag = 0;
      break;
    }
    if( abs(st0->z-st->z) > epsilon )
    {
      flag = 0;
      break;
    }
    if( abs(st0->vx-st->vx) > epsilon )
    {
      flag = 0;
      break;
    }
    if( abs(st0->vy-st->vy) > epsilon )
    {
      flag = 0;
      break;
    }
    if( abs(st0->vz-st->vz) > epsilon )
    {
      flag = 0;
      break;
    }
    st0++;
    st++;
  }
  return flag;
}
int Bs_integrator::converged()
{
  if(flag == 1) return 1;
  else return 0;
}

