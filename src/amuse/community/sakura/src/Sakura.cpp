#include "Sakura.h"

Sakura::Sakura() {
  timestep.set_dt(1e-3);
  Nstep = 0;
  sum_dt = 0;
  sum_inv_dt = 0;
}
Sakura::Sakura(Particles particles) {
  this->particles = particles;
  timestep.set_dt(1e-3);
  Nstep = 0;
  sum_dt = 0;
  sum_inv_dt = 0;
}
Sakura::Sakura(Particles particles, double dt) {
  this->particles = particles;
  timestep.set_dt(dt);
  Nstep = 0;
  sum_dt = 0;
  sum_inv_dt = 0;
}

void Sakura::set_particles(Particles particles) {
  this->particles = particles;

  double t_begin = particles.get_t();
  int N = particles.get_N();
  vector<Particle> particle = particles.get_particles();

}
void Sakura::set_dt(double dt) {
  timestep.set_dt(dt);
}
void Sakura::set_tolerance(double tolerance) {
  twobody.set_tolerance(tolerance);
}

Particles Sakura::get_particles() {
  return this->particles;
}
Particles* Sakura::get_pointer_to_particles() {
  return &particles;
}
double Sakura::get_dt() {
  return timestep.get_dt();
}
double Sakura::get_tolerance() {
  return twobody.get_tolerance();
}
unsigned int Sakura::get_Nstep() {
  return Nstep;
}
double Sakura::get_sum_dt() {
  return sum_dt;
}
double Sakura::get_sum_inv_dt() {
  return sum_inv_dt;
}
Particle* Sakura::get_pointer_to_star(int index) {
  return particles.get_pointer_to_star(index);
}

void Sakura::evolve(double t, Communicator &communicator) {
  double t_begin = particles.get_t();
  int N = particles.get_N();
  vector<Particle> particle = particles.get_particles();

  double t_end = t;
  t = t_begin;	

  // Evolve system	
  while(t < t_end) {

    // Time step size
    double dt = timestep.get_dt(particle);

    //sum_dt += dt;
    //sum_inv_dt += 1.0/dt;
    //Nstep++;

    if(t+dt >= t_end) dt = t_end-t;

    step(particle, dt, communicator);	  
	  
    // Update time
    t += dt;	  
  }

  particles.set_t(t);
  particles.set_particles(particle); 
}
void Sakura::step(vector<Particle> &particle, double dt, Communicator &communicator) {
  int N = particle.size();  

  double Mtot = 0;
  double vxcm = 0;
  double vycm = 0;
  double vzcm = 0;
  for(int i=0; i<N; i++) {
    Mtot += particle[i].get_mass();
    vxcm += particle[i].get_mass() * particle[i].get_vx();
    vycm += particle[i].get_mass() * particle[i].get_vy();
    vzcm += particle[i].get_mass() * particle[i].get_vz();
  }
  vxcm /= Mtot;
  vycm /= Mtot;
  vzcm /= Mtot;

  vector<double> mydp(6*N, 0);
  for(int i=communicator.get_my_begin(); i<communicator.get_my_end(); i++) {
    for(int j=0; j<N; j++) {
      if(i != j) {

	double M = particle[i].get_mass() + particle[j].get_mass();
	double dxn = particle[j].get_x() - particle[i].get_x();
	double dyn = particle[j].get_y() - particle[i].get_y();
	double dzn = particle[j].get_z() - particle[i].get_z();
	double dvxn = particle[j].get_vx() - particle[i].get_vx();
	double dvyn = particle[j].get_vy() - particle[i].get_vy();
	double dvzn = particle[j].get_vz() - particle[i].get_vz();

        //double dr2 = dxn*dxn + dyn*dyn + dzn*dzn;
        //double dv2 = dvxn*dvxn + dvyn*dvyn + dvzn*dvzn;
        //if(dr2 > 10.0*dt*dt*dv2) twobody.solve_by_leapfrog(M, dxn, dyn, dzn, dvxn, dvyn, dvzn, dt);
	//if(dr2 > 0.125*0.125) twobody.solve_by_leapfrog(M, dxn, dyn, dzn, dvxn, dvyn, dvzn, dt);
	//else twobody.solve(M, dxn, dyn, dzn, dvxn, dvyn, dvzn, dt);
        twobody.solve(M, dxn, dyn, dzn, dvxn, dvyn, dvzn, dt);

	double mu = particle[i].get_mass()*particle[j].get_mass() / (particle[i].get_mass()+particle[j].get_mass());

	double dx0 = particle[j].get_x() - particle[i].get_x();
	double dy0 = particle[j].get_y() - particle[i].get_y();
	double dz0 = particle[j].get_z() - particle[i].get_z();
	double dvx0 = particle[j].get_vx() - particle[i].get_vx();
	double dvy0 = particle[j].get_vy() - particle[i].get_vy();
	double dvz0 = particle[j].get_vz() - particle[i].get_vz();

	mydp[i*6+0] -= mu * ( (dxn-dx0) - dvx0*dt );
	mydp[i*6+1] -= mu * ( (dyn-dy0) - dvy0*dt );
	mydp[i*6+2] -= mu * ( (dzn-dz0) - dvz0*dt );
	mydp[i*6+3] -= mu * ( dvxn - dvx0 );
	mydp[i*6+4] -= mu * ( dvyn - dvy0 );
	mydp[i*6+5] -= mu * ( dvzn - dvz0 ); 
      }
    }
  }

  vector<double> dp(6*N, 0);
  communicator.join(mydp, dp);
  communicator.bcast(dp);

  for(int i=0; i<N; i++) {
    dp[i*6+0] /= particle[i].get_mass();
    dp[i*6+1] /= particle[i].get_mass();
    dp[i*6+2] /= particle[i].get_mass();
    dp[i*6+3] /= particle[i].get_mass();
    dp[i*6+4] /= particle[i].get_mass();
    dp[i*6+5] /= particle[i].get_mass();
  }

  for(int i=0; i<N; i++) {
    particle[i].add_to_x(dp[i*6+0] + (particle[i].get_vx())*dt);
    particle[i].add_to_y(dp[i*6+1] + (particle[i].get_vy())*dt);
    particle[i].add_to_z(dp[i*6+2] + (particle[i].get_vz())*dt);

    particle[i].add_to_vx(dp[i*6+3]);
    particle[i].add_to_vy(dp[i*6+4]);
    particle[i].add_to_vz(dp[i*6+5]);
  }
}

void Sakura::initial_calculations() {
  ;
}
double Sakura::get_t() {
  return particles.get_t();
}
void Sakura::set_t(double t) {
  particles.set_t(t);
}
void Sakura::step(double dt, Communicator &communicator) {
  vector<Particle> particle = particles.get_particles();
  step(particle, dt, communicator);
  particles.set_particles(particle);
}
vector<double> Sakura::get_coordinates(vector<Particle> &particle) {
  int N = particle.size();
  vector<double> coordinates(N*7);
  int counter = 0;
  for(int i=0; i<N; i++) {
    coordinates[counter+0] = particle[i].get_mass();
    coordinates[counter+1] = particle[i].get_x();
    coordinates[counter+2] = particle[i].get_y();
    coordinates[counter+3] = particle[i].get_z();
    coordinates[counter+4] = particle[i].get_vx();
    coordinates[counter+5] = particle[i].get_vy();
    coordinates[counter+6] = particle[i].get_vz();
    counter += 7;
  }
  return coordinates;
}
void Sakura::update_particles(vector<Particle> &particle, vector<double> &coordinates) {
  int N = particle.size();
  int counter = 0;
  for(int i=0; i<N; i++) {
    particle[i].set_mass( coordinates[counter+0] );
    particle[i].set_x( coordinates[counter+1] );
    particle[i].set_y( coordinates[counter+2] );
    particle[i].set_z( coordinates[counter+3] );
    particle[i].set_vx( coordinates[counter+4] );
    particle[i].set_vy( coordinates[counter+5] );
    particle[i].set_vz( coordinates[counter+6] );
    counter += 7;
  }
}
void Sakura::update_particles(vector<double> &coordinates) {
  int N = particles.get_N();
  vector<Particle> p(N);
  int counter = 0;
  for(int i=0; i<N; i++) {
    p[i].set_mass( coordinates[counter+0] );
    p[i].set_x( coordinates[counter+1] );
    p[i].set_y( coordinates[counter+2] );
    p[i].set_z( coordinates[counter+3] );
    p[i].set_vx( coordinates[counter+4] );
    p[i].set_vy( coordinates[counter+5] );
    p[i].set_vz( coordinates[counter+6] );
    counter += 7;
  }
}
vector<double> Sakura::get_data() {
  vector<Particle> p = particles.get_particles();
  int N = p.size();
  vector<double> data(7*N, 0);
  for(int i=0; i<N; i++) {    
    data[i*7+0] = p[i].get_mass();
    data[i*7+1] = p[i].get_x();
    data[i*7+2] = p[i].get_y();
    data[i*7+3] = p[i].get_z();
    data[i*7+4] = p[i].get_vx();
    data[i*7+5] = p[i].get_vy();
    data[i*7+6] = p[i].get_vz();
  }
  return data;
}



