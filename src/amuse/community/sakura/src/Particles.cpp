#include "Particles.h"

Particles::Particles() {
  t = 0;
  N = 0;
  particle.clear();
}
Particles::Particles(vector<Particle> particle) {
  t = 0;
  N = particle.size();
  this->particle = particle;
}
Particles::Particles(double t, vector<Particle> particle) {
  this->t = t;
  N = particle.size();
  this->particle = particle;
}
Particles::Particles(double t, int N, vector<Particle> particle) {
  this->t = t;
  this->N = N;
  this->particle = particle;
}

void Particles::set_t(double t) {
  this->t = t;
}
void Particles::set_N(int N) {
  this->N = N;
}
void Particles::set_particles(vector<Particle> particle) {
  this->particle = particle;
}
void Particles::set_particle(int index, Particle particle) {
  this->particle[index] = particle;
}

double Particles::get_t() {
  return t;
}
int Particles::get_N() {
  return N;
}
vector<Particle> Particles::get_particles() {
  return particle;
}
Particle Particles::get_particle(int index) {
  return particle[index];
}
Particle* Particles::get_pointer_to_star(int index) {
  return &particle[index];
}

void Particles::set_data(vector<double> &data) {
  N = data.size()/7;
  particle.resize(N);
  for(int i=0; i<N; i++) {
    particle[i].set_mass(data[i*7]);
    particle[i].set_x(data[i*7+1]);
    particle[i].set_y(data[i*7+2]);
    particle[i].set_z(data[i*7+3]);
    particle[i].set_vx(data[i*7+4]);
    particle[i].set_vy(data[i*7+5]);
    particle[i].set_vz(data[i*7+6]);
  }
}  
vector<double> Particles::get_data() {
  vector<double> data;
  for(int i=0; i<N; i++) {
    data.push_back(particle[i].get_mass());
    data.push_back(particle[i].get_x());
    data.push_back(particle[i].get_y());
    data.push_back(particle[i].get_z());
    data.push_back(particle[i].get_vx());
    data.push_back(particle[i].get_vy());
    data.push_back(particle[i].get_vz());
  }  
  return data;
}

void Particles::add_to_t(double dt) {
  t += dt;
}
void Particles::add_particle(Particle particle) {
  this->particle.push_back(particle);
  N++;
}
void Particles::remove_particle(int index) {
  particle.erase(particle.begin()+index);
  N--;
}

void Particles::print() {
  cout << t << " " << N << endl;
  for(int i=0; i<N; i++) {
    cout << particle[i].get_mass() << " "
         << particle[i].get_x() << " " << particle[i].get_y() << " " << particle[i].get_z() << " "
         << particle[i].get_vx() << " " << particle[i].get_vy() << " " << particle[i].get_vz() << endl;
  }
  cout << endl;
}
void Particles::print(double t_cpu) {
  cout << t << " " << N << " " << t_cpu << endl;
  for(int i=0; i<N; i++) {
    cout << particle[i].get_mass() << " "
         << particle[i].get_x() << " " << particle[i].get_y() << " " << particle[i].get_z() << " "
         << particle[i].get_vx() << " " << particle[i].get_vy() << " " << particle[i].get_vz() << endl;
  }
  cout << endl;
}
void Particles::print(ofstream &odata) {
  odata << t << " " << N << endl;
  for(int i=0; i<N; i++) {
    odata << particle[i].get_mass() << " "
         << particle[i].get_x() << " " << particle[i].get_y() << " " << particle[i].get_z() << " "
         << particle[i].get_vx() << " " << particle[i].get_vy() << " " << particle[i].get_vz() << endl;
  }
  odata << endl;
}
void Particles::print(double t_cpu, ofstream &odata) {
  odata << t << " " << N << " " << t_cpu << endl;
  for(int i=0; i<N; i++) {
    odata << particle[i].get_mass() << " " 
         << particle[i].get_x() << " " << particle[i].get_y() << " " << particle[i].get_z() << " "
         << particle[i].get_vx() << " " << particle[i].get_vy() << " " << particle[i].get_vz() << endl;
  }
  odata << endl;
}



