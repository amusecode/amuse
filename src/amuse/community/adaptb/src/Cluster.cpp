#include "Cluster.h"

//////////////////////////////////////////////
// Constructor
//////////////////////////////////////////////
Cluster::Cluster() {
  t = "0";
  N = 0;
  star.clear();
}

Cluster::Cluster( string file, Force fo ) {
  star.clear();
  ifstream data;
  data.open( file.c_str() );
  if( !data ) {
    cerr << "Could not open " << file << "!" << endl;
    exit(1);
  }
  else {
    data >> t >> N;
    mpreal m, x, y, z, vx, vy, vz;
    for(int i=0; i<N; i++) {
      data >> m >> x >> y >> z >> vx >> vy >> vz;
      Star st(m, x, y, z, vx, vy, vz);
      star.push_back( st );
    }
  }
  data.close();
  force = fo;
}
//////////////////////////////////////////////
// Set
//////////////////////////////////////////////
void Cluster::add_star(int id, mpreal m, mpreal radius, mpreal x, mpreal y, mpreal z, mpreal vx, mpreal vy, mpreal vz) {
      Star st(id, m, radius, x, y, z, vx, vy, vz);
      star.push_back( st );
}
void Cluster::set_t(mpreal T) {
  t = T;
}
void Cluster::set_N(int N) {
  this->N = N;
}
//////////////////////////////////////////////
// Get
//////////////////////////////////////////////
mpreal Cluster::get_t() {
  return t;
}
int Cluster::get_N() {
  return N;
}

Force* Cluster::get_pointer_to_force() {
  Force *p = &force;
  return p;
}

Star* Cluster::get_pointer_to_star() {
  Star *p = &star[0];
  return p;
}
Star* Cluster::get_pointer_to_star(int index) {
  Star *p = &(star.at(index));
  return p;
}
mpreal Cluster::get_E() {
  mpreal EK = "0";
  mpreal EP = "0";
  for(int i=0; i<N; i++) {
    star[i].calc_v2mag();
    EK += star[i].m*star[i].get_v2mag()/2;
  }
  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {
      mpreal dx = star[j].x - star[i].x;
      mpreal dy = star[j].y - star[i].y;
      mpreal dz = star[j].z - star[i].z;
      mpreal dr2 = dx*dx + dy*dy + dz*dz;
      EP -= star[i].m*star[j].m/sqrt(dr2);
    }
  }
  return EK + EP;
}

mpreal Cluster::get_dt() {
  mpreal dt_min = "1e100";
  for(int i=0; i<N; i++) {
    if(star[i].dt < dt_min) dt_min = star[i].dt;
  }
  dt_min = sqrt(dt_min);
  return sqrt(dt_min);
}

//////////////////////////////////////////////
// Calculate
//////////////////////////////////////////////
void Cluster::calc_a()
{
  for(int i=0; i<N; i++) {
    star[i].reset_a();
  }

  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {
      force.gravity(&star[i], &star[j]);
    }
  }  
}
void Cluster::calc_a_dt()
{
  for(int i=0; i<N; i++) {
    star[i].reset_a();
    star[i].reset_dt();
  }

  for(int i=0; i<N-1; i++) {
    for(int j=i+1; j<N; j++) {
      force.gravity_dt(&star[i], &star[j]);
    }
  }  
}
void Cluster::leapfrog(mpreal dt)
{
  for(int i=0; i<N; i++) dynamics.update_r(&star[i], dt);
  calc_a();
  for(int i=0; i<N; i++) dynamics.update_v(&star[i], dt);
}
//////////////////////////////////////////////
// Printers
//////////////////////////////////////////////
void Cluster::print()
{
  vector<Star>::iterator st;
  cout << t << " " << N << endl;
  for(st=star.begin(); st!=star.end(); ++st)
  {
    cout << st->m << " " << st->x << " " << st->y << " " << st->z << " " << st->vx << " " << st->vy << " " << st->vz << endl;
  }
}
void Cluster::print( ofstream &data )
{
  vector<Star>::iterator st;
  data << t << " " << N << endl;
  for(st = star.begin(); st != star.end(); ++st)
  {
    data << st->m << " " << st->x << " " << st->y << " " << st->z << " " << st->vx << " " << st->vy << " " << st->vz << endl;
  }
}

