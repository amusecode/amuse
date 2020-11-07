////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////
#include <iostream>
using namespace std;

#include <fstream>
#include <cstdlib>
#include <map>
#include <vector>

#include "Brutus.h"

////////////////////////////////////////////////////////
// Declare global variables
////////////////////////////////////////////////////////
ofstream odata;
ofstream DEBUG::DEB;

size_t particle_id_counter = 0;

int numBits = 64;
//int numDigits = numBits/4;

string out_directory;
std::map<int, int> local_index_map;

/*
 * We need this result_strings array to ensure that
 * C++ strings are not reclaimed before the function ends
 */
//string result_strings[10];

Brutus *brutus = NULL;
Cluster *cluster = NULL;

mpreal t_begin; //constructor sets to zero by default = "0";
mpreal eta = "0.24";
mpreal t; //constructor sets to zero by default = "0";

mpreal epsilon = "1e-6"; // Bulirsch-Stoer tolerance

std::string arg0;
std::string arg1;
std::string arg2;
std::string arg3;
std::string arg4;
std::string arg5;
std::string arg6;
std::string arg7;
std::string arg8;
std::string arg9;
std::string result;

//vector<mpreal> data, data_radius;

////////////////////////////////////////////////////////
// Amuse interface functions
////////////////////////////////////////////////////////
int initialize_code() {
    odata.open("temp.log");

    mpreal::set_default_prec(numBits);

    brutus = new Brutus();
    cluster = new Cluster();

    particle_id_counter = 0;

    t_begin.setZero();// = "0";
    t = t_begin;

    DEBUG::DEB.open ("debug.log");
    DEBUG::DEB << "Brutus started" <<"\n";
    DEBUG::DEB.flush();
    return 0;
}

// functions with "_string" assign strings to mpreals, and without "_string" assign doubles to mpreals

int new_particle_string(int *particle_identifier, char* mass,
        char* x, char* y, char* z, char* vx, char* vy, char* vz, char* radius) {
    Star s;
    s.m=mass;
    s.x[0]=x;
    s.x[1]=y;
    s.x[2]=z;
    s.v[0]=vx;
    s.v[1]=vy;
    s.v[2]=vz;
    s.r=radius;
    s.id=particle_id_counter;
    cluster->s.push_back(s);

    *particle_identifier = particle_id_counter;
    particle_id_counter++;

    return 0;
}
int new_particle_float64(int *particle_identifier, double mass,
        double x, double y, double z, double vx, double vy, double vz, double radius) {
    DEBUG::DEB <<  "new_particle_float64  "<< "\n"; DEBUG::DEB.flush();
    Star s;
    s.m=mass;
    s.x[0]=x;
    s.x[1]=y;
    s.x[2]=z;
    s.v[0]=vx;
    s.v[1]=vy;
    s.v[2]=vz;
    s.r=radius;
    s.id=particle_id_counter;
    cluster->s.push_back(s);

    *particle_identifier = particle_id_counter;
    particle_id_counter++;
    return 0;
}

int commit_particles() {

    DEBUG::DEB <<  "commit_particles  "<< "\n"; DEBUG::DEB.flush();
//    brutus->set_data(data);

//    int numStar = data.size()/arg_cnt;
    brutus->setup();
    brutus->cl=*cluster;
    return 0;
}

// Begin time
int set_t_begin_string(char* tb) {
    t_begin = tb;
    return 0;
}
int get_t_begin_string(char **tb) {
    result = t_begin.toString();
    *tb = (char*) result.c_str();
    return 0;
}
int set_t_begin(double tb) {
    t_begin = (mpreal)tb;
    return 0;
}
int get_t_begin(double *tb) {
    *tb = t_begin.toDouble();
    return 0;
}

int set_begin_time(double input) {
    t_begin = (mpreal)input;
    return 0;
}
int get_begin_time(double * output) {
    *output = t_begin.toDouble();
    return 0;
}

// Timestep parameter, eta
int set_eta_string(char* myeta) {
    eta = myeta;
    brutus->set_eta(eta);
    return 0;
}
int get_eta_string(char **myeta) {
    result = eta.toString();
    *myeta = (char*) result.c_str();
    return 0;
}
int set_eta(double myeta) {
    eta = (mpreal)myeta;
    brutus->set_eta(eta);
    return 0;
}
int get_eta(double *myeta) {
    *myeta = eta.toDouble();
    return 0;
}

int get_time_step(double* dtt){
    *dtt = eta.toDouble();
    return 0;
}

// Time
int set_t_string(char* tt) {
    t = tt;
    return 0;
}
int get_t_string(char **tt) {
    result = t.toString();
    *tt = (char*) result.c_str();
    return 0;
}
int set_t(double tt) {
    t = (mpreal)tt;
    return 0;
}
int get_t(double *tt) {
    *tt = t.toDouble();
    return 0;
}

int get_time(double* time){
    *time = t.toDouble();
    return 0;
}

// Bulirsch-Stoer tolerance, epsilon
// internally epsilon is mpreal, for the interface (and in practice) a double bs_tolerance is enough
int set_bs_tolerance_string(char *bs_tolerance) {
    epsilon = bs_tolerance;
    brutus->set_tolerance(epsilon);
    return 0;
}
int get_bs_tolerance_string(char **bs_tolerance) {
    result = epsilon.toString();
    *bs_tolerance = (char*) result.c_str();
    return 0;
}
int set_bs_tolerance(double bs_tolerance) {

odata << t << ": changing e from " << epsilon << ", to " << bs_tolerance << endl;

    epsilon = (mpreal) bs_tolerance;
    brutus->set_tolerance(epsilon);

odata << epsilon << " " << brutus->get_tolerance() << endl;

    return 0;
}
int get_bs_tolerance(double *bs_tolerance) {
    *bs_tolerance = epsilon.toDouble();
    return 0;
}
// Word-length, numBits in mantissa
int set_word_length(int mynumBits) {
odata << t << ": changing L from " << numBits << ", to " << mynumBits << endl;
    DEBUG::DEB <<  "set_word_length  "<< cluster->s.size() <<"\n"; DEBUG::DEB.flush();

    numBits = mynumBits;
    mpreal::set_default_prec(numBits);
//    numDigits = (int)abs(log10( pow("2.0", -numBits) )).toLong();
//    brutus->set_numBits(numBits);

//odata << numBits << " " << brutus->get_numBits() << endl;
    for(int i=0;i<cluster->s.size();i++)
    {
        cluster->s[i].x[0].set_prec(numBits);
        cluster->s[i].x[1].set_prec(numBits);
        cluster->s[i].x[2].set_prec(numBits);
        cluster->s[i].v[0].set_prec(numBits);
        cluster->s[i].v[1].set_prec(numBits);
        cluster->s[i].v[2].set_prec(numBits);
        cluster->s[i].a[0].set_prec(numBits);
        cluster->s[i].a[1].set_prec(numBits);
        cluster->s[i].a[2].set_prec(numBits);
        cluster->s[i].a0[0].set_prec(numBits);
        cluster->s[i].a0[1].set_prec(numBits);
        cluster->s[i].a0[2].set_prec(numBits);
        cluster->s[i].a_step[0].set_prec(numBits);
        cluster->s[i].a_step[1].set_prec(numBits);
        cluster->s[i].a_step[2].set_prec(numBits);
        cluster->s[i].m.set_prec(numBits);
        cluster->s[i].r.set_prec(numBits);
    }
    return 0;
}
int get_word_length(int *mynumBits) {
    *mynumBits = numBits;
    return 0;
}

int set_eps2(double eps2) {
    return 0;
}
int get_eps2(double *eps2) {
    return 0;
}

int commit_parameters() {
    //brutus = new Brutus(t_begin, data, epsilon, numBits);

    //Brutus br(t_begin, data, epsilon, numBits);
    //brutus = &br;
    //brutus = new Brutus(t_begin, data, epsilon, numBits);

    //brutus->set_t_begin(t_begin);
    //brutus->set_tolerance(epsilon);
    //brutus->set_numBits(numBits);
    return 0;
}
int recommit_parameters() {
    return commit_parameters();
}

// Get/set particle properties
int get_mass_string(int id, char **mass) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  arg0 = cluster->s[id].m.toString();
  *mass = (char*) arg0.c_str();
  return 0;
}
int set_mass_string(int id, char *mass) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].m = mass;
  return 0;
}
int get_mass(int id, double* mass) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *mass = cluster->s[id].m.toDouble();
  return 0;
}
int set_mass(int id, double mass) {
    DEBUG::DEB <<  "set_mass  "<< "\n"; DEBUG::DEB.flush();
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].m = (mpreal)mass;
  return 0;
}

int get_position_string(int id, char **x, char **y, char **z) {

   DEBUG::DEB <<  "get_position_string  "<< "\n"; DEBUG::DEB.flush();
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  arg1 = cluster->s[id].x[0].toString();
  arg2 = cluster->s[id].x[1].toString();
  arg3 = cluster->s[id].x[2].toString();
  *x = (char*) arg1.c_str();
  *y = (char*) arg2.c_str();
  *z = (char*) arg3.c_str();
  return 0;
}
int set_position_string(int id, char *x, char *y, char *z) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].x[0] = x;
  cluster->s[id].x[1] = y;
  cluster->s[id].x[2] = z;
  return 0;
}
int get_position(int id, double* x, double* y, double* z) {
   DEBUG::DEB <<  "get_position  "<< "\n"; DEBUG::DEB.flush();
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *x = cluster->s[id].x[0].toDouble();
  *y = cluster->s[id].x[1].toDouble();
  *z = cluster->s[id].x[2].toDouble();
  return 0;
}
int set_position(int id, double x, double y, double z) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].x[0] = (mpreal)x;
  cluster->s[id].x[1] = (mpreal)y;
  cluster->s[id].x[2] = (mpreal)z;
  return 0;
}

int get_velocity_string(int id, char **vx, char **vy, char **vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  arg4 = cluster->s[id].v[0].toString();
  arg5 = cluster->s[id].v[1].toString();
  arg6 = cluster->s[id].v[2].toString();
  *vx = (char*) arg4.c_str();
  *vy = (char*) arg5.c_str();
  *vz = (char*) arg6.c_str();
  return 0;
}
int set_velocity_string(int id, char* vx, char* vy, char* vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].v[0] = vx;
  cluster->s[id].v[1] = vy;
  cluster->s[id].v[2] = vz;
  return 0;
}
int get_velocity(int id, double* vx, double* vy, double* vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *vx = cluster->s[id].v[0].toDouble();
  *vy = cluster->s[id].v[1].toDouble();
  *vz = cluster->s[id].v[2].toDouble();
  return 0;
}
int set_velocity(int id, double vx, double vy, double vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].v[0] = (mpreal)vx;
  cluster->s[id].v[1] = (mpreal)vy;
  cluster->s[id].v[2] = (mpreal)vz;
  return 0;
}


int get_state_string(int id, char** m, char** x, char** y, char** z, char** vx, char** vy, char** vz, char** radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
    arg0=cluster->s[id].m.toString();
    arg1=cluster->s[id].x[0].toString();
    arg2=cluster->s[id].x[1].toString();
    arg3=cluster->s[id].x[2].toString();
    arg4=cluster->s[id].v[0].toString();
    arg5=cluster->s[id].v[1].toString();
    arg6=cluster->s[id].v[2].toString();
    arg7= cluster->s[id].r.toString();
  *m = (char*) arg0.c_str();
  *x = (char*) arg1.c_str();
  *y = (char*) arg2.c_str();
  *z = (char*) arg3.c_str();
  *vx = (char*) arg4.c_str();
  *vy = (char*) arg5.c_str();
  *vz = (char*) arg6.c_str();
  *radius = (char*) arg7.c_str();
  return 0;
}
int set_state_string(int id, char* m, char* x, char* y, char* z, char* vx, char* vy, char* vz, char* radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].r = radius;
  cluster->s[id].m = m;
  cluster->s[id].x[0] = x;
  cluster->s[id].x[1] = y;
  cluster->s[id].x[2] = z;
  cluster->s[id].v[0] = vx;
  cluster->s[id].v[1] = vy;
  cluster->s[id].v[2] = vz;
  return 0;
}
int get_state(int id, double* m, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *radius = cluster->s[id].r.toDouble();
  *m = cluster->s[id].m.toDouble();
  *x = cluster->s[id].x[0].toDouble();
  *y = cluster->s[id].x[1].toDouble();
  *z = cluster->s[id].x[2].toDouble();
  *vx = cluster->s[id].v[0].toDouble();
  *vy = cluster->s[id].v[1].toDouble();
  *vz = cluster->s[id].v[2].toDouble();
  return 0;
}
int set_state(int id, double m, double x, double y, double z, double vx, double vy, double vz, double radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].r = (mpreal)radius;
  cluster->s[id].m = (mpreal)m;
  cluster->s[id].x[0] = (mpreal)x;
  cluster->s[id].x[1] = (mpreal)y;
  cluster->s[id].x[2] = (mpreal)z;
  cluster->s[id].v[0] = (mpreal)vx;
  cluster->s[id].v[1] = (mpreal)vy;
  cluster->s[id].v[2] = (mpreal)vz;
  return 0;
}

int get_radius_string(int id, char** radius){
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  arg7 = cluster->s[id].r.toString();
  *radius = (char*) arg7.c_str();
  return 0;
}
int set_radius_string(int id, char* radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].r = radius;
  return 0;
}
int get_radius(int id, double* radius){
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *radius = cluster->s[id].r.toDouble();
  return 0;
}
int set_radius(int id, double radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].r = (mpreal)radius;
  return 0;
}

// Evolve
int evolve_model(double t_end) {
    DEBUG::DEB <<  "evolve_model  "<< "\n"; DEBUG::DEB.flush();
    brutus->evolve((mpreal)t_end);
    t = (mpreal)t_end;
//    data = brutus->get_data();
    return 0;
}

int synchronize_model() {
    DEBUG::DEB <<  "synchronize_model  "<< "\n"; DEBUG::DEB.flush();
//    data = brutus->get_data();
    *cluster=brutus->cl;
    return 0;
}
int cleanup_code() {
    odata.close();
    delete brutus;
    particle_id_counter = 0;
    return 0;
}

int delete_particle(int id) {
  return -2;
}
int recommit_particles() {
    DEBUG::DEB <<  "recommit_particles  "<< "\n"; DEBUG::DEB.flush();
    brutus->cl=*cluster;
    return 0;
//  return -2;
}

int set_brutus_output_directory(char *output_directory){
    out_directory = std::string(output_directory);
    if(out_directory.length() > 0){
        if(*out_directory.rbegin() != '/'){
            out_directory.append("/");
        }
    }
    return 0;
}
int get_brutus_output_directory(char **output_directory){
    *output_directory = (char*) out_directory.c_str();
    return 0;
}

int get_potential(int id, double* pot){return -2;}
int get_gravity_at_point(double m, double x, double y, double z, double* rx, double* ry, double* rz){return -2;}
int get_number_of_particles(int* N){return -2;}
int get_potential_at_point(double m, double x, double y, double z, double* p){return -2;}
int get_total_radius(double* R){return -2;}
int get_index_of_first_particle(int* id){return -2;}
int get_index_of_next_particle(int id, int* idnext){return -2;}

int get_total_mass_string(char **M){
//  int N = data.size()/arg_cnt;
  int N = cluster->s.size();
  mpreal Mtot; //constructor sets to zero by default = "0";
  for(int i=0; i<N; i++) {
//    Mtot += data[i*arg_cnt+arg_m];
    Mtot += cluster->s[i].m;
  }
  arg8=Mtot.toString();
  *M = (char*) arg8.c_str();
  return 0;
}

int get_total_mass(double* M){
//  int N = data.size()/arg_cnt;
  int N = cluster->s.size();
  mpreal Mtot; //constructor sets to zero by default = "0";
  for(int i=0; i<N; i++) {
//    Mtot += data[i*arg_cnt+arg_m];
    Mtot += cluster->s[i].m;
  }
  *M = Mtot.toDouble();
  return 0;
}
int get_potential_energy_m(mpreal* ep) {
//  int N = data.size()/arg_cnt;
  int N = cluster->s.size();
  mpreal eptot; //constructor sets to zero by default = "0";
  for(int i=0; i<N-1; i++) {
    mpreal mi = cluster->s[i].m;
    mpreal xi = cluster->s[i].x[0];
    mpreal yi = cluster->s[i].x[1];
    mpreal zi = cluster->s[i].x[2];
    for(int j=i+1; j<N; j++) {
      mpreal mj = cluster->s[j].m;
      mpreal xj = cluster->s[j].x[0];
      mpreal yj = cluster->s[j].x[1];
      mpreal zj = cluster->s[j].x[2];
      mpreal dx = xj - xi;
      mpreal dy = yj - yi;
      mpreal dz = zj - zi;
      mpreal dr2 = dx*dx + dy*dy + dz*dz;

      eptot -= mi*mj/sqrt(dr2);
    }
  }

  *ep = eptot;
  return 0;
}


int get_potential_energy_string( char **ep) {
    mpreal eptot; //constructor sets to zero by default = "0";
    get_potential_energy_m(&eptot);
    arg9=eptot.toString();
    *ep =(char*) arg9.c_str();
    return 0;
}

int get_kinetic_energy_m(mpreal* ek) {
//  int N = data.size()/arg_cnt;
  int N = cluster->s.size();
  mpreal ektot; //constructor sets to zero by default = "0";
  for(int i=0; i<N; i++) {
    mpreal m  = cluster->s[i].m;
    mpreal vx = cluster->s[i].x[0];
    mpreal vy = cluster->s[i].x[1];
    mpreal vz = cluster->s[i].x[2];
    mpreal v2 = vx*vx + vy*vy + vz*vz;
    ektot += "0.5"*m*v2;
  }
  *ek = ektot;
  return 0;
}


int get_kinetic_energy_string( char **ep) {
    mpreal ektot; //constructor sets to zero by default = "0";
    get_kinetic_energy_m(&ektot);
    arg9=ektot.toString();
    *ep =(char*) arg9.c_str();
    return 0;
}


int get_total_energy_string( char **ep) {
    mpreal ektot; //constructor sets to zero by default = "0";
    mpreal eptot; //constructor sets to zero by default = "0";
    mpreal etot; //constructor sets to zero by default = "0";
    get_potential_energy_m(&eptot);
    get_kinetic_energy_m(&ektot);
    etot = ektot + eptot;
    arg9=etot.toString();
    *ep =(char*) arg9.c_str();
    return 0;
}

int get_potential_energy(double* ep) {
  mpreal eptot; //constructor sets to zero by default = "0";
  get_potential_energy_m(&eptot);
  *ep = eptot.toDouble();
  return 0;
}

int get_kinetic_energy(double* ek) {
  mpreal ektot; //constructor sets to zero by default = "0";
  get_kinetic_energy_m(&ektot);
  *ek = ektot.toDouble();
  return 0;
}

int get_center_of_mass_position(double* x , double* y, double* z){
//  int N = data.size()/arg_cnt;
  int N = cluster->s.size();
  mpreal Mtot; //constructor sets to zero by default = "0";
  for(int i=0; i<N; i++) {
//    Mtot += data[i*arg_cnt+arg_m];
    Mtot += cluster->s[i].m;
  }

  vector<mpreal> rcm(3,"0");
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      rcm[j] += cluster->s[i].m*cluster->s[i].x[j];
    }
  }
  for(int i=0; i<3; i++) rcm[i] /= Mtot;

  *x = rcm[0].toDouble();
  *y = rcm[1].toDouble();
  *z = rcm[2].toDouble();
  return 0;
}
int get_center_of_mass_velocity(double* vx, double* vy, double* vz){
//  int N = data.size()/arg_cnt;
  int N = cluster->s.size();
  mpreal Mtot ; //constructor sets to zero by default ;"0";
  for(int i=0; i<N; i++) {
//    Mtot += data[i*arg_cnt+arg_m];
    Mtot += cluster->s[i].m;
  }

  vector<mpreal> vcm(3,"0");
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      vcm[j] += cluster->s[i].m*cluster->s[i].v[j];
    }
  }
  for(int i=0; i<3; i++) vcm[i] /= Mtot;

  *vx = vcm[0].toDouble();
  *vy = vcm[1].toDouble();
  *vz = vcm[2].toDouble();
  return 0;
}

int get_acceleration(int id, double* ax, double* ay, double* az){return -2;}
int set_acceleration(int id, double ax, double ay, double az){return -2;}

int add_step_acceleration_float64(int id, double a_step_x, double a_step_y, double a_step_z)
#ifdef use_additional_acc
{
    DEBUG::DEB <<  "add_step_acceleration_float64  "<< "\n"; DEBUG::DEB.flush();
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  cluster->s[id].a_step[0]=a_step_x;
  cluster->s[id].a_step[1]=a_step_y;
  cluster->s[id].a_step[2]=a_step_z;

  brutus->cl=*cluster;
  return 0;
}
#else
{return -2;}
#endif // use_additional_acc
