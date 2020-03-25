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

int particle_id_counter = 0;

int numBits = 64;  
int numDigits = numBits/4;

string out_directory;
std::map<int, int> local_index_map;

/*
 * We need this result_strings array to ensure that
 * C++ strings are not reclaimed before the function ends
 */
string result_strings[10];

Brutus *brutus = NULL;

mpreal t_begin = "0";
mpreal eta = "0.24";
mpreal t = "0";

mpreal epsilon = "1e-6"; // Bulirsch-Stoer tolerance

vector<mpreal> data, data_radius;

////////////////////////////////////////////////////////
// Amuse interface functions
////////////////////////////////////////////////////////
int initialize_code() {
    odata.open("temp.log");

    mpreal::set_default_prec(numBits);  

    brutus = new Brutus();

    particle_id_counter = 0;
    data.clear();
    data_radius.clear();

    t_begin = "0";
    t = t_begin;

    return 0;
}

// functions with "_string" assign strings to mpreals, and without "_string" assign doubles to mpreals

int new_particle_string(int *particle_identifier, char* mass, 
        char* x, char* y, char* z, char* vx, char* vy, char* vz, char* radius) {

    data.push_back(mass);
    data.push_back(x);
    data.push_back(y);
    data.push_back(z);
    data.push_back(vx);
    data.push_back(vy);
    data.push_back(vz);

    data_radius.push_back(radius);

    *particle_identifier = particle_id_counter;
    particle_id_counter++;

    return 0;
}
int new_particle_float64(int *particle_identifier, double mass, 
        double x, double y, double z, double vx, double vy, double vz, double radius) {

    data.push_back( (mpreal)mass );
    data.push_back( (mpreal)x );
    data.push_back( (mpreal)y );
    data.push_back( (mpreal)z );
    data.push_back( (mpreal)vx );
    data.push_back( (mpreal)vy );
    data.push_back( (mpreal)vz );

    data_radius.push_back( (mpreal)radius );

    *particle_identifier = particle_id_counter;
    particle_id_counter++;

    return 0;
}

int commit_particles() {
    brutus->set_data(data);

    int numStar = data.size()/7;
    brutus->setup();

    return 0;
}

// Begin time
int set_t_begin_string(char* tb) {
    t_begin = tb;
    return 0;
}
int get_t_begin_string(char **tb) {
    result_strings[0] = t_begin.toString();
    *tb = (char*) result_strings[0].c_str();
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
    result_strings[0] = eta.toString();
    *myeta = (char*) result_strings[0].c_str();
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
    result_strings[0] = t.toString();
    *tt = (char*) result_strings[0].c_str();
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
    result_strings[0] = epsilon.toString();
    *bs_tolerance = (char*) result_strings[0].c_str();
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

    numBits = mynumBits;
    mpreal::set_default_prec(numBits);
    numDigits = (int)abs(log10( pow("2.0", -numBits) )).toLong();
    brutus->set_numBits(numBits);

odata << numBits << " " << brutus->get_numBits() << endl;

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
std::string mass_string;
int get_mass_string(int id, char **mass) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  mass_string = data[id*7+0].toString();
  *mass = (char*) mass_string.c_str();
  return 0;
} 
int set_mass_string(int id, char *mass) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data[id*7+0] = mass;
  return 0;
}
int get_mass(int id, double* mass) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *mass = data[id*7+0].toDouble();
  return 0;
} 
int set_mass(int id, double mass) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data[id*7+0] = (mpreal)mass;
  return 0;
}
std::string position_string_x;
std::string position_string_y;
std::string position_string_z;
int get_position_string(int id, char **x, char **y, char **z) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  position_string_x = data[id*7+1].toString();
  position_string_y = data[id*7+2].toString();
  position_string_z = data[id*7+3].toString();
  *x = (char*) position_string_x.c_str();
  *y = (char*) position_string_y.c_str();
  *z = (char*) position_string_z.c_str();
  return 0;
}
int set_position_string(int id, char *x, char *y, char *z) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data[id*7+1] = x;
  data[id*7+2] = y;
  data[id*7+3] = z;
  return 0;
}
int get_position(int id, double* x, double* y, double* z) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *x = data[id*7+1].toDouble();
  *y = data[id*7+2].toDouble();
  *z = data[id*7+3].toDouble();
  return 0;
}
int set_position(int id, double x, double y, double z) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data[id*7+1] = (mpreal)x;
  data[id*7+2] = (mpreal)y;
  data[id*7+3] = (mpreal)z;
  return 0;
}

std::string velocity_string_vx;
std::string velocity_string_vy;
std::string velocity_string_vz;
int get_velocity_string(int id, char **vx, char **vy, char **vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  velocity_string_vx = data[id*7+4].toString();
  velocity_string_vy = data[id*7+5].toString();
  velocity_string_vz = data[id*7+6].toString();
  *vx = (char*) velocity_string_vx.c_str();
  *vy = (char*) velocity_string_vy.c_str();
  *vz = (char*) velocity_string_vz.c_str();
  return 0;
}
int set_velocity_string(int id, char* vx, char* vy, char* vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data[id*7+4] = vx;
  data[id*7+5] = vy;
  data[id*7+6] = vz;
  return 0;
}
int get_velocity(int id, double* vx, double* vy, double* vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *vx = data[id*7+4].toDouble();
  *vy = data[id*7+5].toDouble();
  *vz = data[id*7+6].toDouble();
  return 0;
}
int set_velocity(int id, double vx, double vy, double vz) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data[id*7+4] = (mpreal)vx;
  data[id*7+5] = (mpreal)vy;
  data[id*7+6] = (mpreal)vz;
  return 0;
}

std::string get_state_strings_m;
std::string get_state_strings_x;
std::string get_state_strings_y;
std::string get_state_strings_z;
std::string get_state_strings_vx;
std::string get_state_strings_vy;
std::string get_state_strings_vz;
std::string get_state_strings_r;
int get_state_string(int id, char** m, char** x, char** y, char** z, char** vx, char** vy, char** vz, char** radius) {  
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
    get_state_strings_m=data[id*7+0].toString();
    get_state_strings_x=data[id*7+1].toString();
    get_state_strings_y=data[id*7+2].toString();
    get_state_strings_z=data[id*7+3].toString();
    get_state_strings_vx=data[id*7+4].toString();
    get_state_strings_vy=data[id*7+5].toString();
    get_state_strings_vz=data[id*7+6].toString();
    get_state_strings_r= data_radius[id].toString();
  *m = (char*) get_state_strings_m.c_str();
  *x = (char*) get_state_strings_x.c_str();
  *y = (char*) get_state_strings_y.c_str();
  *z = (char*) get_state_strings_z.c_str();
  *vx = (char*) get_state_strings_vx.c_str();
  *vy = (char*) get_state_strings_vy.c_str();
  *vz = (char*) get_state_strings_vz.c_str();
  *radius = (char*) get_state_strings_r.c_str();
  return 0;
}
int set_state_string(int id, char* m, char* x, char* y, char* z, char* vx, char* vy, char* vz, char* radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data_radius[id] = radius;
  data[id*7+0] = m;
  data[id*7+1] = x;
  data[id*7+2] = y;
  data[id*7+3] = z;
  data[id*7+4] = vx;
  data[id*7+5] = vy;
  data[id*7+6] = vz;
  return 0;
}
int get_state(int id, double* m, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *radius = data_radius[id].toDouble();
  *m = data[id*7+0].toDouble();
  *x = data[id*7+1].toDouble();
  *y = data[id*7+2].toDouble();
  *z = data[id*7+3].toDouble();
  *vx = data[id*7+4].toDouble();
  *vy = data[id*7+5].toDouble();
  *vz = data[id*7+6].toDouble();
  return 0;
}
int set_state(int id, double m, double x, double y, double z, double vx, double vy, double vz, double radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data_radius[id] = (mpreal)radius;
  data[id*7+0] = (mpreal)m;
  data[id*7+1] = (mpreal)x;
  data[id*7+2] = (mpreal)y;
  data[id*7+3] = (mpreal)z;
  data[id*7+4] = (mpreal)vx;
  data[id*7+5] = (mpreal)vy;
  data[id*7+6] = (mpreal)vz;
  return 0;
}

std::string radius_string;
int get_radius_string(int id, char** radius){ 
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  radius_string = data_radius[id].toString();
  *radius = (char*) radius_string.c_str();
  return 0;
}
int set_radius_string(int id, char* radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data_radius[id] = radius;
  return 0;
}
int get_radius(int id, double* radius){ 
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  *radius = data_radius[id].toDouble();
  return 0;
}
int set_radius(int id, double radius) {
  if (id < 0 || id >= particle_id_counter){
    return -1;
  }
  data_radius[id] = (mpreal)radius;
  return 0;
}

// Evolve
int evolve_model(double t_end) {
    brutus->evolve((mpreal)t_end);
    t = (mpreal)t_end;
    data = brutus->get_data();
    return 0;
}

int synchronize_model() {
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
  return -2;
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

std::string total_mass_string;
int get_total_mass_string(char **M){ 
  int N = data.size()/7;
  mpreal Mtot = "0";
  for(int i=0; i<N; i++) {
    Mtot += data[i*7];
  }
  total_mass_string=Mtot.toString();
  *M = (char*) total_mass_string.c_str();
  return 0;
}

int get_total_mass(double* M){ 
  int N = data.size()/7;
  mpreal Mtot = "0";
  for(int i=0; i<N; i++) {
    Mtot += data[i*7];
  }
  *M = Mtot.toDouble();
  return 0;
}
int get_potential_energy_m(mpreal* ep) {
  int N = data.size()/7;
  mpreal eptot = "0";
  for(int i=0; i<N-1; i++) {
    mpreal mi = data[i*7];
    mpreal xi = data[i*7+1];
    mpreal yi = data[i*7+2];
    mpreal zi = data[i*7+3];
    for(int j=i+1; j<N; j++) {
      mpreal mj = data[j*7];
      mpreal xj = data[j*7+1];
      mpreal yj = data[j*7+2];
      mpreal zj = data[j*7+3];

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

std::string potential_energy_string;

int get_potential_energy_string( char **ep) {
    mpreal eptot = "0";
    get_potential_energy_m(&eptot);
    potential_energy_string=eptot.toString();
    *ep =(char*) potential_energy_string.c_str();
    return 0;
}

int get_kinetic_energy_m(mpreal* ek) {
  int N = data.size()/7;
  mpreal ektot = "0";
  for(int i=0; i<N; i++) {
    mpreal m  = data[i*7];
    mpreal vx = data[i*7+4];
    mpreal vy = data[i*7+5];
    mpreal vz = data[i*7+6];
    mpreal v2 = vx*vx + vy*vy + vz*vz;
    ektot += "0.5"*m*v2;
  }
  *ek = ektot;
  return 0;
}

std::string kinetic_energy_string;

int get_kinetic_energy_string( char **ep) {
    mpreal ektot = "0";
    get_kinetic_energy_m(&ektot);
    kinetic_energy_string=ektot.toString();
    *ep =(char*) kinetic_energy_string.c_str();
    return 0;
}

std::string total_energy_string;

int get_total_energy_string( char **ep) {
    mpreal ektot = "0";   
    mpreal eptot = "0";   
    mpreal etot = "0";
    get_potential_energy_m(&eptot);
    get_kinetic_energy_m(&ektot);
    etot = ektot + eptot;
    total_energy_string=etot.toString();
    *ep =(char*) total_energy_string.c_str();
    return 0;
}

int get_potential_energy(double* ep) {
  mpreal eptot = "0";
  get_potential_energy_m(&eptot);
  *ep = eptot.toDouble();
  return 0;
}

int get_kinetic_energy(double* ek) {
  mpreal ektot = "0";
  get_kinetic_energy_m(&ektot);
  *ek = ektot.toDouble();
  return 0;
}

int get_center_of_mass_position(double* x , double* y, double* z){ 
  int N = data.size()/7;
  mpreal Mtot = "0";
  for(int i=0; i<N; i++) {
    Mtot += data[i*7];
  }

  vector<mpreal> rcm(3,"0");
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      rcm[j] += data[i*7]*data[i*7+(j+1)];
    }
  }
  for(int i=0; i<3; i++) rcm[i] /= Mtot;  

  *x = rcm[0].toDouble();
  *y = rcm[1].toDouble();
  *z = rcm[2].toDouble();
  return 0;
}
int get_center_of_mass_velocity(double* vx, double* vy, double* vz){ 
  int N = data.size()/7;
  mpreal Mtot = "0";
  for(int i=0; i<N; i++) {
    Mtot += data[i*7];
  }

  vector<mpreal> vcm(3,"0");
  for(int i=0; i<N; i++) {
    for(int j=0; j<3; j++) {
      vcm[j] += data[i*7]*data[i*7+(j+4)];
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

