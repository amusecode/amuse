////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////
#include <iostream>
#include <map>
#include "Bs_integrator.h"
#include "Clock.h"

using namespace mpfr;
using namespace std;

////////////////////////////////////////////////////////
// Set precision
////////////////////////////////////////////////////////
int numBits = 64;  
int Lw = numBits/4;
ofstream odata;

////////////////////////////////////////////////////////
// Declare global variables
////////////////////////////////////////////////////////
string out_directory, file_out, file_log;
int sim_state = 0;
int n_max, k_max;
// 'global' mpreals can only be implemented as static members of a class:
class mpreal_globals {
    public:
    static mpreal epsilon;   // Bulirsch-Stoer tolerance
    static mpreal dt_print;
    static mpreal dt_max;
    static mpreal dt_factor;
    static mpreal t_lim;
    static mpreal E0;
    static mpreal E;
    static mpreal dE;
};
mpreal mpreal_globals::epsilon = "1.0e-6"; // Bulirsch-Stoer tolerance
mpreal mpreal_globals::dt_print = "0.1";   // Regular print intervals
mpreal mpreal_globals::dt_max = "0.01";    // Maximum time steps
mpreal mpreal_globals::dt_factor = "0.1";  // time step multiplication factor
mpreal mpreal_globals::t_lim = "3600";     // Maximum CPU time in seconds
mpreal mpreal_globals::E = "0";
mpreal mpreal_globals::E0 = "0";
mpreal mpreal_globals::dE = "0";

////////////////////////////////////////////////////////
// Declare global objects 
////////////////////////////////////////////////////////
std::map<int, int> local_index_map;
int particle_id_counter = 0;
Cluster *cluster = NULL;
Bs_integrator *bs = NULL;
Clock *myclock = NULL;

////////////////////////////////////////////////////////
// Amuse interface functions
////////////////////////////////////////////////////////
int initialize_code() {
    mpreal::set_default_prec(numBits);  
    cout.precision(Lw);
    odata.precision(Lw);
    
    out_directory = "./";
    file_out = "file.out";	// Outfile for phase space coordinates
    file_log = "file.log";	// Outfile for input/output numbers
    
    n_max = 64;			// Bulirsch-Stoer sub step variables
    k_max = 64;
    
    cluster = new Cluster();
    return 0;
}

int new_particle_string(int *particle_identifier, char* mass, 
        char* x, char* y, char* z, char* vx, char* vy, char* vz, char* radius) {
    cluster->add_star(particle_id_counter, mass, radius, x, y, z, vx, vy, vz);
    *particle_identifier = particle_id_counter;
    particle_id_counter++;
    return 0;
}
int new_particle_float64(int *particle_identifier, double mass, 
        double x, double y, double z, double vx, double vy, double vz, double radius) {
    cluster->add_star(particle_id_counter, mass, radius, x, y, z, vx, vy, vz);
    *particle_identifier = particle_id_counter;
    particle_id_counter++;
    return 0;
}
int commit_particles() {
    cluster->set_N(particle_id_counter);
    cluster->print(odata);
    
    cerr << endl;
    cerr << myclock->get_progress() << "%" << endl;
    
    // Initial calculations
    mpreal_globals::E0 = cluster->get_E();
    return 0;
}

// Bulirsch-Stoer tolerance, epsilon
int set_bs_tolerance_string(char *bs_tolerance) {
    mpreal_globals::epsilon = bs_tolerance;
    return 0;
}
int get_bs_tolerance_string(char **bs_tolerance) {
    *bs_tolerance = (char*) mpreal_globals::epsilon.toString().c_str();
    return 0;
}
int set_bs_tolerance_float64(double bs_tolerance) {
    mpreal_globals::epsilon = (mpreal) bs_tolerance;
    return 0;
}
int get_bs_tolerance_float64(double *bs_tolerance) {
    *bs_tolerance = mpreal_globals::epsilon.toDouble();
    return 0;
}
// Word-length, Lw in mantissa
int set_word_length(int myLw) {
    numBits = myLw;
    mpreal::set_default_prec(numBits);  // Do this now, to make sure that ...
    // ... other (mpreal) parameters are stored at this precision!
    cout << "New word-length = " << numBits << endl; 
    return 0;
}
int get_word_length(int *myLw) {
    *myLw = numBits;
    return 0;
}
// Softening squared
int set_eps2(double eps2) {
    cluster->force.softening_sq = (mpreal) eps2;
    return 0;
}
int get_eps2(double *eps2) {
    *eps2 = cluster->force.softening_sq.toDouble();
    return 0;
}
// Regular print intervals
int set_dt_print(double print_interval) {
    mpreal_globals::dt_print = print_interval;
    mpreal_globals::dt_max = mpreal_globals::dt_print / "10";
    return 0;
}
int get_dt_print(double *print_interval) {
    *print_interval = mpreal_globals::dt_print.toDouble();
    return 0;
}
// max cpu time
int set_max_cpu_time(double max_cpu_time) {
    mpreal_globals::t_lim = max_cpu_time;
    return 0;
}
int get_max_cpu_time(double *max_cpu_time) {
    *max_cpu_time = mpreal_globals::t_lim.toDouble();
    return 0;
}
int set_adaptb_output_directory(char *output_directory){
    out_directory = std::string(output_directory);
    if(out_directory.length() > 0){
        if(*out_directory.rbegin() != '/'){
            out_directory.append("/");
        }
    }
    return 0;
}
int get_adaptb_output_directory(char **output_directory){
    *output_directory = (char*) out_directory.c_str();
    return 0;
}
int get_time_step(double* dt){
    if (myclock == NULL) { // Parameters have not been committed yet -> return default
        *dt = mpreal_globals::dt_max.toDouble();
    } else {
        *dt = myclock->get_dt().toDouble();
    }
    return 0;
}

int commit_parameters() {
    cout.precision(numBits/4);
    odata.precision(numBits/4);
    
    bs = new Bs_integrator(mpreal_globals::epsilon, n_max, k_max);
    myclock = new Clock("0", "1", mpreal_globals::dt_print, mpreal_globals::dt_max, mpreal_globals::dt_factor);
    myclock->Start_timer();
    odata.open( (out_directory + file_out).c_str() );
    if( !odata ) {
        cerr << "Could not open " << (out_directory + file_out) << "!" << endl;
        sim_state = 1;
        return -1;
    }
    return 0;
}
int recommit_parameters() {
    return commit_parameters();
}

int get_mass(int id, double*mass) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  *mass = cluster->get_pointer_to_star(id)->m.toDouble();
  return 0;
} 
int set_mass(int id, double mass) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  cluster->get_pointer_to_star(id)->m = mass;
  return 0;
}
int get_position(int id, double* x, double* y, double* z) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  *x = cluster->get_pointer_to_star(id)->x.toDouble();
  *y = cluster->get_pointer_to_star(id)->y.toDouble();
  *z = cluster->get_pointer_to_star(id)->z.toDouble();
  return 0;
}
int set_position(int id, double x, double y, double z) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  cluster->get_pointer_to_star(id)->x = x;
  cluster->get_pointer_to_star(id)->y = y;
  cluster->get_pointer_to_star(id)->z = z;
  return 0;
}
int get_velocity(int id, double* vx, double* vy, double* vz) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  *vx = cluster->get_pointer_to_star(id)->vx.toDouble();
  *vy = cluster->get_pointer_to_star(id)->vy.toDouble();
  *vz = cluster->get_pointer_to_star(id)->vz.toDouble();
  return 0;
}
int set_velocity(int id, double vx, double vy, double vz) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  cluster->get_pointer_to_star(id)->vx = vx;
  cluster->get_pointer_to_star(id)->vy = vy;
  cluster->get_pointer_to_star(id)->vz = vz;
  return 0;
}
int get_state(int id, double* m, double* x, double* y, double* z, double* vx, double* vy, double* vz, double* radius) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  *radius = cluster->get_pointer_to_star(id)->radius.toDouble();
  *m = cluster->get_pointer_to_star(id)->m.toDouble();
  *x = cluster->get_pointer_to_star(id)->x.toDouble();
  *y = cluster->get_pointer_to_star(id)->y.toDouble();
  *z = cluster->get_pointer_to_star(id)->z.toDouble();
  *vx = cluster->get_pointer_to_star(id)->vx.toDouble();
  *vy = cluster->get_pointer_to_star(id)->vy.toDouble();
  *vz = cluster->get_pointer_to_star(id)->vz.toDouble();
  return 0;
}
int set_state(int id, double m, double x, double y, double z, double vx, double vy, double vz, double radius) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  cluster->get_pointer_to_star(id)->radius = radius;
  cluster->get_pointer_to_star(id)->m = m;
  cluster->get_pointer_to_star(id)->x = x;
  cluster->get_pointer_to_star(id)->y = y;
  cluster->get_pointer_to_star(id)->z = z;
  cluster->get_pointer_to_star(id)->vx = vx;
  cluster->get_pointer_to_star(id)->vy = vy;
  cluster->get_pointer_to_star(id)->vz = vz;
  return 0;
}
int get_radius(int id, double* radius){ 
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  *radius = cluster->get_pointer_to_star(id)->radius.toDouble();
  return 0;
}
int set_radius(int id, double radius) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  cluster->get_pointer_to_star(id)->radius = radius;
  return 0;
}

int evolve_model(double t) {
    myclock->set_t_begin_and_end(myclock->get_t(), t);
    while( !myclock->alarm() ) {
        cluster->calc_a_dt();
        myclock->calc_dt( cluster->get_dt() );
        
        mpreal dt = myclock->get_dt();
        bs->integrate(*cluster, dt);
        myclock->set_dt(dt);
        
        if( !bs->converged() ) {
            cerr << "No Convergence Reached, Simulation Aborted!" << endl;
            sim_state = 2;
            myclock->abort();
        } else {
            myclock->tick();
            cluster->set_t( myclock->get_t() );
            
            if( myclock->to_print() ) {
                Cluster cl_exp = *cluster;
                
                cl_exp.calc_a();
                mpreal dt_exp = myclock->get_t_print() - myclock->get_t();
                
                bs->integrate(cl_exp, dt_exp);
                
                cl_exp.set_t( myclock->get_t_print() );
                cl_exp.print(odata);
                
                cout << myclock->get_progress() << "%" << endl;
            }
            
            if( myclock->read() > mpreal_globals::t_lim) {
                sim_state = 3;
                myclock->abort();
            }
            
        }
    
    }
    
    cluster->calc_a();
    mpreal dt = myclock->get_t_end() - myclock->get_t();
    bs->integrate(*cluster, dt);
    cluster->set_t( myclock->get_t_end() );
    myclock->set_t(myclock->get_t_end());
    cluster->print(odata);    
    
    mpreal_globals::E = cluster->get_E();
    mpreal_globals::dE = log10(abs((mpreal_globals::E-mpreal_globals::E0)/mpreal_globals::E0));
    
    return 0;
}
int synchronize_model() {
    return 0;
}
int cleanup_code() {
    mpreal t_cpu;
    if (odata.is_open()){
        odata.close();
        myclock->stop_timer();
        t_cpu = myclock->get_timer();
        delete myclock;
        delete bs;
    } else {
        t_cpu = "0";
    }
    
    odata.open( (out_directory + file_log).c_str() );
    if ( !odata ) {
        cerr << "Could not open " << file_log << "!" << endl;
        return -1;
    } else {
        odata << "sim_state = " << sim_state << endl;
        odata << "N         = " << cluster->get_N() << endl;
        odata << "t_sim     = " << cluster->get_t() << endl;
        odata << "dt_print  = " << mpreal_globals::dt_print << endl;
        odata << "dt_max    = " << mpreal_globals::dt_max << endl;
        odata << "dt_factor = " << mpreal_globals::dt_factor << endl;
        odata << "epsilon   = " << mpreal_globals::epsilon << endl;
        odata << "numBits   = " << numBits << endl;
        odata << "softening = " << cluster->get_pointer_to_force()->softening_sq << endl;
        odata << "t_cpu     = " << t_cpu << endl;
        odata << "dE        = " << mpreal_globals::dE << endl;
    }
    delete cluster;
    return 0;
}

int delete_particle(int id) {
  return -2;
}
int recommit_particles() {
  return -2;
}

int get_potential(int id, double* pot){return -2;}
int get_gravity_at_point(double m, double x, double y, double z, double* rx, double* ry, double* rz){return -2;}
int get_number_of_particles(int* N){return -2;}
int get_time(double* time){return -2;}
int get_potential_at_point(double m, double x, double y, double z, double* p){return -2;}
int get_center_of_mass_position(double* x , double* y, double* z){return -2;}
int get_total_radius(double* R){return -2;}
int get_index_of_first_particle(int* id){return -2;}
int get_index_of_next_particle(int id, int* idnext){return -2;}
int get_total_mass(double* M){return -2;}
int get_kinetic_energy(double* ek){return -2;}
int get_potential_energy(double* ep){return -2;}
int get_center_of_mass_velocity(double* vx, double* vy, double* vz){return -2;}
int get_acceleration(int id, double* ax, double* ay, double* az){return -2;}
int set_acceleration(int id, double ax, double ay, double az){return -2;}

