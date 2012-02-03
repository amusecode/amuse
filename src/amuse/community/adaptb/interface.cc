////////////////////////////////////////////////////////
// Includes
////////////////////////////////////////////////////////
#include <iostream>
#include <map>
using namespace std;

#include "Bs_integrator.h"
#include "Clock.h"

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
string dt_print, dt_max, dt_factor;  
string epsilon;
string softening;
string t_lim;
int sim_state = 0;
int n_max, k_max;
double E0, E, dE;

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

  dt_print = "0.1";		// Regular print intervals
  dt_max = "0.01";		// Maximum time steps
  dt_factor = "0.1";		// time step multiplication factor

  epsilon = "1e-6";		// Bulirsch-Stoer tolerance
  softening = "0";		// Gravitational softening

  t_lim = "3600"; 		// Maximum CPU time in seconds

  n_max = 64;			// Bulirsch-Stoer sub step variables
  k_max = 64;
  
  cluster = new Cluster();
  return 0;
}

int new_particle_string(int *particle_identifier, char* mass, char* radius, char* x, char* y, char* z, char* vx, char* vy, char* vz) {
  cluster->add_star(particle_id_counter, mass, radius, x, y, z, vx, vy, vz);
  *particle_identifier = particle_id_counter;
  particle_id_counter++;
  return 0;
}
int new_particle_float64(int *particle_identifier, double mass, double radius, double x, double y, double z, double vx, double vy, double vz) {
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
    E0 = cluster->get_E();

  return 0;
}

// Begin and simulation time, t_begin and t_sim
//int set_t_begin_and_t_sim(double t_begin, double t_sim) {
//    myclock->set_t_begin_and_end(t_begin, t_sim);
//    return 0;
//}
// Bulirsch-Stoer tolerance, epsilon
int set_bs_tolerance(char *bs_tolerance) {
    epsilon = std::string(bs_tolerance);
    return 0;
}
int get_bs_tolerance(char **bs_tolerance) {
    *bs_tolerance = (char*)epsilon.c_str();
    return 0;
}
// Word-length, Lw in mantissa
int set_word_length(int myLw) {
    numBits = myLw;

    cout << "New word-length = " << numBits << endl; 

    return 0;
}
int get_word_length(int *myLw) {
    *myLw = numBits;
    return 0;
}
// Softening squared
int set_eps2(char *eps2) {
  string eps2_str = std::string(eps2);
  cluster->force.softening_sq = eps2_str.c_str();
  return 0;
}
int get_eps2(char **eps2) {
  string eps2_str = cluster->force.softening_sq;
  *eps2 = (char*)eps2_str.c_str();
  return 0;
}
// Regular print intervals
int set_dt_print(char *print_interval) {
  dt_print = std::string(print_interval);
  mpreal dt = dt_print.c_str();
  dt /= "10";
  dt_max = dt.to_string();
  return 0;
}
int get_dt_print(char **print_interval) {
  *print_interval = (char*)dt_print.c_str();
  return 0;
}
// max cpu time
int set_max_cpu_time(char *max_cpu_time) {
  t_lim = std::string(max_cpu_time);
  return 0;
}
int get_max_cpu_time(char **max_cpu_time) {
  *max_cpu_time = (char*)t_lim.c_str();
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
    *output_directory = (char*)out_directory.c_str();
    return 0;
}

int commit_parameters() {
  mpreal::set_default_prec(numBits);  
  cout.precision(numBits/4);
  odata.precision(numBits/4);

  bs = new Bs_integrator(epsilon.c_str(), n_max, k_max);
  myclock = new Clock("0", "1", dt_print.c_str(), dt_max.c_str(), dt_factor.c_str());
  cluster->get_pointer_to_force()->softening_sq = softening.c_str();

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
  *mass = cluster->get_pointer_to_star(id)->m;
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
  *x = cluster->get_pointer_to_star(id)->x;
  *y = cluster->get_pointer_to_star(id)->y;
  *z = cluster->get_pointer_to_star(id)->z;
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
  *vx = cluster->get_pointer_to_star(id)->vx;
  *vy = cluster->get_pointer_to_star(id)->vy;
  *vz = cluster->get_pointer_to_star(id)->vz;
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
int get_state(int id, double* m, double* radius, double* x, double* y, double* z, double* vx, double* vy, double* vz) {
  if (id < 0 || id >= particle_id_counter){
    return -3;
  }
  *radius = cluster->get_pointer_to_star(id)->radius;
  *m = cluster->get_pointer_to_star(id)->m;
  *x = cluster->get_pointer_to_star(id)->x;
  *y = cluster->get_pointer_to_star(id)->y;
  *z = cluster->get_pointer_to_star(id)->z;
  *vx = cluster->get_pointer_to_star(id)->vx;
  *vy = cluster->get_pointer_to_star(id)->vy;
  *vz = cluster->get_pointer_to_star(id)->vz;
  return 0;
}
int set_state(int id, double m, double radius, double x, double y, double z, double vx, double vy, double vz) {
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
  *radius = cluster->get_pointer_to_star(id)->radius;
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
    cout << myclock->get_t() << endl;
    myclock->set_t_begin_and_end(myclock->get_t(), t);
    while( !myclock->alarm() )
    {
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

        mpreal tcpu_lim = t_lim.c_str();
	if( myclock->read() > tcpu_lim) {
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

    E = cluster->get_E();
    dE = log10(abs((E-E0)/E0));

    return 0;
}
int cleanup_code() {
  if (odata.is_open()){
    odata.close();
  }
  mpreal t_sim = cluster->get_t();
  myclock->stop_timer();

  odata.open( (out_directory + file_log).c_str() );
  if( !odata )
  {
    cerr << "Could not open " << file_log << "!" << endl;
    return -1;
  }
  else
  {
    odata << "sim_state = " << sim_state << endl;

    odata << "N         = " << cluster->get_N() << endl;

    odata << "t_sim     = " << t_sim << endl;
    odata << "dt_print  = " << dt_print << endl;
    odata << "dt_max    = " << dt_max << endl;
    odata << "dt_factor = " << dt_factor << endl;

    odata << "epsilon   = " << epsilon << endl;
    odata << "numBits   = " << numBits << endl;

    odata << "softening = " << softening << endl;

    odata << "t_cpu     = " << myclock->get_timer() << endl;
    odata << "dE        = " << dE << endl;
  }
  delete cluster;
  delete myclock;
  delete bs;
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
int synchronize_model(){return -2;}
int get_time_step(double* dt){return -2;}
int get_kinetic_energy(double* ek){return -2;}
int get_potential_energy(double* ep){return -2;}
int get_center_of_mass_velocity(double* vx, double* vy, double* vz){return -2;}
int get_acceleration(int id, double* ax, double* ay, double* az){return -2;}
int set_acceleration(int id, double ax, double ay, double az){return -2;}

