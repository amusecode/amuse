#include <iostream>
#include <iomanip>
#include <string>
#include <mpi.h>
#include "src/lib/my_errors.h"
#include "src/lib/functions.h"
#include "src/lib/utilis.h"
#include <math.h>

using namespace std;

typedef struct {
    double mass;                                        /// mass
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
} dynamics_state;

double4 *pos_PH;
float4  *vel_PH;
double4 *pos_CH;
double4 *vel_CH;
double4 *a_H0;
double  *step;
double  *local_time;
unsigned int *devices;

map<int, dynamics_state> dm_states;
map<int,dynamics_state>::iterator it;
map<int,dynamics_state>::iterator et;

int size, rank, particle_id_counter = 0;
unsigned int N, NGPU, TPB, FMAX, BFMAX, MAXDIM, GPUMINTHREADS;
double GTIME, GTW, DTMAX, DTMIN, EPS, ETA6, ETA4, DTPRINT, ATIME, plummer_core, plummer_mass;  
bool warm_start;
string GPUNAME;


int echo(int input){
	return 12;
}

int echo_int(int input, int * output){
	*output = echo(input);
	return 0;
}

int initialize_code(){
	HostSafeCall(__MPIstart(&rank, &size, &mpi_float4, &mpi_double4));
	warm_start = 0;  	
	if(rank == 0){
		ofstream generic;
      generic.open("gpu_memory.dat", ios::out);
      generic.close();
      generic.open("cpu_memory.dat", ios::out);
      generic.close();
      generic.open("H6Blog.dat", ios::out);
      generic.close();
      generic.open("energy.dat", ios::out);
      generic<<scientific<<setprecision(5);
      double value = 0.0;
      generic<<value<<"  "<<value<<endl;
      generic.close();
	}
	return 0;
}

int recommit_parameters(){
	return 0;
}

int commit_parameters(){
	GTIME = 0.0;
	ATIME = 1.0e+10;
	GTW = 0.0;

	if(rank == 0){
	   
		ofstream hlog;
      hlog.open("H6Blog.dat", ios::app);
		hlog<<" ================================================== "<<endl;
      hlog<<" N                   : "<<N<<endl;
      hlog<<" Gpus                : "<<NGPU<<endl;
      hlog<<" Threads per block   : "<<TPB<<endl;
      hlog<<" Max time step	    : "<<DTMAX<<endl;
      hlog<<" Min time step	    : "<<DTMIN<<endl;
      hlog<<" Softening           : "<<EPS<<endl;
      hlog<<" eta 6th order       : "<<ETA6<<endl;
      hlog<<" eta 4th order       : "<<ETA4<<endl;
      hlog<<" time for printing   : "<<DTPRINT<<endl;
      hlog<<" Max output files    : "<<FMAX<<endl;
      hlog<<" ================================================== "<<endl;
      hlog.close();

   	HostSafeCall(isDivisible(N, size, NGPU, TPB, &BFMAX));

#ifdef CHECK_ERRORS
		cout<<" Check errors option ENABLED "<<endl;
#else
		cout<<" Check errors option DISABLED "<<endl;
#endif

#ifdef CHECK_TIMES
		ofstream generic;
      generic.open("times.dat", ios::out);
      generic.close();
      cout<<" Check times option ENABLED "<<endl;
#else
      cout<<" Check times option DISABLED "<<endl;
#endif

#ifdef PLUMMER
      cout<<" Plummer Potential option ENABLED "<<endl;
      ofstream hlog;
      hlog.open("H6Blog.dat", ios::app);
      hlog<<" Plummer (core) : "<<plummer_core<<endl;
      hlog<<" Plummer (mass) : "<<plummer_mass<<endl;
      hlog<<" =============================================== "<<endl;
      hlog.close();
#else
      cout<<" Plummer Potential option DISABLED "<<endl;
#endif
	
	   cout<<"NOTE_1: the code works with nbody units ( G = 1 ): please check the parameters, more info are given in the README file"<<endl;
		cout<<"NOTE_2: the evolve method requires an input time divisible for the maximum time step ( 'dt' / 'max_step' == integer ) "<<endl;
	}
   return 0;
}

int new_particle(int *id, double mass, double radius, double x, double y, double z, double vx, double vy, double vz){
	*id = particle_id_counter;
 	dynamics_state state;
	state.mass = mass;
   state.x = x;
   state.y = y;
   state.z = z;
   state.vx = vx;
   state.vy = vy;
   state.vz = vz;
	dm_states[particle_id_counter]= state;
	particle_id_counter++;
	return 0;
}

int commit_particles(){
	N = particle_id_counter;
   pos_PH = new double4 [N];
   vel_PH = new  float4 [N];
   pos_CH = new double4 [N];
   vel_CH = new double4 [N];
   a_H0 = new double4 [3*N];
   step = new double [N];
   local_time = new double [N];
   devices = new unsigned int [NGPU];

	it = dm_states.begin();
   
	for (int i=0; i<N; i++){
	
		pos_PH[i].w = it->second.mass;
      pos_PH[i].x = it->second.x;
      pos_PH[i].y = it->second.y;
      pos_PH[i].z = it->second.z;
		vel_CH[i].x = it->second.vx;
		vel_CH[i].y = it->second.vy;
		vel_CH[i].z = it->second.vz;
		vel_CH[i].w = 0.0;
		pos_CH[i].w = it->second.mass;
      pos_CH[i].x = it->second.x;
      pos_CH[i].y = it->second.y;
      pos_CH[i].z = it->second.z;
		vel_PH[i].x = vel_CH[i].x;
		vel_PH[i].y = vel_CH[i].y;
		vel_PH[i].z = vel_CH[i].z;
		it++;	
	}
	
	HostSafeCall(InitBlocks(pos_PH, vel_PH, TPB, N, BFMAX, ETA4, NGPU, EPS, &MAXDIM, DTMAX, &GPUMINTHREADS, devices, GPUNAME, rank, size, pos_CH, vel_CH, a_H0, step, local_time, &ATIME, plummer_core, plummer_mass));
	return 0;
}

int evolve_model(double dt){
	if(rank==0){
		ofstream hlog;
      hlog.open("H6Blog.dat", ios::app);
      hlog<<endl<<endl;
		hlog<<" ================================================== "<<endl;
      hlog<<" Time of integration : "<<dt<<endl;
      hlog<<" ================================================== "<<endl;
      hlog.close();
	}
   
	MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));

   FMAX = 1000000 + ((GTIME+GTW) / DTPRINT);
	
	HostSafeCall(Hermite6th(dt, &GTIME, &ATIME, local_time, step, N, pos_PH, vel_PH, pos_CH, vel_CH, a_H0, MAXDIM, NGPU, devices, TPB, rank, size, BFMAX, ETA6, ETA4, DTMAX, DTMIN, EPS, DTPRINT, FMAX, warm_start, GTW, GPUMINTHREADS, plummer_core, plummer_mass));
  
	warm_start = 1;
   GTW = 0.0;

	for (int i=0; i<N; i++){
      pos_PH[i].x = pos_CH[i].x;
     	pos_PH[i].y = pos_CH[i].y;
      pos_PH[i].z = pos_CH[i].z;
      vel_PH[i].x = vel_CH[i].x;
      vel_PH[i].y = vel_CH[i].y;
      vel_PH[i].z = vel_CH[i].z;
	}
	return 0;
}

int set_time_begin(double time_begin){
   GTIME = time_begin;
	return 0;
}


int get_time_begin(double *time_begin){
   *time_begin = GTIME;
   return 0;
}


int get_mass(int index_of_the_particle, double * mass){
   *mass = pos_CH[index_of_the_particle].w;
	return 0;
}

int get_time(double * time){
   *time = GTIME;
   return 0;
}

int set_mass(int index_of_the_particle, double mass){
   dm_states[index_of_the_particle].mass = mass;
   return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
   it = dm_states.begin();
   *index_of_the_particle = it->first;
   return 0;
}

int get_total_radius(double * radius){
   *radius = 0.0;
   for(int i=0; i<N; i++){
		double r=sqrt(pos_CH[i].x*pos_CH[i].x+pos_CH[i].y*pos_CH[i].y+pos_CH[i].z*pos_CH[i].z);
		if (r>*radius) *radius = r;
   }
   return 0;
}

int get_potential_at_point(double eps, double x, double y, double z, double * phi){
	*phi = 0.0;
   for(int i=0; i<N ; i++){
		double rij =(x - pos_CH[i].x) * (x - pos_CH[i].x) + (y-pos_CH[i].y) * (y - pos_CH[i].y) + (z - pos_CH[i].z) * (z - pos_CH[i].z);
      *phi += pos_CH[i].w / sqrt(rij + eps * eps); 
	}
   return 0;
}

int get_total_mass(double * mass){
   *mass = 0.0;
   for(int i=0; i<N; i++)  *mass += pos_CH[i].w;
   return 0;
}

int set_eps2(double epsilon){
   EPS = epsilon;
   return 0;
}

int get_eps2(double *epsilon){
   *epsilon = EPS;
   return 0;
}

int set_eps(double epsilon){
   EPS = epsilon;
   return 0;
}

int get_eps(double *epsilon){
   *epsilon = EPS;
   return 0;
}

int set_eta6(double eta6){
   ETA6 = eta6;
   return 0;
}

int get_eta6(double *eta6){
   *eta6 = ETA6;
   return 0;
}

int get_number_of_particles(int * number_of_particles){
   *number_of_particles = N;
   return 0;
}

int set_number_of_particles(int number_of_particles){
   N = number_of_particles;
   return 0;
}

int set_eta4(double eta4){
   ETA4 = eta4;
   return 0;
}

int get_eta4(double *eta4){
   *eta4 = ETA4;
   return 0;
}

int set_Plummer_mass(double Plummer_mass){
   plummer_mass = Plummer_mass;
   return 0;
}

int get_Plummer_mass(double *Plummer_mass){
   *Plummer_mass = plummer_mass;
   return 0;
}

int set_Plummer_core(double Plummer_core){
   plummer_core = Plummer_core;
   return 0;
}

int get_Plummer_core(double *Plummer_core){
   *Plummer_core = plummer_core;
   return 0;
}

int set_number_of_GPU(int number_of_GPU){
   NGPU = number_of_GPU;
   return 0;
}

int get_number_of_GPU(int * number_of_GPU){
   *number_of_GPU = NGPU;
   return 0;
}

int get_number_of_Threads(int * number_of_Threads){
   *number_of_Threads = TPB;
   return 0;
}

int set_number_of_Threads(int number_of_Threads){
   TPB = number_of_Threads;
   return 0;
}

int get_number_of_Print(int * number_of_Print){
   *number_of_Print = FMAX;
   return 0;
}

int set_number_of_Print(int number_of_Print){
   FMAX = number_of_Print;
   return 0;
}

int get_max_time_step(double * max_time_step){
   *max_time_step = DTMAX;
   return 0;
}

int set_max_time_step(double max_time_step){
   DTMAX = pow(2.0,max_time_step);
   return 0;
}

int get_min_time_step(double * min_time_step){
   *min_time_step = DTMIN;
   return 0;
}

int set_min_time_step(double min_time_step){
   DTMIN = pow(2.0,min_time_step);
   return 0;
}

int get_DTPrint(double * DTPrint){
   *DTPrint = DTPRINT;
   return 0;
}

int set_DTPrint(double DTPrint){
   DTPRINT = DTPrint;
   return 0;
}

int get_gpu_name(char ** gpu_name){
   *gpu_name = to_char(GPUNAME);
   return 0;
}

int set_gpu_name(char *gpu_name){
   GPUNAME = to_string(gpu_name);
   return 0;
}

int get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle){
   it = dm_states.begin();
   for(int i=0; i<=index_of_the_particle; i++) it++;
   cout<< it->first<<endl;
   return 0;
}

int delete_particle(int index_of_the_particle){
   dm_states.erase(index_of_the_particle);
	particle_id_counter--;
	return 0;
}

int get_potential(int index_of_the_particle, double * potential){ 
	*potential = 0.0;
   for(int i=0; i<N ; i++){
		if(i != index_of_the_particle){
			double rij =(pos_CH[index_of_the_particle].x - pos_CH[i].x)*(pos_CH[index_of_the_particle].x - pos_CH[i].x) +
				(pos_CH[index_of_the_particle].y - pos_CH[i].y) * (pos_CH[index_of_the_particle].y - pos_CH[i].y) +
				(pos_CH[index_of_the_particle].z - pos_CH[i].z) * (pos_CH[index_of_the_particle].z - pos_CH[i].z);
        	*potential += (pos_CH[i].w + pos_CH[index_of_the_particle].w) / sqrt(rij + EPS * EPS);
		} 
	}
   return 0;
}

int synchronize_model(){
   it = dm_states.begin();

   for (int i=0; i<N; i++){
		it->second.mass = pos_CH[i].w;
      it->second.x = pos_CH[i].x;
      it->second.y = pos_CH[i].y;
      it->second.z = pos_CH[i].z;
      it->second.vx = vel_CH[i].x;
		it->second.vy = vel_CH[i].y;
		it->second.vz = vel_CH[i].z;
		it++;

   }
 
	
	return 0;
}

int set_state(int index_of_the_particle, double mass, double radius, double x, double y, double z, double vx, double vy, double vz){
	dm_states[index_of_the_particle].mass = mass;
   dm_states[index_of_the_particle].x = x;
	dm_states[index_of_the_particle].y = y;
   dm_states[index_of_the_particle].z = z;
   dm_states[index_of_the_particle].vx = vx;
   dm_states[index_of_the_particle].vy = vy;
   dm_states[index_of_the_particle].vz = vz;
   return 0;
}

int get_state(int index_of_the_particle, double * mass, double * radius, double * x, double * y, double * z, double * vx, double * vy, double * vz){
	*mass = pos_CH[index_of_the_particle].w;
   *radius = 0.0;
   *x = pos_CH[index_of_the_particle].x;
   *y = pos_CH[index_of_the_particle].y;
   *z = pos_CH[index_of_the_particle].z;
   *vx = vel_CH[index_of_the_particle].x;
   *vy = vel_CH[index_of_the_particle].y;
   *vz = vel_CH[index_of_the_particle].z;
   return 0;
}

int get_time_step(double * time_step){
   if(rank == 0) cout<<"time steps are printed in the file 'Blocks', see README file to obtain more information"<<endl;  
	return 0;
}

int recommit_particles(){
	
   delete [] pos_PH;
	delete [] vel_PH;
	delete [] pos_CH;
	delete [] vel_CH;
	delete [] a_H0;
	delete [] step;
	delete [] local_time;

	N = particle_id_counter;

   pos_PH = new double4 [N];
   vel_PH = new  float4 [N];
   pos_CH = new double4 [N];
   vel_CH = new double4 [N];
   a_H0 = new double4 [3*N];
   step = new double [N];
   local_time = new double [N];

	it = dm_states.begin();
	
	for (int i=0; i<N; i++){
		pos_PH[i].w = it->second.mass;
      pos_PH[i].x = it->second.x;
      pos_PH[i].y = it->second.y;
      pos_PH[i].z = it->second.z;
      vel_CH[i].x = it->second.vx;
      vel_CH[i].y = it->second.vy;
      vel_CH[i].z = it->second.vz;
      vel_CH[i].w = 0.0;
      pos_CH[i].w = it->second.mass;
		pos_CH[i].x = it->second.x;
		pos_CH[i].y = it->second.y;
		pos_CH[i].z = it->second.z;
		vel_PH[i].x = vel_CH[i].x;
      vel_PH[i].y = vel_CH[i].y;
      vel_PH[i].z = vel_CH[i].z;
		it++;
		
	}
   
	GTW = GTIME;
	GTIME = 0.0;
	ATIME = 1.0e+10;
	warm_start = 1;

	HostSafeCall(InitBlocks(pos_PH, vel_PH, TPB, N, BFMAX, ETA4, NGPU, EPS, &MAXDIM, DTMAX, &GPUMINTHREADS, devices, GPUNAME, rank, size, pos_CH, vel_CH, a_H0, step, local_time, &ATIME, plummer_core, plummer_mass));
   return 0;
}

int get_kinetic_energy(double * kinetic_energy){
	*kinetic_energy = 0.0;
   double K = 0.0;
   unsigned int ppG = N / size;
   for(int i=0; i<ppG; i++)   
		K += 0.5 * pos_CH[i].w * (vel_CH[i].x * vel_CH[i].x + vel_CH[i].y * vel_CH[i].y + vel_CH[i].z * vel_CH[i].z);
	MPISafeCall(MPI_Barrier(MPI_COMM_WORLD));
   MPISafeCall(MPI_Allreduce(&K, kinetic_energy, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));
	return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, double az){
   if (rank == 0) cout<<"accelerations are stored only on GPU, you can't set it"<<endl;
   return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z){
   double m = 0.0;
   *x = *y = *z = 0.0;
   for(int i=0; i<N; i++){
		*x += pos_CH[i].w * pos_CH[i].x;
	   *y += pos_CH[i].w * pos_CH[i].y;
	   *z += pos_CH[i].w * pos_CH[i].z;
	   m += pos_CH[i].w;
   }
   *x /= m;
   *y /= m;
   *z /= m;
   return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){ 
	double m = 0.0;
   *vx = *vy = *vz = 0.0;
   for(int i=0; i<N; i++){
		*vx += pos_CH[i].w * vel_CH[i].x;
      *vy += pos_CH[i].w * vel_CH[i].y;
      *vz += pos_CH[i].w * vel_CH[i].z;
      m += pos_CH[i].w;
   }
   *vx /= m;
   *vy /= m;
   *vz /= m;
   return 0;
}

int get_radius(int index_of_the_particle, double * radius){
   *radius = 0.0;
   return 0;
}

int set_radius(int index_of_the_particle, double radius){
   return 0;
}

int cleanup_code(){
   delete [] pos_PH;
   delete [] vel_PH;
   delete [] pos_CH;
   delete [] vel_CH;
   delete [] a_H0;
   delete [] step;
   delete [] local_time;
	dm_states.clear();
	return 0;
}

int get_potential_energy(double * potential_energy){
   unsigned int ppG = N/(NGPU*size); 
   HostSafeCall(CudaInit(&GPUMINTHREADS, NGPU, rank, devices, GPUNAME));
   HostSafeCall(Calculate_potential_Energy(pos_CH, N, EPS, TPB, NGPU, rank, devices, ppG, potential_energy, plummer_core, plummer_mass));
	return 0;
}

int get_gravity_at_point(double eps, double x, double y, double z, double * forcex, double * forcey, double * forcez){
   *forcex = *forcey = *forcez = 0.0;
   for(int i=0; i<N ; i++){
		double rx = (pos_CH[i].x - x) * (pos_CH[i].x - x);
      double ry = (pos_CH[i].y - y) * (pos_CH[i].y - y);
      double rz = (pos_CH[i].z - z) * (pos_CH[i].z - z);
      double distance = rx * rx + ry * ry + rz * rz + eps * eps; 
      *forcex += pos_CH[i].w * rx / sqrt(distance * distance * distance);
      *forcey += pos_CH[i].w * ry / sqrt(distance * distance * distance);
      *forcez += pos_CH[i].w * rz / sqrt(distance * distance * distance);
   }
   return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, double * vz){ 
   *vx = vel_CH[index_of_the_particle].x;
   *vy = vel_CH[index_of_the_particle].y;
   *vz = vel_CH[index_of_the_particle].z;
   return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, double * z){  
   *x = pos_CH[index_of_the_particle].x;
   *y = pos_CH[index_of_the_particle].y;
   *z = pos_CH[index_of_the_particle].z;
   return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z){
   dm_states[index_of_the_particle].x = x;
   dm_states[index_of_the_particle].y = y;
   dm_states[index_of_the_particle].z = z;
   return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, double * az){
   if( rank == 0) cout<<"accelerations are stored only on GPU, you can't get it"<<endl;
   return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, double vz){
   dm_states[index_of_the_particle].vx = vx;
   dm_states[index_of_the_particle].vy = vy;
   dm_states[index_of_the_particle].vz = vz;
   return 0;
}



