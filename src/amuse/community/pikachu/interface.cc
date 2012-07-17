//~#include <fstream>
//~#include <cmath>
//~//#include <unistd.h>
#include <map>
#include <cstring>
//~
//~#ifdef SAPPORO_GRAPE
//~#include "sapporo.h"
//~#include "../lib/g6/g6lib.h"
//~#else
//~//#include"../lib/g6/g6_dummy.h"
//~#endif //GRAPE6


#include <iostream>

#include "BHtree.h"
#include "particle.h"

#include "const.h"
#include "distribution.h"
#include "system.h"
#include "mpi_interface.h"

#define GLOBAL_VALUE_DEFINE
#include "global.h"

#include "interface.h"
#include "worker_code.h"
// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>

using namespace std;


// Globals
static double begin_time = 0.0;

map<int, dynamics_state> particle_buffer;        // for initialization only
map<int, int> local_index_map;
map<int, int> reverse_index_map;
int particle_id_counter = 0;
bool particles_initialized = false;
bool debug = false;
bool energy_up_to_date = false;

Nbody_System *nbody_system = 0;
char output_dir[STRINGSIZE];
char kernel_dir[STRINGSIZE];
char kernel_dir_bogus[STRINGSIZE];
char *fake_argv[1];
double Ekin, Epot, Etot;
double Egr = 0.0;
double Emerge = 0.0;
double eps2_FS_FS, eps2_FS_BH, eps2_BH_BH;
double rcut_out_FS_FS, rcut_out_FS_BH, rcut_out_BH_BH;
double rsearch_FS_FS, rsearch_FS_BH, rsearch_BH_BH;
double theta2;
int Ncrit, Nleaf, quad_flag;
double eta_s, eta_FS, eta_BH, mass_min;
double search_factor, vel_disp;

// Interface functions

void halt_program(){
    dev_close();
    int error = 0;
    MPI_Abort(MPI_COMM_WORLD, error);
}

void close_program(){
    dev_close();
    MPI_Finalize();
}

int get_kernel_directory(char **kernel_directory){
    *kernel_directory = kernel_dir;
    return 0;
}
int set_kernel_directory(char *kernel_directory){
    int length = strlen(kernel_directory);
    if (length >= STRINGSIZE) {
        return -1;
    }
    strcpy(kernel_dir, kernel_directory);
    strcpy(kernel_dir_bogus, kernel_directory);
    strcat(kernel_dir_bogus, "/bogus");
    return 0;
}

void set_default_parameters() {
    strcpy(output_dir, "./");
    strcpy(kernel_dir, "./");
    strcpy(kernel_dir_bogus, "./bogus");
    fake_argv[0] = kernel_dir_bogus;
    
    NBH_GLB = NBH_LOC = 0;
    NDEAD_GLB = 0;
    NBH_GLB_ORG = NBH_GLB + NDEAD_GLB;
    
    nbody_system->dt_glb = 1.0 / 2048;
    nbody_system->Tsys = 0.0;
    nbody_system->Tend = 0.0;
    nbody_system->snp_id = 0;
    
    rcut_out_FS_FS = 2.0e-3;
    rcut_out_FS_BH = 2.0e-2;
    rcut_out_BH_BH = 1.0e5;
    
    vel_disp = 0.707106781;
    search_factor = 3.0;
    
    rsearch_FS_FS = rcut_out_FS_FS + search_factor * vel_disp * nbody_system->dt_glb;
    rsearch_FS_BH = rcut_out_FS_BH + search_factor * vel_disp * nbody_system->dt_glb;
    rsearch_BH_BH = rcut_out_BH_BH + search_factor * vel_disp * nbody_system->dt_glb;
    
    eps2_FS_FS = 1.0e-8;
    eps2_FS_BH = 1.0e-8;
    eps2_BH_BH = 0.0;
    
    eta_s = 0.005;
    eta_FS = 0.025;
    eta_BH = 0.025;
    theta2 = 0.16;
//~
//~####################################
//~##    sqrt(MFS/eps) * dt = r_cut  ##
//~####################################
    
    quad_flag = 0; // Nathan: "This was missing in the example parameter file, therefore it was used unitialized!"
    Ncrit = 8; // Nathan: "Idem, but I don't think this is actually used anywhere..."
    Nleaf = 8; // Nathan: "Idem, but I don't think this is actually used anywhere..."
}

int initialize_code() {
    //~fout_debug.open("debug.dat");
    cout << setprecision(15);
    cerr << setprecision(15);
    dump_cerr(sizeof(Nbody_System));
    dump_cerr(sizeof(Hard_System));
    dump_cerr(sizeof(Soft_System));
    dump_cerr(sizeof(Particle));
    dump_cerr(sizeof(Particle_Short));
    dump_cerr(sizeof(Cell_Tree));

    MPI_Comm_size(MPI_COMM_WORLD, &NPROC);
    MPI_Comm_rank(MPI_COMM_WORLD, &MYRANK);
    MPI_Barrier(MPI_COMM_WORLD);
    std::cerr << "MPI Initialize: myrank = " << MYRANK
        << "                Nproc = " << NPROC <<std::endl;

    nbody_system = new Nbody_System;
    
    //~set_support_for_condition(COLLISION_DETECTION);
    set_default_parameters();
    return 0;
}

int new_particle(int *particle_identifier, double mass, 
        double x, double y, double z, double vx, double vy, double vz, double radius) {
    dynamics_state new_p;
    
    // Particle is stored in the buffer until (re)commit_particles is called
    particle_id_counter++;
    *particle_identifier = particle_id_counter;
    new_p.mass = mass;
    new_p.radius = radius;
    new_p.x = x;
    new_p.y = y;
    new_p.z = z;
    new_p.vx = vx;
    new_p.vy = vy;
    new_p.vz = vz;
    particle_buffer.insert(pair<int, dynamics_state>(*particle_identifier, new_p));
    return 0;
}

void erase_particle(Particle *part){
    part->have_ngh = 0;
    part->index = 0;
    part->mass = 0.0;
    part->radius = 0.0;
    part->pot = 0.0;
    part->pot_short = 0.0;
    part->pos = 0.0;
    part->vel = 0.0;
    part->acc = 0.0;
    part->acc_long = 0.0;
    part->acc_short = 0.0;
    part->type = star;
    part->idx_ngh_FS = -1;
    part->r2_ngh_FS = -0.0;
    part->idx_ngh_BH = -1;
    part->r2_ngh_BH = -0.0;
    part->Nj = 0;
    part->node_org = -1;
    part->prt_next = NULL;
}

int delete_particle(int particle_identifier) {
    map<int, int>::iterator particle_iter = local_index_map.find(particle_identifier);
    if (particle_iter != local_index_map.end()){
        int index = particle_iter->second;
        local_index_map.erase(particle_iter);
        particle_iter = reverse_index_map.find(index);
        reverse_index_map.erase(particle_iter);
        erase_particle(&(nbody_system->prt_loc[index]));
        return 0;
    }
    
    map<int, dynamics_state>::iterator buffer_iter = particle_buffer.find(particle_identifier);
    if (buffer_iter != particle_buffer.end()){
        particle_buffer.erase(buffer_iter);
        return 0;
    }
    return -3; // Not found!
}


int commit_particles() {
    NFS_LOC = NFS_GLB = particle_buffer.size();
    NFS_GLB_ORG = NFS_GLB + NDEAD_GLB;
    
    int i = 0;
    local_index_map.clear();
    reverse_index_map.clear();
    for (map<int, dynamics_state>::iterator iter = particle_buffer.begin();
            iter != particle_buffer.end(); iter++, i++){
        local_index_map.insert(pair<int, int>((*iter).first, i)); // identifier -> index
        reverse_index_map.insert(pair<int, int>(i, (*iter).first)); // index -> identifier
        nbody_system->prt_loc[i].mass = (*iter).second.mass;
        nbody_system->prt_loc[i].radius = (*iter).second.radius;
        nbody_system->prt_loc[i].pos = Vector3((*iter).second.x, (*iter).second.y, (*iter).second.z);
        nbody_system->prt_loc[i].vel = Vector3((*iter).second.vx, (*iter).second.vy, (*iter).second.vz);
        nbody_system->prt_loc[i].index = i;
        nbody_system->prt_loc[i].type = star;
    }
    particle_buffer.clear();
    
    double mass_min_loc = nbody_system->prt_loc[0].mass;
    for (int i = 1; i < NFS_LOC + NBH_LOC; i++){
        if (nbody_system->prt_loc[i].mass < mass_min_loc){
            mass_min_loc = nbody_system->prt_loc[i].mass;
        }
    }
    mass_min = mpi_min_double(&mass_min_loc);
    nbody_system->hard_system.set(eps2_FS_FS, eps2_FS_BH, eps2_BH_BH,
        rcut_out_FS_FS, rcut_out_FS_BH, rcut_out_BH_BH,
        eta_s, eta_FS, eta_BH, mass_min);
    nbody_system->hard_system.dump();
    
    nbody_system->initialize_division();
    
    // assign particles (prt_loc and BH_glb)
    // determine NBH_LOC and NFS_LOC
    // sum of NBH_LOC + NFS_LOC over all nodes is NFS_GLB+NBH_GLB
    nbody_system->divide_particles();
    cerr<<"finish nbody_system->divide_particles()"<<endl;
    
    nbody_system->calc_soft_forces(fake_argv);
    //nbody_system->dump_tree_strcture(); // for the check of the tree strcture using jeroen's code
    cerr<<"finish nbody_system->calc_soft_forces();"<<endl;
    
    nbody_system->calc_energy(0);
    
    //~nbody_system->write_file(output_dir);
    
    particles_initialized = true;
    return 0;
}

void push_particle_data_back_to_buffer(){
    map<int, int>::iterator iter;
    int i;
    for (iter = local_index_map.begin(); iter != local_index_map.end(); iter++){
        i = iter->second;
        dynamics_state state;
        state.mass = nbody_system->prt_loc[i].mass;
        state.radius = nbody_system->prt_loc[i].radius;
        state.x = nbody_system->prt_loc[i].pos[0];
        state.y = nbody_system->prt_loc[i].pos[1];
        state.z = nbody_system->prt_loc[i].pos[2];
        state.vx = nbody_system->prt_loc[i].vel[0];
        state.vy = nbody_system->prt_loc[i].vel[1];
        state.vz = nbody_system->prt_loc[i].vel[2];
        particle_buffer.insert(pair<int, dynamics_state>(iter->first, state));
        erase_particle(&(nbody_system->prt_loc[i]));
    }
    local_index_map.clear();
    reverse_index_map.clear();
}

int recommit_particles() {
    if (particles_initialized) {
        push_particle_data_back_to_buffer();
        particles_initialized = false;
    }
    return commit_particles();
}

int cleanup_code() {
    if (particles_initialized){
        map<int, int>::iterator iter;
        for (iter = local_index_map.begin(); iter != local_index_map.end(); iter++){
            erase_particle(&(nbody_system->prt_loc[iter->second]));
        }
    }
    local_index_map.clear();
    reverse_index_map.clear();
    particle_buffer.clear();
    particle_id_counter = 0;
    return 0;
}

int evolve_model(double time) {
    // AMUSE STOPPING CONDITIONS SUPPORT
    //~int is_collision_detection_enabled = 0;
    //~is_stopping_condition_enabled(COLLISION_DETECTION, &is_collision_detection_enabled);
    //~reset_stopping_conditions();
    
    int after_stellar_evolution = 0;
    nbody_system->Tend = time;
    long long int loop_end = (long long int) ((time - nbody_system->Tsys) / nbody_system->dt_glb) + 1;
    for(long long int loop = 0; loop < loop_end; loop++){
        nbody_system->evolve(output_dir, fake_argv, after_stellar_evolution);
    }
    energy_up_to_date = false;
    return 0;
}

bool found_particle(int particle_identifier, int *index){
    if (particles_initialized) {
        map<int, int>::iterator iter = local_index_map.find(particle_identifier);
        if (iter != local_index_map.end()){
            *index = (*iter).second;
            return true;
        }
    }
    return false;
}
void get_identifier_of_particle_with_index(int index, int *particle_identifier){
    map<int, int>::iterator iter = reverse_index_map.find(index);
    if (iter == reverse_index_map.end()){
        cerr << "Error: Could not determine the identifier for particle at index: " << index << endl;
        *particle_identifier = -1;
    } else {
        *particle_identifier = (*iter).second;
    }
}

int get_index_of_first_particle(int *particle_identifier) {
    return -2; // Not implemented
}
int get_index_of_next_particle(int particle_identifier, int *next_particle_identifier) {
    return -2; // Not implemented
}


// simulation property getters:
int get_total_mass(double *total_mass) {
    // calculate only on the root mpi process, not on others
    if (MYRANK == 0) {
        *total_mass = 0.0;
        map<int, int>::iterator iter;
        for (iter = local_index_map.begin(); iter != local_index_map.end(); iter++){
            *total_mass += nbody_system->prt_loc[iter->second].mass;
        }
    }
    return 0;
}
int get_total_radius(double *total_radius) {
    return -2; // Not implemented
}
int get_time(double *time) {
    *time = nbody_system->Tsys;
    return 0;
}

int set_begin_time(double input) {
    begin_time = input;
    return 0;
}

int get_begin_time(double * output) {
    *output = begin_time;
    return 0;
}

int get_center_of_mass_position(double *x, double *y, double *z){
    // calculate only on the root mpi process, not on others
    if (MYRANK == 0) {
        *x = *y = *z = 0.0;
        double total_mass = 0.0;
        map<int, int>::iterator iter;
        int i;
        for (iter = local_index_map.begin(); iter != local_index_map.end(); iter++){
            i = iter->second;
            total_mass += nbody_system->prt_loc[i].mass;
            *x += nbody_system->prt_loc[i].mass * nbody_system->prt_loc[i].pos[0];
            *y += nbody_system->prt_loc[i].mass * nbody_system->prt_loc[i].pos[1];
            *z += nbody_system->prt_loc[i].mass * nbody_system->prt_loc[i].pos[2];
        }
        *x /= total_mass;
        *y /= total_mass;
        *z /= total_mass;
    }
    return 0;
}
int get_center_of_mass_velocity(double *vx, double *vy, double *vz) {
    // calculate only on the root mpi process, not on others
    if (MYRANK == 0) {
        *vx = *vy = *vz = 0.0;
        double total_mass = 0.0;
        map<int, int>::iterator iter;
        int i;
        for (iter = local_index_map.begin(); iter != local_index_map.end(); iter++){
            i = iter->second;
            total_mass += nbody_system->prt_loc[i].mass;
            *vx += nbody_system->prt_loc[i].mass * nbody_system->prt_loc[i].vel[0];
            *vy += nbody_system->prt_loc[i].mass * nbody_system->prt_loc[i].vel[1];
            *vz += nbody_system->prt_loc[i].mass * nbody_system->prt_loc[i].vel[2];
        }
        *vx /= total_mass;
        *vy /= total_mass;
        *vz /= total_mass;
    }
    return 0;
}
int get_kinetic_energy(double *kinetic_energy) {
    if (!energy_up_to_date) {
        nbody_system->calc_energy(Ekin, Epot, Etot);
        energy_up_to_date = true;
    }
    *kinetic_energy = Ekin;
    return 0;
}
int get_potential_energy(double *potential_energy) {
    if (!energy_up_to_date) {
        nbody_system->calc_energy(Ekin, Epot, Etot);
        energy_up_to_date = true;
    }
    *potential_energy = Epot;
    return 0;
}
int get_number_of_particles(int *number_of_particles) {
    *number_of_particles = NFS_GLB + NBH_GLB;
    return 0;
}




// particle property getters/setters: (will only work after commit_particles() is called)
int set_mass(int particle_identifier, double mass) {
    int index;
    if (found_particle(particle_identifier, &index)){
        nbody_system->prt_loc[index].mass = mass;
        return 0;
    }
    return -3; // Not found!
}
int get_mass(int particle_identifier, double *mass) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *mass = nbody_system->prt_loc[index].mass;
        return 0;
    }
    return -3; // Not found!
}
int set_radius(int particle_identifier, double radius) {
    int index;
    if (found_particle(particle_identifier, &index)){
        nbody_system->prt_loc[index].radius = radius;
        return 0;
    }
    return -3; // Not found!
}
int get_radius(int particle_identifier, double * radius) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *radius = nbody_system->prt_loc[index].radius;
        return 0;
    }
    return -3; // Not found!
}
int set_position(int particle_identifier, double x, double y, double z) {
    int index;
    if (found_particle(particle_identifier, &index)){
        nbody_system->prt_loc[index].pos = Vector3(x, y, z);
        return 0;
    }
    return -3; // Not found!
}
int get_position(int particle_identifier, double *x, double *y, double *z) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *x = nbody_system->prt_loc[index].pos[0];
        *y = nbody_system->prt_loc[index].pos[1];
        *z = nbody_system->prt_loc[index].pos[2];
        return 0;
    }
    return -3; // Not found!
}
int set_velocity(int particle_identifier, double vx, double vy, double vz) {
    int index;
    if (found_particle(particle_identifier, &index)){
        nbody_system->prt_loc[index].vel = Vector3(vx, vy, vz);
        return 0;
    }
    return -3; // Not found!
}
int get_velocity(int particle_identifier, double *vx, double *vy, double *vz) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *vx = nbody_system->prt_loc[index].vel[0];
        *vy = nbody_system->prt_loc[index].vel[1];
        *vz = nbody_system->prt_loc[index].vel[2];
        return 0;
    }
    return -3; // Not found!
}
int set_state(int particle_identifier, double mass, 
        double x, double y, double z, 
        double vx, double vy, double vz, double radius) {
    int index;
    if (found_particle(particle_identifier, &index)){
        nbody_system->prt_loc[index].mass = mass;
        nbody_system->prt_loc[index].radius = radius;
        nbody_system->prt_loc[index].pos = Vector3(x, y, z);
        nbody_system->prt_loc[index].vel = Vector3(vx, vy, vz);
        return 0;
    }
    return -3; // Not found!
}
int get_state(int particle_identifier, double *mass, 
        double *x, double *y, double *z,
        double *vx, double *vy, double *vz, double *radius) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *mass = nbody_system->prt_loc[index].mass;
        *radius = nbody_system->prt_loc[index].radius;
        *x = nbody_system->prt_loc[index].pos[0];
        *y = nbody_system->prt_loc[index].pos[1];
        *z = nbody_system->prt_loc[index].pos[2];
        *vx = nbody_system->prt_loc[index].vel[0];
        *vy = nbody_system->prt_loc[index].vel[1];
        *vz = nbody_system->prt_loc[index].vel[2];
        return 0;
    }
    return -3; // Not found!
}

int set_acceleration(int particle_identifier, double ax, double ay, double az) {
    return -2; // Not implemented
}

int get_acceleration(int particle_identifier, double *ax, double *ay, double *az) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *ax = nbody_system->prt_loc[index].acc[0];
        *ay = nbody_system->prt_loc[index].acc[1];
        *az = nbody_system->prt_loc[index].acc[2];
        return 0;
    }
    return -3; // Not found!
}
int get_potential(double x, double y, double z, double *V){
    return -2; // Not implemented
}




// parameter getters/setters:
int set_eps2(double epsilon_squared) {
    return -2;
}
int get_eps2(double *epsilon_squared) {
    return -2;
}
int commit_parameters() {
    nbody_system->Tsys = begin_time;
    
    if(rsearch_FS_FS < rcut_out_FS_FS || rsearch_FS_BH < rcut_out_FS_BH || rsearch_BH_BH < rcut_out_BH_BH){
        return -1;
    }
    
    rsearch_FS_FS = rcut_out_FS_FS + search_factor * vel_disp * nbody_system->dt_glb;
    rsearch_FS_BH = rcut_out_FS_BH + search_factor * vel_disp * nbody_system->dt_glb;
    rsearch_BH_BH = rcut_out_BH_BH + search_factor * vel_disp * nbody_system->dt_glb;
    
    nbody_system->soft_system.set(eps2_FS_FS, eps2_FS_BH, eps2_BH_BH,
        rcut_out_FS_FS, rcut_out_FS_BH, rcut_out_BH_BH,
        rsearch_FS_FS, rsearch_FS_BH, rsearch_BH_BH,
        theta2, Ncrit, Nleaf, quad_flag);
    nbody_system->soft_system.dump();
    return 0;
}
int recommit_parameters() {
    return commit_parameters();
}


int get_potential_at_point(double *eps, double *x, double *y, double *z, double *phi, int length) {
    // Create "ghost"-particles to measure potential at their locations
    int* tmp_index = new int[length];
    double* tmp_mass_in = new double[length];
    double* tmp_eps2_in = new double[length];
    double (*tmp_pos_in)[3] = new double[length][3];
    double (*tmp_vel_in)[3] = new double[length][3];
    double (*tmp_acc_in)[3] = new double[length][3];
    double (*tmp_acc_out)[3] = new double[length][3];
    double (*tmp_jrk_out)[3] = new double[length][3];
    double (*tmp_snp_out)[3] = new double[length][3];
    double (*tmp_crk_out)[3] = new double[length][3];
    double *tmp_phi_out = new double[length];
    int *tmp_nnb_out = new int[length];
    double *tmp_nnb_r2_out = new double[length];
    
    for (int i=0; i<length; i++) {
        // Make sure there is no particle in the code with the same index, 
        // since the code would think it is calculating force on itself, and ignore it.
        tmp_index[i] = -1;
        tmp_mass_in[i] = 0.0;
        //~tmp_eps2_in[i] = eps[i]*eps[i] + eps2_fs_fs;
        tmp_pos_in[i][0] = x[i];
        tmp_pos_in[i][1] = y[i];
        tmp_pos_in[i][2] = z[i];
        
        for (int j=0; j<3; j++) {
            tmp_vel_in[i][j] = 0.0;
            tmp_acc_in[i][j] = 0.0;
            tmp_acc_out[i][j] = 0.0;
            tmp_jrk_out[i][j] = 0.0;
            tmp_snp_out[i][j] = 0.0;
            tmp_crk_out[i][j] = 0.0;
        }
        tmp_phi_out[i] = 0.0;
        tmp_nnb_out[i] = 0;
        tmp_nnb_r2_out[i] = 0.0;
    }
    //~calc_force_on_predictors(length, Njp,
        //~tmp_index, tmp_pos_in, 
        //~tmp_vel_in, tmp_acc_in, 
        //~tmp_mass_in, tmp_eps2_in,
        //~tmp_acc_out, tmp_jrk_out,
        //~tmp_snp_out, tmp_crk_out, tmp_phi_out,
        //~tmp_nnb_out, tmp_nnb_r2_out);
    
    MPI_Allreduce(tmp_phi_out, phi,
        length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    delete[] tmp_index;
    delete[] tmp_mass_in;
    delete[] tmp_eps2_in;
    delete[] tmp_pos_in;
    delete[] tmp_vel_in;
    delete[] tmp_acc_in;
    delete[] tmp_acc_out;
    delete[] tmp_jrk_out;
    delete[] tmp_snp_out;
    delete[] tmp_crk_out;
    delete[] tmp_phi_out;
    delete[] tmp_nnb_out;
    delete[] tmp_nnb_r2_out;
    return 0;
}
int get_gravity_at_point(double *eps, double *x, double *y, double *z, 
        double *forcex, double *forcey, double *forcez, int length){
    // Create "ghost"-particles to measure gravity at their locations
    int* tmp_index = new int[length];
    double* tmp_mass_in = new double[length];
    double* tmp_eps2_in = new double[length];
    double (*tmp_pos_in)[3] = new double[length][3];
    double (*tmp_vel_in)[3] = new double[length][3];
    double (*acc)[3] = new double[length][3];
    double (*tmp_acc_out)[3] = new double[length][3];
    double (*tmp_jrk_out)[3] = new double[length][3];
    double (*tmp_snp_out)[3] = new double[length][3];
    double (*tmp_crk_out)[3] = new double[length][3];
    double *tmp_phi_out = new double[length];
    int *tmp_nnb_out = new int[length];
    double *tmp_nnb_r2_out = new double[length];
    
    for (int i=0; i<length; i++) {
        // Make sure there is no particle in the code with the same index, 
        // since the code would think it is calculating force on itself, and ignore it.
        tmp_index[i] = -1;
        tmp_mass_in[i] = 0.0;
        //~tmp_eps2_in[i] = eps[i]*eps[i] + eps2_fs_fs;
        tmp_pos_in[i][0] = x[i];
        tmp_pos_in[i][1] = y[i];
        tmp_pos_in[i][2] = z[i];
        
        for (int j=0; j<3; j++) {
            tmp_vel_in[i][j] = 0.0;
            acc[i][j] = 0.0;
            tmp_acc_out[i][j] = 0.0;
            tmp_jrk_out[i][j] = 0.0;
            tmp_snp_out[i][j] = 0.0;
            tmp_crk_out[i][j] = 0.0;
        }
        tmp_phi_out[i] = 0.0;
        tmp_nnb_out[i] = 0;
        tmp_nnb_r2_out[i] = 0.0;
    }
    //~calc_force_on_predictors(length, Njp,
        //~tmp_index, tmp_pos_in, 
        //~tmp_vel_in, acc, 
        //~tmp_mass_in, tmp_eps2_in,
        //~tmp_acc_out, tmp_jrk_out,
        //~tmp_snp_out, tmp_crk_out, tmp_phi_out,
        //~tmp_nnb_out, tmp_nnb_r2_out);
    
    MPI_Allreduce(tmp_acc_out, acc,
        3*length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    delete[] tmp_index;
    delete[] tmp_mass_in;
    delete[] tmp_eps2_in;
    delete[] tmp_pos_in;
    delete[] tmp_vel_in;
    delete[] tmp_acc_out;
    delete[] tmp_jrk_out;
    delete[] tmp_snp_out;
    delete[] tmp_crk_out;
    delete[] tmp_phi_out;
    delete[] tmp_nnb_out;
    delete[] tmp_nnb_r2_out;
    
    for (int i=0; i<length; i++) {
        forcex[i] = acc[i][0];
        forcey[i] = acc[i][1];
        forcez[i] = acc[i][2];
    }
    delete[] acc;
    return 0;
}


int get_time_step(double *time_step) {
    *time_step = -1;
    return 0; // Not implemented
}


int get_potential(int id, double *phi)
{
  return -1;
}

int synchronize_model()
{
  return 0;
}




