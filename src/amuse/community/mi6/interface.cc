#include <iostream>
#include <fstream>
#include <cmath>
#include <unistd.h>
//~#include <vector>
#include <map>

#ifdef SAPPORO_GRAPE
#include "sapporo.h"
#include "../lib/g6/g6lib.h"
#else
//#include"../lib/g6/g6_dummy.h"
#endif //GRAPE6

//#include "src/IO.h"

//#include<mpi_interface.h>
//#include"drive.h"

#include "interface.h"
#include "worker_code.h"
//#include "../lib/stopcond/stopcond.h"
#include "Vector3.h"
#include "Particle.h"
#include "evolve.h"
#include "energy.h"

using namespace std;



// Globals
double eps2_fs_fs = 0.0;
double eps2_fs_smbh = 0.0;
double eps2_bh_bh = 0.0;
double eta_s = 1.0e-4;
double eta_fs = 0.1;
double eta_smbh = 0.4;
double eta_imbh = 0.4;
double Tsys = 0.0;
double Tmerge = 0.0;
//~double Tend = 0.0;
double Egr = 0.0;
double Emerge = 0.0;
int myrank;
int Nproc;
map<int, dynamics_state> particle_buffer;        // for initialization only
map<int, dynamics_state> black_hole_buffer;      // for initialization only
map<int, int> local_index_map;
int particle_id_counter = 0;
bool particles_initialized = false;
int Ntot = 0;
int Nip_tot = 0; // Ntot - Ndead     (every node have same number)
int Nip = 0; // Nip  (every node have same number)
int NFS = 0;
int NBH = 0;
int Njp_org = 0;
int Njp = 0;
int first_address = 0;
int Ndead = 0; // Nmerge + Naccrete (every node have same number)
int NSMBH = 0;
int NIMBH = 0;
Particle *prt = NULL;
Particle *prt_old = NULL;
//const double dEcrit = 1e-10;
double dEcrit = 5e-5;
//const double dEcrit = 1e30;
int Nstep = 0;
int Nstep_old = 0;
int Nloop = 0;
int Nloop_old = 0;
int state;
//double dt_max = 1.0/8.0;
double dt_max = 1.0/1024.0;
//double dt_max = 1.0/16384.0;
//~double dt_snp;
//~double Tsnp;
int first_loop;
int first_loop_old;
int itr;
int Nmerge = 0;
int Nmerge_loop = 0;
int Naccrete = 0;
int Naccrete_loop = 0;
double Tcal_grav0 = 0.0;
double Tcal_grav1 = 0.0;
double Tcal_comm0 = 0.0;
double Tcal_comm1 = 0.0;
double Tcal_fix0 = 0.0;
double Tcal_fix1 = 0.0;
double Tcal_tot0 = 0.0;
double Tcal_tot1 = 0.0;
double Tcal_all = 0.0;
int EX_FLAG = 0;
int *address = NULL;
int *address_old = NULL;
//~Particle **prt_merged = NULL;
//~Particle **prt_accreted = NULL;
Particle *(prt_merged[10000]);
Particle *(prt_accreted[10000]);
double E0, Ek0, Ep0;
double E1, Ek1, Ep1;
double E1_old, Ek1_old, Ep1_old;
double Tmerge_old;
int Nip_tot_old, Njp_old;
int Ndead_old, Nmerge_old, Naccrete_old;
double Tsys_old;
double Egr_old;
double Emerge_old;


// Interface functions

int initialize_code() {
    cout << setprecision(15);
    cerr << setprecision(15);
    
    MPI_Comm_rank (MPI_COMM_WORLD, &myrank); 
    MPI_Comm_size (MPI_COMM_WORLD, &Nproc); 
    return 0;
}

int new_particle(int *particle_identifier, double mass, double radius, 
        double x, double y, double z, double vx, double vy, double vz) {
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

int new_black_hole(int *particle_identifier, double mass, double radius, 
        double x, double y, double z, double vx, double vy, double vz) {
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
    black_hole_buffer.insert(pair<int, dynamics_state>(*particle_identifier, new_p));
    return 0;
}

int delete_particle(int particle_identifier) {
    // First try deleting it from prt Particle array
    
    // ...
    
    // Then try deleting it from particle_buffer (maybe (re)commit_particles wasn't called yet)
    map<int, dynamics_state>::iterator buffer_iter;
    buffer_iter = particle_buffer.find(particle_identifier);
    if (buffer_iter != particle_buffer.end()){
        particle_buffer.erase(buffer_iter);
        return 0;
    }
    buffer_iter = black_hole_buffer.find(particle_identifier);
    if (buffer_iter != black_hole_buffer.end()){
        black_hole_buffer.erase(buffer_iter);
        return 0;
    }
    
    return -3; // Not found!
}

int commit_particles() {
    NFS = particle_buffer.size();
    NBH = black_hole_buffer.size();
    Ntot = NFS + NBH;
    Nip = Ntot; 
    prt = new Particle[Ntot];
    prt_old = new Particle[Ntot];
    cout << "a" << endl << flush;
    int i = 0;
    for (map<int, dynamics_state>::iterator iter = black_hole_buffer.begin();
            iter != black_hole_buffer.end(); iter++, i++){
        local_index_map.insert(pair<int, int>((*iter).first, i));
        prt[i].mass = (*iter).second.mass;
        prt[i].radius = (*iter).second.radius;
        prt[i].pos = Vector3((*iter).second.x, (*iter).second.y, (*iter).second.z);
        prt[i].vel = Vector3((*iter).second.vx, (*iter).second.vy, (*iter).second.vz);
        prt[i].index = i;
        prt[i].address = i;
        prt[i].type = SMBH;
    }
    black_hole_buffer.clear();
    cout << "a" << endl << flush;
    i = NBH; // Just to be sure...
    for (map<int, dynamics_state>::iterator iter = particle_buffer.begin();
            iter != particle_buffer.end(); iter++, i++){
        local_index_map.insert(pair<int, int>((*iter).first, i));
        prt[i].mass = (*iter).second.mass;
        prt[i].radius = (*iter).second.radius;
        prt[i].pos = Vector3((*iter).second.x, (*iter).second.y, (*iter).second.z);
        prt[i].vel = Vector3((*iter).second.vx, (*iter).second.vy, (*iter).second.vz);
        prt[i].index = i;
        prt[i].address = i;
    }
    particle_buffer.clear();
    cout << "a" << endl << flush;
    
    Nip_tot = Ntot - Ndead;
    divide_proc(Ntot, Njp_org, first_address);
    
    address = new int[Ntot];
    address_old = new int[Ntot];
    for(int i=0; i<Ntot; i++){
        address[i] = i;
    }
    cout << Ntot << " " << Njp_org << " " << first_address << " " << prt[0].address << " " << prt[1].address << endl;
    Njp = 0;
    for(int j=0; j<Ntot; j++){
        if(first_address <= j && j < first_address+Njp_org && prt[j].mass != 0.0){
            prt[j].address = Njp;
            Njp++;
        } else {
            prt[j].address = -1;
        }
    }
    
    cout << address[0] << " " << address[1] << endl;
    cout << prt[0].address << " " << prt[1].address << endl;
    cout << Ntot << " " << NBH << endl;
    cout << Njp << " " << Tsys << endl;
    prt[0].dump();
    prt[1].dump();
    evolve_initialize(prt, address, Ntot, NBH, Njp, Tsys);
    prt[0].dump();
    prt[1].dump();
    
    for(int i=0; i<Ntot; i++){
        prt_old[i] = prt[i];
        address_old[i] = address[i];
    }
    
    cout << address[0] << " " << address[1] << " " << Nip_tot << endl;
    E0, Ek0, Ep0;
    calc_energy(prt, address, Nip_tot, E0, Ek0, Ep0, 0);
    E1 = E0;
    Ek1 = Ek0;
    Ep1 = Ep0;
    E1_old = E0;
    Ek1_old = Ek0;
    Ep1_old = Ep0;
    Tmerge_old = Tmerge;
    Nip_tot_old = Nip_tot;
    Njp_old = Njp;
    
    Ndead_old = Ndead;
    Nmerge_old = Nmerge;
    Naccrete_old = Naccrete;
    
    Tsys_old = Tsys;
    Egr_old = Egr;
    Emerge_old = Emerge;
    
    copy_SMBH_NEW_TO_OLD();
    
    if(myrank == 0){
        cerr<<"E0="<<E0<<endl;
        cerr<<"Ek0="<<Ek0<<endl;
        cerr<<"Ep0="<<Ep0<<endl;
//        write0(prt, Ntot, NBH, Ndead, Tsys, Tmerge, Egr, dirname, snpid, EX_FLAG);
    }
    
    state = 0;
    //~dt_snp = 1.0/32.0;
    //~Tsnp = Tsys + dt_snp;
    first_loop = 1;
    first_loop_old = first_loop;
    //~Particle *(tmp_prt_merged[10000]);
    //~Particle *(tmp_prt_accreted[10000]);
    //~prt_merged = tmp_prt_merged;
    //~prt_accreted = tmp_prt_accreted;
    Tcal_tot0 = MPI_Wtime();
    itr = 0;
    
    particles_initialized = true;
    return 0;
}

int recommit_particles() {
    return -2;
}

int cleanup_code() {
    if (particles_initialized) {
        delete[] prt;
        delete[] prt_old;
        delete[] address;
        delete[] address_old;
    }
    return 0;
}

void iteration(Particle prt[], Particle prt_old[], int address[], int address_old[],
    const int &Ntot, const int &Ndead_old, const int &Nmerge_old, const int &Naccrete_old,
    const int &Nip_tot_old, const int &Njp_old, const double &E1_old, const double &Ep1_old,
    const double &Ek1_old, const double &Egr_old, const double &Emerge_old,
    const double &Tsys_old, const double &Tmerge_old, const int &Nloop_old,
    const int &Nstep_old, const int &first_loop_old, const double &eta_s,
    const double &eta_fs,const double &eta_smbh,const double &eta_imbh, int &Ndead,
    int &Nmerge, int &Naccrete, int &Nip_tot, int &Njp, double &E1, double &Ep1,
    double &Ek1, double &Egr, double &Emerge, double &Tsys, double &Tmerge,
    int &Nloop, int &Nstep, int &first_loop, int &itr, int &Nip){
    
    itr++;
    cerr<<"energy error is too large: iteration="<<itr<<endl;
    for (int i=0; i<Ntot; i++) {
        prt[i] = prt_old[i];
        address[i] = address_old[i];
    }
    
    Ndead = Ndead_old;
    Nip_tot = Nip_tot_old;
    Njp = Njp_old;
    
    set_NJP(Njp);
    setj_to_sapporo(prt, address, Ntot);
    //setj_to_sapporo(prt, address, Njp);
    
    E1 = E1_old;
    Ek1 = Ek1_old;
    Ep1 = Ep1_old;
    Egr = Egr_old;
    Emerge = Emerge_old;
    
    Tsys = Tsys_old;
    Tmerge = Tmerge_old;
    
    Nloop = Nloop_old;
    Nstep = Nstep_old;
    first_loop = first_loop_old;
    //set_eta(eta_s*pow(0.5, itr), eta_fs*pow(0.5, itr), eta_smbh*pow(0.5, itr), eta_imbh*pow(0.5, itr));
    Tcal_grav1 = Tcal_comm1 = Tcal_fix1 = 0.0;
    Tcal_tot0 = MPI_Wtime();
    set_NSTEP(Nstep);
    
    Nmerge = Nmerge_old;
    Naccrete = Naccrete_old;
    set_Ndead(Nmerge, Naccrete);
    Nip = Nip_tot;
    
    copy_SMBH_OLD_TO_NEW();
}

int evolve_model(double time) {
    while (Tsys < time) {
        if (Nloop % 1000 == 0) {
            cerr << "Tsys=" << Tsys << ",  Nloop=" << Nloop << ",  Nstep=" << get_NSTEP() << endl;
        }
        Nloop++;
        if (state == 1) {
            state = 0;
            first_loop = 1;
        }
        
        set_eta(eta_s*pow(0.5, itr), eta_fs*pow(0.5, itr), eta_smbh*pow(0.5, itr), eta_imbh*pow(0.5, itr));
        state = evolve_onestep(prt, address, Nip, Ntot, NBH, Tsys, Tmerge, dt_max, first_loop, Egr, itr);
        
        first_loop = 0;
        if (state == 1) {
            double E1_tmp = 0.0; 
            double Ek1_tmp = 0.0;
            double Ep1_tmp = 0.0;
            calc_energy(prt, Ntot, E1_tmp, Ek1_tmp, Ep1_tmp, 0);
            
            merge_prt();
            Tmerge = Tsys;
            Njp = 0;
            
            for (int j=first_address; j<first_address+Njp_org; j++) {
                if (prt[j].address == -1) {
                    continue;
                }
                prt[j].address = Njp;
                Njp++;
            }
            
            // do somthing to evolve merged stars using SSE
            get_merged_prt(prt_merged, Nmerge_loop);
            Nmerge += Nmerge_loop;
            get_accreted_prt(prt_accreted, Naccrete_loop);
            Naccrete += Naccrete_loop;
            Ndead += Nmerge_loop + Naccrete_loop;
            Nip_tot = Ntot - Ndead;
            
            sort_time_all(prt, address, Ntot);
            
            //~if (myrank == 0) {
                //~fout_merge<<"--------------------------"<<endl;
                //~fout_merge<<"Tsys="<<Tsys<<endl;
                //~fout_merge<<"Ndead="<<Ndead<<endl;
                //~fout_merge<<"merged particles"<<endl;
                //~fout_merge<<"Nmerge_loop="<<Nmerge_loop<<endl;
                //~fout_merge<<"Nmerge="<<Nmerge<<endl;
                //~for(int i=0; i<Nmerge_loop; i++){
                    //~fout_merge<<"index="<<prt_merged[i]->index<<endl;
                    //~fout_merge<<"mass="<<prt_merged[i]->mass<<endl;
                    //~fout_merge<<"pos="<<prt_merged[i]->pos<<endl;
                    //~fout_merge<<"vel="<<prt_merged[i]->vel<<endl;
                    //~fout_merge<<"phi="<<prt_merged[i]->phi<<endl;
                    //~fout_merge<<endl;
                //~}
                //~fout_merge<<endl;
                //~fout_merge<<"accreted particles"<<endl;
                //~fout_merge<<"Naccrete_loop="<<Naccrete_loop<<endl;
                //~fout_merge<<"Naccrete="<<Naccrete<<endl;
                //~for(int i=0; i<Naccrete_loop; i++){
                    //~fout_merge<<"index="<<prt_accreted[i]->index<<endl;
                    //~fout_merge<<"mass="<<prt_accreted[i]->mass<<endl;
                    //~fout_merge<<"pos="<<prt_accreted[i]->pos<<endl;
                    //~fout_merge<<"vel="<<prt_accreted[i]->vel<<endl;
                    //~fout_merge<<"phi="<<prt_accreted[i]->phi<<endl;
                    //~fout_merge<<endl;
                //~}
                //~fout_merge<<endl;
            //~}
            
            for (int i=0; i<Ntot; i++) {
                prt[i].clear();
                prt[i].acc4 = 0.0;
                prt[i].acc5 = 0.0;
            }
            evolve_initialize(prt,  address,  Ntot,  NBH,  Njp,  Tsys);
            calc_energy(prt, Ntot, E1, Ek1, Ep1, 0);
            
            // accumulate dissipation energy through merger
            Emerge += E1 - E1_tmp;
            Nip = Nip_tot;
            Tmerge = Tsys;
        }
        
        //~if(Tsnp <= Tsys){
            //~write0(prt, Ntot, NBH, Ndead, Tsys, Tmerge, Egr, dirname, snpid, EX_FLAG);
            //~Tsnp += dt_snp;
        //~}
        
        //if( fmod(Tsys-Tmerge, dt_max) == 0.0 || state == 1){
        if (fmod(Tsys-Tmerge, dt_max) == 0.0 && state == 0) {
            Tcal_tot1 = MPI_Wtime() - Tcal_tot0;
            //calc_energy(prt, address, Nip_tot, E1, Ek1, Ep1, 0);
            calc_energy(prt, Ntot, E1, Ek1, Ep1, 0);
            if(myrank == 0){
                cerr<<endl;
                cerr<<"Tsys="<<Tsys<<endl;
                cerr<<"E1="<<E1<<endl;
                cerr<<"E1_old="<<E1_old<<endl;
                cerr<<"Ek1="<<Ek1<<endl;
                cerr<<"Ek1_old="<<Ek1_old<<endl;
                cerr<<"Ep1="<<Ep1<<endl;
                cerr<<"Ep1_old="<<Ep1_old<<endl;
                cerr<<"Emerge="<<Emerge<<endl;
                cerr<<"Emerge_old="<<Emerge_old<<endl;
                cerr<<"Egr="<<Egr<<endl;
                cerr<<"Egr_old="<<Egr_old<<endl;
                cerr<<"(E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old/dt_max="<<(E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old/dt_max<<endl;
                cerr<<"(E1-E0-Egr-Emerge)/E0="<<(E1-E0-Egr-Emerge)/E0<<endl;
            }
            if (dEcrit < fabs( (E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old/dt_max ) ) {
                iteration(prt,  prt_old,  address,  address_old,  
                    Ntot,  Ndead_old,  Nmerge_old,  Naccrete_old,
                    Nip_tot_old,  Njp_old,
                    E1_old,  Ep1_old,  Ek1_old,  Egr_old,  Emerge_old,
                    Tsys_old,  Tmerge_old,  
                    Nloop_old,  Nstep_old,  first_loop_old,
                    eta_s, eta_fs, eta_smbh, eta_imbh, 
                    Ndead,   Nmerge,  Naccrete,
                    Nip_tot,  Njp,
                    E1,  Ep1,  Ek1,  Egr,  Emerge,
                    Tsys,  Tmerge,
                    Nloop,  Nstep,  first_loop,
                    itr, 
                    Nip);
                
                goto ITERATION;
            } else {
                if(myrank == 0){
                    cerr<<"no iteration"<<endl;
                    cerr<<endl;
                }
            }
            
            if (myrank == 0) {
                Nstep = get_NSTEP();
                cerr<<"Tsys="<<Tsys<<" ( = "<<Tsys*7.45e3<<"[yr] = "<<Tsys*7.45e3*365<<"[day]"<<endl;
                cerr<<"Tcal_tot="<<Tcal_tot1<<"[sec]"<<endl;
                cerr<<"Tcal_grav="<<Tcal_grav1<<"[sec]"<<endl;
                cerr<<"Tcal_fix="<<Tcal_fix1<<"[sec]"<<endl;
                cerr<<"Tcal_comm="<<Tcal_comm1<<"[sec]"<<endl;
                cerr<<"Nip_tot="<<Nip_tot<<endl;
                cerr<<"Ndead="<<Ndead<<endl;
                cerr<<"Nstep="<<Nstep<<endl;
                cerr<<"Nloop="<<Nloop<<endl;
                cerr<<((double)(Nstep-Nstep_old)*(double)(Nip_tot-Ndead)*97.0)*1e-9/Tcal_tot1<<"[Gflops]"<<endl;
                cerr<<"Nave="<<(Nstep-Nstep_old)/(Nloop-Nloop_old)<<endl;
                cerr<<"Egr="<<Egr<<endl;
                cerr<<"Emerge="<<Emerge<<endl;
                cerr<<"fabs((E1-E0-Egr-Emerge)/E0)="<<fabs((E1-E0-Egr-Emerge)/E0)<<endl;
                cerr<<"fabs((E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old)="<<fabs((E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old)<<endl;
                cerr<<endl;
                
                //~fout_log<<"Tsys="<<Tsys<<endl;
                //~fout_log<<"Tcal_tot1="<<Tcal_tot1<<endl;
                //~fout_log<<((double)(Nstep-Nstep_old)*(double)(Nip_tot-Ndead)*97.0)*1e-9/Tcal_tot1<<"[Gflops]"<<endl;
                //~fout_log<<"Nave="<<(Nstep-Nstep_old)/(Nloop-Nloop_old)<<endl;
                //~fout_log<<"Nloop="<<Nloop<<endl;
                //~fout_log<<"Egr="<<Egr<<endl;
                //~fout_log<<"Emerge="<<Emerge<<endl;
                //~fout_log<<"fabs((E1-E0-Egr-Emerge)/E0)="<<fabs((E1-E0-Egr-Emerge)/E0)<<endl;
                //~fout_log<<"fabs((E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E0)="<<fabs((E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old)<<endl;
                //~fout_log<<"get_NSTEP()="<<get_NSTEP()<<endl;
                //~fout_log<<endl;
                //~fout_BH<<Tsys;
                //~for(int i=0; i<NBH; i++){
                    //~fout_BH<<"   "<<prt[i].pos<<"   "<<prt[i].vel<<"   "<<prt[i].phi;
                //~}
                //~fout_BH<<endl;
            }
            
            if (state == 0) {
                for (int i=0; i<Ntot; i++) {
                    prt_old[i] = prt[i];
                    address_old[i] = address[i];
                }
                
                Ndead_old = Ndead;
                Nmerge_old = Nmerge;
                Naccrete_old = Naccrete;
                
                Nip_tot_old = Nip_tot;
                Njp_old = Njp;
                
                E1_old = E1;
                Ek1_old = Ek1;
                Ep1_old = Ep1;
                
                Egr_old = Egr;
                Emerge_old = Emerge;
                
                Tcal_tot0 = MPI_Wtime();
                Nstep_old = Nstep;
                
                
                Nloop_old = Nloop;
                first_loop_old = first_loop;
                Tcal_grav1 = Tcal_comm1 = Tcal_fix1 = 0.0;
                itr = 0;
                
                Tsys_old = Tsys;
                Tmerge_old = Tmerge;
                
                copy_SMBH_NEW_TO_OLD();
                
            }
        }
        ITERATION: ;
    }
    return 0;
}



bool found_particle(int particle_identifier, int *index){
    if (particles_initialized) {
        //~if (!particle_map_up_to_date)
            //~update_particle_map();
        
        map<int, int>::iterator iter = local_index_map.find(particle_identifier);
        if (iter != local_index_map.end()){
            *index = (*iter).second;
            return true;
        }
    }
    return false;
}

int get_index_of_first_particle(int *particle_identifier) {
    //*index_of_the_particle = prt.front();
    return -2; // Not implemented
}
int get_index_of_next_particle(int particle_identifier, int *next_particle_identifier) {
    return -2; // Not implemented
}

int get_indices_of_colliding_particles(int *index_of_particle1, int *index_of_particle2) {
    return -2; // Not implemented
}



// simulation property getters:
int get_total_mass(double *total_mass) {
    return -2; // Not implemented
}
int get_total_radius(double *total_radius) {
    return -2; // Not implemented
}
int get_time(double *time) {
    return -2; // Not implemented
}
int get_center_of_mass_position(double *x, double *y, double *z){
  //~if (mpi_rank == 0)
    //~{ // calculate only on the root mpi process, not on others
      //~*x=0; *y=0; *z=0;
      //~double m = 0.0;
      //~get_total_mass(&m)
  //~
//~
      //~for (int i = 0; i<Ntot; i++)
      //~{
         //~*x += prt[i].mass*prt[i].pos[0];
         //~*y += prt[i].mass*prt[i].pos[1];
         //~*z += prt[i].mass*prt[i].pos[2];
      //~}
//~
      //~*x /= m;
      //~*y /= m;
      //~*z /= m;
      //~return 0;
      //~}
    return -2; // Not implemented
}
int get_center_of_mass_velocity(double *vx, double *vy, double *vz) {
  //~if (mpi_rank == 0)
    //~{ // calculate only on the root mpi process, not on others
      //~*vx=0; *vy=0; *vz=0;
      //~double m = 0.0;
      //~get_total_mass(&m)
  //~
//~
      //~for (int i = 0; i<Ntot; i++)
      //~{
         //~*x += prt[i].mass*prt[i].veli[0];
         //~*y += prt[i].mass*prt[i].veli[1];
         //~*z += prt[i].mass*prt[i].veli[2];
      //~}
//~
      //~*vx /= m;
      //~*vy /= m;
      //~*vz /= m;
      //~return 0;
      //~}
    return -2; // Not implemented
}
int get_kinetic_energy(double *kinetic_energy) {
    *kinetic_energy = Ek1;
    return 0;
}
int get_potential_energy(double *potential_energy) {
    *potential_energy = Ep1;
    return 0;
}
int get_number_of_particles(int *number_of_particles) {
    *number_of_particles = Ntot - Ndead;
    return 0;
}




// particle property getters/setters: (will only work after commit_particles() is called)
int set_mass(int particle_identifier, double mass) {
    int index;
    if (found_particle(particle_identifier, &index)){
        prt[index].mass = mass;
        return 0;
    }
    return -3; // Not found!
}
int get_mass(int particle_identifier, double *mass) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *mass = prt[index].mass;
        return 0;
    }
    return -3; // Not found!
}
int set_radius(int particle_identifier, double radius) {
    int index;
    if (found_particle(particle_identifier, &index)){
        prt[index].radius = radius;
        return 0;
    }
    return -3; // Not found!
}
int get_radius(int particle_identifier, double * radius) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *radius = prt[index].radius;
        return 0;
    }
    return -3; // Not found!
}
int set_position(int particle_identifier, double x, double y, double z) {
    int index;
    if (found_particle(particle_identifier, &index)){
        prt[index].pos = Vector3(x, y, z);
        return 0;
    }
    return -3; // Not found!
}
int get_position(int particle_identifier, double *x, double *y, double *z) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *x = prt[index].pos[0];
        *y = prt[index].pos[1];
        *z = prt[index].pos[2];
        return 0;
    }
    return -3; // Not found!
}
int set_velocity(int particle_identifier, double vx, double vy, double vz) {
    int index;
    if (found_particle(particle_identifier, &index)){
        prt[index].vel = Vector3(vx, vy, vz);
        return 0;
    }
    return -3; // Not found!
}
int get_velocity(int particle_identifier, double *vx, double *vy, double *vz) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *vx = prt[index].vel[0];
        *vy = prt[index].vel[1];
        *vz = prt[index].vel[2];
        return 0;
    }
    return -3; // Not found!
}
int set_state(int particle_identifier, double mass, double radius, 
        double x, double y, double z, double vx, double vy, double vz) {
    int index;
    if (found_particle(particle_identifier, &index)){
        prt[index].mass = mass;
        prt[index].radius = radius;
        prt[index].pos = Vector3(x, y, z);
        prt[index].vel = Vector3(vx, vy, vz);
        return 0;
    }
    return -3; // Not found!
}
int get_state(int particle_identifier, double *mass, double *radius, 
        double *x, double *y, double *z, double *vx, double *vy, double *vz) {
    int index;
    if (found_particle(particle_identifier, &index)){
        *mass = prt[index].mass;
        *radius = prt[index].radius;
        *x = prt[index].pos[0];
        *y = prt[index].pos[1];
        *z = prt[index].pos[2];
        *vx = prt[index].vel[0];
        *vy = prt[index].vel[1];
        *vz = prt[index].vel[2];
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
        *ax = prt[index].acc[0];
        *ay = prt[index].acc[1];
        *az = prt[index].acc[2];
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
int set_eps2_fs_fs(double epsilon_squared_fs_fs) {
    eps2_fs_fs = epsilon_squared_fs_fs;
    return 0;
}
int get_eps2_fs_fs(double *epsilon_squared_fs_fs) {
    *epsilon_squared_fs_fs = eps2_fs_fs;
    return 0;
}
int set_eps2_fs_bh(double epsilon_squared_fs_smbh) {
    eps2_fs_smbh = epsilon_squared_fs_smbh;
    return 0;
}
int get_eps2_fs_bh(double *epsilon_squared_fs_smbh) {
    *epsilon_squared_fs_smbh = eps2_fs_smbh;
    return 0;
}
int set_eps2_bh_bh(double epsilon_squared_bh_bh) {
    eps2_bh_bh = epsilon_squared_bh_bh;
    return 0;
}
int get_eps2_bh_bh(double *epsilon_squared_bh_bh) {
    *epsilon_squared_bh_bh = eps2_bh_bh;
    return 0;
}
int set_eta_s(double eta_s_in) {
    eta_s = eta_s_in;
    return 0;
}
int get_eta_s(double *eta_s_out) {
    *eta_s_out = eta_s;
    return 0;
}
int set_eta_fs(double eta_fs_in) {
    eta_fs = eta_fs_in;
    return 0;
}
int get_eta_fs(double *eta_fs_out) {
    *eta_fs_out = eta_fs;
    return 0;
}
int set_eta_smbh(double eta_smbh_in) {
    eta_smbh = eta_smbh_in;
    return 0;
}
int get_eta_smbh(double *eta_smbh_out) {
    *eta_smbh_out = eta_smbh;
    return 0;
}
int set_eta_imbh(double eta_imbh_in) {
    eta_imbh = eta_imbh_in;
    return 0;
}
int get_eta_imbh(double *eta_imbh_out) {
    *eta_imbh_out = eta_imbh;
    return 0;
}
int get_time_step(double *time_step) {
    *time_step = -1;
    return 0; // Not implemented
}
int set_max_relative_energy_error(double max_relative_energy_error) {
    dEcrit = max_relative_energy_error;
    return 0;
}
int get_max_relative_energy_error(double *max_relative_energy_error) {
    *max_relative_energy_error = dEcrit;
    return 0;
}
int set_maximum_timestep(double maximum_timestep) {
    dt_max = maximum_timestep;
    return 0;
}
int get_maximum_timestep(double *maximum_timestep) {
    *maximum_timestep = dt_max;
    return 0;
}
int set_smbh_mass(double smbh_mass) {
    Vector3 pos=0.0;
    Vector3 vel=0.0;
    set_SMBH(smbh_mass, pos, vel);
    return 0;
}
int get_smbh_mass(double *smbh_mass) {
    Vector3 pos=0.0;
    Vector3 vel=0.0;
    get_SMBH(*smbh_mass, pos, vel);
    return 0;
}
int set_include_smbh_flag(int value){
    EX_FLAG = value;
    return 0;
}
int get_include_smbh_flag(int *value){
    *value = EX_FLAG;
    return 0;
}
int commit_parameters() {
    set_eps2(eps2_fs_fs, eps2_fs_smbh, eps2_bh_bh);
    set_eta(eta_s, eta_fs, eta_smbh, eta_imbh);
    return 0;
}
int recommit_parameters() {
    return commit_parameters();
}


//energy and gravity values
int get_potential_at_point(double *eps, double *x, double *y, double *z, double *phi, int length) {
    //~if (myrank == 0) { // calculate only on the root mpi process, not on others
        //~double r,rx,ry,rz;
        //~*phi = 0.0;
        //~for (int i = 0; i < Ntot; i++) {
            //~rx = prt[i].pos[0] - x;
            //~ry = prt[i].pos[1] - y;
            //~rz = prt[i].pos[2] - z;
            //~r = sqrt(rx*rx + ry*ry + rz*rz + eps*eps + eps2_fs_fs);
            //~*phi -= prt[i].mass/r;
        //~}
    //~}
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
        tmp_eps2_in[i] = eps[i]*eps[i] + eps2_fs_fs;
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
    calc_force_on_predictors(length, Njp,
        tmp_index, tmp_pos_in, 
        tmp_vel_in, tmp_acc_in, 
        tmp_mass_in, tmp_eps2_in,
        tmp_acc_out, tmp_jrk_out,
        tmp_snp_out, tmp_crk_out, tmp_phi_out,
        tmp_nnb_out, tmp_nnb_r2_out);
    
    MPI_Allreduce(tmp_phi_out, phi,
        length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if (EX_FLAG == 1) {
        // use function from external_field.cc to get the current properties of the SMBH
        double smbh_mass;
        Vector3 smbh_pos, smbh_vel;
        double r2;
        get_SMBH(smbh_mass, smbh_pos, smbh_vel);
        for (int i=0; i<length; i++) {
            Vector3 rij = ((Vector3) (x[i], y[i], z[i])) - smbh_pos[0];
            phi[i] -= smbh_mass * sqrt(rij*rij); //+ eps2_fs_smbh) ???
            
            // PN ???
        }
    }
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
        tmp_eps2_in[i] = eps[i]*eps[i] + eps2_fs_fs;
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
    calc_force_on_predictors(length, Njp,
        tmp_index, tmp_pos_in, 
        tmp_vel_in, acc, 
        tmp_mass_in, tmp_eps2_in,
        tmp_acc_out, tmp_jrk_out,
        tmp_snp_out, tmp_crk_out, tmp_phi_out,
        tmp_nnb_out, tmp_nnb_r2_out);
    
    MPI_Allreduce(tmp_acc_out, acc,
        3*length, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    
    if (EX_FLAG == 1) {
        // use function from external_field.cc to get the current properties of the SMBH
        double smbh_mass;
        Vector3 smbh_pos, smbh_vel;
        double r2;
        get_SMBH(smbh_mass, smbh_pos, smbh_vel);
        for (int i=0; i<length; i++) {
            Vector3 rij = ((Vector3) (x[i], y[i], z[i])) - smbh_pos[0];
            double mjR3 = smbh_mass * pow(rij*rij, 1.5); //+ eps2_fs_smbh) ???
            for (int j=0; j<3; j++) {
                acc[i][j] += -mjR3 * rij[j];
            }
            
            // PN ???
        }
    }
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

//misc
int set_pair_detect_factor(double pair_detect_factor)
{
  return 0;
}

int get_pair_detect_factor(double * pair_detect_factor)
{
  return 0;
}





//get various stopping condition information
int get_stopping_condition_info(int index, int * type, 
  int * number_of_particles)
{
  return 0;
}

int get_stopping_condition_particle_index(int index, 
  int index_of_the_column, int * index_of_particle)
{
  return 0;
}

int set_stopping_condition_out_of_box_parameter(double value)
{
  return 0;
}

int get_stopping_condition_out_of_box_parameter(double * value)
{
  return 0;
}

int get_stopping_condition_number_of_steps_parameter(int * value)
{
  return 0;
}

int set_stopping_condition_number_of_steps_parameter(int value)
{
  return 0;
}

int get_stopping_condition_timeout_parameter(double * value)
{
  return 0;
}

int get_number_of_stopping_conditions_set(int * result)
{
  return 0;
}

int set_stopping_condition_timeout_parameter(double value)
{
  return 0;
}

//~//non-set/get interactions
//~int evolve(double time)
//~{
  //~return 0;
//~}

int get_potential(int id, double *phi)
{
  return -1;
}

int get_smoothing(int id1, int id2, double *eps)
{
  return 0;
}

int set_smoothing(int id1, int id2, double eps)
{
  return 0;
}


int synchronize_model()
{
  return 0;
}

int enable_stopping_condition(int type)
{
  return 0;
}

int is_stopping_condition_set(int type, int * result)
{
  return 0;
}

int is_stopping_condition_enabled(int type, int * result)
{
  return 0;
}

int disable_stopping_condition(int type)
{
  return 0;
}

int has_stopping_condition(int type, int * result)
{
  return 0;
}



