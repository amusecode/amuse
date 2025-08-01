#include<iostream>

#include"BHtree.h"
#include"particle.h"

#include"const.h"
#include"distribution.h"
#include"system.h"
#include"mpi_interface.h"

#define GLOBAL_VALUE_DEFINE
#include"global.h"

void halt_program(){
    dev_close();
    int error = 0;
    MPI_Abort(MPI_COMM_WORLD, error);
};

void close_program(){
    dev_close();
    MPI_Finalize();
};

using namespace std;

int main(int argc, char *argv[]){

    fout_debug.open("debug.dat");

    double Egr = 0.0;
    double Emerge = 0.0;

    cout<<setprecision(15);
    cerr<<setprecision(15);
    dump_cerr(sizeof(Nbody_System));
    dump_cerr(sizeof(Hard_System));
    dump_cerr(sizeof(Soft_System));
    dump_cerr(sizeof(Particle));
    dump_cerr(sizeof(Particle_Short));
    dump_cerr(sizeof(Cell_Tree));

    mpi_initialize(argc, argv, &MYRANK, &NPROC);

    Nbody_System system;

    char param_file[STRINGSIZE];
    char output_dir[STRINGSIZE];
    if(MYRANK == 0){
        sprintf(param_file, argv[1]);
        dump_cerr(param_file);
    }

    system.read_file(param_file,  output_dir);

    system.initialize_division();

    // assign particles (prt_loc and BH_glb)
    // determine NBH_LOC and NFS_LOC
    // sum of NBH_LOC + NFS_LOC over all nodes is NFS_GLB+NBH_GLB
    system.divide_particles();
    cerr<<"finish system.divide_particles()"<<endl;

    system.calc_soft_forces(argv);
    //system.dump_tree_strcture(); // for the check of the tree strcture using jeroen's code
    cerr<<"finish system.calc_soft_forces();"<<endl;

    system.calc_energy(0);

    system.write_file(output_dir);

    int after_stellar_evolution = 0;
    long long int loop_end = (long long int)(system.Tend/system.dt_glb) + 1;
    for(long long int loop = 0; loop<loop_end; loop++){
	system.evolve(output_dir, argv, after_stellar_evolution);
	/*
	if(stellar_evolution){
	    evolve_star();
	    after_stellar_evolution = 1;
	}
	*/
    }

    close_program();

    return 0;
}



