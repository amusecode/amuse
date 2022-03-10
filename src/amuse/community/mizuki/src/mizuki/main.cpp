// Include the standard C++ headers
#include <cmath>
#include <math.h>
#include <cfloat>
#include <cstdio>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <sys/stat.h>
#include <time.h>
// Include the header file of FDPS
#include <particle_simulator.hpp>
// Include the header file of Phantom-GRAPE library
#if defined(ENABLE_PHANTOM_GRAPE_X86)
#include <gp5util.h>
#endif

#include "mathematical_constants.h"
#include "physical_constants.h"
#include "user_defined.hpp"
#include "mizuki.hpp"

static Mizuki* mizuki=NULL;

int main(int argc, char* argv[]){
    // Initialize FDPS
    mizuki = new Mizuki;
    PS::S64 id;
    PS::F64 mass;
    PS::F64 x;
    PS::F64 y;
    PS::F64 z;
    PS::F64 vx;
    PS::F64 vy;
    PS::F64 vz;
    PS::F64 eng;

    mizuki->initialize(argc, argv);
    mizuki->set_defaults();

    // Populate SPH with AMUSE here
    PS::S32 N_sph = 10000;
    for (PS::S32 i = 0; i < N_sph; i++) {
	mass = 1;
	x = i;
	y = 0.;
	z = i / 2.;
	vx = 0.;
	vy = 0.;
	vz = 0.;
	eng = 100.;
        mizuki->add_hydro_particle(
			i, mass, x, y, z, vx, vy, vz, eng
			);
    }
    PS::F64vec pos = mizuki->get_position(50);
    x = pos.x;
    y = pos.y;
    z = pos.z;
    std::cout << id << " " << x << " " << y << " " << z << std::endl;
    
    // Commit particles
    mizuki->commit_particles();

    std::cout << "There are " << mizuki->get_number_of_particles_sph()
	    << " particles!" << std::endl;

    // Finalize FDPS
    mizuki->finalize();
    delete mizuki;
    mizuki = NULL;
    return 0;
}
