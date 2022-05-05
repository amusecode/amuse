#include <vector>
#include <math.h>

using namespace std;


class Particle {

public:
    int index;
    double x, y, z;
    double vx, vy, vz;
    double mass, radius;

    Particle (int index, double mass, 
        double x, double y, double z, 
        double vx, double vy, double vz, 
        double radius):
        index(index), mass(mass), 
        x(x), y(y), z(z), 
        vx(vx), vy(vy), vz(vz), 
        radius(radius) {
    }
    Particle (const Particle& original):index(original.index), 
        mass(original.mass), x(original.x), y(original.y), z(original.z),
        vx(original.vx), vy(original.vy), vz(original.vz), 
    radius(original.radius) {
    }
};

int highest_index = 0;
vector<Particle*> Particles;

double time = 0.;
double dt = 0.01;
double t0 = 0.;

double eps_sq = 0.;


int new_particle (int *index_of_the_particle, double mass, double x, double y, double z, double vx, double vy, double vz, double radius) {

    *index_of_the_particle = highest_index;

    Particle *p = new Particle(highest_index, mass, x, y, z, vx, vy, vz, radius);
    Particles.push_back(p);

    highest_index++;

    return 0;
}

int delete_particle (int index_of_the_particle) {

    if (index_of_the_particle > highest_index) { return -1; }
    Particles.erase(Particles.begin() + index_of_the_particle);
    return 0;
}


int get_number_of_particles (int *value) {

    *value = Particles.size();
    return 0;
}


int cleanup_code () {

    for (int i = 0; i < Particles.size(); i++) {
        delete Particles[i];
    }

    Particles.clear();

    return 0;
}



int get_index_of_first_particle (int *index) {

    if (Particles.size() > 0) {
        *index = (Particles[0])->index; 
        return 0;
    }
    return 1;
}

int get_index_of_next_particle (int index, int *next_index) {

    int N = Particles.size();

    for (int i = 0; i < N; i++) {
        if ((Particles[i])->index == index) {
            if (i == N-1) { 
                return 1;
            }
            else { 
                *next_index = (Particles[i+1])->index; 
                return 0;
            }
        }
    }

    return -1;
}


int evolve_model (double tend) {

    int N = Particles.size();
    int error;

    double *mass = new double[N];

    double *x = new double[N];
    double *y = new double[N];
    double *z = new double[N];

    double *vx = new double[N];
    double *vy = new double[N];
    double *vz = new double[N];

    for (int i = 0; i < N; i++) {
        mass[i] = (Particles[i])->mass;

        x[i] = (Particles[i])->x;
        y[i] = (Particles[i])->y;
        z[i] = (Particles[i])->z;

        vx[i] = (Particles[i])->vx;
        vy[i] = (Particles[i])->vy;
        vz[i] = (Particles[i])->vz;
    }

    while (time + t0 < tend) {
        time += dt;
        gravity_step(mass, x, y, z, vx, vy, vz, N, dt, eps_sq, 1.);
    }

    for (int i = 0; i < N; i++) {
        (Particles[i])->x = x[i];
        (Particles[i])->y = y[i];
        (Particles[i])->z = z[i];

        (Particles[i])->vx = vx[i];
        (Particles[i])->vy = vy[i];
        (Particles[i])->vz = vz[i];
    }

    delete [] mass;
    delete [] x;
    delete [] y;
    delete [] z;
    delete [] vx;
    delete [] vy;
    delete [] vz;

    return 0;
}

int get_time (double *current_time) { *current_time = time + t0; return 0; }



int get_mass (int index_of_the_particle, double *mass) {

    if (index_of_the_particle > highest_index) { return -1; }
    *mass = (Particles[index_of_the_particle])->mass;
    return 0;
}

int set_mass (int index_of_the_particle, double  mass) {

    if (index_of_the_particle > highest_index) { return -1; }
    (Particles[index_of_the_particle])->mass = mass;
    return 0;
}


int get_radius (int index_of_the_particle, double *radius) {

    if (index_of_the_particle > highest_index) { return -1; }
    *radius = (Particles[index_of_the_particle])->radius;
    return 0;
}

int set_radius (int index_of_the_particle, double  radius) {

    if (index_of_the_particle > highest_index) { return -1; }
    (Particles[index_of_the_particle])->radius = radius;
    return 0;
}


int get_position (int index_of_the_particle, double *x, double *y, double *z) {

    if (index_of_the_particle > highest_index) { return -1; }
    *x = (Particles[index_of_the_particle])->x;
    *y = (Particles[index_of_the_particle])->y;
    *z = (Particles[index_of_the_particle])->z;
    return 0;
}

int set_position (int index_of_the_particle, double  x, double  y, double  z) {

    if (index_of_the_particle > highest_index) { return -1; }
    (Particles[index_of_the_particle])->x = x;
    (Particles[index_of_the_particle])->y = y;
    (Particles[index_of_the_particle])->z = z;
    return 0;
}


int get_velocity (int index_of_the_particle, double *vx, double *vy, double *vz) {

    if (index_of_the_particle > highest_index) { return -1; }
    *vx = (Particles[index_of_the_particle])->vx;
    *vy = (Particles[index_of_the_particle])->vy;
    *vz = (Particles[index_of_the_particle])->vz;
    return 0;
}

int set_velocity (int index_of_the_particle, double  vx, double  vy, double  vz) {

    if (index_of_the_particle > highest_index) { return -1; }
    (Particles[index_of_the_particle])->vx = vx;
    (Particles[index_of_the_particle])->vy = vy;
    (Particles[index_of_the_particle])->vz = vz;
    return 0;
}


int get_state (int index_of_the_particle, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *radius) {

    if (index_of_the_particle > highest_index) { return -1; }
    *mass = (Particles[index_of_the_particle])->mass;
    *x = (Particles[index_of_the_particle])->x;
    *y = (Particles[index_of_the_particle])->y;
    *z = (Particles[index_of_the_particle])->z;
    *vx = (Particles[index_of_the_particle])->vx;
    *vy = (Particles[index_of_the_particle])->vy;
    *vz = (Particles[index_of_the_particle])->vz;
    *radius = (Particles[index_of_the_particle])->radius;
    return 0;
}

int set_state (int index_of_the_particle, double  mass, double  x, double  y, double  z, double  vx, double  vy, double  vz, double  radius) {

    if (index_of_the_particle > highest_index) { return -1; }
    (Particles[index_of_the_particle])->mass = mass;
    (Particles[index_of_the_particle])->x = x;
    (Particles[index_of_the_particle])->y = y;
    (Particles[index_of_the_particle])->z = z;
    (Particles[index_of_the_particle])->vx = vx;
    (Particles[index_of_the_particle])->vy = vy;
    (Particles[index_of_the_particle])->vz = vz;
    (Particles[index_of_the_particle])->radius = radius;
    return 0;
}



int get_eps2 (double *eps2) { *eps2 = eps_sq; return 0; }

int set_eps2 (double  eps2) {  eps_sq = eps2; return 0; }


int get_time_step (double *timestep) { *timestep = dt; return 0; }


int get_begin_time (double *begin_time) { *begin_time = t0; return 0; }

int set_begin_time (double  begin_time) {  t0 = begin_time; return 0; }


int commit_parameters () { return 0; }

int recommit_parameters () { return 0; }

int commit_particles () { return 0; }

int recommit_particles () { return 0; }

int synchronize_model () { return 0; }

int initialize_code () { return 0; }


int get_total_mass (double *total_mass) { return -2; }

int get_center_of_mass_position (double *com_x, double *com_y, double *com_z) { return -2; }

int get_center_of_mass_velocity (double *com_vx, double *com_vy, double *com_vz) { return -2; }

int get_total_radius (double *radius) { return -2; }

int get_acceleration (int index_of_the_particle, double *ax, double *ay, double *az) { return -2; }

int set_acceleration (int index_of_the_particle, double  ax, double  ay, double  az) { return -2; }

int get_potential (int index_of_the_particle, double *potential) { return -2; }

int get_kinetic_energy (double *kinetic_energy) { return -1; }

int get_potential_energy (double *potential_energy) { return -1; }
