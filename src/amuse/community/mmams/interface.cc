#include "src/mmas2/src/mmas/mmas.h"
#include "worker_mmams.h"

int number_of_particles = 0;
vector<usm*> usm_models;
vector<mmas*> results;


// Default parameters:
int dump_mixed = 0;
int target_n_shells_mixing = 200;
int target_n_shells = 10000;


int initialize_code(){
    return 0;
}

int cleanup_code(){
    return 0;
}

int commit_parameters(){
    return 0;
}

int recommit_parameters(){
    return commit_parameters();
}

int new_particle(int *index_of_the_particle, double mass){
    usm *new_model = new usm;
    new_model->star_mass = mass;
    usm_models.push_back(new_model);
    *index_of_the_particle = number_of_particles;
    number_of_particles++;
    return 0;
}

int get_number_of_particles(int *number_of_particles_out){
    *number_of_particles_out = number_of_particles;
    return 0;
}

inline void is_file(char *fn) {
    ifstream fin(fn, ifstream::in);
    fin.close();
    if (fin.fail() != 0) {
        cerr << "File \"" << fn << "\" does not seem to exist ! " << endl;
        exit(-1);
    };
}

int read_usm(int *index_of_the_particle, char *usm_file){
    FILE *fmodel = NULL;
    usm *new_model = new usm;
    
    is_file(usm_file);
    fmodel = fopen(usm_file, "r");
    new_model->read(fmodel, 1);
    fclose(fmodel);
    usm_models.push_back(new_model);
    *index_of_the_particle = number_of_particles;
    number_of_particles++;
    return 0;
}

int add_shell(int index_of_the_particle, double cumul_mass, double radius, 
        double density, double pressure, double e_thermal, double entropy, 
        double temperature, double molecular_weight, double H1, double He4, 
        double O16, double N14, double C12, double Ne20, double Mg24, 
        double Si28, double Fe56){
    mass_shell shell;
    shell.mass = cumul_mass;
    shell.radius = radius;
    shell.density = density;
    shell.pressure = pressure;
    shell.e_thermal = e_thermal;
    shell.entropy = entropy;
    shell.temperature = temperature;
    shell.mean_mu = molecular_weight;
    shell.composition.H1 = H1;
    shell.composition.He4 = He4;
    shell.composition.O16 = O16;
    shell.composition.N14 = N14;
    shell.composition.C12 = C12;
    shell.composition.Ne20 = Ne20;
    shell.composition.Mg24 = Mg24;
    shell.composition.Si28 = Si28;
    shell.composition.Fe56 = Fe56;
    usm_models.at(index_of_the_particle)->add_shell(shell);
    return 0;
}

bool particle_exists(int index_of_the_particle){
    cout << (index_of_the_particle < number_of_particles && index_of_the_particle >= 0) << endl << flush;
    return (index_of_the_particle < number_of_particles && index_of_the_particle >= 0);
}

int get_shell(int index_of_the_particle, int index_of_the_shell, double *cumul_mass, 
        double *radius, double *density, double *pressure, double *e_thermal, 
        double *entropy, double *temperature, double *molecular_weight, double *H1, 
        double *He4, double *O16, double *N14, double *C12, double *Ne20, 
        double *Mg24, double *Si28, double *Fe56){
    mass_shell shell;
    
    if (index_of_the_particle < number_of_particles && index_of_the_particle >= 0){
        if (usm_models.at(index_of_the_particle)->get_num_shells()){
            usm_models.at(index_of_the_particle)->build_hashtable();
        } else {
            return -2;
        }
        shell = usm_models.at(index_of_the_particle)->get_shell(index_of_the_shell);
        *cumul_mass = shell.mass;
        *radius = shell.radius;
        *density = shell.density;
        *pressure = shell.pressure;
        *e_thermal = shell.e_thermal;
        *entropy = shell.entropy;
        *temperature = shell.temperature;
        *molecular_weight = shell.mean_mu;
        *H1 = shell.composition.H1;
        *He4 = shell.composition.He4;
        *O16 = shell.composition.O16;
        *N14 = shell.composition.N14;
        *C12 = shell.composition.C12;
        *Ne20 = shell.composition.Ne20;
        *Mg24 = shell.composition.Mg24;
        *Si28 = shell.composition.Si28;
        *Fe56 = shell.composition.Fe56;
    } else {
        return -1;
    }
    return 0;
}

int get_number_of_shells(int index_of_the_particle, int *number_of_shells){
    *number_of_shells = usm_models.at(index_of_the_particle)->get_num_shells();
    return 0;
}

int merge_two_stars(int *id_product, int id_primary, int id_secondary) {
    float r_p   = 0.0;
    float v_inf = 0.0;
    
    mmas *mmams = new mmas(*usm_models.at(id_primary), *usm_models.at(id_secondary), r_p, v_inf);
    results.push_back(mmams);
  
    mmams->merge_stars_consistently(target_n_shells);
    mmams->mixing_product(target_n_shells_mixing);
    if (!dump_mixed) {
//        cerr << "Dumping unmixed product in stdout \n";
        usm_models.push_back(&(mmams->get_product()));
//        mmams->get_product().write(stdout);
    } else {
//        cerr << "Dumping mixed product in stdout \n";
        usm_models.push_back(&(mmams->get_product()));
//        mmams->get_mixed_product().write(stdout);
    }
    *id_product = number_of_particles;
    number_of_particles++;
    return 0;
}

int set_dump_mixed_flag(int dump_mixed_flag){
    dump_mixed = dump_mixed_flag;
    return 0;
}
int get_dump_mixed_flag(int *dump_mixed_flag){
    *dump_mixed_flag = dump_mixed;
    return 0;
}

int set_target_n_shells_mixing(int target_n_shells_mixing_in){
    target_n_shells_mixing = target_n_shells_mixing_in;
    return 0;
}
int get_target_n_shells_mixing(int *target_n_shells_mixing_out){
    *target_n_shells_mixing_out = target_n_shells_mixing;
    return 0;
}

int set_target_n_shells(int target_n_shells_in){
    target_n_shells = target_n_shells_in;
    return 0;
}
int get_target_n_shells(int *target_n_shells_out){
    *target_n_shells_out = target_n_shells;
    return 0;
}
