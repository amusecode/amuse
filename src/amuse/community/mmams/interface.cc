#include "src/mmas2/src/mmas/mmas.h"
#include "src/mmas2/src/eos/eos.h"
#include "worker_code.h"
#include <map>
#include <gsl/gsl_errno.h>

using namespace std;

// Default parameters:
int dump_mixed = 1;
int target_n_shells_mixing = 200;
int target_n_shells = 10000;
int flag_do_shock_heating = 1;


int number_of_particles = 0;
int particle_id_counter = 0;
map<long long, mmas*> results;
map<long long, usm*> usm_models;
long long hashtable_up_to_date_for_particle_with_index = -1;

bool error_occurred = false;
gsl_error_handler_t * previous_error_handler;
void amuse_error_handler (const char * reason, const char * file, int line, int gsl_errno) {
    gsl_stream_printf ("ERROR", file, line, reason);
    fprintf (stderr, "AMUSE GSL error handler invoked.\n");
    error_occurred = true;
}



int initialize_code(){
    previous_error_handler = gsl_set_error_handler(&amuse_error_handler);
    return 0;
}

int cleanup_code(){
    gsl_set_error_handler(previous_error_handler);
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
	usm_models.insert(usm_models.end(), std::pair<long long, usm*>(particle_id_counter, new_model));
    *index_of_the_particle = particle_id_counter;
    number_of_particles++;
    particle_id_counter++;
    return 0;
}

int delete_particle(int index_of_the_particle){
    map<long long, mmas*>::iterator iter1 = results.find(index_of_the_particle);
    map<long long, usm*>::iterator iter2 = usm_models.find(index_of_the_particle);
    
    if (iter2 == usm_models.end())
        return -1;
    if (iter1 != results.end()){
        delete (*iter1).second;
        results.erase(iter1);
    } else {
        delete (*iter2).second;
    }
    usm_models.erase(iter2);
    number_of_particles--;
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
	usm_models.insert(usm_models.end(), std::pair<long long, usm*>(particle_id_counter, new_model));
    *index_of_the_particle = particle_id_counter;
    number_of_particles++;
    particle_id_counter++;
    return 0;
}

int add_shell(int index_of_the_particle, double d_mass, double cumul_mass, 
        double radius, double density, double pressure, 
        double temperature, double luminosity, double molecular_weight, double H1, double He4, 
        double C12, double N14, double O16, double Ne20, double Mg24, 
        double Si28, double Fe56){
    mass_shell shell;
    map<long long, usm*>::iterator it = usm_models.find(index_of_the_particle);
    
    if (it == usm_models.end())
        return -1;
    
    if (hashtable_up_to_date_for_particle_with_index == index_of_the_particle)
        hashtable_up_to_date_for_particle_with_index = -1;
    
    shell.dm = d_mass;
    shell.mass = cumul_mass;
    shell.radius = radius;
    shell.density = density;
    shell.pressure = pressure;
    shell.temperature = temperature;
    shell.luminosity = luminosity;
    shell.mean_mu = molecular_weight;
    shell.entropy = compute_entropy(density, temperature, molecular_weight);
    shell.composition.H1 = H1;
    shell.composition.He4 = He4;
    shell.composition.C12 = C12;
    shell.composition.N14 = N14;
    shell.composition.O16 = O16;
    shell.composition.Ne20 = Ne20;
    shell.composition.Mg24 = Mg24;
    shell.composition.Si28 = Si28;
    shell.composition.Fe56 = Fe56;
    
    if (it->second->get_num_shells() && shell.mass <= it->second->get_last_shell().mass)
        cerr << "Warning: shell ignored, because cumulative mass does not increase" << endl;
    else
        it->second->add_shell(shell);
    return 0;
}

int get_stellar_model_element(int index_of_the_shell, int index_of_the_particle, 
        double *d_mass, double *cumul_mass, double *radius, double *density, 
        double *pressure, double *entropy, double *temperature, double *luminosity, 
        double *molecular_weight, double *H1, double *He4, double *C12, double *N14, 
        double *O16, double *Ne20, double *Mg24, double *Si28, double *Fe56){
    mass_shell shell;
    map<long long, usm*>::iterator it = usm_models.find(index_of_the_particle);
    if (it == usm_models.end())
        return -1;
    
    if (index_of_the_shell >= it->second->get_num_shells())
        return -2;
    
    if (hashtable_up_to_date_for_particle_with_index != index_of_the_particle){
        it->second->build_hashtable();
        hashtable_up_to_date_for_particle_with_index = index_of_the_particle;
    }
    
    shell = it->second->get_shell(index_of_the_shell);
    *d_mass = shell.dm;
    *cumul_mass = shell.mass;
    *radius = shell.radius;
    *density = shell.density;
    *pressure = shell.pressure;
    *entropy = shell.entropy;
    *temperature = shell.temperature;
    *luminosity = shell.luminosity;
    *molecular_weight = shell.mean_mu;
    *H1 = shell.composition.H1;
    *He4 = shell.composition.He4;
    *C12 = shell.composition.C12;
    *N14 = shell.composition.N14;
    *O16 = shell.composition.O16;
    *Ne20 = shell.composition.Ne20;
    *Mg24 = shell.composition.Mg24;
    *Si28 = shell.composition.Si28;
    *Fe56 = shell.composition.Fe56;
    return 0;
}

int get_number_of_zones(int index_of_the_particle, int *number_of_shells){
    map<long long, usm*>::iterator it = usm_models.find(index_of_the_particle);
    if (it == usm_models.end())
        return -1;
    *number_of_shells = it->second->get_num_shells();
    return 0;
}

int merge_two_stars(int *id_product, int id_primary, int id_secondary) {
    float r_p   = 0.0;
    float v_inf = 0.0;
    map<long long, usm*>::iterator it_primary = usm_models.find(id_primary);
    map<long long, usm*>::iterator it_secondary = usm_models.find(id_secondary);
    
    if (it_primary == usm_models.end() || it_secondary == usm_models.end())
        return -1;
    
    it_primary->second->build_hashtable();
    it_secondary->second->build_hashtable();
    
    mmas *mmams = new mmas(*it_primary->second, *it_secondary->second, r_p, v_inf);
    results.insert(results.end(), std::pair<long long, mmas*>(particle_id_counter, mmams));
    
    mmams->merge_stars_consistently(target_n_shells, flag_do_shock_heating);
    mmams->mixing_product(target_n_shells_mixing);
    if (!dump_mixed) {
//        cerr << "Dumping unmixed product in stdout \n";
        usm_models.insert(--usm_models.end(), std::pair<long long, usm*>(particle_id_counter, &(mmams->get_product())));
//        mmams->get_product().write(stdout);
    } else {
//        cerr << "Dumping mixed product in stdout \n";
        usm_models.insert(--usm_models.end(), std::pair<long long, usm*>(particle_id_counter, &(mmams->get_mixed_product())));
//        mmams->get_mixed_product().write(stdout);
    }
    *id_product = particle_id_counter;
    number_of_particles++;
    particle_id_counter++;
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

int set_do_shock_heating_flag(int do_shock_heating_flag_in){
    flag_do_shock_heating = do_shock_heating_flag_in;
    return 0;
}
int get_do_shock_heating_flag(int *do_shock_heating_flag_out){
    *do_shock_heating_flag_out = flag_do_shock_heating;
    return 0;
}
