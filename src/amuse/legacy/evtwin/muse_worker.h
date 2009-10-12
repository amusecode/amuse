extern "C" {
double wrapper_get_age(int id);
double wrapper_get_mass(int id);
double wrapper_get_radius(int id);
void wrapper_set_init_run_name(char * new_init_run_name);
int wrapper_load_zams_star(double mass, double age_tag);
double wrapper_get_temperature(int id);
int wrapper_initialise_twin(char * path, int nstars, double z);
void wrapper_set_init_dat_name(char * new_init_dat_name);
void wrapper_flush_star();
void wrapper_swap_out_star(int id);
void wrapper_select_star(int id);
int wrapper_twin_evolve();
void wrapper_swap_in_star(int id);
double wrapper_get_luminosity(int id);
int wrapper_get_stellar_type(int id);
}
