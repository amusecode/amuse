/* twinmuse.i */
%module twinmuse
%{
extern void wrapper_set_init_dat_name(char *new_init_dat_name);
extern void wrapper_set_init_run_name(char *new_init_run_name);
extern int wrapper_initialise_twin(char *path, int nstars, double z);

extern int wrapper_load_zams_star(double mass, double age_tag);
extern double *wrapper_allocate_stellar_model_data(void);
extern void wrapper_free_stellar_model_data(double *data);
extern int wrapper_import_stellar_merger(int nmesh, int nvar, double *data, double age);
extern void wrapper_export_stellar_model(int jstar, double *data);
extern int wrapper_get_number_of_meshpoints(void);
extern int wrapper_twin_evolve(void);

extern double wrapper_get_luminosity(int id);
extern double wrapper_get_mass(int id);
extern double wrapper_get_radius(int id);
extern double wrapper_get_temperature(int id);
extern double wrapper_get_age(int id);
extern int wrapper_get_stellar_type(int id);

extern void wrapper_swap_in_star(int id);
extern void wrapper_swap_out_star(int id);
extern void wrapper_select_star(int id);
extern void wrapper_flush_star(void);
%}
%rename wrapper_set_init_dat_name set_init_dat_name;
extern void wrapper_set_init_dat_name(char *new_init_dat_name);
%rename wrapper_set_init_run_name set_init_run_name;
extern void wrapper_set_init_run_name(char *new_init_run_name);
%rename wrapper_initialise_twin initialise_twin;
extern int wrapper_initialise_twin(char *path, int nstars, double z);

%rename wrapper_load_zams_star load_zams_star;
extern int wrapper_load_zams_star(double mass, double age_tag);
%rename wrapper_allocate_stellar_model_data allocate_stellar_model_data;
extern double *wrapper_allocate_stellar_model_data(void);
%rename wrapper_free_stellar_model_data free_stellar_model_data;
extern void wrapper_free_stellar_model_data(double *data);
%rename wrapper_import_stellar_merger import_stellar_merger;
extern int wrapper_import_stellar_merger(int nmesh, int nvar, double *data, double age);
%rename wrapper_export_stellar_model export_stellar_model;
extern void wrapper_export_stellar_model(int jstar, double *data);
%rename wrapper_get_number_of_meshpoints get_number_of_meshpoints;
extern int wrapper_get_number_of_meshpoints(void);
%rename wrapper_twin_evolve twin_evolve;
extern int wrapper_twin_evolve(void);

%rename wrapper_get_luminosity get_luminosity;
extern double wrapper_get_luminosity(int id);
%rename wrapper_get_mass get_mass;
extern double wrapper_get_mass(int id);
%rename wrapper_get_radius get_radius;
extern double wrapper_get_radius(int id);
%rename wrapper_get_temperature get_temperature;
extern double wrapper_get_temperature(int id);
%rename wrapper_get_age get_age;
extern double wrapper_get_age(int id);
%rename wrapper_get_stellar_type get_stellar_type;
extern int wrapper_get_stellar_type(int id);

%rename wrapper_swap_in_star swap_in_star;
extern void wrapper_swap_in_star(int id);
%rename wrapper_swap_out_star swap_out_star;
extern void wrapper_swap_out_star(int id);
%rename wrapper_select_star select_star;
extern void wrapper_select_star(int id);
%rename wrapper_flush_star flush_star;
extern void wrapper_flush_star(void);
