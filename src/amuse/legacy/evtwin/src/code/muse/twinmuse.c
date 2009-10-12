#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "config.h"

#ifdef __APPLE__
#define __unix__
#endif

/* external function declaraions for the getter functions.
 * These are needed because the return type is double, not int
 * (really, we should include declarations of all the other functions too, but
 * not including those doesn't cause the code to misbehave).
 */
extern double FORTRAN_NAME(get_luminosity)(int *id);
extern double FORTRAN_NAME(get_mass)(int *id);
extern double FORTRAN_NAME(get_radius)(int *id);
extern double FORTRAN_NAME(get_temperature)(int *id);
extern double FORTRAN_NAME(get_age)(int *id);

extern void FORTRAN_NAME(close_temporary_files)(void);

int wrapper_read_atomic_data_file(char *filename)
{
   /* The string lengths are passed as extra by-value arguments at the end of 
    * the parameter list. This seems to be the standard UNIX way for passing
    * FORTRAN strings around, though VMS at least is different.
    */
#ifdef __unix__
   return FORTRAN_NAME(read_atomic_data_file)(filename, strlen(filename));
#else
   #error "Don't know how to pass strings from C to FORTRAN on this system"
#endif
}

int wrapper_read_reaction_rate_file(char *filename)
{
   /* The string lengths are passed as extra by-value arguments at the end of 
    * the parameter list. This seems to be the standard UNIX way for passing
    * FORTRAN strings around, though VMS at least is different.
    */
#ifdef __unix__
   return FORTRAN_NAME(read_reaction_rate_file)(filename, strlen(filename));
#else
   #error "Don't know how to pass strings from C to FORTRAN on this system"
#endif
}

int wrapper_read_co_opacity_file(char *filename)
{
   /* The string lengths are passed as extra by-value arguments at the end of 
    * the parameter list. This seems to be the standard UNIX way for passing
    * FORTRAN strings around, though VMS at least is different.
    */
#ifdef __unix__
   return FORTRAN_NAME(read_co_opacity_file)(filename, strlen(filename));
#else
   #error "Don't know how to pass strings from C to FORTRAN on this system"
#endif
}

int wrapper_read_opacity_file(char *filename)
{
   /* The string lengths are passed as extra by-value arguments at the end of 
    * the parameter list. This seems to be the standard UNIX way for passing
    * FORTRAN strings around, though VMS at least is different.
    */
#ifdef __unix__
   return FORTRAN_NAME(read_opacity_file)(filename, strlen(filename));
#else
   #error "Don't know how to pass strings from C to FORTRAN on this system"
#endif
}

void wrapper_close_zams_library(void)
{
   FORTRAN_NAME(close_zams_library)();
}

int wrapper_open_zams_library(char *filename)
{
   /* The string lengths are passed as extra by-value arguments at the end of 
    * the parameter list. This seems to be the standard UNIX way for passing
    * FORTRAN strings around, though VMS at least is different.
    */
#ifdef __unix__
   return FORTRAN_NAME(open_zams_library)(filename, strlen(filename));
#else
   #error "Don't know how to pass strings from C to FORTRAN on this system"
#endif
}

int wrapper_read_zams_library_format(char *filename)
{
   /* The string lengths are passed as extra by-value arguments at the end of 
    * the parameter list. This seems to be the standard UNIX way for passing
    * FORTRAN strings around, though VMS at least is different.
    */
#ifdef __unix__
   return FORTRAN_NAME(read_zams_library_format)(filename, strlen(filename));
#else
   #error "Don't know how to pass strings from C to FORTRAN on this system"
#endif
}

void wrapper_set_init_dat_name(char *new_init_dat_name)
{
   /* The string lengths are passed as extra by-value arguments at the end of 
    * the parameter list. This seems to be the standard UNIX way for passing
    * FORTRAN strings around, though VMS at least is different.
    */
#ifdef __unix__
   FORTRAN_NAME(set_init_dat_name)(new_init_dat_name, strlen(new_init_dat_name));
#else
   #error "Don't know how to pass strings from C to FORTRAN on this system"
#endif
}

void wrapper_set_init_run_name(char *new_init_run_name)
{
   /* The string lengths are passed as extra by-value arguments at the end of 
    * the parameter list. This seems to be the standard UNIX way for passing
    * FORTRAN strings around, though VMS at least is different.
    */
#ifdef __unix__
   FORTRAN_NAME(set_init_run_name)(new_init_run_name, strlen(new_init_run_name));
#else
   #error "Don't know how to pass strings from C to FORTRAN on this system"
#endif
}

int wrapper_initialise_twin(char *path, int nstars, double z)
{
   char zstr[32];
   char *s;
   
   snprintf(zstr, sizeof zstr, "%g", z);
   s = zstr+2;
   
   /* Make sure the scratch files are closed and deleted when the library
    * is unloaded
    */
   atexit(FORTRAN_NAME(close_temporary_files));
   
   /* The string lengths are passed as extra by-value arguments at the end of 
    * the parameter list. This seems to be the standard UNIX way for passing
    * FORTRAN strings around, though VMS at least is different.
    */
#ifdef __unix__
   return FORTRAN_NAME(initialise_twin)(path, &nstars, s, strlen(path), strlen(s));
#else
   #error "Don't know how to pass strings from C to FORTRAN on this system"
#endif
}

int wrapper_load_zams_star(double mass, double age_tag)
{
   return FORTRAN_NAME(load_zams_star)(&mass, &age_tag);
}

/* Allocate storage space for a stellar model, for use with
 * export_stellar_model() and import_stellar_merger()
 */
double *wrapper_allocate_stellar_model_data(void)
{
   return malloc(13 * wrapper_get_number_of_meshpoints() * sizeof(double));
}

/* Free storage space allocated by allocate_stellar_model_data() */
void wrapper_free_stellar_model_data(double *data)
{
   free(data);
}

int wrapper_import_stellar_merger(int nmesh, int nvar, double *data, double age)
{
   return GLOBAL_FORTRAN_NAME(import_stellar_merger)(&nmesh, &nvar, data, &age);
}

void wrapper_export_stellar_model(int jstar, double *data)
{
   GLOBAL_FORTRAN_NAME(export_stellar_model)(&jstar, data);
}

int wrapper_get_number_of_meshpoints(void)
{
   return GLOBAL_FORTRAN_NAME(get_number_of_meshpoints)();
}

int wrapper_twin_evolve(void)
{
   return FORTRAN_NAME(twin_evolve)();
}

double wrapper_get_luminosity(int id)
{
   return FORTRAN_NAME(get_luminosity)(&id);
}

double wrapper_get_mass(int id)
{
   return FORTRAN_NAME(get_mass)(&id);
}

double wrapper_get_radius(int id)
{
   return FORTRAN_NAME(get_radius)(&id);
}

double wrapper_get_temperature(int id)
{
   return FORTRAN_NAME(get_temperature)(&id);
}

double wrapper_get_age(int id)
{
   return FORTRAN_NAME(get_age)(&id);
}

int wrapper_get_stellar_type(int id)
{
   return FORTRAN_NAME(get_stellar_type)(&id);
}

void wrapper_swap_in_star(int id)
{
   FORTRAN_NAME(swap_in_star)(&id);
}

void wrapper_swap_out_star(int id)
{
   FORTRAN_NAME(swap_out_star)(&id);
}

void wrapper_select_star(int id)
{
   FORTRAN_NAME(select_star)(&id);
}

void wrapper_flush_star(void)
{
   FORTRAN_NAME(flush_star)();
}
