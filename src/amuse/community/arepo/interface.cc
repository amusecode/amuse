#include "worker_code.h"

#include "src/main/allvars.h"
#include "src/main/proto.h"

int initialize_code(){
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  /* output a welcome message */
  hello();

  /* initialize CPU-time/Wallclock-time measurement */
  init_cpu_log();

  determine_compute_nodes();

  for(PTask = 0; NTask > (1 << PTask); PTask++)
    ;

  begrun0();

  strcpy(ParameterFile, "param.txt");  /* Removing command line parsing. argv[1] replaced with "param.txt". */
  RestartFlag = 0;

  begrun1(); /* set-up run  */

  char fname[MAXLEN_PATH];
  strcpy(fname, All.InitCondFile);

  /* now we can load the file */

#ifdef READ_DM_AS_GAS
      read_ic(fname, (RestartFlag == 14) ? 0x02 : LOAD_TYPES);
#else  /* #ifdef READ_DM_AS_GAS */
      read_ic(fname, (RestartFlag == 14) ? 0x01 : LOAD_TYPES);
#endif /* #ifdef READ_DM_AS_GAS #else */

  /* init returns a status code, where a value of >=0 means that endrun() should be called. */
  int status = init();

  if(status >= 0)
    {
      if(status > 0)
        mpi_printf("init() returned with %d\n", status);

      cleanup_code();
    }

  begrun2();
  return 0;
}

int cleanup_code(){
  mpi_printf("Code run for %f seconds!\n", timediff(StartOfRun, second()));
  mpi_printf("endrun called, calling MPI_Finalize()\nbye!\n\n");
  fflush(stdout);

#ifdef HAVE_HDF5
  /*The hdf5 library will sometimes register an atexit() handler that calls its
   * error handler. In AREPO this is set to my_hdf_error_handler, which calls
   * MPI_Abort. Calling MPI_Abort after MPI_Finalize is not allowed.
   * Hence unset the HDF error handler here
   */
  H5Eset_auto(NULL, NULL);
#endif /* #ifdef HAVE_HDF5 */

  MPI_Finalize();
  exit(0);
  return 0;
}

int get_mass(int index_of_the_particle, double * mass){
  return 0;
}

int commit_particles(){
  return 0;
}

int get_time(double * time){
  return 0;
}

int set_mass(int index_of_the_particle, double mass){
  return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
  return 0;
}

int get_total_radius(double * radius){
  return 0;
}

int new_particle(int * index_of_the_particle, double mass, double x, 
  double y, double z, double vx, double vy, double vz, double radius){
  return 0;
}

int get_total_mass(double * mass){
  return 0;
}

int evolve_model(double time){
  return 0;
}

int set_eps2(double epsilon_squared){
  return 0;
}

int get_begin_time(double * time){
  return 0;
}

int get_eps2(double * epsilon_squared){
  return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
  int * index_of_the_next_particle){
  return 0;
}

int delete_particle(int index_of_the_particle){
  return 0;
}

int get_potential(int index_of_the_particle, double * potential){
  return 0;
}

int synchronize_model(){
  return 0;
}

int set_state(int index_of_the_particle, double mass, double x, double y, 
  double z, double vx, double vy, double vz, double radius){
  return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
  double * y, double * z, double * vx, double * vy, double * vz, 
  double * radius){
  return 0;
}

int get_time_step(double * time_step){
  return 0;
}

int recommit_particles(){
  return 0;
}

int get_kinetic_energy(double * kinetic_energy){
  return 0;
}

int get_number_of_particles(int * number_of_particles){
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
  double az){
  return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z){
  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
  return 0;
}

int get_radius(int index_of_the_particle, double * radius){
  return 0;
}

int set_begin_time(double time){
  return 0;
}

int set_radius(int index_of_the_particle, double radius){
  return 0;
}

int recommit_parameters(){
  return 0;
}

int get_potential_energy(double * potential_energy){
  return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, 
  double * vz){
  return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, 
  double * z){
  return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z){
  return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
  double * az){
  return 0;
}

int commit_parameters(){
  return 0;
}

int set_parameters(char * param_file){
  return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, 
  double vz){
  return 0;
}

