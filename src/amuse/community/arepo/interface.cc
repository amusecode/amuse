#include <cstdio>
#include <cstring>
#include <map>

#ifndef NOMPI
#include <mpi.h>
#endif

#include "worker_code.h"
#include "interface.h"


// general interface functions:

using namespace std;


// Global ID_RLOOKUP will be initalized by create_ID_reverse_lookup in
// initialize_code()
map<MyIDType, size_t> ID_RLOOKUP;
static void create_ID_reverse_lookup();


void set_default_parameters(){
  // Relevant files
  strcpy(All.InitCondFile, "./ICs");
  strcpy(All.OutputDir,   "./output");
  strcpy(All.SnapshotFileBase, "snap");
  strcpy(All.OutputListFilename, "./output_list.txt");

  // File formats
  All.ICFormat = 1;
  All.SnapFormat = 1;

  // CPU-time LimitUBelowThisDensity
  All.TimeLimitCPU = 93000;
  All.CpuTimeBetRestartFile = 12000;
  All.ResubmitOn = 0;
  strcpy(All.ResubmitCommand, "my-scriptfile");

  // Memory allocation
  All.MaxMemSize = 2500;

  // Characteristics of run
  All.TimeBegin = 0.0;
  All.TimeMax = 2.7;
  All.TimeStep = 0.00314159;

  // Basic code options that set simulation type
  All.ComovingIntegrationOn = 0;
  All.PeriodicBoundariesOn = 0;
  All.CoolingOn = 0;
  All.StarformationOn = 0;

  // Cosmological parameters
  All.Omega0 = 0.0;
  All.OmegaLambda = 0.0;
  All.OmegaBaryon = 0.0;
  All.HubbleParam = 1.0;
  All.BoxSize = 100000.0;

  // Output frequency and output parameters
  All.OutputListOn = 1;
  All.TimeBetSnapshot = 0.1;
  All.TimeOfFirstSnapshot = 0.0;
  All.TimeBetStatistics = 0.01;
  All.NumFilesPerSnapshot = 1;
  All.NumFilesWrittenInParallel = 1;

  // Integration timing accuracy
  All.TypeOfTimestepCriterion = 0;
  All.ErrTolIntAccuracy = 0.012;
  All.CourantFac = 0.3;
  All.MaxSizeTimestep = 0.05;
  All.MinSizeTimestep = 2.0e-9;

  // Treatment of empty space and temp limits
  All.InitGasTemp = 244.8095;
  All.MinGasTemp = 5.0;
  All.MinimumDensityOnStartUp = 1.0e-20;
  All.LimitUBelowThisDensity = 0.0;
  All.LimitUBelowCertainDensityToThisValue = 0.0;
  All.MinEgySpec = 0.0;

  // Tree algorithm, force accuracy, domain update frequency
  All.TypeOfOpeningCriterion = 1;
  All.ErrTolTheta = 0.7;
  All.ErrTolForceAcc = 0.0025;
  All.MultipleDomains = 8;
  All.TopNodeFactor = 2.5;
  All.ActivePartFracForNewDomainDecomp = 0.01;

  // Initial density estimates
  All.DesNumNgb = 64;
  All.MaxNumNgbDeviation = 4;

  // System of Units
  All.UnitLength_in_cm = 3.085678e21;
  All.UnitMass_in_g = 1.989e43;
  All.UnitVelocity_in_cm_per_s = 1e5;

  // Gravitational softening lengths
  All.SofteningComoving[0] = 1.0;
  All.SofteningComoving[1] = 1.0;
  All.SofteningMaxPhys[0] = 1.0;
  All.SofteningMaxPhys[1] = 1.0;
  All.GasSoftFactor = 2.5;


  All.SofteningTypeOfPartType[0] = 0;
  All.SofteningTypeOfPartType[1] = 1;
  All.SofteningTypeOfPartType[2] = 1;
  All.SofteningTypeOfPartType[3] = 1;
  All.SofteningTypeOfPartType[4] = 1;
  All.SofteningTypeOfPartType[5] = 1;
  #ifdef ADAPTIVE_HYDRO_SOFTENING
    All.MinimumComovingHydroSoftening = 1.0;
    All.AdaptiveHydroSofteningSpacing = 1.2;
  #endif

  // Mesh regularization options
  All.CellShapingSpeed = 0.5;
  All.CellShapingFactor = 1.0;

  // parameters that are fixed for AMUSE:
  All.TreeAllocFactor = 0.8; // Memory allocation parameter
  All.ResubmitOn = 0;              // Keep this turned off!
  All.OutputListOn = 0;            // Keep this turned off
  All.GravityConstantInternal = 0; // Keep this turned off
}

int initialize_code(){

  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  /* output a welcome message */
  hello();

  /* initialize CPU-time/Wallclock-time measurement */
  init_cpu_log();

  determine_compute_nodes();
  // Needed to check available memory
  mpi_report_committable_memory();

  set_default_parameters();

  RestartFlag = 0;

  // May not need to do this (we want AMUSE to manage this)
  // MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
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
        printf("init() returned with %d\n", status);

      cleanup_code();
    }

  begrun2();

  return 0;
}

int run_sim() {
  /* This run command is for the Arepo simulation */
  run();
  return 0;
}

int cleanup_code(){
  printf("Code run for %f seconds!\n", timediff(StartOfRun, second()));
  printf("endrun called, calling MPI_Finalize()\nbye!\n\n");
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

static void create_ID_reverse_lookup() {
  map<MyIDType, size_t> id_rlookup_local;
  for (size_t i = 0; i < NumPart; i++) {
    MyIDType id = P[i].ID;
    id_rlookup_local[id] = i;
  }
  ID_RLOOKUP = id_rlookup_local;
}

static int find_particle_with_ID(int particle_id) {
  for (int j = 0; j < 2; ++j) {
    // If there's a problem, rebuild ID_RLOOKUP and retry once
    // These cautious checks may add unnecessary overhead

    auto it = ID_RLOOKUP.find(particle_id);

    if (it == ID_RLOOKUP.end()) {
      // particle_id wasn't in the map - rebuild ID_RLOOKUP and try again
      cout << "AMUSE: Rebuilding particle_ID lookup (ID not found).\n";
      create_ID_reverse_lookup();
      continue;
    }

    size_t particle_pos = (*it).second;

    if (P[particle_pos].ID != particle_id) {
      // particle_id had the wrong value - rebuild ID_RLOOKUP and try again
      cout << "AMUSE: Rebuilding particle ID lookup table (ID index changed).\n";
      create_ID_reverse_lookup();
      continue;
    }

    return particle_pos;
  }
  return -1;
}

int get_mass(int index_of_the_particle, double * mass){
  int p = find_particle_with_ID(index_of_the_particle);
  if (p >= 0) {
    *mass = P[p].Mass;
    return 0;
  }
  return -3;
}

int commit_particles(){
  return 0;
}

int get_time(double * time){
  // Return error code if calling from non-zero task
  if (ThisTask) {return 0;}
  *time = All.Time;
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
  printf("AMUSE interface: setting TimeMax from %g to %g\n", All.TimeMax, time);
  All.TimeMax = time;
  //All.TimeStep = time - All.Time;
  run();
  return 0;
}

int set_eps2(double epsilon_squared){
  // This looks bizarre
  if (ThisTask) {return 0;}
  return -2;
}

int get_begin_time(double * time){
  if (ThisTask) {return 0;}
  *time = All.TimeBegin;
  return 0;
}

int get_eps2(double * epsilon_squared){
  if (ThisTask) {return 0;}
  return -2;
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
    // Arepo has delaunay cell radii.
  return 0;
}

int get_time_step(double * time_step){
  if (ThisTask) {return 0;}
  *time_step = All.TimeStep;
  return 0;
}

int recommit_particles(){
  return 0;
}

int get_kinetic_energy(double * kinetic_energy){
  return 0;
}

int get_number_of_particles(int * number_of_particles){
  if (ThisTask) {return 0;}
  *number_of_particles = All.TotNumPart;
  return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay,
  double az){
  return -2;
}

int get_center_of_mass_position(double * x, double * y, double * z){
  return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
  return 0;
}

int get_radius(int index_of_the_particle, double * radius){
  return -2;
}

int set_begin_time(double time){
  All.TimeBegin = time;
  return 0;
}

int set_radius(int index_of_the_particle, double radius){
  return -2;
}

int recommit_parameters(){
  return 0;
}

int get_potential_energy(double * potential_energy){
  return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy,
  double * vz){
  int p = find_particle_with_ID(index_of_the_particle);
  if (p >= 0) {
    *vx = P[p].Vel[0];
    *vy = P[p].Vel[1];
    *vz = P[p].Vel[2];
    return 0;
  }
  return -3;
}

int get_position(int index_of_the_particle, double * x, double * y,
  double * z){
    int p = find_particle_with_ID(index_of_the_particle);
    if (p >= 0) {
      *x = P[p].Pos[0];
      *y = P[p].Pos[1];
      *z = P[p].Pos[2];
      return 0;
    }
    return -3;
}

int set_position(int index_of_the_particle, double x, double y, double z){
  return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay,
  double * az){
    // TODO: allvars.h defines a GravAccel vector if this is what we want?
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


