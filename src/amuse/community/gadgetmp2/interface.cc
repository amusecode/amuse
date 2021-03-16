/* ################################################################################## */
/* ###                                                                            ### */
/* ###                                 Gadgetmp2                                  ### */
/* ###                                                                            ### */
/* ###   Original: Gadget2 in the version used in Amuse                           ### */
/* ###   Author: Gadget2 and Amuse contributors                                   ### */
/* ###                                                                            ### */
/* ###   Modified: July 2020                                                      ### */
/* ###   Author: Thomas Schano                                                    ### */
/* ###                                                                            ### */
/* ###   Changes are intended to enable precise calculations in                   ### */
/* ###   non periodic small domain simulations in which comoving parts            ### */
/* ###   are simulated in std precision                                           ### */
/* ###                                                                            ### */
/* ################################################################################## */
#ifndef NOMPI
#include <mpi.h>
#include <amuse_mpi.h>
#endif
#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <math.h>
#include "interface.h"
#include "worker_code.h"
//AMUSE STOPPING CONDITIONS
#include "stopcond.h"

//using namespace std;

const bool debug = false;

bool particles_initialized = false;
bool outfiles_opened = false;
bool global_quantities_of_system_up_to_date = false;
bool potential_energy_also_up_to_date = false;
bool density_up_to_date = false;
bool particle_map_up_to_date = false;
bool interpret_kicks_as_feedback = false;
bool interpret_heat_as_feedback = true;
long long particle_id_counter = 0;
long long dm_particles_in_buffer = 0;
long long sph_particles_in_buffer = 0;
long long index_of_highest_mapped_particle = 0;
map<long long, dynamics_state> dm_states;
map<long long, sph_state> sph_states;
map<long long, int> local_index_map;

double redshift_begin_parameter = 20.0;
double redshift_max_parameter = 0.0;

int numBits = my_float::get_default_prec();
int numDigits = numBits/4;

gadgetmp2 a;
std::string arg0;
std::string arg1;
std::string arg2;
std::string arg3;
std::string arg4;
std::string arg5;
std::string arg6;
std::string arg7;
std::string arg8;
std::string arg9;

// general interface functions:

void set_default_parameters(){
    // parameters that can be changed from AMUSE
    gadgetmp2::All.TimeLimitCPU = 36000;
    gadgetmp2::All.ComovingIntegrationOn = 0;
    gadgetmp2::All.TypeOfTimestepCriterion = 0;
#ifdef PERIODIC
    gadgetmp2::All.PeriodicBoundariesOn = 1;
#else
    gadgetmp2::All.PeriodicBoundariesOn = 0;
#endif
    gadgetmp2::All.Time = 0.0;
    gadgetmp2::All.TimeBegin = 0.0;
    gadgetmp2::All.TimeMax = 100.0;
    gadgetmp2::All.Omega0 = 0;
    gadgetmp2::All.OmegaLambda = 0;
    gadgetmp2::All.OmegaBaryon = 0;
    gadgetmp2::All.HubbleParam = 0.7;
    gadgetmp2::All.BoxSize = 1.0;
    gadgetmp2::All.TimeBetStatistics = 0.1;
    gadgetmp2::All.ErrTolIntAccuracy = 0.025;
    gadgetmp2::All.CourantFac = 0.15;
    gadgetmp2::All.MaxSizeTimestep = 0.01;
    gadgetmp2::All.MinSizeTimestep = 0.0;
    gadgetmp2::All.ErrTolTheta = 0.5;
    gadgetmp2::All.TypeOfOpeningCriterion = 1;
    gadgetmp2::All.ErrTolForceAcc = 0.005;
    gadgetmp2::All.TreeDomainUpdateFrequency = 0.05;
    gadgetmp2::All.DesNumNgb = 50;
    gadgetmp2::All.MaxNumNgbDeviation = 5.;
    gadgetmp2::All.ArtBulkViscConst = 0.5;
    gadgetmp2::All.ArtBulkViscBeta = 1.;
    gadgetmp2::All.MinGasTemp = 0;
    gadgetmp2::All.UnitLength_in_cm = 3.085678e21;
    gadgetmp2::All.UnitMass_in_g = 1.989e43;
    gadgetmp2::All.UnitVelocity_in_cm_per_s = 1e5;
    gadgetmp2::All.UnitTime_in_s = gadgetmp2::All.UnitLength_in_cm / gadgetmp2::All.UnitVelocity_in_cm_per_s;
    gadgetmp2::All.MinGasHsmlFractional = 0.0;
    gadgetmp2::All.SofteningGas = 0.01;
    gadgetmp2::All.SofteningHalo = 0.01;
    gadgetmp2::All.SofteningGasMaxPhys = 0.0;
    gadgetmp2::All.SofteningHaloMaxPhys = 0.0;
    gadgetmp2::set_softenings();
    strcpy(gadgetmp2::All.OutputDir,   ".");
    strcpy(gadgetmp2::All.EnergyFile,  "energy.txt");
    strcpy(gadgetmp2::All.InfoFile,    "info.txt");
    strcpy(gadgetmp2::All.TimingsFile, "timings.txt");
    strcpy(gadgetmp2::All.CpuFile,     "cpu.txt");

    // parameters that are fixed for AMUSE:
    gadgetmp2::All.PartAllocFactor = 1.5; // Memory allocation parameter
    gadgetmp2::All.TreeAllocFactor = 0.8; // Memory allocation parameter
    gadgetmp2::All.BufferSize = 25;       // Memory allocation parameter
    gadgetmp2::All.ResubmitOn = 0;              // Keep this turned off!
    gadgetmp2::All.OutputListOn = 0;            // Keep this turned off!
    gadgetmp2::All.GravityConstantInternal = 0; // Keep this turned off!

    // parameters that are unused for AMUSE:
    strcpy(gadgetmp2::All.InitCondFile, "");
    strcpy(gadgetmp2::All.RestartFile, "");
    strcpy(gadgetmp2::All.SnapshotFileBase, "");
    strcpy(gadgetmp2::All.OutputListFilename, "");
    strcpy(gadgetmp2::All.ResubmitCommand, "");
    gadgetmp2::All.ICFormat = 1;
    gadgetmp2::All.SnapFormat = 1;
    gadgetmp2::All.TimeBetSnapshot = 100.0;
    gadgetmp2::All.TimeOfFirstSnapshot = 100.0;
    gadgetmp2::All.CpuTimeBetRestartFile = 36000.0;
    gadgetmp2::All.NumFilesPerSnapshot = 1;
    gadgetmp2::All.NumFilesWrittenInParallel = 1;
    gadgetmp2::All.InitGasTemp = 0;
    gadgetmp2::All.MaxRMSDisplacementFac = 0.2; // parameter for PM; PM is currently not supported
}

int initialize_code(){
    double t0, t1;
#ifndef NOMPI
    get_comm_world(&gadgetmp2::GADGET_WORLD);
    MPI_Comm_rank(gadgetmp2::GADGET_WORLD, &gadgetmp2::ThisTask);
    MPI_Comm_size(gadgetmp2::GADGET_WORLD, &gadgetmp2::NTask);
#else
    gadgetmp2::ThisTask = 0;
    gadgetmp2::NTask = 1;
#endif
    for(gadgetmp2::PTask = 0; gadgetmp2::NTask > (1 << gadgetmp2::PTask); gadgetmp2::PTask++);
    gadgetmp2::RestartFlag = gadgetmp2::All.TotNumPart = gadgetmp2::All.TotN_gas = 0;
    gadgetmp2::All.CPU_TreeConstruction = gadgetmp2::All.CPU_TreeWalk = gadgetmp2::All.CPU_Gravity = gadgetmp2::All.CPU_Potential = gadgetmp2::All.CPU_Domain =
        gadgetmp2::All.CPU_Snapshot = gadgetmp2::All.CPU_Total = gadgetmp2::All.CPU_CommSum = gadgetmp2::All.CPU_Imbalance = gadgetmp2::All.CPU_Hydro =
        gadgetmp2::All.CPU_HydCompWalk = gadgetmp2::All.CPU_HydCommSumm = gadgetmp2::All.CPU_HydImbalance =
        gadgetmp2::All.CPU_EnsureNgb = gadgetmp2::All.CPU_Predict = gadgetmp2::All.CPU_TimeLine = gadgetmp2::All.CPU_PM = gadgetmp2::All.CPU_Peano = 0;
    gadgetmp2::CPUThisRun = 0;
    t0 = gadgetmp2::second();
    if(gadgetmp2::ThisTask == 0){
        printf("\nThis is Gadget, version `%s'.\n", GADGETVERSION);
        printf("\nRunning on %d processors.\n", gadgetmp2::NTask);
    }

    t1 = gadgetmp2::second();
    gadgetmp2::CPUThisRun += gadgetmp2::timediff(t0, t1);
    gadgetmp2::All.CPU_Total += gadgetmp2::timediff(t0, t1);

    //AMUSE STOPPING CONDITIONS SUPPORT
    set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
    set_support_for_condition(DENSITY_LIMIT_DETECTION);
    set_support_for_condition(INTERNAL_ENERGY_LIMIT_DETECTION);
    mpi_setup_stopping_conditions();

    set_default_parameters();

//    gadgetmp2::open_outputfiles();
//    outfiles_opened = true;
//    gadgetmp2::allocate_commbuffers();        /* ... allocate buffer-memory for  exchange during force computation */
    return 0;
}

int cleanup_code(){
    if (outfiles_opened)
        gadgetmp2::close_outputfiles();
    if (particles_initialized){
        gadgetmp2::free_memory();
        gadgetmp2::ngb_treefree();
        gadgetmp2::force_treefree();
    }
    return 0;
}


int check_parameters(){
#ifndef NOMPI
    MPI_Bcast(&gadgetmp2::All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, gadgetmp2::GADGET_WORLD);
#endif
    if (gadgetmp2::ThisTask){
        return 0;
    }
    if(sizeof(long long) != 8){
        printf("\nType `long long' is not 64 bit on this platform. Stopping.\n\n");
        return -4;
    }
    if(sizeof(int) != 4){
        printf("\nType `int' is not 32 bit on this platform. Stopping.\n\n");
        return -4;
    }
    if(sizeof(float) != 4){
        printf("\nType `float' is not 32 bit on this platform. Stopping.\n\n");
        return -4;
    }
    if(sizeof(double) != 8){
        printf("\nType `double' is not 64 bit on this platform. Stopping.\n\n");
        return -4;
    }
    if(gadgetmp2::All.NumFilesWrittenInParallel < 1){
        printf("NumFilesWrittenInParallel MUST be at least 1\n");
        return -4;
    }
    if(gadgetmp2::All.NumFilesWrittenInParallel > gadgetmp2::NTask){
        printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
        return -4;
    }
#ifdef PERIODIC
    if(gadgetmp2::All.PeriodicBoundariesOn == 0){
        printf("Code was compiled with periodic boundary conditions switched on.\n");
        printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
        return -4;
    }
#else
    if(gadgetmp2::All.PeriodicBoundariesOn == 1){
        printf("Code was compiled with periodic boundary conditions switched off.\n");
        printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
        return -4;
    }
#endif
    if(gadgetmp2::All.TypeOfTimestepCriterion >= 1){
        printf("The specified timestep criterion\n");
        printf("is not valid\n");
        return -4;
    }
#if defined(LONG_X) ||  defined(LONG_Y) || defined(LONG_Z)
#ifndef NOGRAVITY
    printf("Code was compiled with LONG_X/Y/Z, but not with NOGRAVITY.\n");
    printf("Stretched periodic boxes are not implemented for gravity yet.\n");
    return -4;
#endif
#endif
    cout << "Parameters successfully committed." << endl << flush;
    return 0;
}

int commit_parameters(){
    gadgetmp2::open_outputfiles();
    outfiles_opened = true;
    gadgetmp2::allocate_commbuffers();        /* ... allocate buffer-memory for particle exchange during force computation */
    gadgetmp2::set_units();
#if defined(PERIODIC) && (!defined(PMGRID) || defined(FORCETEST))
    ewald_init();
#endif
    gadgetmp2::random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(gadgetmp2::random_generator, 42);        /* start-up seed */
#ifdef PMGRID
    long_range_init();
#endif
    gadgetmp2::set_random_numbers();

    for(int i = 0; i < 6; i++)
        gadgetmp2::All.MassTable[i] = 0;

    if(gadgetmp2::All.ComovingIntegrationOn){
        gadgetmp2::All.TimeBegin = 1.0 / (1.0 + redshift_begin_parameter);
        gadgetmp2::All.TimeMax = 1.0 / (1.0 + redshift_max_parameter);
        gadgetmp2::All.Timebase_interval = (log(gadgetmp2::All.TimeMax) - log(gadgetmp2::All.TimeBegin)) / TIMEBASE;
    }else{
        gadgetmp2::All.Timebase_interval = (gadgetmp2::All.TimeMax - gadgetmp2::All.TimeBegin) / TIMEBASE;
    }
    gadgetmp2::All.Time = gadgetmp2::All.TimeBegin;
    gadgetmp2::All.Ti_Current = 0;
    gadgetmp2::All.NumCurrentTiStep = 0;        /* setup some counters */
    gadgetmp2::All.SnapshotFileCount = 0;
    gadgetmp2::All.TotNumOfForces = gadgetmp2::All.NumForcesSinceLastDomainDecomp = 0;
    gadgetmp2::All.TimeLastStatistics = gadgetmp2::All.TimeBegin - gadgetmp2::All.TimeBetStatistics;
#ifdef PMGRID
    gadgetmp2::All.PM_Ti_endstep = gadgetmp2::All.PM_Ti_begstep = 0;
#endif
#ifdef FLEXSTEPS
    gadgetmp2::All.PresentMinStep = TIMEBASE;
#endif
    gadgetmp2::All.Ti_nextoutput = -1; // find_next_outputtime(gadgetmp2::All.Ti_Current);
    return check_parameters();
}
int recommit_parameters(){
    return check_parameters();
}

int commit_particles(){
    double t0, t1;
    int i, j;
#ifndef ISOTHERM_EQS
    my_float a3;
#endif
    my_float a, a_inv;

    if (gadgetmp2::All.ComovingIntegrationOn) {
        a = gadgetmp2::All.Time;
        a_inv = (my_float)"1.0" / gadgetmp2::All.Time;
    } else {
        a = a_inv = "1";
    }
    t0 = gadgetmp2::second();
    gadgetmp2::All.TotNumPart = dm_particles_in_buffer + sph_particles_in_buffer;
    gadgetmp2::All.TotN_gas = sph_particles_in_buffer;
    gadgetmp2::All.MaxPart = gadgetmp2::All.PartAllocFactor * (gadgetmp2::All.TotNumPart / gadgetmp2::NTask);        /* sets the maximum number of particles that may */
    gadgetmp2::All.MaxPartSph = gadgetmp2::All.PartAllocFactor * (gadgetmp2::All.TotN_gas / gadgetmp2::NTask);        /* sets the maximum number of particles that may
                                                                        reside on a processor */
    gadgetmp2::NumPart = dm_states.size()+sph_states.size();
    gadgetmp2::N_gas = sph_states.size();
    gadgetmp2::allocate_memory();
    // initialize sph particles
    i = 0;
    for (map<long long, sph_state>::iterator state_iter = sph_states.begin();
            state_iter != sph_states.end(); state_iter++, i++){
        gadgetmp2::P[i].ID = (*state_iter).first;
        gadgetmp2::P[i].Mass = (*state_iter).second.mass;
        gadgetmp2::P[i].Pos[0] = (*state_iter).second.x * a_inv;
        gadgetmp2::P[i].Pos[1] = (*state_iter).second.y * a_inv;
        gadgetmp2::P[i].Pos[2] = (*state_iter).second.z * a_inv;
        gadgetmp2::P[i].Vel[0] = (*state_iter).second.vx * a;
        gadgetmp2::P[i].Vel[1] = (*state_iter).second.vy * a;
        gadgetmp2::P[i].Vel[2] = (*state_iter).second.vz * a;
        gadgetmp2::P[i].Type = 0; // SPH particles (dark matter particles have type 1)
        gadgetmp2::SphP[i].Entropy = (*state_iter).second.u;
        gadgetmp2::SphP[i].Density = "-1";
        gadgetmp2::SphP[i].Hsml.setZero(); // = 0;
#ifdef MORRIS97VISC
        SphP[i].Alpha = (*state_iter).second.alpha;
        SphP[i].DAlphaDt = (*state_iter).second.dalphadt;
#endif
    }
    sph_states.clear();

    // initialize dark matter particles
    i = gadgetmp2::N_gas;
    for (map<long long, dynamics_state>::iterator state_iter = dm_states.begin();
            state_iter != dm_states.end(); state_iter++, i++){
        gadgetmp2::P[i].ID = (*state_iter).first;
        gadgetmp2::P[i].Mass = (*state_iter).second.mass;
        gadgetmp2::P[i].Pos[0] = (*state_iter).second.x * a_inv;
        gadgetmp2::P[i].Pos[1] = (*state_iter).second.y * a_inv;
        gadgetmp2::P[i].Pos[2] = (*state_iter).second.z * a_inv;
        gadgetmp2::P[i].Vel[0] = (*state_iter).second.vx * a;
        gadgetmp2::P[i].Vel[1] = (*state_iter).second.vy * a;
        gadgetmp2::P[i].Vel[2] = (*state_iter).second.vz * a;
        gadgetmp2::P[i].radius = (*state_iter).second.radius;
        gadgetmp2::P[i].Type = 1; // dark matter particles (SPH particles have type 0)
    }
    dm_states.clear();
    gadgetmp2::All.TimeBegin += gadgetmp2::All.Ti_Current * gadgetmp2::All.Timebase_interval;
    gadgetmp2::All.Ti_Current = 0;
    gadgetmp2::All.Time = gadgetmp2::All.TimeBegin;
    gadgetmp2::set_softenings();
    for(i = 0; i < gadgetmp2::NumPart; i++){        /*  start-up initialization */
        for(j = 0; j < 3; j++)
            gadgetmp2::P[i].GravAccel[j].setZero(); // = 0;
#ifdef PMGRID
        for(j = 0; j < 3; j++)
            gadgetmp2::P[i].GravPM[j].setZero(); // = 0;
#endif
        gadgetmp2::P[i].Ti_endstep = 0;
        gadgetmp2::P[i].Ti_begstep = 0;
        gadgetmp2::P[i].OldAcc.setZero(); // = 0;
        gadgetmp2::P[i].GravCost = "1";
        gadgetmp2::P[i].Potential.setZero(); // = 0;
    }
#ifdef FLEXSTEPS
    for(i = 0; i < gadgetmp2::NumPart; i++)        /*  start-up initialization */
        P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
#endif
    for(i = 0; i < gadgetmp2::N_gas; i++){        /* initialize sph_properties */
        for(j = 0; j < 3; j++){
            gadgetmp2::SphP[i].VelPred[j] = gadgetmp2::P[i].Vel[j];
            gadgetmp2::SphP[i].HydroAccel[j].setZero(); // = 0;
        }
        gadgetmp2::SphP[i].DtEntropy.setZero(); // = 0;
#ifdef TIMESTEP_UPDATE
        gadgetmp2::SphP[i].FeedbackFlag = 0;
        for(j = 0; j < 3; j++)
            gadgetmp2::SphP[i].FeedAccel[j].setZero(); // = 0;
#endif
    }

    if((gadgetmp2::All.MaxPart < 1000)){
        if (gadgetmp2::ThisTask == 0){
            cout << "Gadget takes "<< gadgetmp2::All.PartAllocFactor  << " times the number of particles on a processors as a maximum."<<endl;
            cout << "For large numbers of particles some room is always available for storing nodes from other processors." << endl;
            cout << "For smaller numbers, this assumption is incorrect."<<endl;
            cout << "Changed allocation of tree to include more nodes."<<endl;
        }
        //gadgetmp2::All.MaxPart = MAXTOPNODES;
        gadgetmp2::All.TreeAllocFactor = 4000.0 / gadgetmp2::All.MaxPart;
        gadgetmp2::ngb_treeallocate(MAX_NGB);
        gadgetmp2::force_treeallocate(10 * gadgetmp2::All.TreeAllocFactor * gadgetmp2::All.MaxPart, 10 * gadgetmp2::All.MaxPart);
    } else {
        gadgetmp2::ngb_treeallocate(MAX_NGB);
        gadgetmp2::force_treeallocate(gadgetmp2::All.TreeAllocFactor * gadgetmp2::All.MaxPart, gadgetmp2::All.MaxPart);
    }

    gadgetmp2::All.NumForcesSinceLastDomainDecomp = 1 + gadgetmp2::All.TotNumPart * gadgetmp2::All.TreeDomainUpdateFrequency;
    gadgetmp2::Flag_FullStep = 1;                /* to ensure that Peano-Hilber order is done */
    gadgetmp2::domain_Decomposition();        /* do initial domain decomposition (gives equal numbers of particles) */
    update_particle_map();

    index_of_highest_mapped_particle = local_index_map.rbegin()->first;
#ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &index_of_highest_mapped_particle, 1, MPI_LONG_LONG_INT, MPI_MAX, gadgetmp2::GADGET_WORLD);
#endif
    gadgetmp2::ngb_treebuild();                /* will build tree */
    gadgetmp2::setup_smoothinglengths();
    gadgetmp2::TreeReconstructFlag = 1;
  /* at this point, the entropy variable normally contains the
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly.
   * Once the density has been computed, we can convert thermal energy to entropy.
   */
#ifndef ISOTHERM_EQS
    if(gadgetmp2::All.ComovingIntegrationOn){a3 = gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time;}else{a3 = 1;}
    for(i = 0; i < gadgetmp2::N_gas; i++)
        gadgetmp2::SphP[i].Entropy = gadgetmp2::const_GAMMA_MINUS1 * gadgetmp2::SphP[i].Entropy / pow(gadgetmp2::SphP[i].Density / a3, gadgetmp2::const_GAMMA_MINUS1);
#endif
#ifdef PMGRID
    long_range_init_regionsize();
#endif
    if(gadgetmp2::All.ComovingIntegrationOn)
        gadgetmp2::init_drift_table();
    t1 = gadgetmp2::second();
    gadgetmp2::CPUThisRun += gadgetmp2::timediff(t0, t1);
    gadgetmp2::All.CPU_Total += gadgetmp2::timediff(t0, t1);
    particles_initialized = true;
    if (gadgetmp2::ThisTask == 0){
        cout << flush;
    }
    return 0;
}

void push_particle_data_on_state_vectors(){
    map<long long, int>::iterator iter;
    int i;
    my_float a_inv, a;
#ifndef ISOTHERM_EQS
    my_float a3;
    if(gadgetmp2::All.ComovingIntegrationOn){a3 = gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }
#endif
    if (gadgetmp2::All.ComovingIntegrationOn) {
        a = gadgetmp2::All.Time;
        a_inv = (my_float)"1.0" / gadgetmp2::All.Time;
    } else {
        a = a_inv = "1";
    }
    for (iter = local_index_map.begin(); iter != local_index_map.end(); iter++){
        i = (*iter).second;
        if (gadgetmp2::P[i].Type == 0){
            // store sph particle data
            sph_state state;
            state.mass = gadgetmp2::P[i].Mass;
            state.x =    gadgetmp2::P[i].Pos[0] * a;
            state.y =    gadgetmp2::P[i].Pos[1] * a;
            state.z =    gadgetmp2::P[i].Pos[2] * a;
            state.vx =   gadgetmp2::P[i].Vel[0] * a_inv;
            state.vy =   gadgetmp2::P[i].Vel[1] * a_inv;
            state.vz =   gadgetmp2::P[i].Vel[2] * a_inv;
#ifdef ISOTHERM_EQS
            state.u = SphP[i].Entropy;
#else
            state.u = gadgetmp2::SphP[i].Entropy * pow(gadgetmp2::SphP[i].Density / a3, gadgetmp2::const_GAMMA_MINUS1) / gadgetmp2::const_GAMMA_MINUS1;
#endif

#ifdef MORRIS97VISC
            state.alpha = SphP[i].Alpha;
            state.dalphadt = SphP[i].DAlphaDt;
#endif


            sph_states.insert(std::pair<long long, sph_state>(gadgetmp2::P[i].ID, state));
        } else {
            // store dark matter particle data
            dynamics_state state;
            state.mass = gadgetmp2::P[i].Mass;
            state.x =    gadgetmp2::P[i].Pos[0] * a;
            state.y =    gadgetmp2::P[i].Pos[1] * a;
            state.z =    gadgetmp2::P[i].Pos[2] * a;
            state.vx =   gadgetmp2::P[i].Vel[0] * a_inv;
            state.vy =   gadgetmp2::P[i].Vel[1] * a_inv;
            state.vz =   gadgetmp2::P[i].Vel[2] * a_inv;
            state.radius =   gadgetmp2::P[i].radius;
            dm_states.insert(std::pair<long long, dynamics_state>(gadgetmp2::P[i].ID, state));
        }
    }
}

int recommit_particles(){
    if (particles_initialized){
        push_particle_data_on_state_vectors();
        gadgetmp2::free_memory();
        gadgetmp2::ngb_treefree();
        gadgetmp2::force_treefree();
    }
    return commit_particles();
}

bool drift_to_t_end(int ti_end){
    bool done;
    int n, min, min_glob, flag, *temp;
    my_float timeold;
    double t0, t1;
    t0 = gadgetmp2::second();
    timeold = gadgetmp2::All.Time;
    for(n = 1, min = gadgetmp2::P[0].Ti_endstep; n < gadgetmp2::NumPart; n++){
        if(min > gadgetmp2::P[n].Ti_endstep){
            min = gadgetmp2::P[n].Ti_endstep;
        }
    }
#ifndef NOMPI
    MPI_Allreduce(&min, &min_glob, 1, MPI_INT, MPI_MIN, gadgetmp2::GADGET_WORLD);
#else
    min_glob = min;
#endif
    /* We check whether this is a full step where all particles are synchronized */
    flag = 1;
    for(n = 0; n < gadgetmp2::NumPart; n++)
        if(gadgetmp2::P[n].Ti_endstep > min_glob)
            flag = 0;
#ifndef NOMPI
    MPI_Allreduce(&flag, &gadgetmp2::Flag_FullStep, 1, MPI_INT, MPI_MIN, gadgetmp2::GADGET_WORLD);
#else
    gadgetmp2::Flag_FullStep = flag;
#endif
#ifdef PMGRID
    if(min_glob >= gadgetmp2::All.PM_Ti_endstep){
        min_glob = gadgetmp2::All.PM_Ti_endstep;
        gadgetmp2::Flag_FullStep = 1;
    }
#endif
    /* Determine 'NumForceUpdate', i.e. the number of particles on this processor that are going to be active */
    for(n = 0, gadgetmp2::NumForceUpdate = 0; n < gadgetmp2::NumPart; n++){
        if(gadgetmp2::P[n].Ti_endstep == min_glob)
#ifdef SELECTIVE_NO_GRAVITY
          if(!((1 << P[n].Type) & (SELECTIVE_NO_GRAVITY)))
#endif
            gadgetmp2::NumForceUpdate++;
    }
    /* note: NumForcesSinceLastDomainDecomp has type "long long" */
    temp = new int[gadgetmp2::NTask];

#ifndef NOMPI
    MPI_Allgather(&gadgetmp2::NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, gadgetmp2::GADGET_WORLD);
#else
    temp[0] = gadgetmp2::NumForceUpdate;
#endif
    for(n = 0; n < gadgetmp2::NTask; n++)
        gadgetmp2::All.NumForcesSinceLastDomainDecomp += temp[n];
    free(temp);
    t1 = gadgetmp2::second();
    gadgetmp2::All.CPU_Predict += gadgetmp2::timediff(t0, t1);
    if (min_glob >= ti_end){
        min_glob = ti_end;
        done = true;
    } else {
        done = false;
    }
    gadgetmp2::move_particles(gadgetmp2::All.Ti_Current, min_glob);
    gadgetmp2::All.Ti_Current = min_glob;
    if(gadgetmp2::All.ComovingIntegrationOn)
        gadgetmp2::All.Time = gadgetmp2::All.TimeBegin * exp(gadgetmp2::All.Ti_Current * gadgetmp2::All.Timebase_interval);
    else
        gadgetmp2::All.Time = gadgetmp2::All.TimeBegin + gadgetmp2::All.Ti_Current * gadgetmp2::All.Timebase_interval;
    gadgetmp2::All.TimeStep = gadgetmp2::All.Time - timeold.toDouble();
    return done;
}

bool check_density_stopping_condition(){
    int stopping_condition_is_set;
    double minimum_density_parameter, maximum_density_parameter;
    get_stopping_condition_minimum_density_parameter(&minimum_density_parameter);
    get_stopping_condition_maximum_density_parameter(&maximum_density_parameter);
    for (int i=0; i<gadgetmp2::N_gas; i++) {
        if ( (gadgetmp2::SphP[i].Density < minimum_density_parameter) ||
             (gadgetmp2::SphP[i].Density > maximum_density_parameter)) {
            int stopping_index  = next_index_for_stopping_condition();
            if (stopping_index >= 0) {
                cout << "set_stopping_condition_info returned: " <<
                    set_stopping_condition_info(stopping_index, DENSITY_LIMIT_DETECTION) << endl;
                cout << "set_stopping_condition_particle_index returned: " <<
                    set_stopping_condition_particle_index(stopping_index, 0, gadgetmp2::P[i].ID) << endl;
            }
        }
    }

    mpi_collect_stopping_conditions();
    is_stopping_condition_set(DENSITY_LIMIT_DETECTION, &stopping_condition_is_set);
    return stopping_condition_is_set;
}
bool check_internal_energy_stopping_condition(){
    int stopping_condition_is_set;
    my_float internal_energy;
    double minimum_internal_energy_parameter, maximum_internal_energy_parameter;
    get_stopping_condition_minimum_internal_energy_parameter(&minimum_internal_energy_parameter);
    get_stopping_condition_maximum_internal_energy_parameter(&maximum_internal_energy_parameter);

#ifndef ISOTHERM_EQS
    my_float a3;
    if(gadgetmp2::All.ComovingIntegrationOn){a3 = gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time;}else{a3 = 1;}
#endif

    for (int i=0; i<gadgetmp2::N_gas; i++) {
#ifdef ISOTHERM_EQS
        internal_energy = SphP[i].Entropy;
#else
        internal_energy = gadgetmp2::SphP[i].Entropy *
                pow(gadgetmp2::SphP[i].Density / a3, gadgetmp2::const_GAMMA_MINUS1) / gadgetmp2::const_GAMMA_MINUS1;
#endif
        if ( (internal_energy < minimum_internal_energy_parameter) ||
             (internal_energy > maximum_internal_energy_parameter)) {
            int stopping_index  = next_index_for_stopping_condition();
            if (stopping_index > 0) {
                cout << "set_stopping_condition_info returned: " <<
                    set_stopping_condition_info(stopping_index, INTERNAL_ENERGY_LIMIT_DETECTION) << endl;
                cout << "set_stopping_condition_particle_index returned: " <<
                    set_stopping_condition_particle_index(stopping_index, 0, gadgetmp2::P[i].ID) << endl;
            }
        }
    }

    mpi_collect_stopping_conditions();
    is_stopping_condition_set(INTERNAL_ENERGY_LIMIT_DETECTION, &stopping_condition_is_set);
    return stopping_condition_is_set;
}


int evolve_model_generic(double t_end){
    bool done;
    double t0, t1;
    int Ti_end, stopflag = 0;
    // AMUSE STOPPING CONDITIONS
    int is_number_of_steps_detection_enabled;
    int is_density_limit_detection_enabled;
    int is_internal_energy_limit_detection_enabled;
    int max_number_of_steps;
    int number_of_steps_innerloop = 0;

    reset_stopping_conditions();
    is_stopping_condition_enabled(NUMBER_OF_STEPS_DETECTION,
        &is_number_of_steps_detection_enabled);
    is_stopping_condition_enabled(DENSITY_LIMIT_DETECTION,
        &is_density_limit_detection_enabled);
    is_stopping_condition_enabled(INTERNAL_ENERGY_LIMIT_DETECTION,
        &is_internal_energy_limit_detection_enabled);
    get_stopping_condition_number_of_steps_parameter(&max_number_of_steps);

    // .......

    if (t_end > gadgetmp2::All.TimeMax)
        return -7;
    gadgetmp2::ZeroTimestepEncountered = 0;
    if(gadgetmp2::All.ComovingIntegrationOn){
        Ti_end = (log(t_end / gadgetmp2::All.TimeBegin) / gadgetmp2::All.Timebase_interval);
    } else {
        Ti_end = ((t_end - gadgetmp2::All.TimeBegin) / gadgetmp2::All.Timebase_interval);
    }
    if (Ti_end >= gadgetmp2::All.Ti_Current){
        global_quantities_of_system_up_to_date = density_up_to_date = false;
        done = drift_to_t_end(Ti_end); /* find next synchronization point and drift particles to MIN(this time, t_end). */
        while (!done && gadgetmp2::All.Ti_Current < TIMEBASE && gadgetmp2::All.Time <= gadgetmp2::All.TimeMax) {
            t0 = gadgetmp2::second();
            gadgetmp2::every_timestep_stuff();        /* write some info to log-files */
            gadgetmp2::domain_Decomposition();        /* do domain decomposition if needed */
            particle_map_up_to_date = false;
            gadgetmp2::compute_accelerations(0);        /* compute accelerations for
                * the particles that are to be advanced  */

            // AMUSE stopping conditions: density and internal energy check
            // After compute_accelerations(), SPH particle densities are up to date
            if (is_density_limit_detection_enabled) {
                if (check_density_stopping_condition()) break;
            }
            if (is_internal_energy_limit_detection_enabled) {
                if (check_internal_energy_stopping_condition()) break;
            }

            /* check whether we want a full energy statistics */
            if((gadgetmp2::All.Time - gadgetmp2::All.TimeLastStatistics) >= gadgetmp2::All.TimeBetStatistics) {
#ifdef COMPUTE_POTENTIAL_ENERGY
                compute_potential();
#endif
                gadgetmp2::energy_statistics();        /* compute and output energy statistics */
                gadgetmp2::All.TimeLastStatistics += gadgetmp2::All.TimeBetStatistics;
            }
            gadgetmp2::advance_and_find_timesteps();        /* 'kick' active particles in
                            * momentum space and compute new
                            * timesteps for them  */
            done = drift_to_t_end(Ti_end);
            gadgetmp2::All.NumCurrentTiStep++;

            /* Check whether we need to interrupt the run */
#ifndef NOMPI
            MPI_Allreduce(MPI_IN_PLACE, &gadgetmp2::ZeroTimestepEncountered, 1, MPI_INT, MPI_MAX, gadgetmp2::GADGET_WORLD);
#endif
            if(gadgetmp2::ZeroTimestepEncountered)
                return -8;

            if(gadgetmp2::ThisTask == 0) {
                /* are we running out of CPU-time ? If yes, interrupt run. */
                if(gadgetmp2::CPUThisRun > 0.85 * gadgetmp2::All.TimeLimitCPU){
                    printf("reaching time-limit. stopping.\n");
                    stopflag = 2;
                }
            }
#ifndef NOMPI
            MPI_Bcast(&stopflag, 1, MPI_INT, 0, gadgetmp2::GADGET_WORLD);
#endif
            if(stopflag)
                return -5;
            t1 = gadgetmp2::second();
            gadgetmp2::All.CPU_Total += gadgetmp2::timediff(t0, t1);
            gadgetmp2::CPUThisRun += gadgetmp2::timediff(t0, t1);
            //
            if(is_number_of_steps_detection_enabled) {
                number_of_steps_innerloop++;
                if(number_of_steps_innerloop > max_number_of_steps) {
                    int stopping_index  = next_index_for_stopping_condition();
                    set_stopping_condition_info(stopping_index, NUMBER_OF_STEPS_DETECTION);
                }
            }
            if (set_conditions & enabled_conditions) {
                break;
            }
        }
    } else {
        return -6;
    }
    if (gadgetmp2::ThisTask == 0)
        cout << flush;
    if (gadgetmp2::All.Ti_Current > TIMEBASE || gadgetmp2::All.Time > gadgetmp2::All.TimeMax)
        return -7;
    return 0;
}
int evolve_to_redshift(double redshift){
    if (!gadgetmp2::All.ComovingIntegrationOn) {return -9;}
    return evolve_model_generic(1.0 / (1.0 + redshift));
}
int evolve_model(double t_end){
    if (gadgetmp2::All.ComovingIntegrationOn) {return -9;}
    return evolve_model_generic(t_end);
}

int synchronize_model() {
    return 0;
}

int construct_tree_if_needed(void){
    double tstart, tend;
    if (!particles_initialized)
        return -1;
    tstart = gadgetmp2::second();
    if (gadgetmp2::TreeReconstructFlag){
        if(gadgetmp2::ThisTask == 0)
            printf("Tree construction.\n");
        gadgetmp2::force_treebuild(gadgetmp2::NumPart);
        gadgetmp2::TreeReconstructFlag = 0;
        if(gadgetmp2::ThisTask == 0)
            printf("Tree construction done.\n");
    }
    tend = gadgetmp2::second();
    gadgetmp2::All.CPU_TreeConstruction += gadgetmp2::timediff(tstart, tend);
    return 0;
}

int new_dm_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double radius){
    particle_id_counter++;
    if (gadgetmp2::ThisTask == 0)
        *id = particle_id_counter;
    // Divide the particles equally over all Tasks, Gadget will redistribute them later.
    if (gadgetmp2::ThisTask == (dm_particles_in_buffer % gadgetmp2::NTask)){
        dynamics_state state;
        state.mass = mass;
        state.x = x;
        state.y = y;
        state.z = z;
        state.vx = vx;
        state.vy = vy;
        state.vz = vz;
        dm_states.insert(std::pair<long long, dynamics_state>(particle_id_counter, state));
    }
    dm_particles_in_buffer++;
    return 0;
}

int new_sph_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double u){
    exit(0);
    particle_id_counter++;
    if (gadgetmp2::ThisTask == 0)
        *id = particle_id_counter;

    // Divide the sph particles equally over all Tasks, Gadget will redistribute them later.
    if (gadgetmp2::ThisTask == (sph_particles_in_buffer % gadgetmp2::NTask)){
        sph_state state;
        state.mass = mass;
        state.x = x;
        state.y = y;
        state.z = z;
        state.vx = vx;
        state.vy = vy;
        state.vz = vz;
        state.u = u;
#ifdef MORRIS97VISC
        state.alpha = gadgetmp2::All.ArtBulkViscConst;
        state.dalphadt = 0;
#endif
        sph_states.insert(std::pair<long long, sph_state>(particle_id_counter, state));
    }
    sph_particles_in_buffer++;
    return 0;
}

int delete_particle(int id){
    int found = 0;
    map<long long, int>::iterator it;
    map<long long, dynamics_state>::iterator dyn_it;
    map<long long, sph_state>::iterator sph_it;

    if (!particle_map_up_to_date)
        update_particle_map();

    it = local_index_map.find(id);
    if (it != local_index_map.end()){
        local_index_map.erase(it);
        found = 1 + gadgetmp2::P[(*it).second].Type; // 1 for sph; 2 for dm
    }
#ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, gadgetmp2::GADGET_WORLD);
#endif
    if (found){
        if (found == 2)
            dm_particles_in_buffer--;
        else
            sph_particles_in_buffer--;
        return 0;
    }

    dyn_it = dm_states.find(id);
    if (dyn_it != dm_states.end()){
        dm_states.erase(dyn_it);
        found = 1;
    }
#ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_LOR, gadgetmp2::GADGET_WORLD);
#endif
    if (found){
        dm_particles_in_buffer--;
        return 0;
    }

    sph_it = sph_states.find(id);
    if (sph_it != sph_states.end()){
        sph_states.erase(sph_it);
        found = 1;
    }
#ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_LOR, gadgetmp2::GADGET_WORLD);
#endif
    if (found){
        sph_particles_in_buffer--;
        return 0;
    }
    return -3;
}

// parameter getters/setters:

int get_time_step(double *timestep){
    if (gadgetmp2::ThisTask) {return 0;}
    *timestep = gadgetmp2::All.TimeStep;
    return 0;
}
int set_time_step(double timestep){
    if (gadgetmp2::ThisTask) {return 0;}
    return -2;
}
int get_epsilon(double *epsilon){
    if (gadgetmp2::ThisTask) {return 0;}
    gadgetmp2::set_softenings();
    *epsilon = gadgetmp2::All.SofteningTable[1];
    return 0;
}
int set_epsilon(double epsilon){
    gadgetmp2::All.SofteningHalo = epsilon;
    gadgetmp2::set_softenings();
    return 0;
}
int get_eps2(double *epsilon_squared){
    if (gadgetmp2::ThisTask) {return 0;}
    return -2;
}
int set_eps2(double epsilon_squared){
    if (gadgetmp2::ThisTask) {return 0;}
    return -2;
}
int get_epsgas(double *gas_epsilon){
    if (gadgetmp2::ThisTask) {return 0;}
    gadgetmp2::set_softenings();
    *gas_epsilon = gadgetmp2::All.SofteningTable[0];
    return 0;
}
int set_epsgas(double gas_epsilon){
    gadgetmp2::All.SofteningGas = gas_epsilon;
    gadgetmp2::set_softenings();
    return 0;
}
int get_unit_mass(double *code_mass_unit){
    if (gadgetmp2::ThisTask) {return 0;}
    *code_mass_unit = gadgetmp2::All.UnitMass_in_g;
    return 0;
}
int set_unit_mass(double code_mass_unit){
    gadgetmp2::All.UnitMass_in_g = code_mass_unit;
    return 0;
}
int get_unit_length(double *code_length_unit){
    if (gadgetmp2::ThisTask) {return 0;}
    *code_length_unit = gadgetmp2::All.UnitLength_in_cm;
    return 0;
}
int set_unit_length(double code_length_unit){
    gadgetmp2::All.UnitLength_in_cm = code_length_unit;
    return 0;
}
int get_unit_time(double *code_time_unit){
    if (gadgetmp2::ThisTask) {return 0;}
    *code_time_unit = gadgetmp2::All.UnitTime_in_s;
    return 0;
}
int set_unit_time(double code_time_unit){
    gadgetmp2::All.UnitVelocity_in_cm_per_s = gadgetmp2::All.UnitLength_in_cm / code_time_unit;
    gadgetmp2::All.UnitTime_in_s = code_time_unit;
    return 0;
}
int get_unit_velocity(double *code_velocity_unit){
    if (gadgetmp2::ThisTask) {return 0;}
    *code_velocity_unit = gadgetmp2::All.UnitVelocity_in_cm_per_s;
    return 0;
}
int set_unit_velocity(double code_velocity_unit){
    gadgetmp2::All.UnitVelocity_in_cm_per_s = code_velocity_unit;
    gadgetmp2::All.UnitTime_in_s = gadgetmp2::All.UnitLength_in_cm / gadgetmp2::All.UnitVelocity_in_cm_per_s;
    return 0;
}

int get_gadget_output_directory(char **output_directory){
    if (gadgetmp2::ThisTask) {return 0;}
    *output_directory = gadgetmp2::All.OutputDir;
    return 0;
}
int get_viscosity_switch(char **viscosity_switch){
    if (gadgetmp2::ThisTask) {return 0;}
#ifdef MONAGHAN83VISC
    *viscosity_switch = (char *)"Monaghan & Gingold 1983";
#elif MORRIS97VISC
    *viscosity_switch = (char *)"Morris & Monaghan 1997";
#else
    *viscosity_switch = (char *)"standard Gadget2 viscosity";
#endif
    return 0;
}

int set_gadget_output_directory(char *output_directory){
    int length = strlen(output_directory);
    if (length >= MAXLEN_FILENAME)
        return -1;
    strcpy(gadgetmp2::All.OutputDir, output_directory);
#ifdef WIN32
    const char sep[] = "\\";
#else
    const char sep[] = "/";
#endif // WIN32
    if(length > 0) {
        if(gadgetmp2::All.OutputDir[length - 1] != sep[0]) {
            strcat(gadgetmp2::All.OutputDir, sep);
        }
    }
    return 0;
}
int get_nogravity(int *no_gravity_flag){
    if (gadgetmp2::ThisTask) {return 0;}
#ifdef NOGRAVITY
    *no_gravity_flag = 1;
#else
    *no_gravity_flag = 0;
#endif
    return 0;
}
int get_gdgop(int *gadget_cell_opening_flag){
    if (gadgetmp2::ThisTask) {return 0;}
    *gadget_cell_opening_flag = gadgetmp2::All.TypeOfOpeningCriterion;
    return 0;
}
int set_gdgop(int gadget_cell_opening_flag){
    gadgetmp2::All.TypeOfOpeningCriterion = gadget_cell_opening_flag;
    return 0;
}
int get_isotherm(int *isothermal_flag){
    if (gadgetmp2::ThisTask) {return 0;}
#ifdef ISOTHERM_EQS
    *isothermal_flag = 1;
#else
    *isothermal_flag = 0;
#endif
    return 0;
}
int get_eps_is_h(int *eps_is_h_flag){
    if (gadgetmp2::ThisTask) {return 0;}
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) &&  defined(UNEQUALSOFTENINGS)
    *eps_is_h_flag = 1;
#else
    *eps_is_h_flag = 0;
#endif
    return 0;
}
int get_nsmooth(int *nsmooth){
    if (gadgetmp2::ThisTask) {return 0;}
    *nsmooth = gadgetmp2::All.DesNumNgb;
    return 0;
}
int set_nsmooth(int nsmooth){
    gadgetmp2::All.DesNumNgb = nsmooth;
    return 0;
}
int get_bh_tol(double *opening_angle){
    if (gadgetmp2::ThisTask) {return 0;}
    *opening_angle = gadgetmp2::All.ErrTolTheta;
    return 0;
}
int set_bh_tol(double opening_angle){
    gadgetmp2::All.ErrTolTheta = opening_angle;
    return 0;
}
int get_gdgtol(double *gadget_cell_opening_constant){
    if (gadgetmp2::ThisTask) {return 0;}
    *gadget_cell_opening_constant = gadgetmp2::All.ErrTolForceAcc;
    return 0;
}
int set_gdgtol(double gadget_cell_opening_constant){
    gadgetmp2::All.ErrTolForceAcc = gadget_cell_opening_constant;
    return 0;
}
int get_gamma(double *gamma){
    if (gadgetmp2::ThisTask) {return 0;}
    *gamma = gadgetmp2::const_GAMMA.toDouble();
    return 0;
}
int get_alpha(double *artificial_viscosity_alpha){
    if (gadgetmp2::ThisTask) {return 0;}
    *artificial_viscosity_alpha = gadgetmp2::All.ArtBulkViscConst;
    return 0;
}
int set_alpha(double artificial_viscosity_alpha){
    gadgetmp2::All.ArtBulkViscConst = artificial_viscosity_alpha;
    return 0;
}
int get_beta(double *artificial_viscosity_beta){
    if (gadgetmp2::ThisTask) {return 0;}
    *artificial_viscosity_beta = gadgetmp2::All.ArtBulkViscBeta;
    return 0;
}
int set_beta(double artificial_viscosity_beta){
    gadgetmp2::All.ArtBulkViscBeta = artificial_viscosity_beta;
    return 0;
}
int get_courant(double *courant){
    if (gadgetmp2::ThisTask) {return 0;}
    *courant = gadgetmp2::All.CourantFac*2.0;
    return 0;
}
int set_courant(double courant){
    gadgetmp2::All.CourantFac = courant/2.0;
    return 0;
}
int get_nsmtol(double *n_neighbour_tol){
    if (gadgetmp2::ThisTask) {return 0;}
    *n_neighbour_tol = (gadgetmp2::All.MaxNumNgbDeviation / gadgetmp2::All.DesNumNgb);
    return 0;
}
int set_nsmtol(double n_neighbour_tol){
    gadgetmp2::All.MaxNumNgbDeviation = n_neighbour_tol * gadgetmp2::All.DesNumNgb;
    return 0;
}

int get_energy_file(char **energy_file){
    if (gadgetmp2::ThisTask) {return 0;}
    *energy_file = gadgetmp2::All.EnergyFile;
    return 0;
}
int set_energy_file(char *energy_file){
    strcpy(gadgetmp2::All.EnergyFile, energy_file);
    return 0;
}
int get_info_file(char **info_file){
    if (gadgetmp2::ThisTask) {return 0;}
    *info_file = gadgetmp2::All.InfoFile;
    return 0;
}
int set_info_file(char *info_file){
    strcpy(gadgetmp2::All.InfoFile, info_file);
    return 0;
}
int get_timings_file(char **timings_file){
    if (gadgetmp2::ThisTask) {return 0;}
    *timings_file = gadgetmp2::All.TimingsFile;
    return 0;
}
int set_timings_file(char *timings_file){
    strcpy(gadgetmp2::All.TimingsFile, timings_file);
    return 0;
}
int get_cpu_file(char **cpu_file){
    if (gadgetmp2::ThisTask) {return 0;}
    *cpu_file = gadgetmp2::All.CpuFile;
    return 0;
}
int set_cpu_file(char *cpu_file){
    strcpy(gadgetmp2::All.CpuFile, cpu_file);
    return 0;
}

int get_time_limit_cpu(double *time_limit_cpu){
    if (gadgetmp2::ThisTask) {return 0;}
    *time_limit_cpu = gadgetmp2::All.TimeLimitCPU;
    return 0;
}
int set_time_limit_cpu(double time_limit_cpu){
    gadgetmp2::All.TimeLimitCPU = time_limit_cpu;
    return 0;
}
int get_comoving_integration_flag(bool *comoving_integration_flag){
    if (gadgetmp2::ThisTask) {return 0;}
    *comoving_integration_flag = gadgetmp2::All.ComovingIntegrationOn;
    return 0;
}
int set_comoving_integration_flag(bool comoving_integration_flag){
    gadgetmp2::All.ComovingIntegrationOn = comoving_integration_flag;
    return 0;
}
int get_type_of_timestep_criterion(int *type_of_timestep_criterion){
    if (gadgetmp2::ThisTask) {return 0;}
    *type_of_timestep_criterion = gadgetmp2::All.TypeOfTimestepCriterion;
    return 0;
}
int set_type_of_timestep_criterion(int type_of_timestep_criterion){
    gadgetmp2::All.TypeOfTimestepCriterion = type_of_timestep_criterion;
    return 0;
}
int get_begin_time(double *time_begin){
    if (gadgetmp2::ThisTask) {return 0;}
    *time_begin = gadgetmp2::All.TimeBegin;
    return 0;
}
int set_begin_time(double time_begin){
    gadgetmp2::All.TimeBegin = time_begin;
    return 0;
}
int get_time_max(double *time_max){
    if (gadgetmp2::ThisTask) {return 0;}
    *time_max = gadgetmp2::All.TimeMax;
    return 0;
}
int set_time_max(double time_max){
    gadgetmp2::All.TimeMax = time_max;
    return 0;
}
int get_redshift_begin(double *redshift_begin){
    if (gadgetmp2::ThisTask) {return 0;}
    *redshift_begin = redshift_begin_parameter;
    return 0;
}
int set_redshift_begin(double redshift_begin){
    redshift_begin_parameter = redshift_begin;
    return 0;
}
int get_redshift_max(double *redshift_max){
    if (gadgetmp2::ThisTask) {return 0;}
    *redshift_max = redshift_max_parameter;
    return 0;
}
int set_redshift_max(double redshift_max){
    redshift_max_parameter = redshift_max;
    return 0;
}
int get_omega_zero(double *omega_zero){
    if (gadgetmp2::ThisTask) {return 0;}
    *omega_zero = gadgetmp2::All.Omega0;
    return 0;
}
int set_omega_zero(double omega_zero){
    gadgetmp2::All.Omega0 = omega_zero;
    return 0;
}
int get_omega_lambda(double *omega_lambda){
    if (gadgetmp2::ThisTask) {return 0;}
    *omega_lambda = gadgetmp2::All.OmegaLambda;
    return 0;
}
int set_omega_lambda(double omega_lambda){
    gadgetmp2::All.OmegaLambda = omega_lambda;
    return 0;
}
int get_omega_baryon(double *omega_baryon){
    if (gadgetmp2::ThisTask) {return 0;}
    *omega_baryon = gadgetmp2::All.OmegaBaryon;
    return 0;
}
int set_omega_baryon(double omega_baryon){
    gadgetmp2::All.OmegaBaryon = omega_baryon;
    return 0;
}
int get_hubble_param(double *hubble_param){
    if (gadgetmp2::ThisTask) {return 0;}
    *hubble_param = gadgetmp2::All.HubbleParam;
    return 0;
}
int set_hubble_param(double hubble_param){
    gadgetmp2::All.HubbleParam = hubble_param;
    return 0;
}
int get_err_tol_int_accuracy(double *err_tol_int_accuracy){
    if (gadgetmp2::ThisTask) {return 0;}
    *err_tol_int_accuracy = gadgetmp2::All.ErrTolIntAccuracy;
    return 0;
}
int set_err_tol_int_accuracy(double err_tol_int_accuracy){
    gadgetmp2::All.ErrTolIntAccuracy = err_tol_int_accuracy;
    return 0;
}
int get_max_size_timestep(double *max_size_timestep){
    if (gadgetmp2::ThisTask) {return 0;}
    *max_size_timestep = gadgetmp2::All.MaxSizeTimestep;
    return 0;
}
int set_max_size_timestep(double max_size_timestep){
    gadgetmp2::All.MaxSizeTimestep = max_size_timestep;
    return 0;
}
int get_min_size_timestep(double *min_size_timestep){
    if (gadgetmp2::ThisTask) {return 0;}
    *min_size_timestep = gadgetmp2::All.MinSizeTimestep;
    return 0;
}
int set_min_size_timestep(double min_size_timestep){
    gadgetmp2::All.MinSizeTimestep = min_size_timestep;
    return 0;
}
int get_tree_domain_update_frequency(double *tree_domain_update_frequency){
    if (gadgetmp2::ThisTask) {return 0;}
    *tree_domain_update_frequency = gadgetmp2::All.TreeDomainUpdateFrequency;
    return 0;
}
int set_tree_domain_update_frequency(double tree_domain_update_frequency){
    gadgetmp2::All.TreeDomainUpdateFrequency = tree_domain_update_frequency;
    return 0;
}
int get_time_between_statistics(double *time_between_statistics){
    if (gadgetmp2::ThisTask) {return 0;}
    *time_between_statistics = gadgetmp2::All.TimeBetStatistics;
    return 0;
}
int set_time_between_statistics(double time_between_statistics){
    gadgetmp2::All.TimeBetStatistics = time_between_statistics;
    return 0;
}
int get_min_gas_temp(double *min_gas_temp){
    if (gadgetmp2::ThisTask) {return 0;}
    *min_gas_temp = gadgetmp2::All.MinGasTemp;
    return 0;
}
int set_min_gas_temp(double min_gas_temp){
    gadgetmp2::All.MinGasTemp = min_gas_temp;
    return 0;
}
int get_min_gas_hsmooth_fractional(double *min_gas_hsmooth_fractional){
    if (gadgetmp2::ThisTask) {return 0;}
    *min_gas_hsmooth_fractional = gadgetmp2::All.MinGasHsmlFractional;
    return 0;
}
int set_min_gas_hsmooth_fractional(double min_gas_hsmooth_fractional){
    gadgetmp2::All.MinGasHsmlFractional = min_gas_hsmooth_fractional;
    return 0;
}
int get_softening_gas_max_phys(double *softening_gas_max_phys){
    if (gadgetmp2::ThisTask) {return 0;}
    *softening_gas_max_phys = gadgetmp2::All.SofteningGasMaxPhys;
    return 0;
}
int set_softening_gas_max_phys(double softening_gas_max_phys){
    gadgetmp2::All.SofteningGasMaxPhys = softening_gas_max_phys;
    return 0;
}
int get_softening_halo_max_phys(double *softening_halo_max_phys){
    if (gadgetmp2::ThisTask) {return 0;}
    *softening_halo_max_phys = gadgetmp2::All.SofteningHaloMaxPhys;
    return 0;
}
int set_softening_halo_max_phys(double softening_halo_max_phys){
    gadgetmp2::All.SofteningHaloMaxPhys = softening_halo_max_phys;
    return 0;
}

int get_box_size(double *value)
{
    if (gadgetmp2::ThisTask) {return 0;}
    *value = gadgetmp2::All.BoxSize;
    return 0;
}

int set_box_size(double value)
{
    gadgetmp2::All.BoxSize = value;
    return 0;
}

int get_periodic_boundaries_flag(bool *value)
{
    if (gadgetmp2::ThisTask) {return 0;}
    *value = gadgetmp2::All.PeriodicBoundariesOn;
    return 0;
}

int set_periodic_boundaries_flag(bool value)
{
// gadgetmp2::All.PeriodicBoundariesOn is read only because compile time determined
    return -2;
}

int get_interpret_kicks_as_feedback_flag(bool *value)
{
    if (gadgetmp2::ThisTask) {return 0;}
    *value = interpret_kicks_as_feedback;
    return 0;
}

int set_interpret_kicks_as_feedback_flag(bool value)
{
    interpret_kicks_as_feedback = value;
    return 0;
}

int get_interpret_heat_as_feedback_flag(bool *value) {
    if (gadgetmp2::ThisTask) {return 0;}
    *value = interpret_heat_as_feedback;
    return 0;
}

int set_interpret_heat_as_feedback_flag(bool value) {
    interpret_heat_as_feedback = value;
    return 0;
}


// particle property getters/setters: (will only work after commit_particles() is called)

int get_index_of_first_particle(int *index_of_the_particle){
    return get_index_of_next_particle(0, index_of_the_particle);
}
int get_index_of_next_particle(int index_of_the_particle, int *index_of_the_next_particle){
    map<long long, int>::iterator it;
    long long next_local_index = 0;

    if (!particles_initialized)
        return -1;

    if (!particle_map_up_to_date)
        update_particle_map();

    it = local_index_map.lower_bound(index_of_the_particle + 1);
    if (it != local_index_map.end()){
        next_local_index = (*it).first;
    } else {
        next_local_index = index_of_highest_mapped_particle + 1;
    }

    if (gadgetmp2::ThisTask == 0){
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, &next_local_index, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, gadgetmp2::GADGET_WORLD);
#endif
        *index_of_the_next_particle = next_local_index;
        if (next_local_index < index_of_highest_mapped_particle){
            return 0;
        } else if (next_local_index == index_of_highest_mapped_particle){
            return 1;
        } else {
            return -1;
        }
    } else {
#ifndef NOMPI
        MPI_Reduce(&next_local_index, NULL, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, gadgetmp2::GADGET_WORLD);
#endif
        return 0;
    }
}

void update_particle_map(void){
    local_index_map.clear();
    for(int i = 0; i < gadgetmp2::NumPart; i++) {
        local_index_map.insert(std::pair<long long, int>(gadgetmp2::P[i].ID, i));
    }
    particle_map_up_to_date = true;
}
int found_particle(int index_of_the_particle, int *local_index){
    map<long long, int>::iterator it;
    if (!particles_initialized || index_of_the_particle < 1 ||
            index_of_the_particle > index_of_highest_mapped_particle)
        return 0;

    if (!particle_map_up_to_date)
        update_particle_map();

    it = local_index_map.find(index_of_the_particle);
    if (it != local_index_map.end()){
        *local_index = (*it).second;
        return 1;
    }
    return 0;
}

int get_mass(int *index, double *mass, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = gadgetmp2::P[local_index].Mass.toDouble();
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                mass[i] = 0;
            } else
                mass[i] = buffer[i];
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}

int check_counts_and_free(int *count, int length){
    int errors = 0;
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        return 0;
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1)
                errors++;
        }
    }
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}

int set_mass(int *index, double *mass, int length){
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            gadgetmp2::P[local_index].Mass = mass[i];
            count[i] = 1;
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}

int get_radius(int index, double *radius){
    return -2;
}

int get_position_comoving(int *index, double *x, double *y, double *z, int length){
    int errors = 0;
    double *buffer = new double[length*3];
    int *count = new int[length];
    int local_index;
#ifdef PERIODIC
    double boxSize = gadgetmp2::All.BoxSize;
    double boxHalf = 0.5 * gadgetmp2::All.BoxSize;
#endif


    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = gadgetmp2::P[local_index].Pos[0].toDouble();
            buffer[i+length] = gadgetmp2::P[local_index].Pos[1].toDouble();
            buffer[i+2*length] = gadgetmp2::P[local_index].Pos[2].toDouble();
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
#ifdef PERIODIC
        for (int i = 0; i < 3*length; i++){
            if (buffer[i] > boxHalf){
                buffer[i] -= boxSize;
            }
        }
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                x[i] = 0;
                y[i] = 0;
                z[i] = 0;
            } else {
                x[i] = buffer[i];
                y[i] = buffer[i+length];
                z[i] = buffer[i+2*length];
            }
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}
int get_position(int *index, double *x, double *y, double *z, int length){
    int result = get_position_comoving(index, x, y, z, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        for (int i = 0; i < length; i++){
            x[i] *= gadgetmp2::All.Time;
            y[i] *= gadgetmp2::All.Time;
            z[i] *= gadgetmp2::All.Time;
        }
    }
    return result;
}

int set_position_comoving(int *index, double *x, double *y, double *z, int length){
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            gadgetmp2::P[local_index].Pos[0] = x[i];
            gadgetmp2::P[local_index].Pos[1] = y[i];
            gadgetmp2::P[local_index].Pos[2] = z[i];
            count[i] = 1;
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}
int set_position(int *index, double *x, double *y, double *z, int length){
    if(gadgetmp2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmp2::All.Time;
        for (int i = 0; i < length; i++){
            x[i] *= a_inv;
            y[i] *= a_inv;
            z[i] *= a_inv;
        }
    }
    return set_position_comoving(index, x, y, z, length);
}

int get_velocity_gadget_u(int *index, double *vx, double *vy, double *vz, int length){
    int errors = 0;
    double *buffer = new double[length*3];
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = gadgetmp2::P[local_index].Vel[0].toDouble();
            buffer[i+length] = gadgetmp2::P[local_index].Vel[1].toDouble();
            buffer[i+2*length] = gadgetmp2::P[local_index].Vel[2].toDouble();
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                vx[i] = 0;
                vy[i] = 0;
                vz[i] = 0;
            } else {
                vx[i] = buffer[i];
                vy[i] = buffer[i+length];
                vz[i] = buffer[i+2*length];
            }
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}
int get_velocity_comoving(int *index, double *vx, double *vy, double *vz, int length){
    int result = get_velocity_gadget_u(index, vx, vy, vz, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        double a2_inv = 1.0 / (gadgetmp2::All.Time * gadgetmp2::All.Time);
        for (int i = 0; i < length; i++){
            vx[i] *= a2_inv;
            vy[i] *= a2_inv;
            vz[i] *= a2_inv;
        }
    }
    return result;
}
int get_velocity(int *index, double *vx, double *vy, double *vz, int length){
    int result = get_velocity_gadget_u(index, vx, vy, vz, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmp2::All.Time;
        for (int i = 0; i < length; i++){
            vx[i] *= a_inv;
            vy[i] *= a_inv;
            vz[i] *= a_inv;
        }
    }
    return result;
}

int set_velocity_gadget_u(int *index, double *vx, double *vy, double *vz, int length){
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            gadgetmp2::P[local_index].Vel[0] = vx[i];
            gadgetmp2::P[local_index].Vel[1] = vy[i];
            gadgetmp2::P[local_index].Vel[2] = vz[i];
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_kicks_as_feedback && gadgetmp2::P[local_index].Type == 0) {
                gadgetmp2::SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if(interpret_kicks_as_feedback && gadgetmp2::P[local_index].Type == 0 && gadgetmp2::P[local_index].Ti_endstep != gadgetmp2::All.Ti_Current) {
                gadgetmp2::make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}
int set_velocity_comoving(int *index, double *vx, double *vy, double *vz, int length){
    if(gadgetmp2::All.ComovingIntegrationOn) {
        double a2 = (gadgetmp2::All.Time * gadgetmp2::All.Time);
        for (int i = 0; i < length; i++){
            vx[i] *= a2;
            vy[i] *= a2;
            vz[i] *= a2;
        }
    }
    return set_velocity_gadget_u(index, vx, vy, vz, length);
}
int set_velocity(int *index, double *vx, double *vy, double *vz, int length){
    if(gadgetmp2::All.ComovingIntegrationOn) {
        for (int i = 0; i < length; i++){
            vx[i] *= gadgetmp2::All.Time;
            vy[i] *= gadgetmp2::All.Time;
            vz[i] *= gadgetmp2::All.Time;
        }
    }
    return set_velocity_gadget_u(index, vx, vy, vz, length);
}

int get_state_gadget(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *radius, int length) {
    int errors = 0;
    double *buffer = new double[length*8];
    int *count = new int[length];
    int local_index;
#ifdef PERIODIC
    double boxSize = gadgetmp2::All.BoxSize;
    double boxHalf = 0.5 * gadgetmp2::All.BoxSize;
#endif

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = gadgetmp2::P[local_index].Mass.toDouble();
            buffer[i+length] = gadgetmp2::P[local_index].Pos[0].toDouble();
            buffer[i+2*length] = gadgetmp2::P[local_index].Pos[1].toDouble();
            buffer[i+3*length] = gadgetmp2::P[local_index].Pos[2].toDouble();
            buffer[i+4*length] = gadgetmp2::P[local_index].Vel[0].toDouble();
            buffer[i+5*length] = gadgetmp2::P[local_index].Vel[1].toDouble();
            buffer[i+6*length] = gadgetmp2::P[local_index].Vel[2].toDouble();
            buffer[i+7*length] = gadgetmp2::P[local_index].radius.toDouble();
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
            buffer[i+3*length] = 0;
            buffer[i+4*length] = 0;
            buffer[i+5*length] = 0;
            buffer[i+6*length] = 0;
            buffer[i+7*length] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*8, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*8, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
#ifdef PERIODIC
        for (int i = length; i < 4*length; i++){
            if (buffer[i] > boxHalf){
                buffer[i] -= boxSize;
            }
        }
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                mass[i] = 0;
                x[i] = 0;
                y[i] = 0;
                z[i] = 0;
                vx[i] = 0;
                vy[i] = 0;
                vz[i] = 0;
                radius[i] = 0;
            } else {
                mass[i] = buffer[i];
                x[i] = buffer[i+length];
                y[i] = buffer[i+2*length];
                z[i] = buffer[i+3*length];
                vx[i] = buffer[i+4*length];
                vy[i] = buffer[i+5*length];
                vz[i] = buffer[i+6*length];
                radius[i] = buffer[i+7*length];
            }
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}
int get_state_comoving(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length) {
    double* radius;
    int result = get_state_gadget(index, mass, x, y, z, vx, vy, vz, radius, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        double a2_inv = 1.0 / (gadgetmp2::All.Time * gadgetmp2::All.Time);
        for (int i = 0; i < length; i++){
            vx[i] *= a2_inv;
            vy[i] *= a2_inv;
            vz[i] *= a2_inv;
        }
    }
    return result;
}
int get_state(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *radius, int length) {
    int result = get_state_gadget(index, mass, x, y, z, vx, vy, vz, radius, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmp2::All.Time;
        for (int i = 0; i < length; i++){
            x[i] *= gadgetmp2::All.Time;
            y[i] *= gadgetmp2::All.Time;
            z[i] *= gadgetmp2::All.Time;
            vx[i] *= a_inv;
            vy[i] *= a_inv;
            vz[i] *= a_inv;
            radius[1] *= 1;
        }
    }
    return result;
}

int set_state_gadget(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length){
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            gadgetmp2::P[local_index].Mass = mass[i];
            gadgetmp2::P[local_index].Pos[0] = x[i];
            gadgetmp2::P[local_index].Pos[1] = y[i];
            gadgetmp2::P[local_index].Pos[2] = z[i];
            gadgetmp2::P[local_index].Vel[0] = vx[i];
            gadgetmp2::P[local_index].Vel[1] = vy[i];
            gadgetmp2::P[local_index].Vel[2] = vz[i];
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_kicks_as_feedback && gadgetmp2::P[local_index].Type == 0) {
                gadgetmp2::SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if(interpret_kicks_as_feedback && gadgetmp2::P[local_index].Type == 0 && gadgetmp2::P[local_index].Ti_endstep != gadgetmp2::All.Ti_Current) {
                gadgetmp2::make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}
int set_state_comoving(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length){
    if(gadgetmp2::All.ComovingIntegrationOn) {
        double a2 = (gadgetmp2::All.Time * gadgetmp2::All.Time);
        for (int i = 0; i < length; i++){
            vx[i] *= a2;
            vy[i] *= a2;
            vz[i] *= a2;
        }
    }
    return set_state_gadget(index, mass, x, y, z, vx, vy, vz, length);
}
int set_state(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length){
    if(gadgetmp2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmp2::All.Time;
        for (int i = 0; i < length; i++){
            x[i] *= a_inv;
            y[i] *= a_inv;
            z[i] *= a_inv;
            vx[i] *= gadgetmp2::All.Time;
            vy[i] *= gadgetmp2::All.Time;
            vz[i] *= gadgetmp2::All.Time;
        }
    }
    return set_state_gadget(index, mass, x, y, z, vx, vy, vz, length);
}

int get_state_sph_gadget(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *internal_energy, int length) {
    int errors = 0;
    double *buffer = new double[length*8];
    int *count = new int[length];
    int local_index;
#ifdef PERIODIC
    double boxSize = gadgetmp2::All.BoxSize;
    double boxHalf = 0.5 * gadgetmp2::All.BoxSize;
#endif
#ifndef ISOTHERM_EQS
    double a3;

    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }
    if(gadgetmp2::All.ComovingIntegrationOn){a3 = (gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time);}else{a3 = 1;}
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmp2::P[local_index].Mass.toDouble();
            buffer[i+length] = gadgetmp2::P[local_index].Pos[0].toDouble();
            buffer[i+2*length] = gadgetmp2::P[local_index].Pos[1].toDouble();
            buffer[i+3*length] = gadgetmp2::P[local_index].Pos[2].toDouble();
            buffer[i+4*length] = gadgetmp2::P[local_index].Vel[0].toDouble();
            buffer[i+5*length] = gadgetmp2::P[local_index].Vel[1].toDouble();
            buffer[i+6*length] = gadgetmp2::P[local_index].Vel[2].toDouble();
#ifdef ISOTHERM_EQS
            buffer[i+7*length] = gadgetmp2::SphP[local_index].Entropy.toDouble();
#else
            buffer[i+7*length] = (gadgetmp2::SphP[local_index].Entropy *
                pow(gadgetmp2::SphP[local_index].Density / a3, gadgetmp2::const_GAMMA_MINUS1) / gadgetmp2::const_GAMMA_MINUS1).toDouble();
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
            buffer[i+3*length] = 0;
            buffer[i+4*length] = 0;
            buffer[i+5*length] = 0;
            buffer[i+6*length] = 0;
            buffer[i+7*length] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*8, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*8, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
#ifdef PERIODIC
        for (int i = length; i < 4*length; i++){
            if (buffer[i] > boxHalf){
                buffer[i] -= boxSize;
            }
        }
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                mass[i] = 0;
                x[i] = 0;
                y[i] = 0;
                z[i] = 0;
                vx[i] = 0;
                vy[i] = 0;
                vz[i] = 0;
                internal_energy[i] = 0;
            } else {
                mass[i] = buffer[i];
                x[i] = buffer[i+length];
                y[i] = buffer[i+2*length];
                z[i] = buffer[i+3*length];
                vx[i] = buffer[i+4*length];
                vy[i] = buffer[i+5*length];
                vz[i] = buffer[i+6*length];
                internal_energy[i] = buffer[i+7*length];
            }
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}
int get_state_sph(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *internal_energy, int length) {
    int result = get_state_sph_gadget(index, mass, x, y, z, vx, vy, vz, internal_energy, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmp2::All.Time;
        for (int i = 0; i < length; i++){
            x[i] *= gadgetmp2::All.Time;
            y[i] *= gadgetmp2::All.Time;
            z[i] *= gadgetmp2::All.Time;
            vx[i] *= a_inv;
            vy[i] *= a_inv;
            vz[i] *= a_inv;
        }
    }
    return result;
}

int set_state_sph_gadget(int *index, double *mass, double *x, double *y, double *z,
        double *vx, double *vy, double *vz, double *internal_energy, int length){
    int *count = new int[length];
    int local_index;
#ifndef ISOTHERM_EQS
    double a3;

    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }
    if(gadgetmp2::All.ComovingIntegrationOn){a3 = (gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time);}else{a3 = 1;}
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            gadgetmp2::P[local_index].Mass = mass[i];
            gadgetmp2::P[local_index].Pos[0] = x[i];
            gadgetmp2::P[local_index].Pos[1] = y[i];
            gadgetmp2::P[local_index].Pos[2] = z[i];
            gadgetmp2::P[local_index].Vel[0] = vx[i];
            gadgetmp2::P[local_index].Vel[1] = vy[i];
            gadgetmp2::P[local_index].Vel[2] = vz[i];
#ifdef ISOTHERM_EQS
            gadgetmp2::SphP[local_index].Entropy = internal_energy[i];
#else
            gadgetmp2::SphP[local_index].Entropy = gadgetmp2::const_GAMMA_MINUS1 * internal_energy[i] /
                pow(gadgetmp2::SphP[local_index].Density / a3, gadgetmp2::const_GAMMA_MINUS1);
#endif
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_heat_as_feedback || interpret_kicks_as_feedback) {
                gadgetmp2::SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if ((interpret_heat_as_feedback || interpret_kicks_as_feedback) &&
                    gadgetmp2::P[local_index].Ti_endstep != gadgetmp2::All.Ti_Current) {
                gadgetmp2::make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}
int set_state_sph(int *index, double *mass, double *x, double *y, double *z,
        double *vx, double *vy, double *vz, double *internal_energy, int length){
    if(gadgetmp2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmp2::All.Time;
        for (int i = 0; i < length; i++){
            x[i] *= a_inv;
            y[i] *= a_inv;
            z[i] *= a_inv;
            vx[i] *= gadgetmp2::All.Time;
            vy[i] *= gadgetmp2::All.Time;
            vz[i] *= gadgetmp2::All.Time;
        }
    }
    return set_state_sph_gadget(index, mass, x, y, z, vx, vy, vz, internal_energy, length);
}

int get_acceleration_comoving(int *index, double * ax, double * ay, double * az, int length){
    int errors = 0;
    double *buffer = new double[length*3];
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = gadgetmp2::P[local_index].GravAccel[0].toDouble();
            buffer[i+length] = gadgetmp2::P[local_index].GravAccel[1].toDouble();
            buffer[i+2*length] = gadgetmp2::P[local_index].GravAccel[2].toDouble();
            if(gadgetmp2::P[local_index].Type == 0){
                buffer[i] += gadgetmp2::SphP[local_index].HydroAccel[0].toDouble();
                buffer[i+length] += gadgetmp2::SphP[local_index].HydroAccel[1].toDouble();
                buffer[i+2*length] += gadgetmp2::SphP[local_index].HydroAccel[2].toDouble();
            }
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                ax[i] = 0;
                ay[i] = 0;
                az[i] = 0;
            } else {
                ax[i] = buffer[i];
                ay[i] = buffer[i+length];
                az[i] = buffer[i+2*length];
            }
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}
int get_acceleration(int *index, double * ax, double * ay, double * az, int length){
    int result = get_acceleration_comoving(index, ax, ay, az, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        for (int i = 0; i < length; i++){
            ax[i] *= gadgetmp2::All.Time;
            ay[i] *= gadgetmp2::All.Time;
            az[i] *= gadgetmp2::All.Time;
        }
    }
    return result;
}

int set_acceleration(int index, double ax, double ay, double az){
    return -2;
}

int get_internal_energy(int *index, double *internal_energy, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;
#ifndef ISOTHERM_EQS
    double a3;

    if(gadgetmp2::All.ComovingIntegrationOn){a3 = (gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time);}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            count[i] = 1;
#ifdef ISOTHERM_EQS
            buffer[i] = gadgetmp2::SphP[local_index].Entropy.toDouble();
#else
            buffer[i] = (gadgetmp2::SphP[local_index].Entropy *
                pow(gadgetmp2::SphP[local_index].Density / a3, gadgetmp2::const_GAMMA_MINUS1) / gadgetmp2::const_GAMMA_MINUS1).toDouble();
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                internal_energy[i] = 0;
            } else
                internal_energy[i] = buffer[i];
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}

int set_internal_energy(int *index, double *internal_energy, int length){
    int *count = new int[length];
    int local_index;
#ifndef ISOTHERM_EQS
    double a3;

    if(gadgetmp2::All.ComovingIntegrationOn){a3 = (gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time);}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
#ifdef ISOTHERM_EQS
            gadgetmp2::SphP[local_index].Entropy = internal_energy[i];
#else
            gadgetmp2::SphP[local_index].Entropy = gadgetmp2::const_GAMMA_MINUS1 * internal_energy[i] /
                pow(gadgetmp2::SphP[local_index].Density / a3, gadgetmp2::const_GAMMA_MINUS1);
#endif
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_heat_as_feedback) {
                gadgetmp2::SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if(interpret_heat_as_feedback && gadgetmp2::P[local_index].Ti_endstep != gadgetmp2::All.Ti_Current) {
                gadgetmp2::make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}

int get_smoothing_length_comoving(int *index, double *smoothing_length, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmp2::SphP[local_index].Hsml.toDouble();
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                smoothing_length[i] = 0;
            } else
                smoothing_length[i] = buffer[i];
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}
int get_smoothing_length(int *index, double *smoothing_length, int length){
    int result = get_smoothing_length_comoving(index, smoothing_length, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        for (int i = 0; i < length; i++){
            smoothing_length[i] *= gadgetmp2::All.Time;
        }
    }
    return result;
}


int get_alpha_visc(int *index, double *alpha_visc, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            count[i] = 1;
#ifdef MORRIS97VISC
            buffer[i] = gadgetmp2::SphP[local_index].Alpha.toDouble();
#else
	    buffer[i] = gadgetmp2::All.ArtBulkViscConst;
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                alpha_visc[i] = 0;
            } else
                alpha_visc[i] = buffer[i];
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}

int get_dalphadt_visc(int *index, double *dalphadt_visc, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            count[i] = 1;
#ifdef MORRIS97VISC
            buffer[i] = gadgetmp2::SphP[local_index].DAlphaDt;
#else
            buffer[i] = 0;
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                dalphadt_visc[i] = 0;
            } else
                dalphadt_visc[i] = buffer[i];
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}




int get_density_comoving(int *index, double *density_out, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;
    double a3;

    if(gadgetmp2::All.ComovingIntegrationOn){a3 = (gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time);}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmp2::SphP[local_index].Density.toDouble() / a3;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                density_out[i] = 0;
            } else
                density_out[i] = buffer[i];
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}
int get_density(int *index, double *density_out, int length){
    int result = get_density_comoving(index, density_out, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        double a3_inv = 1.0 / (gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time);
        for (int i = 0; i < length; i++){
            density_out[i] *= a3_inv;
        }
    }
    return result;
}

int get_pressure_comoving(int *index, double *pressure_out, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;
    double a;

    if(gadgetmp2::All.ComovingIntegrationOn){a = gadgetmp2::All.Time;}else{a = 1;}
    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmp2::SphP[local_index].Pressure.toDouble() / a;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                pressure_out[i] = 0;
            } else
                pressure_out[i] = buffer[i];
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}
int get_pressure(int *index, double *pressure_out, int length){
    int result = get_pressure_comoving(index, pressure_out, length);
    if(gadgetmp2::ThisTask == 0 && gadgetmp2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmp2::All.Time;
        for (int i = 0; i < length; i++){
            pressure_out[i] *= a_inv;
        }
    }
    return result;
}

int get_d_internal_energy_dt(int *index, double *d_internal_energy_dt_out, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    double hubble;
    if (gadgetmp2::All.ComovingIntegrationOn){
        hubble = (gadgetmp2::All.Hubble * std::sqrt(gadgetmp2::All.Omega0 / (gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time)
            + (1 - gadgetmp2::All.Omega0 - gadgetmp2::All.OmegaLambda) / (gadgetmp2::All.Time * gadgetmp2::All.Time) + gadgetmp2::All.OmegaLambda));
    } else {
        hubble = 1;
    }

#ifndef ISOTHERM_EQS
    //~double a3;
    //~if(gadgetmp2::All.ComovingIntegrationOn){a3 = gadgetmp2::All.Time * gadgetmp2::All.Time * gadgetmp2::All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            count[i] = 1;
#ifdef ISOTHERM_EQS
            buffer[i] = gadgetmp2::SphP[local_index].DtEntropy * hubble;
#else
            buffer[i] = (- gadgetmp2::SphP[local_index].Pressure * gadgetmp2::SphP[local_index].DivVel /
                gadgetmp2::SphP[local_index].Density * hubble).toDouble();
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                d_internal_energy_dt_out[i] = 0;
            } else
                d_internal_energy_dt_out[i] = buffer[i];
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}

int get_n_neighbours(int *index, double *n_neighbours, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    if (!density_up_to_date){
        gadgetmp2::density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmp2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmp2::SphP[local_index].NumNgb.toDouble();
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        for (int i = 0; i < length; i++){
            if (count[i] != 1){
                errors++;
                n_neighbours[i] = 0;
            } else
                n_neighbours[i] = buffer[i];
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors){
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}
int get_epsilon_dm_part(int *index, double *epsilon, int length){
    gadgetmp2::set_softenings();
    if (gadgetmp2::ThisTask) {return 0;}
    double a;
    if (gadgetmp2::All.ComovingIntegrationOn) {a = gadgetmp2::All.Time;} else {a = 1;}
    for (int i = 0; i < length; i++){
        epsilon[i] = a * gadgetmp2::All.SofteningTable[1];
    }
    return 0;
}
int get_epsilon_gas_part(int *index, double *epsilon, int length){
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) &&  defined(UNEQUALSOFTENINGS)
    return get_smoothing_length(index, epsilon, length);
#else
    gadgetmp2::set_softenings();
    if (gadgetmp2::ThisTask) {return 0;}
    double a;
    if (gadgetmp2::All.ComovingIntegrationOn) {a = gadgetmp2::All.Time;} else {a = 1;}
    for (int i = 0; i < length; i++){
        epsilon[i] = a * gadgetmp2::All.SofteningTable[0];
    }
    return 0;
#endif
}



// simulation property getters:

void update_global_quantities(bool do_potential){
    if (do_potential) {
        gadgetmp2::compute_potential();
        potential_energy_also_up_to_date = true;
    } else {potential_energy_also_up_to_date = false;}
    gadgetmp2::compute_global_quantities_of_system();
    global_quantities_of_system_up_to_date = true;
}
int get_time(double *time){
    if (gadgetmp2::ThisTask) {return 0;}
    if (gadgetmp2::All.ComovingIntegrationOn) {return -9;}
    *time = gadgetmp2::All.Time;
    return 0;
}


int get_redshift(double *redshift){
    if (gadgetmp2::ThisTask) {return 0;}
    if (!gadgetmp2::All.ComovingIntegrationOn) {return -9;}
    *redshift = 1.0 / gadgetmp2::All.Time - 1.0;
    return 0;
}
int get_total_radius(double *radius){
    double r_squared, local_max = 0;
    int i, j;

    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    for (i = 0; i < gadgetmp2::NumPart; i++){
        for (r_squared = 0, j = 0; j < 3; j++)
            r_squared += ((gadgetmp2::SysState.CenterOfMass[j]-gadgetmp2::P[i].Pos[j])*(gadgetmp2::SysState.CenterOfMass[j]-gadgetmp2::P[i].Pos[j])).toDouble();
        if (r_squared > local_max)
            local_max = r_squared;
    }

    if(gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(&local_max, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, &local_max, 1, MPI_DOUBLE, MPI_MAX, 0, gadgetmp2::GADGET_WORLD);
#endif
        if (gadgetmp2::All.ComovingIntegrationOn){
            *radius = gadgetmp2::All.Time * std::sqrt(local_max);
        } else {
            *radius = std::sqrt(local_max);
        }
    }
    return 0;
}

int get_total_mass(double *mass){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (gadgetmp2::ThisTask) {return 0;}
    *mass = gadgetmp2::SysState.Mass.toDouble();
    return 0;
}

int get_potential(int *index, double *potential, int length) {
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    if (!potential_energy_also_up_to_date) {
        gadgetmp2::compute_potential();
        potential_energy_also_up_to_date = true;
    }

    for (int i = 0; i < length; i++) {
        if (found_particle(index[i], &local_index)) {
            count[i] = 1;
            buffer[i] = gadgetmp2::P[local_index].Potential.toDouble();
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }

    if (gadgetmp2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmp2::GADGET_WORLD);
#endif
        double a2;
        if (gadgetmp2::All.ComovingIntegrationOn) {a2 = (gadgetmp2::All.Time * gadgetmp2::All.Time);} else {a2 = 1;}
        for (int i = 0; i < length; i++) {
            if (count[i] != 1){
                errors++;
                potential[i] = 0;
            } else {
                potential[i] = a2 * buffer[i];
            }
        }
    }
    delete[] buffer;
    delete[] count;
    if (errors) {
        cout << "Number of particles not found: " << errors << endl;
        return -3;
    }
    return 0;
}

int get_kinetic_energy(double *kinetic_energy){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (gadgetmp2::ThisTask) {return 0;}
    *kinetic_energy = gadgetmp2::SysState.EnergyKin.toDouble();
    return 0;
}
int get_potential_energy(double *potential_energy){
    if (!(global_quantities_of_system_up_to_date && potential_energy_also_up_to_date))
        update_global_quantities(true);
    if (gadgetmp2::ThisTask) {return 0;}
    *potential_energy = gadgetmp2::SysState.EnergyPot.toDouble();
    return 0;
}
int get_thermal_energy(double *thermal_energy){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (gadgetmp2::ThisTask) {return 0;}
    *thermal_energy = gadgetmp2::SysState.EnergyInt.toDouble();
    return 0;
}
int get_number_of_particles(int *number_of_particles){
    if (gadgetmp2::ThisTask) {return 0;}
    *number_of_particles = gadgetmp2::All.TotNumPart;
    return 0;
}
int get_center_of_mass_position(double *x, double *y, double *z){
#ifdef PERIODIC
    return -2;
#endif
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (gadgetmp2::ThisTask) {return 0;}
    double a;
    if (gadgetmp2::All.ComovingIntegrationOn) {a = gadgetmp2::All.Time;} else {a = 1;}
    *x = a * gadgetmp2::SysState.CenterOfMass[0].toDouble();
    *y = a * gadgetmp2::SysState.CenterOfMass[1].toDouble();
    *z = a * gadgetmp2::SysState.CenterOfMass[2].toDouble();
    return 0;
}
int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (gadgetmp2::ThisTask) {return 0;}
    double a_inv;
    if (gadgetmp2::All.ComovingIntegrationOn) {a_inv = 1.0 / gadgetmp2::All.Time;} else {a_inv = 1;}
    *vx = a_inv * (gadgetmp2::SysState.Momentum[0]/gadgetmp2::SysState.Mass).toDouble();
    *vy = a_inv * (gadgetmp2::SysState.Momentum[1]/gadgetmp2::SysState.Mass).toDouble();
    *vz = a_inv * (gadgetmp2::SysState.Momentum[2]/gadgetmp2::SysState.Mass).toDouble();
    return 0;
}
int get_gravity_at_point(double eps, double x, double y, double z, double *forcex, double *forcey, double *forcez){
    return -2;
}
int get_potential_at_point(double eps, double x, double y, double z, double * phi){
    return -2;
}
int get_hydro_state_at_point(double x, double y, double z, double vx, double vy, double vz,
        double * rho, double * rhovx, double * rhovy, double * rhovz, double * rhoe){
    my_float pos[3], vel[3];
    my_float h_out, ngb_out, dhsml_out, rho_out, rhov_out[3], rhov2_out, rhoe_out;
    int error;
    double a, a_inv, a3_inv, a4_inv, a5_inv;
    if (gadgetmp2::All.ComovingIntegrationOn) {
        a = gadgetmp2::All.Time;
        a_inv = 1.0 / gadgetmp2::All.Time;
        a3_inv = a_inv * a_inv * a_inv;
        a4_inv = a3_inv * a_inv;
        a5_inv = a4_inv * a_inv;
    } else {
        a = a_inv = a3_inv = a4_inv = a5_inv = 1;
    }

    error = construct_tree_if_needed();
    if (error) {return error;}

    pos[0] = a_inv * x;
    pos[1] = a_inv * y;
    pos[2] = a_inv * z;
#ifdef PERIODIC
    for (int i = 0; i < 3; i++){
        if (pos[i] < 0.0){
            pos[i] += gadgetmp2::All.BoxSize;
        }
    }
#endif
    vel[0] = a * vx;
    vel[1] = a * vy;
    vel[2] = a * vz;
    gadgetmp2::hydro_state_at_point(pos, vel, &h_out, &ngb_out, &dhsml_out, &rho_out, rhov_out, &rhov2_out, &rhoe_out);
    if (gadgetmp2::ThisTask) {return 0;}
    *rho   = rho_out.toDouble() * a3_inv;
    *rhovx = rhov_out[0].toDouble() * a4_inv;
    *rhovy = rhov_out[1].toDouble() * a4_inv;
    *rhovz = rhov_out[2].toDouble() * a4_inv;
#ifdef ISOTHERM_EQS
    *rhoe = a3_inv * rhoe_out + a5_inv * 0.5*(rhov_out[0]*rhov_out[0] + rhov_out[1]*rhov_out[1] + rhov_out[2]*rhov_out[2]) / rho_out;
#else
    *rhoe = (a3_inv * rhoe_out * (pow(rho_out * a3_inv, gadgetmp2::const_GAMMA_MINUS1) / gadgetmp2::const_GAMMA_MINUS1) + a5_inv * 0.5*rhov2_out).toDouble();
#endif
    return 0;
}
// Word-length, numBits in mantissa
int set_word_length(int mynumBits) {
    numBits = mynumBits;
    size_t i;
    my_float::set_default_prec(numBits);
    numDigits = (int)abs(log10( pow("2.0", -numBits) )).toLong();
    if (gadgetmp2::CommBuffer!= nullptr) gadgetmp2::allocate_commbuffers();
    if (gadgetmp2::P!=nullptr)
    {
    for (i=0; i<gadgetmp2::All.MaxPart; i++)
    gadgetmp2::P[i].change_prec();
    }
    if (gadgetmp2::SphP!=nullptr)
    {
    for (i=0; i<gadgetmp2::All.MaxPartSph; i++)
    gadgetmp2::SphP[i].change_prec();
    }
    return 0;
}
int get_word_length(int *mynumBits) {
    numBits=my_float::get_default_prec();
    *mynumBits = numBits;
    return 0;
}

int get_total_energy_string( char **ep) {
    if (gadgetmp2::ThisTask == 0)
    {
    arg1 = gadgetmp2::SysState.EnergyTot.toString();
    *ep =(char*)arg1.c_str();
    }
    return 0;
}

int get_state_string(int id, char** m, char** x, char** y, char** z, char** vx, char** vy, char** vz, char** radius) {
    int retval=0;
    bool flag=false;
    int local_index;
    my_float_buff* arg=my_float_buff::place_buffer(8, gadgetmp2::CommBuffer);
    #ifndef NOMPI
    MPI_Status status;
    #endif
    if(found_particle(id, &local_index)){
        retval++;
        flag=true;
        arg[0]=gadgetmp2::P[local_index].Mass;
        arg[1]=gadgetmp2::P[local_index].Pos[0];
        arg[2]=gadgetmp2::P[local_index].Pos[1];
        arg[3]=gadgetmp2::P[local_index].Pos[2];
        arg[4]=gadgetmp2::P[local_index].Vel[0];
        arg[5]=gadgetmp2::P[local_index].Vel[1];
        arg[6]=gadgetmp2::P[local_index].Vel[2];
        arg[7]=gadgetmp2::P[local_index].radius;
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    if(retval == 1)
    {
        if (gadgetmp2::ThisTask == 0)
        {
            if(flag == false)
            {
                #ifndef NOMPI
                MPI_Recv(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,MPI_ANY_SOURCE,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD,&status);
                my_float_buff::re_org_buff(arg);
                #endif
            }
            arg0=arg[0].toString();
            arg1=arg[1].toString();
            arg2=arg[2].toString();
            arg3=arg[3].toString();
            arg4=arg[4].toString();
            arg5=arg[5].toString();
            arg6=arg[6].toString();
            arg7=arg[7].toString();
            *m = (char*) arg0.c_str();
            *x = (char*) arg1.c_str();
            *y = (char*) arg2.c_str();
            *z = (char*) arg3.c_str();
            *vx = (char*) arg4.c_str();
            *vy = (char*) arg5.c_str();
            *vz = (char*) arg6.c_str();
            *radius = (char*) arg7.c_str();
        }else{
            #ifndef NOMPI
            if(flag == true)
            {
                MPI_Send(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,0,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD);
            }
            #endif
        }
    }
    return retval-1;
}

int get_position_string(int id, char **x, char **y, char **z) {
    int retval=0;
    bool flag=false;
    int local_index;
    my_float_buff* arg=my_float_buff::place_buffer(8, gadgetmp2::CommBuffer);
    #ifndef NOMPI
    MPI_Status status;
    #endif
    if(found_particle(id, &local_index)){
        retval++;
        flag=true;
        arg[0]=gadgetmp2::P[local_index].Pos[0];
        arg[1]=gadgetmp2::P[local_index].Pos[1];
        arg[2]=gadgetmp2::P[local_index].Pos[2];
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    if(retval == 1)
    {
        if (gadgetmp2::ThisTask == 0)
        {
            if(flag == false)
            {
                #ifndef NOMPI
                MPI_Recv(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,MPI_ANY_SOURCE,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD,&status);
                my_float_buff::re_org_buff(arg);
                #endif
            }
            arg0=arg[0].toString();
            arg1=arg[1].toString();
            arg2=arg[2].toString();
            *x = (char*) arg0.c_str();
            *y = (char*) arg1.c_str();
            *z = (char*) arg2.c_str();
        }else{
            #ifndef NOMPI
            if(flag == true)
            {
                MPI_Send(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,0,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD);
            }
            #endif
        }
    }
    return retval-1;
}

int get_velocity_string(int id, char **vx, char **vy, char **vz) {
    int retval=0;
    bool flag=false;
    int local_index;
    my_float_buff* arg=my_float_buff::place_buffer(8, gadgetmp2::CommBuffer);
    #ifndef NOMPI
    MPI_Status status;
    #endif
    if(found_particle(id, &local_index)){
        retval++;
        flag=true;
        arg[0]=gadgetmp2::P[local_index].Vel[0];
        arg[1]=gadgetmp2::P[local_index].Vel[1];
        arg[2]=gadgetmp2::P[local_index].Vel[2];
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    if(retval == 1)
    {
        if (gadgetmp2::ThisTask == 0)
        {
            if(flag == false)
            {
                #ifndef NOMPI
                MPI_Recv(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,MPI_ANY_SOURCE,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD,&status);
                my_float_buff::re_org_buff(arg);
                #endif
            }
            arg0=arg[0].toString();
            arg1=arg[1].toString();
            arg2=arg[2].toString();
            *vx = (char*) arg0.c_str();
            *vy = (char*) arg1.c_str();
            *vz = (char*) arg2.c_str();
        }else{
            #ifndef NOMPI
            if(flag == true)
            {
                MPI_Send(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,0,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD);
            }
            #endif
        }
    }
    return retval-1;
}

int get_mass_string(int id, char **mass) {
    int retval=0;
    bool flag=false;
    int local_index;
    my_float_buff* arg=my_float_buff::place_buffer(1, gadgetmp2::CommBuffer);
    #ifndef NOMPI
    MPI_Status status;
    #endif
    if(found_particle(id, &local_index)){
        retval++;
        flag=true;
        arg[0]=gadgetmp2::P[local_index].Mass;
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    if(retval == 1)
    {
        if (gadgetmp2::ThisTask == 0)
        {
            if(flag == false)
            {
                #ifndef NOMPI
                MPI_Recv(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,MPI_ANY_SOURCE,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD,&status);
                my_float_buff::re_org_buff(arg);
                #endif
            }
            arg0=arg[0].toString();
            *mass = (char*) arg0.c_str();
        }else{
            #ifndef NOMPI
            if(flag == true)
            {
                MPI_Send(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,0,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD);
            }
            #endif
        }
    }
    return retval-1;
}

int get_radius_string(int id, char** radius){
    int retval=0;
    bool flag=false;
    int local_index;
    my_float_buff* arg=my_float_buff::place_buffer(1, gadgetmp2::CommBuffer);
    #ifndef NOMPI
    MPI_Status status;
    #endif
    if(found_particle(id, &local_index)){
        retval++;
        flag=true;
        arg[0]=gadgetmp2::P[local_index].radius;
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    if(retval == 1)
    {
        if (gadgetmp2::ThisTask == 0)
        {
            if(flag == false)
            {
                #ifndef NOMPI
                MPI_Recv(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,MPI_ANY_SOURCE,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD,&status);
                my_float_buff::re_org_buff(arg);
                #endif
            }
            arg0=arg[0].toString();
            *radius = (char*) arg0.c_str();
        }else{
            #ifndef NOMPI
            if(flag == true)
            {
                MPI_Send(gadgetmp2::CommBuffer,my_float_buff::get_buff_size(arg),MPI_BYTE,0,TAG_AMUSE_SYNC, gadgetmp2::GADGET_WORLD);
            }
            #endif
        }
    }
    return retval-1;
}

int set_state_string(int id, char* m, char* x, char* y, char* z, char* vx, char* vy, char* vz, char* radius) {
    int retval=0;
    int local_index;
    if(found_particle(id, &local_index)){
        retval++;
        gadgetmp2::P[local_index].Mass=m;
        gadgetmp2::P[local_index].Pos[0]=x;
        gadgetmp2::P[local_index].Pos[1]=y;
        gadgetmp2::P[local_index].Pos[2]=z;
        gadgetmp2::P[local_index].Vel[0]=vx;
        gadgetmp2::P[local_index].Vel[1]=vy;
        gadgetmp2::P[local_index].Vel[2]=vz;
        gadgetmp2::P[local_index].radius=radius;
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    return retval-1;
}

int set_mass_string(int id, char* m) {
    int retval=0;
    int local_index;
    if(found_particle(id, &local_index)){
        retval++;
        gadgetmp2::P[local_index].Mass=m;
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    return retval-1;
}

int set_radius_string(int id, char* radius) {
    int retval=0;
    int local_index;
    if(found_particle(id, &local_index)){
        retval++;
        gadgetmp2::P[local_index].radius=radius;
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    return retval-1;
}

int set_radius(int id, double radius) {
    int retval=0;
    int local_index;
    if(found_particle(id, &local_index)){
        retval++;
        gadgetmp2::P[local_index].radius=radius;
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    return retval-1;
}

int set_position_string(int id,char* x, char* y, char* z) {
    int retval=0;
    int local_index;
    if(found_particle(id, &local_index)){
        retval++;
        gadgetmp2::P[local_index].Pos[0]=x;
        gadgetmp2::P[local_index].Pos[1]=y;
        gadgetmp2::P[local_index].Pos[2]=z;
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    return retval-1;
}

int set_velocity_string(int id, char* vx, char* vy, char* vz) {
    int retval=0;
    int local_index;
    if(found_particle(id, &local_index)){
        retval++;
        gadgetmp2::P[local_index].Vel[0]=vx;
        gadgetmp2::P[local_index].Vel[1]=vy;
        gadgetmp2::P[local_index].Vel[2]=vz;
    };
    #ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &retval, 1, MPI_INT, MPI_SUM, gadgetmp2::GADGET_WORLD);
    #endif
    return retval-1;
}

int new_dm_particle_string(int *id, char* mass, char* x, char* y, char* z, char* vx, char* vy, char* vz, char* radius){
    particle_id_counter++;
    if (gadgetmp2::ThisTask == 0)
        *id = particle_id_counter;
    // Divide the particles equally over all Tasks, Gadget will redistribute them later.
    if (gadgetmp2::ThisTask == (dm_particles_in_buffer % gadgetmp2::NTask)){
        dynamics_state state;
        state.mass = mass;
        state.x = x;
        state.y = y;
        state.z = z;
        state.vx = vx;
        state.vy = vy;
        state.vz = vz;
        dm_states.insert(std::pair<long long, dynamics_state>(particle_id_counter, state));
    }
    dm_particles_in_buffer++;
    return 0;
}
