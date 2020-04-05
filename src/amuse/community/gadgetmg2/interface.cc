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

using namespace std;

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
gadgetmg2 a;

// general interface functions:

void set_default_parameters(){
    // parameters that can be changed from AMUSE
    gadgetmg2::All.TimeLimitCPU = 36000;
    gadgetmg2::All.ComovingIntegrationOn = 0;
    gadgetmg2::All.TypeOfTimestepCriterion = 0;
#ifdef PERIODIC
    gadgetmg2::All.PeriodicBoundariesOn = 1;
#else
    gadgetmg2::All.PeriodicBoundariesOn = 0;
#endif
    gadgetmg2::All.Time = 0.0;
    gadgetmg2::All.TimeBegin = 0.0;
    gadgetmg2::All.TimeMax = 100.0;
    gadgetmg2::All.Omega0 = 0;
    gadgetmg2::All.OmegaLambda = 0;
    gadgetmg2::All.OmegaBaryon = 0;
    gadgetmg2::All.HubbleParam = 0.7;
    gadgetmg2::All.BoxSize = 1.0;
    gadgetmg2::All.TimeBetStatistics = 0.1;
    gadgetmg2::All.ErrTolIntAccuracy = 0.025;
    gadgetmg2::All.CourantFac = 0.15;
    gadgetmg2::All.MaxSizeTimestep = 0.01;
    gadgetmg2::All.MinSizeTimestep = 0.0;
    gadgetmg2::All.ErrTolTheta = 0.5;
    gadgetmg2::All.TypeOfOpeningCriterion = 1;
    gadgetmg2::All.ErrTolForceAcc = 0.005;
    gadgetmg2::All.TreeDomainUpdateFrequency = 0.05;
    gadgetmg2::All.DesNumNgb = 50;
    gadgetmg2::All.MaxNumNgbDeviation = 5.;
    gadgetmg2::All.ArtBulkViscConst = 0.5;
    gadgetmg2::All.ArtBulkViscBeta = 1.;
    gadgetmg2::All.MinGasTemp = 0;
    gadgetmg2::All.UnitLength_in_cm = 3.085678e21;
    gadgetmg2::All.UnitMass_in_g = 1.989e43;
    gadgetmg2::All.UnitVelocity_in_cm_per_s = 1e5;
    gadgetmg2::All.UnitTime_in_s = gadgetmg2::All.UnitLength_in_cm / gadgetmg2::All.UnitVelocity_in_cm_per_s;
    gadgetmg2::All.MinGasHsmlFractional = 0.0;
    gadgetmg2::All.SofteningGas = 0.01;
    gadgetmg2::All.SofteningHalo = 0.01;
    gadgetmg2::All.SofteningGasMaxPhys = 0.0;
    gadgetmg2::All.SofteningHaloMaxPhys = 0.0;
    gadgetmg2::set_softenings();
    strcpy(gadgetmg2::All.OutputDir,   ".");
    strcpy(gadgetmg2::All.EnergyFile,  "energy.txt");
    strcpy(gadgetmg2::All.InfoFile,    "info.txt");
    strcpy(gadgetmg2::All.TimingsFile, "timings.txt");
    strcpy(gadgetmg2::All.CpuFile,     "cpu.txt");

    // parameters that are fixed for AMUSE:
    gadgetmg2::All.PartAllocFactor = 1.5; // Memory allocation parameter
    gadgetmg2::All.TreeAllocFactor = 0.8; // Memory allocation parameter
    gadgetmg2::All.BufferSize = 25;       // Memory allocation parameter
    gadgetmg2::All.ResubmitOn = 0;              // Keep this turned off!
    gadgetmg2::All.OutputListOn = 0;            // Keep this turned off!
    gadgetmg2::All.GravityConstantInternal = 0; // Keep this turned off!

    // parameters that are unused for AMUSE:
    strcpy(gadgetmg2::All.InitCondFile, "");
    strcpy(gadgetmg2::All.RestartFile, "");
    strcpy(gadgetmg2::All.SnapshotFileBase, "");
    strcpy(gadgetmg2::All.OutputListFilename, "");
    strcpy(gadgetmg2::All.ResubmitCommand, "");
    gadgetmg2::All.ICFormat = 1;
    gadgetmg2::All.SnapFormat = 1;
    gadgetmg2::All.TimeBetSnapshot = 100.0;
    gadgetmg2::All.TimeOfFirstSnapshot = 100.0;
    gadgetmg2::All.CpuTimeBetRestartFile = 36000.0;
    gadgetmg2::All.NumFilesPerSnapshot = 1;
    gadgetmg2::All.NumFilesWrittenInParallel = 1;
    gadgetmg2::All.InitGasTemp = 0;
    gadgetmg2::All.MaxRMSDisplacementFac = 0.2; // parameter for PM; PM is currently not supported
}

int initialize_code(){
    double t0, t1;
#ifndef NOMPI
    get_comm_world(&gadgetmg2::GADGET_WORLD);
    MPI_Comm_rank(gadgetmg2::GADGET_WORLD, &gadgetmg2::ThisTask);
    MPI_Comm_size(gadgetmg2::GADGET_WORLD, &gadgetmg2::NTask);
#else
    gadgetmg2::ThisTask = 0;
    gadgetmg2::NTask = 1;
#endif
    for(gadgetmg2::PTask = 0; gadgetmg2::NTask > (1 << gadgetmg2::PTask); gadgetmg2::PTask++);
    gadgetmg2::RestartFlag = gadgetmg2::All.TotNumPart = gadgetmg2::All.TotN_gas = 0;
    gadgetmg2::All.CPU_TreeConstruction = gadgetmg2::All.CPU_TreeWalk = gadgetmg2::All.CPU_Gravity = gadgetmg2::All.CPU_Potential = gadgetmg2::All.CPU_Domain =
        gadgetmg2::All.CPU_Snapshot = gadgetmg2::All.CPU_Total = gadgetmg2::All.CPU_CommSum = gadgetmg2::All.CPU_Imbalance = gadgetmg2::All.CPU_Hydro =
        gadgetmg2::All.CPU_HydCompWalk = gadgetmg2::All.CPU_HydCommSumm = gadgetmg2::All.CPU_HydImbalance =
        gadgetmg2::All.CPU_EnsureNgb = gadgetmg2::All.CPU_Predict = gadgetmg2::All.CPU_TimeLine = gadgetmg2::All.CPU_PM = gadgetmg2::All.CPU_Peano = 0;
    gadgetmg2::CPUThisRun = 0;
    t0 = gadgetmg2::second();
    if(gadgetmg2::ThisTask == 0){
        printf("\nThis is Gadget, version `%s'.\n", GADGETVERSION);
        printf("\nRunning on %d processors.\n", gadgetmg2::NTask);
    }

    t1 = gadgetmg2::second();
    gadgetmg2::CPUThisRun += gadgetmg2::timediff(t0, t1);
    gadgetmg2::All.CPU_Total += gadgetmg2::timediff(t0, t1);

    //AMUSE STOPPING CONDITIONS SUPPORT
    set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
    set_support_for_condition(DENSITY_LIMIT_DETECTION);
    set_support_for_condition(INTERNAL_ENERGY_LIMIT_DETECTION);
    mpi_setup_stopping_conditions();

    set_default_parameters();
    return 0;
}

int cleanup_code(){
    if (outfiles_opened)
        gadgetmg2::close_outputfiles();
    if (particles_initialized){
        gadgetmg2::free_memory();
        gadgetmg2::ngb_treefree();
        gadgetmg2::force_treefree();
    }
    return 0;
}


int check_parameters(){
#ifndef NOMPI
    MPI_Bcast(&gadgetmg2::All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, gadgetmg2::GADGET_WORLD);
#endif
    if (gadgetmg2::ThisTask){
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
    if(gadgetmg2::All.NumFilesWrittenInParallel < 1){
        printf("NumFilesWrittenInParallel MUST be at least 1\n");
        return -4;
    }
    if(gadgetmg2::All.NumFilesWrittenInParallel > gadgetmg2::NTask){
        printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
        return -4;
    }
#ifdef PERIODIC
    if(gadgetmg2::All.PeriodicBoundariesOn == 0){
        printf("Code was compiled with periodic boundary conditions switched on.\n");
        printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
        return -4;
    }
#else
    if(gadgetmg2::All.PeriodicBoundariesOn == 1){
        printf("Code was compiled with periodic boundary conditions switched off.\n");
        printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
        return -4;
    }
#endif
    if(gadgetmg2::All.TypeOfTimestepCriterion >= 1){
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
    gadgetmg2::allocate_commbuffers();        /* ... allocate buffer-memory for particle
                                   exchange during force computation */
    gadgetmg2::set_units();
#if defined(PERIODIC) && (!defined(PMGRID) || defined(FORCETEST))
    ewald_init();
#endif
    gadgetmg2::random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(gadgetmg2::random_generator, 42);        /* start-up seed */
#ifdef PMGRID
    long_range_init();
#endif
    gadgetmg2::set_random_numbers();

    for(int i = 0; i < 6; i++)
        gadgetmg2::All.MassTable[i] = 0;

    if(gadgetmg2::All.ComovingIntegrationOn){
        gadgetmg2::All.TimeBegin = 1.0 / (1.0 + redshift_begin_parameter);
        gadgetmg2::All.TimeMax = 1.0 / (1.0 + redshift_max_parameter);
        gadgetmg2::All.Timebase_interval = (log(gadgetmg2::All.TimeMax) - log(gadgetmg2::All.TimeBegin)) / TIMEBASE;
    }else{
        gadgetmg2::All.Timebase_interval = (gadgetmg2::All.TimeMax - gadgetmg2::All.TimeBegin) / TIMEBASE;
    }
    gadgetmg2::All.Time = gadgetmg2::All.TimeBegin;
    gadgetmg2::All.Ti_Current = 0;
    gadgetmg2::All.NumCurrentTiStep = 0;        /* setup some counters */
    gadgetmg2::All.SnapshotFileCount = 0;
    gadgetmg2::All.TotNumOfForces = gadgetmg2::All.NumForcesSinceLastDomainDecomp = 0;
    gadgetmg2::All.TimeLastStatistics = gadgetmg2::All.TimeBegin - gadgetmg2::All.TimeBetStatistics;
#ifdef PMGRID
    gadgetmg2::All.PM_Ti_endstep = gadgetmg2::All.PM_Ti_begstep = 0;
#endif
#ifdef FLEXSTEPS
    gadgetmg2::All.PresentMinStep = TIMEBASE;
#endif
    gadgetmg2::All.Ti_nextoutput = -1; // find_next_outputtime(gadgetmg2::All.Ti_Current);
    gadgetmg2::open_outputfiles();
    outfiles_opened = true;
    return check_parameters();
}
int recommit_parameters(){
    return check_parameters();
}

int commit_particles(){
    double t0, t1;
    int i, j;
#ifndef ISOTHERM_EQS
    double a3;
#endif
    double a, a_inv;

    if (gadgetmg2::All.ComovingIntegrationOn) {
        a = gadgetmg2::All.Time;
        a_inv = 1.0 / gadgetmg2::All.Time;
    } else {
        a = a_inv = 1;
    }
    t0 = gadgetmg2::second();
    gadgetmg2::All.TotNumPart = dm_particles_in_buffer + sph_particles_in_buffer;
    gadgetmg2::All.TotN_gas = sph_particles_in_buffer;
    gadgetmg2::All.MaxPart = gadgetmg2::All.PartAllocFactor * (gadgetmg2::All.TotNumPart / gadgetmg2::NTask);        /* sets the maximum number of particles that may */
    gadgetmg2::All.MaxPartSph = gadgetmg2::All.PartAllocFactor * (gadgetmg2::All.TotN_gas / gadgetmg2::NTask);        /* sets the maximum number of particles that may
                                                                        reside on a processor */
    gadgetmg2::NumPart = dm_states.size()+sph_states.size();
    gadgetmg2::N_gas = sph_states.size();
    gadgetmg2::allocate_memory();
    // initialize sph particles
    i = 0;
    for (map<long long, sph_state>::iterator state_iter = sph_states.begin();
            state_iter != sph_states.end(); state_iter++, i++){
        gadgetmg2::P[i].ID = (*state_iter).first;
        gadgetmg2::P[i].Mass = (*state_iter).second.mass;
        gadgetmg2::P[i].Pos[0] = (*state_iter).second.x * a_inv;
        gadgetmg2::P[i].Pos[1] = (*state_iter).second.y * a_inv;
        gadgetmg2::P[i].Pos[2] = (*state_iter).second.z * a_inv;
        gadgetmg2::P[i].Vel[0] = (*state_iter).second.vx * a;
        gadgetmg2::P[i].Vel[1] = (*state_iter).second.vy * a;
        gadgetmg2::P[i].Vel[2] = (*state_iter).second.vz * a;
        gadgetmg2::P[i].Type = 0; // SPH particles (dark matter particles have type 1)
        gadgetmg2::SphP[i].Entropy = (*state_iter).second.u;
        gadgetmg2::SphP[i].Density = -1;
        gadgetmg2::SphP[i].Hsml = 0;
#ifdef MORRIS97VISC
        SphP[i].Alpha = (*state_iter).second.alpha;
        SphP[i].DAlphaDt = (*state_iter).second.dalphadt;
#endif
    }
    sph_states.clear();

    // initialize dark matter particles
    i = gadgetmg2::N_gas;
    for (map<long long, dynamics_state>::iterator state_iter = dm_states.begin();
            state_iter != dm_states.end(); state_iter++, i++){
        gadgetmg2::P[i].ID = (*state_iter).first;
        gadgetmg2::P[i].Mass = (*state_iter).second.mass;
        gadgetmg2::P[i].Pos[0] = (*state_iter).second.x * a_inv;
        gadgetmg2::P[i].Pos[1] = (*state_iter).second.y * a_inv;
        gadgetmg2::P[i].Pos[2] = (*state_iter).second.z * a_inv;
        gadgetmg2::P[i].Vel[0] = (*state_iter).second.vx * a;
        gadgetmg2::P[i].Vel[1] = (*state_iter).second.vy * a;
        gadgetmg2::P[i].Vel[2] = (*state_iter).second.vz * a;
        gadgetmg2::P[i].Type = 1; // dark matter particles (SPH particles have type 0)
    }
    dm_states.clear();
    gadgetmg2::All.TimeBegin += gadgetmg2::All.Ti_Current * gadgetmg2::All.Timebase_interval;
    gadgetmg2::All.Ti_Current = 0;
    gadgetmg2::All.Time = gadgetmg2::All.TimeBegin;
    gadgetmg2::set_softenings();
    for(i = 0; i < gadgetmg2::NumPart; i++){        /*  start-up initialization */
        for(j = 0; j < 3; j++)
            gadgetmg2::P[i].GravAccel[j] = 0;
#ifdef PMGRID
        for(j = 0; j < 3; j++)
            gadgetmg2::P[i].GravPM[j] = 0;
#endif
        gadgetmg2::P[i].Ti_endstep = 0;
        gadgetmg2::P[i].Ti_begstep = 0;
        gadgetmg2::P[i].OldAcc = 0;
        gadgetmg2::P[i].GravCost = 1;
        gadgetmg2::P[i].Potential = 0;
    }
#ifdef FLEXSTEPS
    for(i = 0; i < gadgetmg2::NumPart; i++)        /*  start-up initialization */
        P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
#endif
    for(i = 0; i < gadgetmg2::N_gas; i++){        /* initialize sph_properties */
        for(j = 0; j < 3; j++){
            gadgetmg2::SphP[i].VelPred[j] = gadgetmg2::P[i].Vel[j];
            gadgetmg2::SphP[i].HydroAccel[j] = 0;
        }
        gadgetmg2::SphP[i].DtEntropy = 0;
#ifdef TIMESTEP_UPDATE
        gadgetmg2::SphP[i].FeedbackFlag = 0;
        for(j = 0; j < 3; j++)
            gadgetmg2::SphP[i].FeedAccel[j] = 0;
#endif
    }

    if((gadgetmg2::All.MaxPart < 1000)){
        if (gadgetmg2::ThisTask == 0){
            cout << "Gadget takes "<< gadgetmg2::All.PartAllocFactor  << " times the number of particles on a processors as a maximum."<<endl;
            cout << "For large numbers of particles some room is always available for storing nodes from other processors." << endl;
            cout << "For smaller numbers, this assumption is incorrect."<<endl;
            cout << "Changed allocation of tree to include more nodes."<<endl;
        }
        //gadgetmg2::All.MaxPart = MAXTOPNODES;
        gadgetmg2::All.TreeAllocFactor = 4000.0 / gadgetmg2::All.MaxPart;
        gadgetmg2::ngb_treeallocate(MAX_NGB);
        gadgetmg2::force_treeallocate(10 * gadgetmg2::All.TreeAllocFactor * gadgetmg2::All.MaxPart, 10 * gadgetmg2::All.MaxPart);
    } else {
        gadgetmg2::ngb_treeallocate(MAX_NGB);
        gadgetmg2::force_treeallocate(gadgetmg2::All.TreeAllocFactor * gadgetmg2::All.MaxPart, gadgetmg2::All.MaxPart);
    }

    gadgetmg2::All.NumForcesSinceLastDomainDecomp = 1 + gadgetmg2::All.TotNumPart * gadgetmg2::All.TreeDomainUpdateFrequency;
    gadgetmg2::Flag_FullStep = 1;                /* to ensure that Peano-Hilber order is done */

    gadgetmg2::domain_Decomposition();        /* do initial domain decomposition (gives equal numbers of particles) */
    update_particle_map();
    index_of_highest_mapped_particle = local_index_map.rbegin()->first;
#ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &index_of_highest_mapped_particle, 1, MPI_LONG_LONG_INT, MPI_MAX, gadgetmg2::GADGET_WORLD);
#endif
    gadgetmg2::ngb_treebuild();                /* will build tree */
    gadgetmg2::setup_smoothinglengths();
    gadgetmg2::TreeReconstructFlag = 1;
  /* at this point, the entropy variable normally contains the
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly.
   * Once the density has been computed, we can convert thermal energy to entropy.
   */
#ifndef ISOTHERM_EQS
    if(gadgetmg2::All.ComovingIntegrationOn){a3 = gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time;}else{a3 = 1;}
    for(i = 0; i < gadgetmg2::N_gas; i++)
        gadgetmg2::SphP[i].Entropy = GAMMA_MINUS1 * gadgetmg2::SphP[i].Entropy / pow(gadgetmg2::SphP[i].Density / a3, GAMMA_MINUS1);
#endif
#ifdef PMGRID
    long_range_init_regionsize();
#endif
    if(gadgetmg2::All.ComovingIntegrationOn)
        gadgetmg2::init_drift_table();
    t1 = gadgetmg2::second();
    gadgetmg2::CPUThisRun += gadgetmg2::timediff(t0, t1);
    gadgetmg2::All.CPU_Total += gadgetmg2::timediff(t0, t1);
    particles_initialized = true;
    if (gadgetmg2::ThisTask == 0){
        cout << flush;
    }
    return 0;
}

void push_particle_data_on_state_vectors(){
    map<long long, int>::iterator iter;
    int i;
    double a_inv, a;
#ifndef ISOTHERM_EQS
    double a3;
    if(gadgetmg2::All.ComovingIntegrationOn){a3 = gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmg2::density();
        density_up_to_date = true;
    }
#endif
    if (gadgetmg2::All.ComovingIntegrationOn) {
        a = gadgetmg2::All.Time;
        a_inv = 1.0 / gadgetmg2::All.Time;
    } else {
        a = a_inv = 1;
    }
    for (iter = local_index_map.begin(); iter != local_index_map.end(); iter++){
        i = (*iter).second;
        if (gadgetmg2::P[i].Type == 0){
            // store sph particle data
            sph_state state;
            state.mass = gadgetmg2::P[i].Mass;
            state.x =    gadgetmg2::P[i].Pos[0] * a;
            state.y =    gadgetmg2::P[i].Pos[1] * a;
            state.z =    gadgetmg2::P[i].Pos[2] * a;
            state.vx =   gadgetmg2::P[i].Vel[0] * a_inv;
            state.vy =   gadgetmg2::P[i].Vel[1] * a_inv;
            state.vz =   gadgetmg2::P[i].Vel[2] * a_inv;
#ifdef ISOTHERM_EQS
            state.u = SphP[i].Entropy;
#else
            state.u = gadgetmg2::SphP[i].Entropy * pow(gadgetmg2::SphP[i].Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif

#ifdef MORRIS97VISC
            state.alpha = SphP[i].Alpha;
            state.dalphadt = SphP[i].DAlphaDt;
#endif


            sph_states.insert(std::pair<long long, sph_state>(gadgetmg2::P[i].ID, state));
        } else {
            // store dark matter particle data
            dynamics_state state;
            state.mass = gadgetmg2::P[i].Mass;
            state.x =    gadgetmg2::P[i].Pos[0] * a;
            state.y =    gadgetmg2::P[i].Pos[1] * a;
            state.z =    gadgetmg2::P[i].Pos[2] * a;
            state.vx =   gadgetmg2::P[i].Vel[0] * a_inv;
            state.vy =   gadgetmg2::P[i].Vel[1] * a_inv;
            state.vz =   gadgetmg2::P[i].Vel[2] * a_inv;
            dm_states.insert(std::pair<long long, dynamics_state>(gadgetmg2::P[i].ID, state));
        }
    }
}

int recommit_particles(){
    if (particles_initialized){
        push_particle_data_on_state_vectors();
        gadgetmg2::free_memory();
        gadgetmg2::ngb_treefree();
        gadgetmg2::force_treefree();
    }

    return commit_particles();
}

bool drift_to_t_end(int ti_end){
    bool done;
    int n, min, min_glob, flag, *temp;
    double timeold;
    double t0, t1;
    t0 = gadgetmg2::second();
    timeold = gadgetmg2::All.Time;
    for(n = 1, min = gadgetmg2::P[0].Ti_endstep; n < gadgetmg2::NumPart; n++){
        if(min > gadgetmg2::P[n].Ti_endstep){
            min = gadgetmg2::P[n].Ti_endstep;
        }
    }
#ifndef NOMPI
    MPI_Allreduce(&min, &min_glob, 1, MPI_INT, MPI_MIN, gadgetmg2::GADGET_WORLD);
#else
    min_glob = min;
#endif
    /* We check whether this is a full step where all particles are synchronized */
    flag = 1;
    for(n = 0; n < gadgetmg2::NumPart; n++)
        if(gadgetmg2::P[n].Ti_endstep > min_glob)
            flag = 0;
#ifndef NOMPI
    MPI_Allreduce(&flag, &gadgetmg2::Flag_FullStep, 1, MPI_INT, MPI_MIN, gadgetmg2::GADGET_WORLD);
#else
    Flag_FullStep = flag;
#endif
#ifdef PMGRID
    if(min_glob >= gadgetmg2::All.PM_Ti_endstep){
        min_glob = gadgetmg2::All.PM_Ti_endstep;
        Flag_FullStep = 1;
    }
#endif
    /* Determine 'NumForceUpdate', i.e. the number of particles on this processor that are going to be active */
    for(n = 0, gadgetmg2::NumForceUpdate = 0; n < gadgetmg2::NumPart; n++){
        if(gadgetmg2::P[n].Ti_endstep == min_glob)
#ifdef SELECTIVE_NO_GRAVITY
          if(!((1 << P[n].Type) & (SELECTIVE_NO_GRAVITY)))
#endif
            gadgetmg2::NumForceUpdate++;
    }
    /* note: NumForcesSinceLastDomainDecomp has type "long long" */
    temp = new int[gadgetmg2::NTask];

#ifndef NOMPI
    MPI_Allgather(&gadgetmg2::NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, gadgetmg2::GADGET_WORLD);
#else
    temp[0] = NumForceUpdate;
#endif
    for(n = 0; n < gadgetmg2::NTask; n++)
        gadgetmg2::All.NumForcesSinceLastDomainDecomp += temp[n];
    free(temp);
    t1 = gadgetmg2::second();
    gadgetmg2::All.CPU_Predict += gadgetmg2::timediff(t0, t1);
    if (min_glob >= ti_end){
        min_glob = ti_end;
        done = true;
    } else {
        done = false;
    }
    gadgetmg2::move_particles(gadgetmg2::All.Ti_Current, min_glob);
    gadgetmg2::All.Ti_Current = min_glob;
    if(gadgetmg2::All.ComovingIntegrationOn)
        gadgetmg2::All.Time = gadgetmg2::All.TimeBegin * exp(gadgetmg2::All.Ti_Current * gadgetmg2::All.Timebase_interval);
    else
        gadgetmg2::All.Time = gadgetmg2::All.TimeBegin + gadgetmg2::All.Ti_Current * gadgetmg2::All.Timebase_interval;
    gadgetmg2::All.TimeStep = gadgetmg2::All.Time - timeold;
    return done;
}

bool check_density_stopping_condition(){
    int stopping_condition_is_set;
    double minimum_density_parameter, maximum_density_parameter;
    get_stopping_condition_minimum_density_parameter(&minimum_density_parameter);
    get_stopping_condition_maximum_density_parameter(&maximum_density_parameter);
    for (int i=0; i<gadgetmg2::N_gas; i++) {
        if ( (gadgetmg2::SphP[i].Density < minimum_density_parameter) ||
             (gadgetmg2::SphP[i].Density > maximum_density_parameter)) {
            int stopping_index  = next_index_for_stopping_condition();
            if (stopping_index >= 0) {
                cout << "set_stopping_condition_info returned: " <<
                    set_stopping_condition_info(stopping_index, DENSITY_LIMIT_DETECTION) << endl;
                cout << "set_stopping_condition_particle_index returned: " <<
                    set_stopping_condition_particle_index(stopping_index, 0, gadgetmg2::P[i].ID) << endl;
            }
        }
    }

    mpi_collect_stopping_conditions();
    is_stopping_condition_set(DENSITY_LIMIT_DETECTION, &stopping_condition_is_set);
    return stopping_condition_is_set;
}
bool check_internal_energy_stopping_condition(){
    int stopping_condition_is_set;
    double internal_energy;
    double minimum_internal_energy_parameter, maximum_internal_energy_parameter;
    get_stopping_condition_minimum_internal_energy_parameter(&minimum_internal_energy_parameter);
    get_stopping_condition_maximum_internal_energy_parameter(&maximum_internal_energy_parameter);

#ifndef ISOTHERM_EQS
    double a3;
    if(gadgetmg2::All.ComovingIntegrationOn){a3 = gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time;}else{a3 = 1;}
#endif

    for (int i=0; i<gadgetmg2::N_gas; i++) {
#ifdef ISOTHERM_EQS
        internal_energy = SphP[i].Entropy;
#else
        internal_energy = gadgetmg2::SphP[i].Entropy *
                pow(gadgetmg2::SphP[i].Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif
        if ( (internal_energy < minimum_internal_energy_parameter) ||
             (internal_energy > maximum_internal_energy_parameter)) {
            int stopping_index  = next_index_for_stopping_condition();
            if (stopping_index > 0) {
                cout << "set_stopping_condition_info returned: " <<
                    set_stopping_condition_info(stopping_index, INTERNAL_ENERGY_LIMIT_DETECTION) << endl;
                cout << "set_stopping_condition_particle_index returned: " <<
                    set_stopping_condition_particle_index(stopping_index, 0, gadgetmg2::P[i].ID) << endl;
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

    if (t_end > gadgetmg2::All.TimeMax)
        return -7;
    gadgetmg2::ZeroTimestepEncountered = 0;
    if(gadgetmg2::All.ComovingIntegrationOn){
        Ti_end = log(t_end / gadgetmg2::All.TimeBegin) / gadgetmg2::All.Timebase_interval;
    } else {
        Ti_end = (t_end - gadgetmg2::All.TimeBegin) / gadgetmg2::All.Timebase_interval;
    }
    if (Ti_end >= gadgetmg2::All.Ti_Current){
        global_quantities_of_system_up_to_date = density_up_to_date = false;
        done = drift_to_t_end(Ti_end); /* find next synchronization point and drift particles to MIN(this time, t_end). */
        while (!done && gadgetmg2::All.Ti_Current < TIMEBASE && gadgetmg2::All.Time <= gadgetmg2::All.TimeMax) {
            t0 = gadgetmg2::second();
            gadgetmg2::every_timestep_stuff();        /* write some info to log-files */
            gadgetmg2::domain_Decomposition();        /* do domain decomposition if needed */
            particle_map_up_to_date = false;
            gadgetmg2::compute_accelerations(0);        /* compute accelerations for
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
            if((gadgetmg2::All.Time - gadgetmg2::All.TimeLastStatistics) >= gadgetmg2::All.TimeBetStatistics) {
#ifdef COMPUTE_POTENTIAL_ENERGY
                compute_potential();
#endif
                gadgetmg2::energy_statistics();        /* compute and output energy statistics */
                gadgetmg2::All.TimeLastStatistics += gadgetmg2::All.TimeBetStatistics;
            }
            gadgetmg2::advance_and_find_timesteps();        /* 'kick' active particles in
                            * momentum space and compute new
                            * timesteps for them  */
            done = drift_to_t_end(Ti_end);
            gadgetmg2::All.NumCurrentTiStep++;

            /* Check whether we need to interrupt the run */
#ifndef NOMPI
            MPI_Allreduce(MPI_IN_PLACE, &gadgetmg2::ZeroTimestepEncountered, 1, MPI_INT, MPI_MAX, gadgetmg2::GADGET_WORLD);
#endif
            if(gadgetmg2::ZeroTimestepEncountered)
                return -8;

            if(gadgetmg2::ThisTask == 0) {
                /* are we running out of CPU-time ? If yes, interrupt run. */
                if(gadgetmg2::CPUThisRun > 0.85 * gadgetmg2::All.TimeLimitCPU){
                    printf("reaching time-limit. stopping.\n");
                    stopflag = 2;
                }
            }
#ifndef NOMPI
            MPI_Bcast(&stopflag, 1, MPI_INT, 0, gadgetmg2::GADGET_WORLD);
#endif
            if(stopflag)
                return -5;

            t1 = gadgetmg2::second();
            gadgetmg2::All.CPU_Total += gadgetmg2::timediff(t0, t1);
            gadgetmg2::CPUThisRun += gadgetmg2::timediff(t0, t1);
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
    if (gadgetmg2::ThisTask == 0)
        cout << flush;
    if (gadgetmg2::All.Ti_Current > TIMEBASE || gadgetmg2::All.Time > gadgetmg2::All.TimeMax)
        return -7;
    return 0;
}
int evolve_to_redshift(double redshift){
    if (!gadgetmg2::All.ComovingIntegrationOn) {return -9;}
    return evolve_model_generic(1.0 / (1.0 + redshift));
}
int evolve_model(double t_end){
    if (gadgetmg2::All.ComovingIntegrationOn) {return -9;}
    return evolve_model_generic(t_end);
}

int synchronize_model() {
    return 0;
}

int construct_tree_if_needed(void){
    double tstart, tend;
    if (!particles_initialized)
        return -1;
    tstart = gadgetmg2::second();
    if (gadgetmg2::TreeReconstructFlag){
        if(gadgetmg2::ThisTask == 0)
            printf("Tree construction.\n");
        gadgetmg2::force_treebuild(gadgetmg2::NumPart);
        gadgetmg2::TreeReconstructFlag = 0;
        if(gadgetmg2::ThisTask == 0)
            printf("Tree construction done.\n");
    }
    tend = gadgetmg2::second();
    gadgetmg2::All.CPU_TreeConstruction += gadgetmg2::timediff(tstart, tend);
    return 0;
}

int new_dm_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz){
    particle_id_counter++;
    if (gadgetmg2::ThisTask == 0)
        *id = particle_id_counter;
    // Divide the particles equally over all Tasks, Gadget will redistribute them later.
    if (gadgetmg2::ThisTask == (dm_particles_in_buffer % gadgetmg2::NTask)){
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
    if (gadgetmg2::ThisTask == 0)
        *id = particle_id_counter;

    // Divide the sph particles equally over all Tasks, Gadget will redistribute them later.
    if (gadgetmg2::ThisTask == (sph_particles_in_buffer % gadgetmg2::NTask)){
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
        state.alpha = gadgetmg2::All.ArtBulkViscConst;
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
        found = 1 + gadgetmg2::P[(*it).second].Type; // 1 for sph; 2 for dm
    }
#ifndef NOMPI
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, gadgetmg2::GADGET_WORLD);
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
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_LOR, gadgetmg2::GADGET_WORLD);
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
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_LOR, gadgetmg2::GADGET_WORLD);
#endif
    if (found){
        sph_particles_in_buffer--;
        return 0;
    }
    return -3;
}

// parameter getters/setters:

int get_time_step(double *timestep){
    if (gadgetmg2::ThisTask) {return 0;}
    *timestep = gadgetmg2::All.TimeStep;
    return 0;
}
int set_time_step(double timestep){
    if (gadgetmg2::ThisTask) {return 0;}
    return -2;
}
int get_epsilon(double *epsilon){
    if (gadgetmg2::ThisTask) {return 0;}
    gadgetmg2::set_softenings();
    *epsilon = gadgetmg2::All.SofteningTable[1];
    return 0;
}
int set_epsilon(double epsilon){
    gadgetmg2::All.SofteningHalo = epsilon;
    gadgetmg2::set_softenings();
    return 0;
}
int get_eps2(double *epsilon_squared){
    if (gadgetmg2::ThisTask) {return 0;}
    return -2;
}
int set_eps2(double epsilon_squared){
    if (gadgetmg2::ThisTask) {return 0;}
    return -2;
}
int get_epsgas(double *gas_epsilon){
    if (gadgetmg2::ThisTask) {return 0;}
    gadgetmg2::set_softenings();
    *gas_epsilon = gadgetmg2::All.SofteningTable[0];
    return 0;
}
int set_epsgas(double gas_epsilon){
    gadgetmg2::All.SofteningGas = gas_epsilon;
    gadgetmg2::set_softenings();
    return 0;
}
int get_unit_mass(double *code_mass_unit){
    if (gadgetmg2::ThisTask) {return 0;}
    *code_mass_unit = gadgetmg2::All.UnitMass_in_g;
    return 0;
}
int set_unit_mass(double code_mass_unit){
    gadgetmg2::All.UnitMass_in_g = code_mass_unit;
    return 0;
}
int get_unit_length(double *code_length_unit){
    if (gadgetmg2::ThisTask) {return 0;}
    *code_length_unit = gadgetmg2::All.UnitLength_in_cm;
    return 0;
}
int set_unit_length(double code_length_unit){
    gadgetmg2::All.UnitLength_in_cm = code_length_unit;
    return 0;
}
int get_unit_time(double *code_time_unit){
    if (gadgetmg2::ThisTask) {return 0;}
    *code_time_unit = gadgetmg2::All.UnitTime_in_s;
    return 0;
}
int set_unit_time(double code_time_unit){
    gadgetmg2::All.UnitVelocity_in_cm_per_s = gadgetmg2::All.UnitLength_in_cm / code_time_unit;
    gadgetmg2::All.UnitTime_in_s = code_time_unit;
    return 0;
}
int get_unit_velocity(double *code_velocity_unit){
    if (gadgetmg2::ThisTask) {return 0;}
    *code_velocity_unit = gadgetmg2::All.UnitVelocity_in_cm_per_s;
    return 0;
}
int set_unit_velocity(double code_velocity_unit){
    gadgetmg2::All.UnitVelocity_in_cm_per_s = code_velocity_unit;
    gadgetmg2::All.UnitTime_in_s = gadgetmg2::All.UnitLength_in_cm / gadgetmg2::All.UnitVelocity_in_cm_per_s;
    return 0;
}

int get_gadget_output_directory(char **output_directory){
    if (gadgetmg2::ThisTask) {return 0;}
    *output_directory = gadgetmg2::All.OutputDir;
    return 0;
}
int get_viscosity_switch(char **viscosity_switch){
    if (gadgetmg2::ThisTask) {return 0;}
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
    strcpy(gadgetmg2::All.OutputDir, output_directory);
#ifdef WIN32
    const char sep[] = "\\";
#else
    const char sep[] = "/";
#endif // WIN32
    if(length > 0) {
        if(gadgetmg2::All.OutputDir[length - 1] != sep[0]) {
            strcat(gadgetmg2::All.OutputDir, sep);
        }
    }
    return 0;
}
int get_nogravity(int *no_gravity_flag){
    if (gadgetmg2::ThisTask) {return 0;}
#ifdef NOGRAVITY
    *no_gravity_flag = 1;
#else
    *no_gravity_flag = 0;
#endif
    return 0;
}
int get_gdgop(int *gadget_cell_opening_flag){
    if (gadgetmg2::ThisTask) {return 0;}
    *gadget_cell_opening_flag = gadgetmg2::All.TypeOfOpeningCriterion;
    return 0;
}
int set_gdgop(int gadget_cell_opening_flag){
    gadgetmg2::All.TypeOfOpeningCriterion = gadget_cell_opening_flag;
    return 0;
}
int get_isotherm(int *isothermal_flag){
    if (gadgetmg2::ThisTask) {return 0;}
#ifdef ISOTHERM_EQS
    *isothermal_flag = 1;
#else
    *isothermal_flag = 0;
#endif
    return 0;
}
int get_eps_is_h(int *eps_is_h_flag){
    if (gadgetmg2::ThisTask) {return 0;}
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) &&  defined(UNEQUALSOFTENINGS)
    *eps_is_h_flag = 1;
#else
    *eps_is_h_flag = 0;
#endif
    return 0;
}
int get_nsmooth(int *nsmooth){
    if (gadgetmg2::ThisTask) {return 0;}
    *nsmooth = gadgetmg2::All.DesNumNgb;
    return 0;
}
int set_nsmooth(int nsmooth){
    gadgetmg2::All.DesNumNgb = nsmooth;
    return 0;
}
int get_bh_tol(double *opening_angle){
    if (gadgetmg2::ThisTask) {return 0;}
    *opening_angle = gadgetmg2::All.ErrTolTheta;
    return 0;
}
int set_bh_tol(double opening_angle){
    gadgetmg2::All.ErrTolTheta = opening_angle;
    return 0;
}
int get_gdgtol(double *gadget_cell_opening_constant){
    if (gadgetmg2::ThisTask) {return 0;}
    *gadget_cell_opening_constant = gadgetmg2::All.ErrTolForceAcc;
    return 0;
}
int set_gdgtol(double gadget_cell_opening_constant){
    gadgetmg2::All.ErrTolForceAcc = gadget_cell_opening_constant;
    return 0;
}
int get_gamma(double *gamma){
    if (gadgetmg2::ThisTask) {return 0;}
    *gamma = GAMMA;
    return 0;
}
int get_alpha(double *artificial_viscosity_alpha){
    if (gadgetmg2::ThisTask) {return 0;}
    *artificial_viscosity_alpha = gadgetmg2::All.ArtBulkViscConst;
    return 0;
}
int set_alpha(double artificial_viscosity_alpha){
    gadgetmg2::All.ArtBulkViscConst = artificial_viscosity_alpha;
    return 0;
}
int get_beta(double *artificial_viscosity_beta){
    if (gadgetmg2::ThisTask) {return 0;}
    *artificial_viscosity_beta = gadgetmg2::All.ArtBulkViscBeta;
    return 0;
}
int set_beta(double artificial_viscosity_beta){
    gadgetmg2::All.ArtBulkViscBeta = artificial_viscosity_beta;
    return 0;
}
int get_courant(double *courant){
    if (gadgetmg2::ThisTask) {return 0;}
    *courant = gadgetmg2::All.CourantFac*2.0;
    return 0;
}
int set_courant(double courant){
    gadgetmg2::All.CourantFac = courant/2.0;
    return 0;
}
int get_nsmtol(double *n_neighbour_tol){
    if (gadgetmg2::ThisTask) {return 0;}
    *n_neighbour_tol = gadgetmg2::All.MaxNumNgbDeviation / gadgetmg2::All.DesNumNgb;
    return 0;
}
int set_nsmtol(double n_neighbour_tol){
    gadgetmg2::All.MaxNumNgbDeviation = n_neighbour_tol * gadgetmg2::All.DesNumNgb;
    return 0;
}

int get_energy_file(char **energy_file){
    if (gadgetmg2::ThisTask) {return 0;}
    *energy_file = gadgetmg2::All.EnergyFile;
    return 0;
}
int set_energy_file(char *energy_file){
    strcpy(gadgetmg2::All.EnergyFile, energy_file);
    return 0;
}
int get_info_file(char **info_file){
    if (gadgetmg2::ThisTask) {return 0;}
    *info_file = gadgetmg2::All.InfoFile;
    return 0;
}
int set_info_file(char *info_file){
    strcpy(gadgetmg2::All.InfoFile, info_file);
    return 0;
}
int get_timings_file(char **timings_file){
    if (gadgetmg2::ThisTask) {return 0;}
    *timings_file = gadgetmg2::All.TimingsFile;
    return 0;
}
int set_timings_file(char *timings_file){
    strcpy(gadgetmg2::All.TimingsFile, timings_file);
    return 0;
}
int get_cpu_file(char **cpu_file){
    if (gadgetmg2::ThisTask) {return 0;}
    *cpu_file = gadgetmg2::All.CpuFile;
    return 0;
}
int set_cpu_file(char *cpu_file){
    strcpy(gadgetmg2::All.CpuFile, cpu_file);
    return 0;
}

int get_time_limit_cpu(double *time_limit_cpu){
    if (gadgetmg2::ThisTask) {return 0;}
    *time_limit_cpu = gadgetmg2::All.TimeLimitCPU;
    return 0;
}
int set_time_limit_cpu(double time_limit_cpu){
    gadgetmg2::All.TimeLimitCPU = time_limit_cpu;
    return 0;
}
int get_comoving_integration_flag(bool *comoving_integration_flag){
    if (gadgetmg2::ThisTask) {return 0;}
    *comoving_integration_flag = gadgetmg2::All.ComovingIntegrationOn;
    return 0;
}
int set_comoving_integration_flag(bool comoving_integration_flag){
    gadgetmg2::All.ComovingIntegrationOn = comoving_integration_flag;
    return 0;
}
int get_type_of_timestep_criterion(int *type_of_timestep_criterion){
    if (gadgetmg2::ThisTask) {return 0;}
    *type_of_timestep_criterion = gadgetmg2::All.TypeOfTimestepCriterion;
    return 0;
}
int set_type_of_timestep_criterion(int type_of_timestep_criterion){
    gadgetmg2::All.TypeOfTimestepCriterion = type_of_timestep_criterion;
    return 0;
}
int get_begin_time(double *time_begin){
    if (gadgetmg2::ThisTask) {return 0;}
    *time_begin = gadgetmg2::All.TimeBegin;
    return 0;
}
int set_begin_time(double time_begin){
    gadgetmg2::All.TimeBegin = time_begin;
    return 0;
}
int get_time_max(double *time_max){
    if (gadgetmg2::ThisTask) {return 0;}
    *time_max = gadgetmg2::All.TimeMax;
    return 0;
}
int set_time_max(double time_max){
    gadgetmg2::All.TimeMax = time_max;
    return 0;
}
int get_redshift_begin(double *redshift_begin){
    if (gadgetmg2::ThisTask) {return 0;}
    *redshift_begin = redshift_begin_parameter;
    return 0;
}
int set_redshift_begin(double redshift_begin){
    redshift_begin_parameter = redshift_begin;
    return 0;
}
int get_redshift_max(double *redshift_max){
    if (gadgetmg2::ThisTask) {return 0;}
    *redshift_max = redshift_max_parameter;
    return 0;
}
int set_redshift_max(double redshift_max){
    redshift_max_parameter = redshift_max;
    return 0;
}
int get_omega_zero(double *omega_zero){
    if (gadgetmg2::ThisTask) {return 0;}
    *omega_zero = gadgetmg2::All.Omega0;
    return 0;
}
int set_omega_zero(double omega_zero){
    gadgetmg2::All.Omega0 = omega_zero;
    return 0;
}
int get_omega_lambda(double *omega_lambda){
    if (gadgetmg2::ThisTask) {return 0;}
    *omega_lambda = gadgetmg2::All.OmegaLambda;
    return 0;
}
int set_omega_lambda(double omega_lambda){
    gadgetmg2::All.OmegaLambda = omega_lambda;
    return 0;
}
int get_omega_baryon(double *omega_baryon){
    if (gadgetmg2::ThisTask) {return 0;}
    *omega_baryon = gadgetmg2::All.OmegaBaryon;
    return 0;
}
int set_omega_baryon(double omega_baryon){
    gadgetmg2::All.OmegaBaryon = omega_baryon;
    return 0;
}
int get_hubble_param(double *hubble_param){
    if (gadgetmg2::ThisTask) {return 0;}
    *hubble_param = gadgetmg2::All.HubbleParam;
    return 0;
}
int set_hubble_param(double hubble_param){
    gadgetmg2::All.HubbleParam = hubble_param;
    return 0;
}
int get_err_tol_int_accuracy(double *err_tol_int_accuracy){
    if (gadgetmg2::ThisTask) {return 0;}
    *err_tol_int_accuracy = gadgetmg2::All.ErrTolIntAccuracy;
    return 0;
}
int set_err_tol_int_accuracy(double err_tol_int_accuracy){
    gadgetmg2::All.ErrTolIntAccuracy = err_tol_int_accuracy;
    return 0;
}
int get_max_size_timestep(double *max_size_timestep){
    if (gadgetmg2::ThisTask) {return 0;}
    *max_size_timestep = gadgetmg2::All.MaxSizeTimestep;
    return 0;
}
int set_max_size_timestep(double max_size_timestep){
    gadgetmg2::All.MaxSizeTimestep = max_size_timestep;
    return 0;
}
int get_min_size_timestep(double *min_size_timestep){
    if (gadgetmg2::ThisTask) {return 0;}
    *min_size_timestep = gadgetmg2::All.MinSizeTimestep;
    return 0;
}
int set_min_size_timestep(double min_size_timestep){
    gadgetmg2::All.MinSizeTimestep = min_size_timestep;
    return 0;
}
int get_tree_domain_update_frequency(double *tree_domain_update_frequency){
    if (gadgetmg2::ThisTask) {return 0;}
    *tree_domain_update_frequency = gadgetmg2::All.TreeDomainUpdateFrequency;
    return 0;
}
int set_tree_domain_update_frequency(double tree_domain_update_frequency){
    gadgetmg2::All.TreeDomainUpdateFrequency = tree_domain_update_frequency;
    return 0;
}
int get_time_between_statistics(double *time_between_statistics){
    if (gadgetmg2::ThisTask) {return 0;}
    *time_between_statistics = gadgetmg2::All.TimeBetStatistics;
    return 0;
}
int set_time_between_statistics(double time_between_statistics){
    gadgetmg2::All.TimeBetStatistics = time_between_statistics;
    return 0;
}
int get_min_gas_temp(double *min_gas_temp){
    if (gadgetmg2::ThisTask) {return 0;}
    *min_gas_temp = gadgetmg2::All.MinGasTemp;
    return 0;
}
int set_min_gas_temp(double min_gas_temp){
    gadgetmg2::All.MinGasTemp = min_gas_temp;
    return 0;
}
int get_min_gas_hsmooth_fractional(double *min_gas_hsmooth_fractional){
    if (gadgetmg2::ThisTask) {return 0;}
    *min_gas_hsmooth_fractional = gadgetmg2::All.MinGasHsmlFractional;
    return 0;
}
int set_min_gas_hsmooth_fractional(double min_gas_hsmooth_fractional){
    gadgetmg2::All.MinGasHsmlFractional = min_gas_hsmooth_fractional;
    return 0;
}
int get_softening_gas_max_phys(double *softening_gas_max_phys){
    if (gadgetmg2::ThisTask) {return 0;}
    *softening_gas_max_phys = gadgetmg2::All.SofteningGasMaxPhys;
    return 0;
}
int set_softening_gas_max_phys(double softening_gas_max_phys){
    gadgetmg2::All.SofteningGasMaxPhys = softening_gas_max_phys;
    return 0;
}
int get_softening_halo_max_phys(double *softening_halo_max_phys){
    if (gadgetmg2::ThisTask) {return 0;}
    *softening_halo_max_phys = gadgetmg2::All.SofteningHaloMaxPhys;
    return 0;
}
int set_softening_halo_max_phys(double softening_halo_max_phys){
    gadgetmg2::All.SofteningHaloMaxPhys = softening_halo_max_phys;
    return 0;
}

int get_box_size(double *value)
{
    if (gadgetmg2::ThisTask) {return 0;}
    *value = gadgetmg2::All.BoxSize;
    return 0;
}

int set_box_size(double value)
{
    gadgetmg2::All.BoxSize = value;
    return 0;
}

int get_periodic_boundaries_flag(bool *value)
{
    if (gadgetmg2::ThisTask) {return 0;}
    *value = gadgetmg2::All.PeriodicBoundariesOn;
    return 0;
}

int set_periodic_boundaries_flag(bool value)
{
// gadgetmg2::All.PeriodicBoundariesOn is read only because compile time determined
    return -2;
}

int get_interpret_kicks_as_feedback_flag(bool *value)
{
    if (gadgetmg2::ThisTask) {return 0;}
    *value = interpret_kicks_as_feedback;
    return 0;
}

int set_interpret_kicks_as_feedback_flag(bool value)
{
    interpret_kicks_as_feedback = value;
    return 0;
}

int get_interpret_heat_as_feedback_flag(bool *value) {
    if (gadgetmg2::ThisTask) {return 0;}
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

    if (gadgetmg2::ThisTask == 0){
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, &next_local_index, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, gadgetmg2::GADGET_WORLD);
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
        MPI_Reduce(&next_local_index, NULL, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, gadgetmg2::GADGET_WORLD);
#endif
        return 0;
    }
}

void update_particle_map(void){
    local_index_map.clear();
    for(int i = 0; i < gadgetmg2::NumPart; i++) {
        local_index_map.insert(std::pair<long long, int>(gadgetmg2::P[i].ID, i));
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
            buffer[i] = gadgetmg2::P[local_index].Mass;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
        return 0;
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
            gadgetmg2::P[local_index].Mass = mass[i];
            count[i] = 1;
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}

int get_radius(int index, double *radius){
    return -2;
}

int set_radius(int index, double radius){
    return -2;
}

int get_position_comoving(int *index, double *x, double *y, double *z, int length){
    int errors = 0;
    double *buffer = new double[length*3];
    int *count = new int[length];
    int local_index;
#ifdef PERIODIC
    double boxSize = gadgetmg2::All.BoxSize;
    double boxHalf = 0.5 * gadgetmg2::All.BoxSize;
#endif

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = gadgetmg2::P[local_index].Pos[0];
            buffer[i+length] = gadgetmg2::P[local_index].Pos[1];
            buffer[i+2*length] = gadgetmg2::P[local_index].Pos[2];
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        for (int i = 0; i < length; i++){
            x[i] *= gadgetmg2::All.Time;
            y[i] *= gadgetmg2::All.Time;
            z[i] *= gadgetmg2::All.Time;
        }
    }
    return result;
}

int set_position_comoving(int *index, double *x, double *y, double *z, int length){
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            gadgetmg2::P[local_index].Pos[0] = x[i];
            gadgetmg2::P[local_index].Pos[1] = y[i];
            gadgetmg2::P[local_index].Pos[2] = z[i];
            count[i] = 1;
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}
int set_position(int *index, double *x, double *y, double *z, int length){
    if(gadgetmg2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmg2::All.Time;
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
            buffer[i] = gadgetmg2::P[local_index].Vel[0];
            buffer[i+length] = gadgetmg2::P[local_index].Vel[1];
            buffer[i+2*length] = gadgetmg2::P[local_index].Vel[2];
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        double a2_inv = 1.0 / (gadgetmg2::All.Time * gadgetmg2::All.Time);
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
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmg2::All.Time;
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
            gadgetmg2::P[local_index].Vel[0] = vx[i];
            gadgetmg2::P[local_index].Vel[1] = vy[i];
            gadgetmg2::P[local_index].Vel[2] = vz[i];
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_kicks_as_feedback && gadgetmg2::P[local_index].Type == 0) {
                gadgetmg2::SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if(interpret_kicks_as_feedback && gadgetmg2::P[local_index].Type == 0 && gadgetmg2::P[local_index].Ti_endstep != gadgetmg2::All.Ti_Current) {
                gadgetmg2::make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}
int set_velocity_comoving(int *index, double *vx, double *vy, double *vz, int length){
    if(gadgetmg2::All.ComovingIntegrationOn) {
        double a2 = gadgetmg2::All.Time * gadgetmg2::All.Time;
        for (int i = 0; i < length; i++){
            vx[i] *= a2;
            vy[i] *= a2;
            vz[i] *= a2;
        }
    }
    return set_velocity_gadget_u(index, vx, vy, vz, length);
}
int set_velocity(int *index, double *vx, double *vy, double *vz, int length){
    if(gadgetmg2::All.ComovingIntegrationOn) {
        for (int i = 0; i < length; i++){
            vx[i] *= gadgetmg2::All.Time;
            vy[i] *= gadgetmg2::All.Time;
            vz[i] *= gadgetmg2::All.Time;
        }
    }
    return set_velocity_gadget_u(index, vx, vy, vz, length);
}

int get_state_gadget(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length) {
    int errors = 0;
    double *buffer = new double[length*7];
    int *count = new int[length];
    int local_index;
#ifdef PERIODIC
    double boxSize = gadgetmg2::All.BoxSize;
    double boxHalf = 0.5 * gadgetmg2::All.BoxSize;
#endif

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = gadgetmg2::P[local_index].Mass;
            buffer[i+length] = gadgetmg2::P[local_index].Pos[0];
            buffer[i+2*length] = gadgetmg2::P[local_index].Pos[1];
            buffer[i+3*length] = gadgetmg2::P[local_index].Pos[2];
            buffer[i+4*length] = gadgetmg2::P[local_index].Vel[0];
            buffer[i+5*length] = gadgetmg2::P[local_index].Vel[1];
            buffer[i+6*length] = gadgetmg2::P[local_index].Vel[2];
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
            buffer[i+3*length] = 0;
            buffer[i+4*length] = 0;
            buffer[i+5*length] = 0;
            buffer[i+6*length] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*7, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*7, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
            } else {
                mass[i] = buffer[i];
                x[i] = buffer[i+length];
                y[i] = buffer[i+2*length];
                z[i] = buffer[i+3*length];
                vx[i] = buffer[i+4*length];
                vy[i] = buffer[i+5*length];
                vz[i] = buffer[i+6*length];
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
    int result = get_state_gadget(index, mass, x, y, z, vx, vy, vz, length);
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        double a2_inv = 1.0 / (gadgetmg2::All.Time * gadgetmg2::All.Time);
        for (int i = 0; i < length; i++){
            vx[i] *= a2_inv;
            vy[i] *= a2_inv;
            vz[i] *= a2_inv;
        }
    }
    return result;
}
int get_state(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length) {
    int result = get_state_gadget(index, mass, x, y, z, vx, vy, vz, length);
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmg2::All.Time;
        for (int i = 0; i < length; i++){
            x[i] *= gadgetmg2::All.Time;
            y[i] *= gadgetmg2::All.Time;
            z[i] *= gadgetmg2::All.Time;
            vx[i] *= a_inv;
            vy[i] *= a_inv;
            vz[i] *= a_inv;
        }
    }
    return result;
}

int set_state_gadget(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length){
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            gadgetmg2::P[local_index].Mass = mass[i];
            gadgetmg2::P[local_index].Pos[0] = x[i];
            gadgetmg2::P[local_index].Pos[1] = y[i];
            gadgetmg2::P[local_index].Pos[2] = z[i];
            gadgetmg2::P[local_index].Vel[0] = vx[i];
            gadgetmg2::P[local_index].Vel[1] = vy[i];
            gadgetmg2::P[local_index].Vel[2] = vz[i];
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_kicks_as_feedback && gadgetmg2::P[local_index].Type == 0) {
                gadgetmg2::SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if(interpret_kicks_as_feedback && gadgetmg2::P[local_index].Type == 0 && gadgetmg2::P[local_index].Ti_endstep != gadgetmg2::All.Ti_Current) {
                gadgetmg2::make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}
int set_state_comoving(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length){
    if(gadgetmg2::All.ComovingIntegrationOn) {
        double a2 = gadgetmg2::All.Time * gadgetmg2::All.Time;
        for (int i = 0; i < length; i++){
            vx[i] *= a2;
            vy[i] *= a2;
            vz[i] *= a2;
        }
    }
    return set_state_gadget(index, mass, x, y, z, vx, vy, vz, length);
}
int set_state(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length){
    if(gadgetmg2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmg2::All.Time;
        for (int i = 0; i < length; i++){
            x[i] *= a_inv;
            y[i] *= a_inv;
            z[i] *= a_inv;
            vx[i] *= gadgetmg2::All.Time;
            vy[i] *= gadgetmg2::All.Time;
            vz[i] *= gadgetmg2::All.Time;
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
    double boxSize = gadgetmg2::All.BoxSize;
    double boxHalf = 0.5 * gadgetmg2::All.BoxSize;
#endif
#ifndef ISOTHERM_EQS
    double a3;

    if (!density_up_to_date){
        gadgetmg2::density();
        density_up_to_date = true;
    }
    if(gadgetmg2::All.ComovingIntegrationOn){a3 = gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time;}else{a3 = 1;}
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmg2::P[local_index].Mass;
            buffer[i+length] = gadgetmg2::P[local_index].Pos[0];
            buffer[i+2*length] = gadgetmg2::P[local_index].Pos[1];
            buffer[i+3*length] = gadgetmg2::P[local_index].Pos[2];
            buffer[i+4*length] = gadgetmg2::P[local_index].Vel[0];
            buffer[i+5*length] = gadgetmg2::P[local_index].Vel[1];
            buffer[i+6*length] = gadgetmg2::P[local_index].Vel[2];
#ifdef ISOTHERM_EQS
            buffer[i+7*length] = gadgetmg2::SphP[local_index].Entropy;
#else
            buffer[i+7*length] = gadgetmg2::SphP[local_index].Entropy *
                pow(gadgetmg2::SphP[local_index].Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
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
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*8, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*8, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmg2::All.Time;
        for (int i = 0; i < length; i++){
            x[i] *= gadgetmg2::All.Time;
            y[i] *= gadgetmg2::All.Time;
            z[i] *= gadgetmg2::All.Time;
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
        gadgetmg2::density();
        density_up_to_date = true;
    }
    if(gadgetmg2::All.ComovingIntegrationOn){a3 = gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time;}else{a3 = 1;}
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            gadgetmg2::P[local_index].Mass = mass[i];
            gadgetmg2::P[local_index].Pos[0] = x[i];
            gadgetmg2::P[local_index].Pos[1] = y[i];
            gadgetmg2::P[local_index].Pos[2] = z[i];
            gadgetmg2::P[local_index].Vel[0] = vx[i];
            gadgetmg2::P[local_index].Vel[1] = vy[i];
            gadgetmg2::P[local_index].Vel[2] = vz[i];
#ifdef ISOTHERM_EQS
            gadgetmg2::SphP[local_index].Entropy = internal_energy[i];
#else
            gadgetmg2::SphP[local_index].Entropy = GAMMA_MINUS1 * internal_energy[i] /
                pow(gadgetmg2::SphP[local_index].Density / a3, GAMMA_MINUS1);
#endif
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_heat_as_feedback || interpret_kicks_as_feedback) {
                gadgetmg2::SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if ((interpret_heat_as_feedback || interpret_kicks_as_feedback) &&
                    gadgetmg2::P[local_index].Ti_endstep != gadgetmg2::All.Ti_Current) {
                gadgetmg2::make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}
int set_state_sph(int *index, double *mass, double *x, double *y, double *z,
        double *vx, double *vy, double *vz, double *internal_energy, int length){
    if(gadgetmg2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmg2::All.Time;
        for (int i = 0; i < length; i++){
            x[i] *= a_inv;
            y[i] *= a_inv;
            z[i] *= a_inv;
            vx[i] *= gadgetmg2::All.Time;
            vy[i] *= gadgetmg2::All.Time;
            vz[i] *= gadgetmg2::All.Time;
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
            buffer[i] = gadgetmg2::P[local_index].GravAccel[0];
            buffer[i+length] = gadgetmg2::P[local_index].GravAccel[1];
            buffer[i+2*length] = gadgetmg2::P[local_index].GravAccel[2];
            if(gadgetmg2::P[local_index].Type == 0){
                buffer[i] += gadgetmg2::SphP[local_index].HydroAccel[0];
                buffer[i+length] += gadgetmg2::SphP[local_index].HydroAccel[1];
                buffer[i+2*length] += gadgetmg2::SphP[local_index].HydroAccel[2];
            }
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length*3, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        for (int i = 0; i < length; i++){
            ax[i] *= gadgetmg2::All.Time;
            ay[i] *= gadgetmg2::All.Time;
            az[i] *= gadgetmg2::All.Time;
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

    if(gadgetmg2::All.ComovingIntegrationOn){a3 = gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmg2::density();
        density_up_to_date = true;
    }
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            count[i] = 1;
#ifdef ISOTHERM_EQS
            buffer[i] = gadgetmg2::SphP[local_index].Entropy;
#else
            buffer[i] = gadgetmg2::SphP[local_index].Entropy *
                pow(gadgetmg2::SphP[local_index].Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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

    if(gadgetmg2::All.ComovingIntegrationOn){a3 = gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmg2::density();
        density_up_to_date = true;
    }
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
#ifdef ISOTHERM_EQS
            gadgetmg2::SphP[local_index].Entropy = internal_energy[i];
#else
            gadgetmg2::SphP[local_index].Entropy = GAMMA_MINUS1 * internal_energy[i] /
                pow(gadgetmg2::SphP[local_index].Density / a3, GAMMA_MINUS1);
#endif
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_heat_as_feedback) {
                gadgetmg2::SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if(interpret_heat_as_feedback && gadgetmg2::P[local_index].Ti_endstep != gadgetmg2::All.Ti_Current) {
                gadgetmg2::make_it_active(local_index);
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
        gadgetmg2::density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmg2::SphP[local_index].Hsml;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        for (int i = 0; i < length; i++){
            smoothing_length[i] *= gadgetmg2::All.Time;
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
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            count[i] = 1;
#ifdef MORRIS97VISC
            buffer[i] = gadgetmg2::SphP[local_index].Alpha;
#else
	    buffer[i] = gadgetmg2::All.ArtBulkViscConst;
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            count[i] = 1;
#ifdef MORRIS97VISC
            buffer[i] = gadgetmg2::SphP[local_index].DAlphaDt;
#else
            buffer[i] = 0;
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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

    if(gadgetmg2::All.ComovingIntegrationOn){a3 = gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmg2::density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmg2::SphP[local_index].Density / a3;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        double a3_inv = 1.0 / (gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time);
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

    if(gadgetmg2::All.ComovingIntegrationOn){a = gadgetmg2::All.Time;}else{a = 1;}
    if (!density_up_to_date){
        gadgetmg2::density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmg2::SphP[local_index].Pressure / a;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
    if(gadgetmg2::ThisTask == 0 && gadgetmg2::All.ComovingIntegrationOn) {
        double a_inv = 1.0 / gadgetmg2::All.Time;
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
    if (gadgetmg2::All.ComovingIntegrationOn){
        hubble = gadgetmg2::All.Hubble * sqrt(gadgetmg2::All.Omega0 / (gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time)
            + (1 - gadgetmg2::All.Omega0 - gadgetmg2::All.OmegaLambda) / (gadgetmg2::All.Time * gadgetmg2::All.Time) + gadgetmg2::All.OmegaLambda);
    } else {
        hubble = 1;
    }

#ifndef ISOTHERM_EQS
    //~double a3;
    //~if(gadgetmg2::All.ComovingIntegrationOn){a3 = gadgetmg2::All.Time * gadgetmg2::All.Time * gadgetmg2::All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        gadgetmg2::density();
        density_up_to_date = true;
    }
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            count[i] = 1;
#ifdef ISOTHERM_EQS
            buffer[i] = gadgetmg2::SphP[local_index].DtEntropy * hubble;
#else
            buffer[i] = - gadgetmg2::SphP[local_index].Pressure * gadgetmg2::SphP[local_index].DivVel /
                gadgetmg2::SphP[local_index].Density * hubble;
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
        gadgetmg2::density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && gadgetmg2::P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = gadgetmg2::SphP[local_index].NumNgb;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
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
    gadgetmg2::set_softenings();
    if (gadgetmg2::ThisTask) {return 0;}
    double a;
    if (gadgetmg2::All.ComovingIntegrationOn) {a = gadgetmg2::All.Time;} else {a = 1;}
    for (int i = 0; i < length; i++){
        epsilon[i] = a * gadgetmg2::All.SofteningTable[1];
    }
    return 0;
}
int get_epsilon_gas_part(int *index, double *epsilon, int length){
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) &&  defined(UNEQUALSOFTENINGS)
    return get_smoothing_length(index, epsilon, length);
#else
    gadgetmg2::set_softenings();
    if (gadgetmg2::ThisTask) {return 0;}
    double a;
    if (gadgetmg2::All.ComovingIntegrationOn) {a = gadgetmg2::All.Time;} else {a = 1;}
    for (int i = 0; i < length; i++){
        epsilon[i] = a * gadgetmg2::All.SofteningTable[0];
    }
    return 0;
#endif
}



// simulation property getters:

void update_global_quantities(bool do_potential){
    if (do_potential) {
        gadgetmg2::compute_potential();
        potential_energy_also_up_to_date = true;
    } else {potential_energy_also_up_to_date = false;}
    gadgetmg2::compute_global_quantities_of_system();
    global_quantities_of_system_up_to_date = true;
}
int get_time(double *time){
    if (gadgetmg2::ThisTask) {return 0;}
    if (gadgetmg2::All.ComovingIntegrationOn) {return -9;}
    *time = gadgetmg2::All.Time;
    return 0;
}


int get_redshift(double *redshift){
    if (gadgetmg2::ThisTask) {return 0;}
    if (!gadgetmg2::All.ComovingIntegrationOn) {return -9;}
    *redshift = 1.0 / gadgetmg2::All.Time - 1.0;
    return 0;
}
int get_total_radius(double *radius){
    double r_squared, local_max = 0;
    int i, j;

    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    for (i = 0; i < gadgetmg2::NumPart; i++){
        for (r_squared = 0, j = 0; j < 3; j++)
            r_squared += (gadgetmg2::SysState.CenterOfMass[j]-gadgetmg2::P[i].Pos[j])*(gadgetmg2::SysState.CenterOfMass[j]-gadgetmg2::P[i].Pos[j]);
        if (r_squared > local_max)
            local_max = r_squared;
    }

    if(gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(&local_max, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, &local_max, 1, MPI_DOUBLE, MPI_MAX, 0, gadgetmg2::GADGET_WORLD);
#endif
        if (gadgetmg2::All.ComovingIntegrationOn){
            *radius = gadgetmg2::All.Time * sqrt(local_max);
        } else {
            *radius = sqrt(local_max);
        }
    }
    return 0;
}

int get_total_mass(double *mass){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (gadgetmg2::ThisTask) {return 0;}
    *mass = gadgetmg2::SysState.Mass;
    return 0;
}

int get_potential(int *index, double *potential, int length) {
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    if (!potential_energy_also_up_to_date) {
        gadgetmg2::compute_potential();
        potential_energy_also_up_to_date = true;
    }

    for (int i = 0; i < length; i++) {
        if (found_particle(index[i], &local_index)) {
            count[i] = 1;
            buffer[i] = gadgetmg2::P[local_index].Potential;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }

    if (gadgetmg2::ThisTask) {
#ifndef NOMPI
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
    } else {
#ifndef NOMPI
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, gadgetmg2::GADGET_WORLD);
#endif
        double a2;
        if (gadgetmg2::All.ComovingIntegrationOn) {a2 = gadgetmg2::All.Time * gadgetmg2::All.Time;} else {a2 = 1;}
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
    if (gadgetmg2::ThisTask) {return 0;}
    *kinetic_energy = gadgetmg2::SysState.EnergyKin;
    return 0;
}
int get_potential_energy(double *potential_energy){
    if (!(global_quantities_of_system_up_to_date && potential_energy_also_up_to_date))
        update_global_quantities(true);
    if (gadgetmg2::ThisTask) {return 0;}
    *potential_energy = gadgetmg2::SysState.EnergyPot;
    return 0;
}
int get_thermal_energy(double *thermal_energy){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (gadgetmg2::ThisTask) {return 0;}
    *thermal_energy = gadgetmg2::SysState.EnergyInt;
    return 0;
}
int get_number_of_particles(int *number_of_particles){
    if (gadgetmg2::ThisTask) {return 0;}
    *number_of_particles = gadgetmg2::All.TotNumPart;
    return 0;
}
int get_center_of_mass_position(double *x, double *y, double *z){
#ifdef PERIODIC
    return -2;
#endif
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (gadgetmg2::ThisTask) {return 0;}
    double a;
    if (gadgetmg2::All.ComovingIntegrationOn) {a = gadgetmg2::All.Time;} else {a = 1;}
    *x = a * gadgetmg2::SysState.CenterOfMass[0];
    *y = a * gadgetmg2::SysState.CenterOfMass[1];
    *z = a * gadgetmg2::SysState.CenterOfMass[2];
    return 0;
}
int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (gadgetmg2::ThisTask) {return 0;}
    double a_inv;
    if (gadgetmg2::All.ComovingIntegrationOn) {a_inv = 1.0 / gadgetmg2::All.Time;} else {a_inv = 1;}
    *vx = a_inv * gadgetmg2::SysState.Momentum[0]/gadgetmg2::SysState.Mass;
    *vy = a_inv * gadgetmg2::SysState.Momentum[1]/gadgetmg2::SysState.Mass;
    *vz = a_inv * gadgetmg2::SysState.Momentum[2]/gadgetmg2::SysState.Mass;
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
    double pos[3], vel[3];
    double h_out, ngb_out, dhsml_out, rho_out, rhov_out[3], rhov2_out, rhoe_out;
    int error;
    double a, a_inv, a3_inv, a4_inv, a5_inv;
    if (gadgetmg2::All.ComovingIntegrationOn) {
        a = gadgetmg2::All.Time;
        a_inv = 1.0 / gadgetmg2::All.Time;
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
            pos[i] += gadgetmg2::All.BoxSize;
        }
    }
#endif
    vel[0] = a * vx;
    vel[1] = a * vy;
    vel[2] = a * vz;
    gadgetmg2::hydro_state_at_point(pos, vel, &h_out, &ngb_out, &dhsml_out, &rho_out, rhov_out, &rhov2_out, &rhoe_out);
    if (gadgetmg2::ThisTask) {return 0;}
    *rho   = rho_out * a3_inv;
    *rhovx = rhov_out[0] * a4_inv;
    *rhovy = rhov_out[1] * a4_inv;
    *rhovz = rhov_out[2] * a4_inv;
#ifdef ISOTHERM_EQS
    *rhoe = a3_inv * rhoe_out + a5_inv * 0.5*(rhov_out[0]*rhov_out[0] + rhov_out[1]*rhov_out[1] + rhov_out[2]*rhov_out[2]) / rho_out;
#else
    *rhoe = a3_inv * rhoe_out * (pow(rho_out * a3_inv, GAMMA_MINUS1) / GAMMA_MINUS1) + a5_inv * 0.5*rhov2_out;
#endif
    return 0;
}
