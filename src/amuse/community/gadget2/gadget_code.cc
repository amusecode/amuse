#include <mpi.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <map>
#include <math.h>
#include "gadget_code.h"
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

// general interface functions:

void set_default_parameters(){
    // parameters that can be changed from AMUSE
    All.TimeLimitCPU = 36000;
    All.ComovingIntegrationOn = 0;
    All.TypeOfTimestepCriterion = 0;
    All.PeriodicBoundariesOn = 0;
    All.Time = 0.0;
    All.TimeBegin = 0.0;
    All.TimeMax = 100.0;
    All.Omega0 = 0;
    All.OmegaLambda = 0;
    All.OmegaBaryon = 0;
    All.HubbleParam = 0.7;
    All.BoxSize = 1.0;
    All.TimeBetStatistics = 0.1;
    All.ErrTolIntAccuracy = 0.025;
    All.CourantFac = 0.15;
    All.MaxSizeTimestep = 0.01;
    All.MinSizeTimestep = 0.0;
    All.ErrTolTheta = 0.5;
    All.TypeOfOpeningCriterion = 1;
    All.ErrTolForceAcc = 0.005;
    All.TreeDomainUpdateFrequency = 0.05;
    All.DesNumNgb = 50;
    All.MaxNumNgbDeviation = 5.;
    All.ArtBulkViscConst = 0.5;
    All.MinGasTemp = 0;
    All.UnitLength_in_cm = 3.085678e21;
    All.UnitMass_in_g = 1.989e43;
    All.UnitVelocity_in_cm_per_s = 1e5;
    All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
    All.MinGasHsmlFractional = 0.0;
    All.SofteningGas = 0.01;
    All.SofteningHalo = 0.01;
    All.SofteningGasMaxPhys = 0.0;
    All.SofteningHaloMaxPhys = 0.0;
    set_softenings();
    strcpy(All.OutputDir,   ".");
    strcpy(All.EnergyFile,  "energy.txt");
    strcpy(All.InfoFile,    "info.txt");
    strcpy(All.TimingsFile, "timings.txt");
    strcpy(All.CpuFile,     "cpu.txt");
    
    // parameters that are fixed for AMUSE:
    All.PartAllocFactor = 1.5; // Memory allocation parameter
    All.TreeAllocFactor = 0.8; // Memory allocation parameter
    All.BufferSize = 25;       // Memory allocation parameter
    All.ResubmitOn = 0;              // Keep this turned off!
    All.OutputListOn = 0;            // Keep this turned off!
    All.GravityConstantInternal = 0; // Keep this turned off!
    
    // parameters that are unused for AMUSE:
    strcpy(All.InitCondFile, "");
    strcpy(All.RestartFile, "");
    strcpy(All.SnapshotFileBase, "");
    strcpy(All.OutputListFilename, "");
    strcpy(All.ResubmitCommand, "");
    All.ICFormat = 1;
    All.SnapFormat = 1;
    All.TimeBetSnapshot = 100.0;
    All.TimeOfFirstSnapshot = 100.0;
    All.CpuTimeBetRestartFile = 36000.0;
    All.NumFilesPerSnapshot = 1;
    All.NumFilesWrittenInParallel = 1;
    All.InitGasTemp = 0;
    All.MaxRMSDisplacementFac = 0.2; // parameter for PM; PM is currently not supported
}

int initialize_code(){
    double t0, t1;
    MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
    MPI_Comm_size(MPI_COMM_WORLD, &NTask);
    for(PTask = 0; NTask > (1 << PTask); PTask++);
    RestartFlag = All.TotNumPart = All.TotN_gas = 0;
    All.CPU_TreeConstruction = All.CPU_TreeWalk = All.CPU_Gravity = All.CPU_Potential = All.CPU_Domain =
        All.CPU_Snapshot = All.CPU_Total = All.CPU_CommSum = All.CPU_Imbalance = All.CPU_Hydro =
        All.CPU_HydCompWalk = All.CPU_HydCommSumm = All.CPU_HydImbalance =
        All.CPU_EnsureNgb = All.CPU_Predict = All.CPU_TimeLine = All.CPU_PM = All.CPU_Peano = 0;
    CPUThisRun = 0;
    t0 = second();
    if(ThisTask == 0){
        printf("\nThis is Gadget, version `%s'.\n", GADGETVERSION);
        printf("\nRunning on %d processors.\n", NTask);
    }

    t1 = second();
    CPUThisRun += timediff(t0, t1);
    All.CPU_Total += timediff(t0, t1);

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
        close_outputfiles();
    if (particles_initialized){
        free_memory();
        ngb_treefree();
        force_treefree();
    }
    return 0;
}

int check_parameters(){
    MPI_Bcast(&All, sizeof(struct global_data_all_processes), MPI_BYTE, 0, MPI_COMM_WORLD);
    if (ThisTask)
        return 0;
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
    if(All.NumFilesWrittenInParallel < 1){
        printf("NumFilesWrittenInParallel MUST be at least 1\n");
        return -4;
    }
    if(All.NumFilesWrittenInParallel > NTask){
        printf("NumFilesWrittenInParallel MUST be smaller than number of processors\n");
        return -4;
    }
#ifdef PERIODIC
    if(All.PeriodicBoundariesOn == 0){
        printf("Code was compiled with periodic boundary conditions switched on.\n");
        printf("You must set `PeriodicBoundariesOn=1', or recompile the code.\n");
        return -4;
    }
#else
    if(All.PeriodicBoundariesOn == 1){
        printf("Code was compiled with periodic boundary conditions switched off.\n");
        printf("You must set `PeriodicBoundariesOn=0', or recompile the code.\n");
        return -4;
    }
#endif
    if(All.TypeOfTimestepCriterion >= 1){
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
    allocate_commbuffers();        /* ... allocate buffer-memory for particle
                                   exchange during force computation */
    set_units();
#if defined(PERIODIC) && (!defined(PMGRID) || defined(FORCETEST))
    ewald_init();
#endif
    random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(random_generator, 42);        /* start-up seed */
#ifdef PMGRID
    long_range_init();
#endif
    set_random_numbers();

    for(int i = 0; i < 6; i++)
        All.MassTable[i] = 0;
    All.Time = All.TimeBegin;
    if(All.ComovingIntegrationOn){
        All.Timebase_interval = (log(All.TimeMax) - log(All.TimeBegin)) / TIMEBASE;
    }else{
        All.Timebase_interval = (All.TimeMax - All.TimeBegin) / TIMEBASE;
    }
    All.Ti_Current = 0;
    All.NumCurrentTiStep = 0;        /* setup some counters */
    All.SnapshotFileCount = 0;
    All.TotNumOfForces = All.NumForcesSinceLastDomainDecomp = 0;
    All.TimeLastStatistics = All.TimeBegin - All.TimeBetStatistics;
#ifdef PMGRID
    All.PM_Ti_endstep = All.PM_Ti_begstep = 0;
#endif
#ifdef FLEXSTEPS
    All.PresentMinStep = TIMEBASE;
#endif
    All.Ti_nextoutput = -1; // find_next_outputtime(All.Ti_Current);
    open_outputfiles();
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

    t0 = second();
    All.TotNumPart = dm_particles_in_buffer + sph_particles_in_buffer;
    All.TotN_gas = sph_particles_in_buffer;
    All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);        /* sets the maximum number of particles that may */
    All.MaxPartSph = All.PartAllocFactor * (All.TotN_gas / NTask);        /* sets the maximum number of particles that may
                                                                        reside on a processor */
    NumPart = dm_states.size()+sph_states.size();
    N_gas = sph_states.size();
    allocate_memory();

    // initialize sph particles
    i = 0;
    for (map<long long, sph_state>::iterator state_iter = sph_states.begin();
            state_iter != sph_states.end(); state_iter++, i++){
        P[i].ID = (*state_iter).first;
        P[i].Mass = (*state_iter).second.mass;
        P[i].Pos[0] = (*state_iter).second.x;
        P[i].Pos[1] = (*state_iter).second.y;
        P[i].Pos[2] = (*state_iter).second.z;
        P[i].Vel[0] = (*state_iter).second.vx;
        P[i].Vel[1] = (*state_iter).second.vy;
        P[i].Vel[2] = (*state_iter).second.vz;
        P[i].Type = 0; // SPH particles (dark matter particles have type 1)
        SphP[i].Entropy = (*state_iter).second.u;
        SphP[i].Density = -1;
        SphP[i].Hsml = 0;
    }
    sph_states.clear();

    // initialize dark matter particles
    i = N_gas;
    for (map<long long, dynamics_state>::iterator state_iter = dm_states.begin();
            state_iter != dm_states.end(); state_iter++, i++){
        P[i].ID = (*state_iter).first;
        P[i].Mass = (*state_iter).second.mass;
        P[i].Pos[0] = (*state_iter).second.x;
        P[i].Pos[1] = (*state_iter).second.y;
        P[i].Pos[2] = (*state_iter).second.z;
        P[i].Vel[0] = (*state_iter).second.vx;
        P[i].Vel[1] = (*state_iter).second.vy;
        P[i].Vel[2] = (*state_iter).second.vz;
        P[i].Type = 1; // dark matter particles (SPH particles have type 0)
    }
    dm_states.clear();
    All.TimeBegin += All.Ti_Current * All.Timebase_interval;
    All.Ti_Current = 0;
    All.Time = All.TimeBegin;
    set_softenings();
    for(i = 0; i < NumPart; i++){        /*  start-up initialization */
        for(j = 0; j < 3; j++)
            P[i].GravAccel[j] = 0;
#ifdef PMGRID
        for(j = 0; j < 3; j++)
            P[i].GravPM[j] = 0;
#endif
        P[i].Ti_endstep = 0;
        P[i].Ti_begstep = 0;
        P[i].OldAcc = 0;
        P[i].GravCost = 1;
        P[i].Potential = 0;
    }
#ifdef FLEXSTEPS
    for(i = 0; i < NumPart; i++)        /*  start-up initialization */
        P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
#endif
    for(i = 0; i < N_gas; i++){        /* initialize sph_properties */
        for(j = 0; j < 3; j++){
            SphP[i].VelPred[j] = P[i].Vel[j];
            SphP[i].HydroAccel[j] = 0;
        }
        SphP[i].DtEntropy = 0;
    }

    ngb_treeallocate(MAX_NGB);
    if((All.MaxPart < 1000) && (All.TreeAllocFactor <= 1.0)){
        All.TreeAllocFactor = 4000.0/All.MaxPart;
        if (ThisTask == 0){
            cout << "Gadget assumes large numbers of particles while allocating memory. " << endl << "Changed "
                "TreeAllocFactor to " << All.TreeAllocFactor << " to allocate enough memory" << endl <<
                "for this run with " << All.TotNumPart << " particles only." << endl;
        }
    }
    force_treeallocate(All.TreeAllocFactor * 10*All.MaxPart, 10*All.MaxPart);
    All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
    Flag_FullStep = 1;                /* to ensure that Peano-Hilber order is done */
    domain_Decomposition();        /* do initial domain decomposition (gives equal numbers of particles) */
    update_particle_map();
    index_of_highest_mapped_particle = local_index_map.rbegin()->first;
    MPI_Allreduce(MPI_IN_PLACE, &index_of_highest_mapped_particle, 1, MPI_LONG_LONG_INT, MPI_MAX, MPI_COMM_WORLD);
    ngb_treebuild();                /* will build tree */
    setup_smoothinglengths();
    TreeReconstructFlag = 1;
  /* at this point, the entropy variable normally contains the
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly.
   * Once the density has been computed, we can convert thermal energy to entropy.
   */
#ifndef ISOTHERM_EQS
    if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
    for(i = 0; i < N_gas; i++)
        SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
#endif
#ifdef PMGRID
    long_range_init_regionsize();
#endif
    if(All.ComovingIntegrationOn)
        init_drift_table();
    t1 = second();
    CPUThisRun += timediff(t0, t1);
    All.CPU_Total += timediff(t0, t1);

    particles_initialized = true;
    if (ThisTask == 0)
        cout << flush;
    return 0;
}

void push_particle_data_on_state_vectors(){
    map<long long, int>::iterator iter;
    int i;
#ifndef ISOTHERM_EQS
    double a3;

    if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        density();
        density_up_to_date = true;
    }
#endif
    for (iter = local_index_map.begin(); iter != local_index_map.end(); iter++){
        i = (*iter).second;
        if (P[i].Type == 0){
            // store sph particle data
            sph_state state;
            state.mass = P[i].Mass;
            state.x =    P[i].Pos[0];
            state.y =    P[i].Pos[1];
            state.z =    P[i].Pos[2];
            state.vx =   P[i].Vel[0];
            state.vy =   P[i].Vel[1];
            state.vz =   P[i].Vel[2];
#ifdef ISOTHERM_EQS
            state.u = SphP[i].Entropy;
#else
            state.u = SphP[i].Entropy * pow(SphP[i].Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif
            sph_states.insert(std::pair<long long, sph_state>(P[i].ID, state));
        } else {
            // store dark matter particle data
            dynamics_state state;
            state.mass = P[i].Mass;
            state.x =    P[i].Pos[0];
            state.y =    P[i].Pos[1];
            state.z =    P[i].Pos[2];
            state.vx =   P[i].Vel[0];
            state.vy =   P[i].Vel[1];
            state.vz =   P[i].Vel[2];
            dm_states.insert(std::pair<long long, dynamics_state>(P[i].ID, state));
        }
    }
}

int recommit_particles(){
    if (particles_initialized){
        push_particle_data_on_state_vectors();
        free_memory();
        ngb_treefree();
        force_treefree();
    }
    return commit_particles();
}

bool drift_to_t_end(int ti_end){
    bool done;
    int n, min, min_glob, flag, *temp;
    double timeold;
    double t0, t1;
    t0 = second();
    timeold = All.Time;
    for(n = 1, min = P[0].Ti_endstep; n < NumPart; n++)
        if(min > P[n].Ti_endstep)
            min = P[n].Ti_endstep;
    MPI_Allreduce(&min, &min_glob, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
    /* We check whether this is a full step where all particles are synchronized */
    flag = 1;
    for(n = 0; n < NumPart; n++)
        if(P[n].Ti_endstep > min_glob)
            flag = 0;
    MPI_Allreduce(&flag, &Flag_FullStep, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
#ifdef PMGRID
    if(min_glob >= All.PM_Ti_endstep){
        min_glob = All.PM_Ti_endstep;
        Flag_FullStep = 1;
    }
#endif
    /* Determine 'NumForceUpdate', i.e. the number of particles on this processor that are going to be active */
    for(n = 0, NumForceUpdate = 0; n < NumPart; n++){
        if(P[n].Ti_endstep == min_glob)
#ifdef SELECTIVE_NO_GRAVITY
          if(!((1 << P[n].Type) & (SELECTIVE_NO_GRAVITY)))
#endif
            NumForceUpdate++;
    }
    /* note: NumForcesSinceLastDomainDecomp has type "long long" */
    temp = (int*) malloc(NTask * sizeof(int));
    MPI_Allgather(&NumForceUpdate, 1, MPI_INT, temp, 1, MPI_INT, MPI_COMM_WORLD);
    for(n = 0; n < NTask; n++)
        All.NumForcesSinceLastDomainDecomp += temp[n];
    free(temp);
    t1 = second();
    All.CPU_Predict += timediff(t0, t1);
    if (min_glob >= ti_end){
        min_glob = ti_end;
        done = true;
    } else {
        done = false;
    }
    move_particles(All.Ti_Current, min_glob);
    All.Ti_Current = min_glob;
    if(All.ComovingIntegrationOn)
        All.Time = All.TimeBegin * exp(All.Ti_Current * All.Timebase_interval);
    else
        All.Time = All.TimeBegin + All.Ti_Current * All.Timebase_interval;
    All.TimeStep = All.Time - timeold;
    return done;
}

bool check_density_stopping_condition(){
    int stopping_condition_is_set;
    double minimum_density_parameter, maximum_density_parameter;
    get_stopping_condition_minimum_density_parameter(&minimum_density_parameter);
    get_stopping_condition_maximum_density_parameter(&maximum_density_parameter);
    for (int i=0; i<N_gas; i++) {
        if ( (SphP[i].Density < minimum_density_parameter) || 
             (SphP[i].Density > maximum_density_parameter)) {
            int stopping_index  = next_index_for_stopping_condition();
            if (stopping_index >= 0) {
                cout << "set_stopping_condition_info returned: " << 
                    set_stopping_condition_info(stopping_index, DENSITY_LIMIT_DETECTION) << endl;
                cout << "set_stopping_condition_particle_index returned: " << 
                    set_stopping_condition_particle_index(stopping_index, 0, P[i].ID) << endl;
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
    if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
#endif
    
    for (int i=0; i<N_gas; i++) {
#ifdef ISOTHERM_EQS
        internal_energy = SphP[i].Entropy;
#else
        internal_energy = SphP[i].Entropy *
                pow(SphP[i].Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif
        if ( (internal_energy < minimum_internal_energy_parameter) || 
             (internal_energy > maximum_internal_energy_parameter)) {
            int stopping_index  = next_index_for_stopping_condition();
            if (stopping_index > 0) {
                cout << "set_stopping_condition_info returned: " << 
                    set_stopping_condition_info(stopping_index, INTERNAL_ENERGY_LIMIT_DETECTION) << endl;
                cout << "set_stopping_condition_particle_index returned: " << 
                    set_stopping_condition_particle_index(stopping_index, 0, P[i].ID) << endl;
            }
        }
    }
    
    mpi_collect_stopping_conditions();
    is_stopping_condition_set(INTERNAL_ENERGY_LIMIT_DETECTION, &stopping_condition_is_set);
    return stopping_condition_is_set;
}


int evolve_model(double t_end){
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
    
    if (t_end > All.TimeMax)
        return -7;
    ZeroTimestepEncountered = 0;
    Ti_end = (t_end - All.TimeBegin) / All.Timebase_interval;
    if (Ti_end >= All.Ti_Current){
        global_quantities_of_system_up_to_date = density_up_to_date = false;
        done = drift_to_t_end(Ti_end); /* find next synchronization point and drift particles to MIN(this time, t_end). */
        while (!done && All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax) {
            t0 = second();
            every_timestep_stuff();        /* write some info to log-files */
            domain_Decomposition();        /* do domain decomposition if needed */
            particle_map_up_to_date = false;
            compute_accelerations(0);        /* compute accelerations for
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
            if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics) {
#ifdef COMPUTE_POTENTIAL_ENERGY
                compute_potential();
#endif
                energy_statistics();        /* compute and output energy statistics */
                All.TimeLastStatistics += All.TimeBetStatistics;
            }
            advance_and_find_timesteps();        /* 'kick' active particles in
                            * momentum space and compute new
                            * timesteps for them  */
            done = drift_to_t_end(Ti_end);
            All.NumCurrentTiStep++;
            
            /* Check whether we need to interrupt the run */
            MPI_Allreduce(MPI_IN_PLACE, &ZeroTimestepEncountered, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
            if(ZeroTimestepEncountered)
                return -8;
            
            if(ThisTask == 0) {
                /* are we running out of CPU-time ? If yes, interrupt run. */
                if(CPUThisRun > 0.85 * All.TimeLimitCPU){
                    printf("reaching time-limit. stopping.\n");
                    stopflag = 2;
                }
            }
            MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(stopflag)
                return -5;
            
            t1 = second();
            All.CPU_Total += timediff(t0, t1);
            CPUThisRun += timediff(t0, t1);
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
    if (ThisTask == 0)
        cout << flush;
    if (All.Ti_Current > TIMEBASE || All.Time > All.TimeMax)
        return -7;
    return 0;
}

int synchronize_model() {
    return 0;
}

int contruct_tree_if_needed(void){
    double tstart, tend;
    if (!particles_initialized)
        return -1;
    tstart = second();
    if (TreeReconstructFlag){
        if(ThisTask == 0)
            printf("Tree construction.\n");
        force_treebuild(NumPart);
        TreeReconstructFlag = 0;
        if(ThisTask == 0)
            printf("Tree construction done.\n");
    }
    tend = second();
    All.CPU_TreeConstruction += timediff(tstart, tend);
    return 0;
}

int new_dm_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz){
    particle_id_counter++;
    if (ThisTask == 0)
        *id = particle_id_counter;

    // Divide the particles equally over all Tasks, Gadget will redistribute them later.
    if (ThisTask == (dm_particles_in_buffer % NTask)){
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
    particle_id_counter++;
    if (ThisTask == 0)
        *id = particle_id_counter;

    // Divide the sph particles equally over all Tasks, Gadget will redistribute them later.
    if (ThisTask == (sph_particles_in_buffer % NTask)){
        sph_state state;
        state.mass = mass;
        state.x = x;
        state.y = y;
        state.z = z;
        state.vx = vx;
        state.vy = vy;
        state.vz = vz;
        state.u = u;
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
        found = 1 + P[(*it).second].Type; // 1 for sph; 2 for dm
    }
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
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
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    if (found){
        dm_particles_in_buffer--;
        return 0;
    }

    sph_it = sph_states.find(id);
    if (sph_it != sph_states.end()){
        sph_states.erase(sph_it);
        found = 1;
    }
    MPI_Allreduce(MPI_IN_PLACE, &found, 1, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
    if (found){
        sph_particles_in_buffer--;
        return 0;
    }
    return -3;
}

// parameter getters/setters:

int get_time_step(double *timestep){
    if (ThisTask) {return 0;}
    *timestep = All.TimeStep;
    return 0;
}
int set_time_step(double timestep){
    if (ThisTask) {return 0;}
    return -2;
}
int get_epsilon(double *epsilon){
    if (ThisTask) {return 0;}
    set_softenings();
    *epsilon = All.SofteningTable[1];
    return 0;
}
int set_epsilon(double epsilon){
    All.SofteningHalo = epsilon;
    set_softenings();
    return 0;
}
int get_eps2(double *epsilon_squared){
    if (ThisTask) {return 0;}
    return -2;
}
int set_eps2(double epsilon_squared){
    if (ThisTask) {return 0;}
    return -2;
}
int get_epsgas(double *gas_epsilon){
    if (ThisTask) {return 0;}
    set_softenings();
    *gas_epsilon = All.SofteningTable[0];
    return 0;
}
int set_epsgas(double gas_epsilon){
    All.SofteningGas = gas_epsilon;
    set_softenings();
    return 0;
}
int get_unit_mass(double *code_mass_unit){
    if (ThisTask) {return 0;}
    *code_mass_unit = All.UnitMass_in_g;
    return 0;
}
int set_unit_mass(double code_mass_unit){
    All.UnitMass_in_g = code_mass_unit;
    return 0;
}
int get_unit_length(double *code_length_unit){
    if (ThisTask) {return 0;}
    *code_length_unit = All.UnitLength_in_cm;
    return 0;
}
int set_unit_length(double code_length_unit){
    All.UnitLength_in_cm = code_length_unit;
    return 0;
}
int get_unit_time(double *code_time_unit){
    if (ThisTask) {return 0;}
    *code_time_unit = All.UnitTime_in_s;
    return 0;
}
int set_unit_time(double code_time_unit){
    All.UnitVelocity_in_cm_per_s = All.UnitLength_in_cm / code_time_unit;
    All.UnitTime_in_s = code_time_unit;
    return 0;
}
int get_unit_velocity(double *code_velocity_unit){
    if (ThisTask) {return 0;}
    *code_velocity_unit = All.UnitVelocity_in_cm_per_s;
    return 0;
}
int set_unit_velocity(double code_velocity_unit){
    All.UnitVelocity_in_cm_per_s = code_velocity_unit;
    All.UnitTime_in_s = All.UnitLength_in_cm / All.UnitVelocity_in_cm_per_s;
    return 0;
}

int get_gadget_output_directory(char **output_directory){
    if (ThisTask) {return 0;}
    *output_directory = All.OutputDir;
    return 0;
}
int set_gadget_output_directory(char *output_directory){
    int length = strlen(output_directory);
    if (length >= MAXLEN_FILENAME)
        return -1;
    strcpy(All.OutputDir, output_directory);
    if(length > 0)
                if(All.OutputDir[length - 1] != '/')
            strcat(All.OutputDir, "/");
    return 0;
}
int get_nogravity(int *no_gravity_flag){
    if (ThisTask) {return 0;}
#ifdef NOGRAVITY
    *no_gravity_flag = 1;
#else
    *no_gravity_flag = 0;
#endif
    return 0;
}
int get_gdgop(int *gadget_cell_opening_flag){
    if (ThisTask) {return 0;}
    *gadget_cell_opening_flag = All.TypeOfOpeningCriterion;
    return 0;
}
int set_gdgop(int gadget_cell_opening_flag){
    All.TypeOfOpeningCriterion = gadget_cell_opening_flag;
    return 0;
}
int get_isotherm(int *isothermal_flag){
    if (ThisTask) {return 0;}
#ifdef ISOTHERM_EQS
    *isothermal_flag = 1;
#else
    *isothermal_flag = 0;
#endif
    return 0;
}
int get_eps_is_h(int *eps_is_h_flag){
    if (ThisTask) {return 0;}
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) &&  defined(UNEQUALSOFTENINGS)
    *eps_is_h_flag = 1;
#else
    *eps_is_h_flag = 0;
#endif
    return 0;
}
int get_nsmooth(int *nsmooth){
    if (ThisTask) {return 0;}
    *nsmooth = All.DesNumNgb;
    return 0;
}
int set_nsmooth(int nsmooth){
    All.DesNumNgb = nsmooth;
    return 0;
}
int get_bh_tol(double *opening_angle){
    if (ThisTask) {return 0;}
    *opening_angle = All.ErrTolTheta;
    return 0;
}
int set_bh_tol(double opening_angle){
    All.ErrTolTheta = opening_angle;
    return 0;
}
int get_gdgtol(double *gadget_cell_opening_constant){
    if (ThisTask) {return 0;}
    *gadget_cell_opening_constant = All.ErrTolForceAcc;
    return 0;
}
int set_gdgtol(double gadget_cell_opening_constant){
    All.ErrTolForceAcc = gadget_cell_opening_constant;
    return 0;
}
int get_gamma(double *gamma){
    if (ThisTask) {return 0;}
    *gamma = GAMMA;
    return 0;
}
int get_alpha(double *artificial_viscosity_alpha){
    if (ThisTask) {return 0;}
    *artificial_viscosity_alpha = All.ArtBulkViscConst;
    return 0;
}
int set_alpha(double artificial_viscosity_alpha){
    All.ArtBulkViscConst = artificial_viscosity_alpha;
    return 0;
}
int get_courant(double *courant){
    if (ThisTask) {return 0;}
    *courant = All.CourantFac*2.0;
    return 0;
}
int set_courant(double courant){
    All.CourantFac = courant/2.0;
    return 0;
}
int get_nsmtol(double *n_neighbour_tol){
    if (ThisTask) {return 0;}
    *n_neighbour_tol = All.MaxNumNgbDeviation / All.DesNumNgb;
    return 0;
}
int set_nsmtol(double n_neighbour_tol){
    All.MaxNumNgbDeviation = n_neighbour_tol * All.DesNumNgb;
    return 0;
}

int get_energy_file(char **energy_file){
    if (ThisTask) {return 0;}
    *energy_file = All.EnergyFile;
    return 0;
}
int set_energy_file(char *energy_file){
    strcpy(All.EnergyFile, energy_file);
    return 0;
}
int get_info_file(char **info_file){
    if (ThisTask) {return 0;}
    *info_file = All.InfoFile;
    return 0;
}
int set_info_file(char *info_file){
    strcpy(All.InfoFile, info_file);
    return 0;
}
int get_timings_file(char **timings_file){
    if (ThisTask) {return 0;}
    *timings_file = All.TimingsFile;
    return 0;
}
int set_timings_file(char *timings_file){
    strcpy(All.TimingsFile, timings_file);
    return 0;
}
int get_cpu_file(char **cpu_file){
    if (ThisTask) {return 0;}
    *cpu_file = All.CpuFile;
    return 0;
}
int set_cpu_file(char *cpu_file){
    strcpy(All.CpuFile, cpu_file);
    return 0;
}

int get_time_limit_cpu(double *time_limit_cpu){
    if (ThisTask) {return 0;}
    *time_limit_cpu = All.TimeLimitCPU;
    return 0;
}
int set_time_limit_cpu(double time_limit_cpu){
    All.TimeLimitCPU = time_limit_cpu;
    return 0;
}
int get_comoving_integration_flag(int *comoving_integration_flag){
    if (ThisTask) {return 0;}
    *comoving_integration_flag = All.ComovingIntegrationOn;
    return 0;
}
int set_comoving_integration_flag(int comoving_integration_flag){
    All.ComovingIntegrationOn = comoving_integration_flag;
    return 0;
}
int get_type_of_timestep_criterion(int *type_of_timestep_criterion){
    if (ThisTask) {return 0;}
    *type_of_timestep_criterion = All.TypeOfTimestepCriterion;
    return 0;
}
int set_type_of_timestep_criterion(int type_of_timestep_criterion){
    All.TypeOfTimestepCriterion = type_of_timestep_criterion;
    return 0;
}
int get_time_begin(double *time_begin){
    if (ThisTask) {return 0;}
    *time_begin = All.TimeBegin;
    return 0;
}
int set_time_begin(double time_begin){
    All.TimeBegin = time_begin;
    return 0;
}
int get_time_max(double *time_max){
    if (ThisTask) {return 0;}
    *time_max = All.TimeMax;
    return 0;
}
int set_time_max(double time_max){
    All.TimeMax = time_max;
    return 0;
}
int get_omega_zero(double *omega_zero){
    if (ThisTask) {return 0;}
    *omega_zero = All.Omega0;
    return 0;
}
int set_omega_zero(double omega_zero){
    All.Omega0 = omega_zero;
    return 0;
}
int get_omega_lambda(double *omega_lambda){
    if (ThisTask) {return 0;}
    *omega_lambda = All.OmegaLambda;
    return 0;
}
int set_omega_lambda(double omega_lambda){
    All.OmegaLambda = omega_lambda;
    return 0;
}
int get_omega_baryon(double *omega_baryon){
    if (ThisTask) {return 0;}
    *omega_baryon = All.OmegaBaryon;
    return 0;
}
int set_omega_baryon(double omega_baryon){
    All.OmegaBaryon = omega_baryon;
    return 0;
}
int get_hubble_param(double *hubble_param){
    if (ThisTask) {return 0;}
    *hubble_param = All.HubbleParam;
    return 0;
}
int set_hubble_param(double hubble_param){
    All.HubbleParam = hubble_param;
    return 0;
}
int get_err_tol_int_accuracy(double *err_tol_int_accuracy){
    if (ThisTask) {return 0;}
    *err_tol_int_accuracy = All.ErrTolIntAccuracy;
    return 0;
}
int set_err_tol_int_accuracy(double err_tol_int_accuracy){
    All.ErrTolIntAccuracy = err_tol_int_accuracy;
    return 0;
}
int get_max_size_timestep(double *max_size_timestep){
    if (ThisTask) {return 0;}
    *max_size_timestep = All.MaxSizeTimestep;
    return 0;
}
int set_max_size_timestep(double max_size_timestep){
    All.MaxSizeTimestep = max_size_timestep;
    return 0;
}
int get_min_size_timestep(double *min_size_timestep){
    if (ThisTask) {return 0;}
    *min_size_timestep = All.MinSizeTimestep;
    return 0;
}
int set_min_size_timestep(double min_size_timestep){
    All.MinSizeTimestep = min_size_timestep;
    return 0;
}
int get_tree_domain_update_frequency(double *tree_domain_update_frequency){
    if (ThisTask) {return 0;}
    *tree_domain_update_frequency = All.TreeDomainUpdateFrequency;
    return 0;
}
int set_tree_domain_update_frequency(double tree_domain_update_frequency){
    All.TreeDomainUpdateFrequency = tree_domain_update_frequency;
    return 0;
}
int get_time_between_statistics(double *time_between_statistics){
    if (ThisTask) {return 0;}
    *time_between_statistics = All.TimeBetStatistics;
    return 0;
}
int set_time_between_statistics(double time_between_statistics){
    All.TimeBetStatistics = time_between_statistics;
    return 0;
}
int get_min_gas_temp(double *min_gas_temp){
    if (ThisTask) {return 0;}
    *min_gas_temp = All.MinGasTemp;
    return 0;
}
int set_min_gas_temp(double min_gas_temp){
    All.MinGasTemp = min_gas_temp;
    return 0;
}
int get_min_gas_hsmooth_fractional(double *min_gas_hsmooth_fractional){
    if (ThisTask) {return 0;}
    *min_gas_hsmooth_fractional = All.MinGasHsmlFractional;
    return 0;
}
int set_min_gas_hsmooth_fractional(double min_gas_hsmooth_fractional){
    All.MinGasHsmlFractional = min_gas_hsmooth_fractional;
    return 0;
}
int get_softening_gas_max_phys(double *softening_gas_max_phys){
    if (ThisTask) {return 0;}
    *softening_gas_max_phys = All.SofteningGasMaxPhys;
    return 0;
}
int set_softening_gas_max_phys(double softening_gas_max_phys){
    All.SofteningGasMaxPhys = softening_gas_max_phys;
    return 0;
}
int get_softening_halo_max_phys(double *softening_halo_max_phys){
    if (ThisTask) {return 0;}
    *softening_halo_max_phys = All.SofteningHaloMaxPhys;
    return 0;
}
int set_softening_halo_max_phys(double softening_halo_max_phys){
    All.SofteningHaloMaxPhys = softening_halo_max_phys;
    return 0;
}

int get_box_size(double *value)
{
    if (ThisTask) {return 0;}
    *value = All.BoxSize;
    return 0;
}

int set_box_size(double value)
{
    All.BoxSize = value;
    return 0;
}

int get_periodic_boundaries_flag(int *value)
{
    if (ThisTask) {return 0;}
    *value = All.PeriodicBoundariesOn;
    return 0;
}

int set_periodic_boundaries_flag(int value)
{
    All.PeriodicBoundariesOn = value;
    return 0;
}

int get_interpret_kicks_as_feedback_flag(int *value)
{
    if (ThisTask) {return 0;}
    *value = interpret_kicks_as_feedback;
    return 0;
}

int set_interpret_kicks_as_feedback_flag(int value)
{
    interpret_kicks_as_feedback = value;
    return 0;
}

int get_interpret_heat_as_feedback_flag(int *value) {
    if (ThisTask) {return 0;}
    *value = interpret_heat_as_feedback;
    return 0;
}

int set_interpret_heat_as_feedback_flag(int value) {
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

    if (ThisTask == 0){
        MPI_Reduce(MPI_IN_PLACE, &next_local_index, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        *index_of_the_next_particle = next_local_index;
        if (next_local_index < index_of_highest_mapped_particle){
            return 0;
        } else if (next_local_index == index_of_highest_mapped_particle){
            return 1;
        } else {
            return -1;
        }
    } else {
        MPI_Reduce(&next_local_index, NULL, 1, MPI_LONG_LONG_INT, MPI_MIN, 0, MPI_COMM_WORLD);
        return 0;
    }
}

void update_particle_map(void){
    local_index_map.clear();
    for(int i = 0; i < NumPart; i++) {
        local_index_map.insert(std::pair<long long, int>(P[i].ID, i));
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
            buffer[i] = P[local_index].Mass;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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
    if(ThisTask) {
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        return 0;
    } else {
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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
            P[local_index].Mass = mass[i];
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

int get_position(int *index, double *x, double *y, double *z, int length){
    int errors = 0;
    double *buffer = new double[length*3];
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = P[local_index].Pos[0];
            buffer[i+length] = P[local_index].Pos[1];
            buffer[i+2*length] = P[local_index].Pos[2];
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length*3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length*3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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

int set_position(int *index, double *x, double *y, double *z, int length){
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            P[local_index].Pos[0] = x[i];
            P[local_index].Pos[1] = y[i];
            P[local_index].Pos[2] = z[i];
            count[i] = 1;
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}

int get_velocity(int *index, double *vx, double *vy, double *vz, int length){
    int errors = 0;
    double *buffer = new double[length*3];
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = P[local_index].Vel[0];
            buffer[i+length] = P[local_index].Vel[1];
            buffer[i+2*length] = P[local_index].Vel[2];
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length*3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length*3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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

int set_velocity(int *index, double *vx, double *vy, double *vz, int length){
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            P[local_index].Vel[0] = vx[i];
            P[local_index].Vel[1] = vy[i];
            P[local_index].Vel[2] = vz[i];
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_kicks_as_feedback && P[local_index].Type == 0) {
                SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if(interpret_kicks_as_feedback && P[local_index].Type == 0 && P[local_index].Ti_endstep != All.Ti_Current) {
                make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}

int get_state(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length) {
    int errors = 0;
    double *buffer = new double[length*7];
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = P[local_index].Mass;
            buffer[i+length] = P[local_index].Pos[0];
            buffer[i+2*length] = P[local_index].Pos[1];
            buffer[i+3*length] = P[local_index].Pos[2];
            buffer[i+4*length] = P[local_index].Vel[0];
            buffer[i+5*length] = P[local_index].Vel[1];
            buffer[i+6*length] = P[local_index].Vel[2];
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
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length*7, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length*7, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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

int set_state(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, int length){
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            P[local_index].Mass = mass[i];
            P[local_index].Pos[0] = x[i];
            P[local_index].Pos[1] = y[i];
            P[local_index].Pos[2] = z[i];
            P[local_index].Vel[0] = vx[i];
            P[local_index].Vel[1] = vy[i];
            P[local_index].Vel[2] = vz[i];
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_kicks_as_feedback && P[local_index].Type == 0) {
                SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if(interpret_kicks_as_feedback && P[local_index].Type == 0 && P[local_index].Ti_endstep != All.Ti_Current) {
                make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}

int get_state_sph(int *index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *internal_energy, int length) {
    int errors = 0;
    double *buffer = new double[length*8];
    int *count = new int[length];
    int local_index;
#ifndef ISOTHERM_EQS
    double a3;

    if (!density_up_to_date){
        density();
        density_up_to_date = true;
    }
    if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = P[local_index].Mass;
            buffer[i+length] = P[local_index].Pos[0];
            buffer[i+2*length] = P[local_index].Pos[1];
            buffer[i+3*length] = P[local_index].Pos[2];
            buffer[i+4*length] = P[local_index].Vel[0];
            buffer[i+5*length] = P[local_index].Vel[1];
            buffer[i+6*length] = P[local_index].Vel[2];
#ifdef ISOTHERM_EQS
            buffer[i+7*length] = SphP[local_index].Entropy;
#else
            buffer[i+7*length] = SphP[local_index].Entropy *
                pow(SphP[local_index].Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
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
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length*8, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length*8, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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

int set_state_sph(int *index, double *mass, double *x, double *y, double *z,
        double *vx, double *vy, double *vz, double *internal_energy, int length){
    int *count = new int[length];
    int local_index;
#ifndef ISOTHERM_EQS
    double a3;

    if (!density_up_to_date){
        density();
        density_up_to_date = true;
    }
    if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && P[local_index].Type == 0){
            P[local_index].Mass = mass[i];
            P[local_index].Pos[0] = x[i];
            P[local_index].Pos[1] = y[i];
            P[local_index].Pos[2] = z[i];
            P[local_index].Vel[0] = vx[i];
            P[local_index].Vel[1] = vy[i];
            P[local_index].Vel[2] = vz[i];
#ifdef ISOTHERM_EQS
            SphP[local_index].Entropy = internal_energy[i];
#else
            SphP[local_index].Entropy = GAMMA_MINUS1 * internal_energy[i] /
                pow(SphP[local_index].Density / a3, GAMMA_MINUS1);
#endif
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_heat_as_feedback || interpret_kicks_as_feedback) {
                SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if ((interpret_heat_as_feedback || interpret_kicks_as_feedback) && 
                    P[local_index].Ti_endstep != All.Ti_Current) {
                make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}

int get_acceleration(int *index, double * ax, double * ay, double * az, int length){
    int errors = 0;
    double *buffer = new double[length*3];
    int *count = new int[length];
    int local_index;

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index)){
            count[i] = 1;
            buffer[i] = P[local_index].GravAccel[0];
            buffer[i+length] = P[local_index].GravAccel[1];
            buffer[i+2*length] = P[local_index].GravAccel[2];
            if(P[local_index].Type == 0){
                buffer[i] += SphP[local_index].HydroAccel[0];
                buffer[i+length] += SphP[local_index].HydroAccel[1];
                buffer[i+2*length] += SphP[local_index].HydroAccel[2];
            }
        } else {
            count[i] = 0;
            buffer[i] = 0;
            buffer[i+length] = 0;
            buffer[i+2*length] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length*3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length*3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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
    
    if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        density();
        density_up_to_date = true;
    }
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && P[local_index].Type == 0){
            count[i] = 1;
#ifdef ISOTHERM_EQS
            buffer[i] = SphP[local_index].Entropy;
#else
            buffer[i] = SphP[local_index].Entropy *
                pow(SphP[local_index].Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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
    
    if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        density();
        density_up_to_date = true;
    }
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && P[local_index].Type == 0){
#ifdef ISOTHERM_EQS
            SphP[local_index].Entropy = internal_energy[i];
#else
            SphP[local_index].Entropy = GAMMA_MINUS1 * internal_energy[i] /
                pow(SphP[local_index].Density / a3, GAMMA_MINUS1);
#endif
            count[i] = 1;
#ifdef TIMESTEP_UPDATE
            if (interpret_heat_as_feedback) {
                SphP[local_index].FeedbackFlag = 2;
            }
#endif
#ifdef TIMESTEP_LIMITER
            if(interpret_heat_as_feedback && P[local_index].Ti_endstep != All.Ti_Current) {
                make_it_active(local_index);
            }
#endif
        } else count[i] = 0;
    }
    global_quantities_of_system_up_to_date = false;
    return check_counts_and_free(count, length);
}

int get_smoothing_length(int *index, double *smoothing_length, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    if (!density_up_to_date){
        density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = SphP[local_index].Hsml;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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

int get_density(int *index, double *density_out, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    if (!density_up_to_date){
        density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = SphP[local_index].Density;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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
int get_pressure(int *index, double *pressure_out, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    if (!density_up_to_date){
        density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = SphP[local_index].Pressure;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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
int get_d_internal_energy_dt(int *index, double *d_internal_energy_dt_out, int length){
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;
#ifndef ISOTHERM_EQS
    double a3;
    
    if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
    if (!density_up_to_date){
        density();
        density_up_to_date = true;
    }
#endif
    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && P[local_index].Type == 0){
            count[i] = 1;
#ifdef ISOTHERM_EQS
            buffer[i] = SphP[local_index].DtEntropy;
#else
            buffer[i] = - SphP[local_index].Pressure * SphP[local_index].DivVel /
                SphP[local_index].Density;
#endif
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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
        density();
        density_up_to_date = true;
    }

    for (int i = 0; i < length; i++){
        if(found_particle(index[i], &local_index) && P[local_index].Type == 0){
            count[i] = 1;
            buffer[i] = SphP[local_index].NumNgb;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }
    if(ThisTask) {
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
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
    set_softenings();
    if (ThisTask) {return 0;}
    for (int i = 0; i < length; i++)
        epsilon[i] = All.SofteningTable[1];
    return 0;
}
int get_epsilon_gas_part(int *index, double *epsilon, int length){
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) &&  defined(UNEQUALSOFTENINGS)
    return get_smoothing_length(index, epsilon, length);
#else
    set_softenings();
    if (ThisTask) {return 0;}
    for (int i = 0; i < length; i++)
        epsilon[i] = All.SofteningTable[0];
    return 0;
#endif
}



// simulation property getters:

void update_global_quantities(bool do_potential){
    if (do_potential) {
        compute_potential();
        potential_energy_also_up_to_date = true;
    } else {potential_energy_also_up_to_date = false;}
    compute_global_quantities_of_system();
    global_quantities_of_system_up_to_date = true;
}
int get_time(double *time){
    if (ThisTask) {return 0;}
    *time = All.Time;
    return 0;
}
int get_total_radius(double *radius){
    double r_squared, local_max = 0;
    int i, j;

    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    for (i = 0; i < NumPart; i++){
        for (r_squared = 0, j = 0; j < 3; j++)
            r_squared += (SysState.CenterOfMass[j]-P[i].Pos[j])*(SysState.CenterOfMass[j]-P[i].Pos[j]);
        if (r_squared > local_max)
            local_max = r_squared;
    }

    if(ThisTask) {
        MPI_Reduce(&local_max, NULL, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, &local_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        *radius = sqrt(local_max);
    }
    return 0;
}

int get_total_mass(double *mass){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (ThisTask) {return 0;}
    *mass = SysState.Mass;
    return 0;
}

int get_potential(int *index, double *potential, int length) {
    int errors = 0;
    double *buffer = new double[length];
    int *count = new int[length];
    int local_index;

    if (!potential_energy_also_up_to_date) {
        compute_potential();
        potential_energy_also_up_to_date = true;
    }

    for (int i = 0; i < length; i++) {
        if (found_particle(index[i], &local_index)) {
            count[i] = 1;
            buffer[i] = P[local_index].Potential;
        } else {
            count[i] = 0;
            buffer[i] = 0;
        }
    }

    if (ThisTask) {
        MPI_Reduce(buffer, NULL, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(count, NULL, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, buffer, length, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, count, length, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        for (int i = 0; i < length; i++) {
            if (count[i] != 1){
                errors++;
                potential[i] = 0;
            } else {
                potential[i] = buffer[i];
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
    if (ThisTask) {return 0;}
    *kinetic_energy = SysState.EnergyKin;
    return 0;
}
int get_potential_energy(double *potential_energy){
    if (!(global_quantities_of_system_up_to_date && potential_energy_also_up_to_date))
        update_global_quantities(true);
    if (ThisTask) {return 0;}
    *potential_energy = SysState.EnergyPot;
    return 0;
}
int get_thermal_energy(double *thermal_energy){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (ThisTask) {return 0;}
    *thermal_energy = SysState.EnergyInt;
    return 0;
}
int get_number_of_particles(int *number_of_particles){
    if (ThisTask) {return 0;}
    *number_of_particles = All.TotNumPart;
    return 0;
}
int get_center_of_mass_position(double *x, double *y, double *z){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (ThisTask) {return 0;}
    *x = SysState.CenterOfMass[0];
    *y = SysState.CenterOfMass[1];
    *z = SysState.CenterOfMass[2];
    return 0;
}
int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    if (ThisTask) {return 0;}
    *vx = SysState.Momentum[0]/SysState.Mass;
    *vy = SysState.Momentum[1]/SysState.Mass;
    *vz = SysState.Momentum[2]/SysState.Mass;
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
#ifndef ISOTHERM_EQS
    double a3;
#endif
    error = contruct_tree_if_needed();
    if (error)
        return error;
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    vel[0] = vx;
    vel[1] = vy;
    vel[2] = vz;
    hydro_state_at_point(pos, vel, &h_out, &ngb_out, &dhsml_out, &rho_out, rhov_out, &rhov2_out, &rhoe_out);
    if (ThisTask) {return 0;}
    *rho   = rho_out;
    *rhovx = rhov_out[0];
    *rhovy = rhov_out[1];
    *rhovz = rhov_out[2];
#ifdef ISOTHERM_EQS
    *rhoe = rhoe_out + (rhov_out[0]*rhov_out[0] + rhov_out[1]*rhov_out[1] + rhov_out[2]*rhov_out[2]) / rho_out;
#else
    if (All.ComovingIntegrationOn) {a3 = All.Time * All.Time * All.Time;} else {a3 = 1;}
    *rhoe = rhoe_out * (pow(rho_out / a3, GAMMA_MINUS1) / GAMMA_MINUS1) + rhov2_out;
#endif
    return 0;
}
