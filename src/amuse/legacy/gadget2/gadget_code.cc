#include <mpi.h>
#include <iostream>
#include <string.h>
#include <vector>
#include <math.h>
#include "gadget_code.h"
#include "worker_code.h"
//AMUSE STOPPING CONDITIONS
#include "stopcond.h"

//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>


using namespace std;

const bool debug = true;

bool particles_initialized = false;
bool outfiles_opened = false;
bool global_quantities_of_system_up_to_date = false;
bool potential_energy_also_up_to_date = false;
bool density_up_to_date = false;
int particles_counter = 0;
vector<dynamics_state> ds;        // for initialization only
vector<sph_state> sph_ds;         // for initialization only

// global static parameters




// Interface functions:


// general interface functions:

int set_parameterfile_path(char *parameterfile_path){
    if (strlen(parameterfile_path) >= MAXLEN_FILENAME)
        return -1;
    strcpy(ParameterFile, parameterfile_path);
    if (debug)
        cout << "Parameters will be read from: " << ParameterFile << endl;
    return 0;
}

int initialize_code(){
    double t0, t1;
    char command[200];
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
    read_parameter_file(ParameterFile);	/* ... read in parameters for this run */
    sprintf(command, "rm -f %s%s", All.OutputDir, "parameters-usedvalues");
    system(command);
    allocate_commbuffers();	/* ... allocate buffer-memory for particle 
				   exchange during force computation */
    set_units();
#if defined(PERIODIC) && (!defined(PMGRID) || defined(FORCETEST))
    ewald_init();
#endif
    random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
    gsl_rng_set(random_generator, 42);	/* start-up seed */
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
    All.NumCurrentTiStep = 0;	/* setup some counters */
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
    
    t1 = second();
    CPUThisRun += timediff(t0, t1);
    All.CPU_Total += timediff(t0, t1);

    //AMUSE STOPPING CONDITIONS SUPPORT
    set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
    
    return 0;
}

int cleanup_code(){
    if (outfiles_opened)
        close_outputfiles();
    if (particles_initialized)
        free_memory();
    return 0;
}

int commit_parameters(){
    char command[200];
    sprintf(command, "mv -f %s%s %s%s", ParameterFile, "-usedvalues", All.OutputDir, "parameters-usedvalues");
    if (strlen(command) > 198){
        cerr << "Error: ParameterFile and/or OutputDir names too long. " << command << endl;
        return -1;
    }
    system(command);
    open_outputfiles();
    outfiles_opened = true;
    return 0;
}
int recommit_parameters(){
    return 0;
}

int commit_particles(){
    double t0, t1, a3;
    int i, j;
    t0 = second();
    All.TotNumPart = ds.size()+sph_ds.size();
    All.TotN_gas = sph_ds.size();
    NumPart = All.TotNumPart;
    N_gas = All.TotN_gas;
    All.MaxPart = All.PartAllocFactor * (All.TotNumPart / NTask);	/* sets the maximum number of particles that may */
    All.MaxPartSph = All.PartAllocFactor * (All.TotN_gas / NTask);	/* sets the maximum number of particles that may 
                                                                        reside on a processor */
    allocate_memory();
    for(i = 0; i < N_gas; i++){	/* initialize sph particles */
        P[i].ID = sph_ds.at(i).id;
        P[i].Mass = sph_ds.at(i).mass;
        P[i].Pos[0] = sph_ds.at(i).x;
        P[i].Pos[1] = sph_ds.at(i).y;
        P[i].Pos[2] = sph_ds.at(i).z;
        P[i].Vel[0] = sph_ds.at(i).vx;
        P[i].Vel[1] = sph_ds.at(i).vy;
        P[i].Vel[2] = sph_ds.at(i).vz;
        P[i].Type = 0; // SPH particles (dark matter particles have type 1)
        SphP[i].Entropy = sph_ds.at(i).u;
        SphP[i].Density = -1;
        SphP[i].Hsml = 0;
    }
    for (i = N_gas; i < NumPart; i++){	/* initialize dark matter particles */
        j=i-N_gas;
        P[i].ID = ds.at(j).id;
        P[i].Mass = ds.at(j).mass;
        P[i].Pos[0] = ds.at(j).x;
        P[i].Pos[1] = ds.at(j).y;
        P[i].Pos[2] = ds.at(j).z;
        P[i].Vel[0] = ds.at(j).vx;
        P[i].Vel[1] = ds.at(j).vy;
        P[i].Vel[2] = ds.at(j).vz;
        P[i].Type = 1; // dark matter particles (SPH particles have type 0)
    }
    All.Ti_Current = 0;
    set_softenings();
    for(i = 0; i < NumPart; i++){	/*  start-up initialization */
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
    for(i = 0; i < NumPart; i++)	/*  start-up initialization */
        P[i].FlexStepGrp = (int) (TIMEBASE * get_random_number(P[i].ID));
#endif
    for(i = 0; i < N_gas; i++){	/* initialize sph_properties */
        for(j = 0; j < 3; j++){
            SphP[i].VelPred[j] = P[i].Vel[j];
            SphP[i].HydroAccel[j] = 0;
        }
        SphP[i].DtEntropy = 0;
    }
    
    ngb_treeallocate(MAX_NGB);
    if((All.MaxPart < 1000) && (All.TreeAllocFactor <= 1.0)){
        All.TreeAllocFactor = 4000.0/All.MaxPart;
        cout << "Gadget assumes large numbers of particles while allocating memory. " << endl << "Changed "
            "TreeAllocFactor to " << All.TreeAllocFactor << " to allocate enough memory" << endl << 
            "for this run with " << All.TotNumPart << " particles only." << endl;
    }
    force_treeallocate(All.TreeAllocFactor * 10*All.MaxPart, 10*All.MaxPart);
    All.NumForcesSinceLastDomainDecomp = 1 + All.TotNumPart * All.TreeDomainUpdateFrequency;
    Flag_FullStep = 1;		/* to ensure that Peano-Hilber order is done */
    domain_Decomposition();	/* do initial domain decomposition (gives equal numbers of particles) */
    ngb_treebuild();		/* will build tree */
    setup_smoothinglengths();
    TreeReconstructFlag = 1;
  /* at this point, the entropy variable normally contains the 
   * internal energy, read in from the initial conditions file, unless the file
   * explicitly signals that the initial conditions contain the entropy directly. 
   * Once the density has been computed, we can convert thermal energy to entropy.
   */
#ifndef ISOTHERM_EQS
//  if(header.flag_entropy_instead_u == 0){
    if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
    for(i = 0; i < N_gas; i++)
        SphP[i].Entropy = GAMMA_MINUS1 * SphP[i].Entropy / pow(SphP[i].Density / a3, GAMMA_MINUS1);
//  }
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
    return 0;
}
int recommit_particles(){
    struct particle_data Pcurrent;
    struct sph_particle_data SphPcurrent;
    double a3;
    if (particles_initialized){
        for (vector<dynamics_state>::iterator state_iter = ds.begin(); state_iter != ds.end(); state_iter++){
            if(!(find_particle((*state_iter).id, &Pcurrent))){
                (*state_iter).mass = Pcurrent.Mass;
                (*state_iter).x = Pcurrent.Pos[0];
                (*state_iter).y = Pcurrent.Pos[1];
                (*state_iter).z = Pcurrent.Pos[2];
                (*state_iter).vx = Pcurrent.Vel[0];
                (*state_iter).vy = Pcurrent.Vel[1];
                (*state_iter).vz = Pcurrent.Vel[2];
            }
        }
        if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
        for (vector<sph_state>::iterator state_iter = sph_ds.begin(); state_iter != sph_ds.end(); state_iter++){
            if(!(find_particle((*state_iter).id, &Pcurrent)) && !(find_sph_particle((*state_iter).id, &SphPcurrent))){
                (*state_iter).mass = Pcurrent.Mass;
                (*state_iter).x = Pcurrent.Pos[0];
                (*state_iter).y = Pcurrent.Pos[1];
                (*state_iter).z = Pcurrent.Pos[2];
                (*state_iter).vx = Pcurrent.Vel[0];
                (*state_iter).vy = Pcurrent.Vel[1];
                (*state_iter).vz = Pcurrent.Vel[2];
                (*state_iter).u = SphPcurrent.Entropy * pow(SphPcurrent.Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
            }
        }
        free_memory();
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
int evolve(double t_end){
    bool done;
    double t0, t1;
    int Ti_end, stopflag = 0;
    // AMUSE STOPPING CONDITIONS
    int is_number_of_steps_detection_enabled;
    int max_number_of_steps;
    int number_of_steps_innerloop = 0;
    int error;

    error = is_stopping_condition_enabled(NUMBER_OF_STEPS_DETECTION, 
					  &is_number_of_steps_detection_enabled);
    get_stopping_condition_number_of_steps_parameter(&max_number_of_steps);

    // .......

    Ti_end = (t_end - All.TimeBegin) / All.Timebase_interval;
    if (Ti_end >= All.Ti_Current){
        global_quantities_of_system_up_to_date = density_up_to_date = false;
        done = drift_to_t_end(Ti_end); /* find next synchronization point and drift particles to MIN(this time, t_end). */
        while (!done && All.Ti_Current < TIMEBASE && All.Time <= All.TimeMax) {
            t0 = second();
            every_timestep_stuff();	/* write some info to log-files */
            domain_Decomposition();	/* do domain decomposition if needed */
            compute_accelerations(0);	/* compute accelerations for 
                        * the particles that are to be advanced  */
            /* check whether we want a full energy statistics */
            if((All.Time - All.TimeLastStatistics) >= All.TimeBetStatistics) {
#ifdef COMPUTE_POTENTIAL_ENERGY
                compute_potential();
#endif
                energy_statistics();	/* compute and output energy statistics */
                All.TimeLastStatistics += All.TimeBetStatistics;
            }
            advance_and_find_timesteps();	/* 'kick' active particles in
                            * momentum space and compute new
                            * timesteps for them  */
            done = drift_to_t_end(Ti_end);
            All.NumCurrentTiStep++;
            
            /* Check whether we need to interrupt the run */
            if(ThisTask == 0) {
                /* are we running out of CPU-time ? If yes, interrupt run. */
                if(CPUThisRun > 0.85 * All.TimeLimitCPU){printf("reaching time-limit. stopping.\n"); stopflag = 2;}
            }
            MPI_Bcast(&stopflag, 1, MPI_INT, 0, MPI_COMM_WORLD);
            if(stopflag) return -1;
            
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
    } else {return -1;}
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
    dynamics_state state;
    state.id = ++particles_counter;
    state.mass = mass;
    state.x = x;
    state.y = y;
    state.z = z;
    state.vx = vx;
    state.vy = vy;
    state.vz = vz;
    ds.push_back(state);
    *id = state.id;
    return 0;
}
int new_sph_particle(int *id, double mass, double x, double y, double z, double vx, double vy, double vz, double u){
    sph_state state;
    state.id = ++particles_counter;
    state.mass = mass;
    state.x = x;
    state.y = y;
    state.z = z;
    state.vx = vx;
    state.vy = vy;
    state.vz = vz;
    state.u = u;
    sph_ds.push_back(state);
    *id = state.id;
    return 0;
}
int delete_particle(int id){
    for (vector<dynamics_state>::iterator state_iter = ds.begin(); state_iter != ds.end(); state_iter++){
        if ((*state_iter).id == id){
            ds.erase(state_iter);
            return 0;
        }
    }
    for (vector<sph_state>::iterator state_iter = sph_ds.begin(); state_iter != sph_ds.end(); state_iter++){
        if ((*state_iter).id == id){
            sph_ds.erase(state_iter);
            return 0;
        }
    }
    return -1;
}


// parameter getters/setters:

int get_time_step(double *timestep){
    *timestep = All.TimeStep;
    return 0;
}
int set_time_step(double timestep){
    All.TimeStep = timestep;
    return 0;
}
int get_epsilon(double *epsilon){
    *epsilon = All.SofteningHalo;
    return 0;
}
int set_epsilon(double epsilon){
    All.SofteningHalo = epsilon;
    return 0;
}
int get_eps2(double *epsilon_squared){
    return -2;
}
int set_eps2(double epsilon_squared){
    return -2;
}
int get_epsgas(double *gas_epsilon){
    *gas_epsilon = All.SofteningGas;
    return 0;
}
int set_epsgas(double gas_epsilon){
    All.SofteningGas = gas_epsilon;
    return 0;
}
int get_unit_mass(double *code_mass_unit){
    *code_mass_unit = All.UnitMass_in_g;
    return 0;
}
int set_unit_mass(double code_mass_unit){
    All.UnitMass_in_g = code_mass_unit;
    set_units();
    return 0;
}
int get_unit_length(double *code_length_unit){
    *code_length_unit = All.UnitLength_in_cm;
    return 0;
}
int set_unit_length(double code_length_unit){
    All.UnitLength_in_cm = code_length_unit;
    set_units();
    return 0;
}
int get_unit_time(double *code_time_unit){
    *code_time_unit = All.UnitTime_in_s;
    return 0;
}
int set_unit_time(double code_time_unit){
    All.UnitVelocity_in_cm_per_s = All.UnitLength_in_cm / code_time_unit;
    set_units();
    return 0;
}
int get_unit_velocity(double *code_velocity_unit){
    *code_velocity_unit = All.UnitVelocity_in_cm_per_s;
    return 0;
}
int set_unit_velocity(double code_velocity_unit){
    All.UnitVelocity_in_cm_per_s = code_velocity_unit;
    set_units();
    return 0;
}

int get_gadget_output_directory(char **output_directory){
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
#ifdef NOGRAVITY
    *no_gravity_flag = 1;
#else
    *no_gravity_flag = 0;
#endif
    return 0;
}
int get_gdgop(int *gadget_cell_opening_flag){
    *gadget_cell_opening_flag = All.TypeOfOpeningCriterion;
    return 0;
}
int set_gdgop(int gadget_cell_opening_flag){
    All.TypeOfOpeningCriterion = gadget_cell_opening_flag;
    return 0;
}
int get_isotherm(int *isothermal_flag){
#ifdef ISOTHERM_EQS
    *isothermal_flag = 1;
#else
    *isothermal_flag = 0;
#endif
    return 0;
}
int get_eps_is_h(int *eps_is_h_flag){
#if defined(ADAPTIVE_GRAVSOFT_FORGAS) &&  defined(UNEQUALSOFTENINGS)
    *eps_is_h_flag = 1;
#else
    *eps_is_h_flag = 0;
#endif
    return 0;
}
int get_nsmooth(int *nsmooth){
    *nsmooth = All.DesNumNgb;
    return 0;
}
int set_nsmooth(int nsmooth){
    All.DesNumNgb = nsmooth;
    return 0;
}
int get_bh_tol(double *opening_angle){
    *opening_angle = All.ErrTolTheta;
    return 0;
}
int set_bh_tol(double opening_angle){
    All.ErrTolTheta = opening_angle;
    return 0;
}
int get_gdgtol(double *gadget_cell_opening_constant){
    *gadget_cell_opening_constant = All.ErrTolForceAcc;
    return 0;
}
int set_gdgtol(double gadget_cell_opening_constant){
    All.ErrTolForceAcc = gadget_cell_opening_constant;
    return 0;
}
int get_gamma(double *gamma){
    *gamma = GAMMA;
    return 0;
}
int get_alpha(double *artificial_viscosity_alpha){
    *artificial_viscosity_alpha = All.ArtBulkViscConst;
    return 0;
}
int set_alpha(double artificial_viscosity_alpha){
    All.ArtBulkViscConst = artificial_viscosity_alpha;
    return 0;
}
int get_courant(double *courant){
    *courant = All.CourantFac*2.0;
    return 0;
}
int set_courant(double courant){
    All.CourantFac = courant/2.0;
    return 0;
}
int get_nsmtol(double *n_neighbour_tol){
    *n_neighbour_tol = All.MaxNumNgbDeviation / All.DesNumNgb;
    return 0;
}
int set_nsmtol(double n_neighbour_tol){
    All.MaxNumNgbDeviation = round(n_neighbour_tol * All.DesNumNgb);
    return 0;
}


// particle property getters/setters: (will only work after commit_particles() is called)

int get_index_of_first_particle(int *index_of_the_particle){
    return get_index_of_next_particle(0, index_of_the_particle);
}
int get_index_of_next_particle(int index_of_the_particle, int *index_of_the_next_particle){
    struct particle_data Pcurrent;
    bool found = false;
    for (int i = index_of_the_particle+1; i <= particles_counter; i++){
        if (!(find_particle(i, &Pcurrent))) {
            if (found){
                return 0;
            } else {
                *index_of_the_next_particle = Pcurrent.ID;
                found = true;
            }
        }
    }
    if (found) return 1; // This was the last particle.
    return -1; // No particle found.
}
int find_particle(int index_of_the_particle, struct particle_data *Pfound){
    for(int i = 0; i < All.TotNumPart; i++){
        if(P[i].ID == (unsigned int) index_of_the_particle){
            *Pfound = P[i];
            return 0;
        }
    }
    return -1;
}
int find_sph_particle(int index_of_the_particle, struct sph_particle_data *SphPfound){
    for(int i = 0; i < All.TotN_gas; i++){
        if(P[i].ID == (unsigned int) index_of_the_particle){
            *SphPfound = SphP[i];
            return 0;
        }
    }
    return -1;
}
int get_mass(int index, double *mass){
    struct particle_data Pcurrent;
    if(!(find_particle(index, &Pcurrent))){
        *mass = Pcurrent.Mass;
        return 0;
    }
    return -1;
}
int set_mass(int index, double mass){
    struct particle_data Pcurrent;
    if(!(find_particle(index, &Pcurrent))){
        Pcurrent.Mass = mass;
        return 0;
    }
    return -1;
}
int get_radius(int index, double *radius){
    return -1;
}
int set_radius(int index, double radius){
    return -1;
}
int get_position(int index, double *x, double *y, double *z){
    struct particle_data Pcurrent;
    if(!(find_particle(index, &Pcurrent))){
        *x = Pcurrent.Pos[0];
        *y = Pcurrent.Pos[1];
        *z = Pcurrent.Pos[2];
        return 0;
    }
    return -1;
}
int set_position(int index, double x, double y, double z){
    struct particle_data Pcurrent;
    if(!(find_particle(index, &Pcurrent))){
        Pcurrent.Pos[0] = x;
        Pcurrent.Pos[1] = y;
        Pcurrent.Pos[2] = z;
        return 0;
    }
    return -1;
}
int get_velocity(int index, double *vx, double *vy, double *vz){
    struct particle_data Pcurrent;
    if(!(find_particle(index, &Pcurrent))){
        *vx = Pcurrent.Vel[0];
        *vy = Pcurrent.Vel[1];
        *vz = Pcurrent.Vel[2];
        return 0;
    }
  return -1;
}
int set_velocity(int index, double vx, double vy, double vz){
    struct particle_data Pcurrent;
    if(!(find_particle(index, &Pcurrent))){
        Pcurrent.Vel[0] = vx;
        Pcurrent.Vel[1] = vy;
        Pcurrent.Vel[2] = vz;
        return 0;
    }
    return -1;
}
int get_state(int index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz) {
    struct particle_data Pcurrent;
    if(!(find_particle(index, &Pcurrent))){
        *mass = Pcurrent.Mass;
        *x = Pcurrent.Pos[0];
        *y = Pcurrent.Pos[1];
        *z = Pcurrent.Pos[2];
        *vx = Pcurrent.Vel[0];
        *vy = Pcurrent.Vel[1];
        *vz = Pcurrent.Vel[2];
        return 0;
    }
    return -1;
}
int set_state(int index, double mass, double x, double y, double z, double vx, double vy, double vz){
    struct particle_data Pcurrent;
    if(!(find_particle(index, &Pcurrent))){
        Pcurrent.Mass = mass;
        Pcurrent.Pos[0] = x;
        Pcurrent.Pos[1] = y;
        Pcurrent.Pos[2] = z;
        Pcurrent.Vel[0] = vx;
        Pcurrent.Vel[1] = vy;
        Pcurrent.Vel[2] = vz;
        return 0;
    }
    return -1;
}
int get_state_sph(int index, double *mass, double *x, double *y, double *z, double *vx, double *vy, double *vz, double *internal_energy) {
    struct particle_data Pcurrent;
    struct sph_particle_data SphPcurrent;
    double a3;
    if(!(find_sph_particle(index, &SphPcurrent))){
        if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
        *internal_energy = SphPcurrent.Entropy * pow(SphPcurrent.Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
        if(!(find_particle(index, &Pcurrent))){
            *mass = Pcurrent.Mass;
            *x = Pcurrent.Pos[0];
            *y = Pcurrent.Pos[1];
            *z = Pcurrent.Pos[2];
            *vx = Pcurrent.Vel[0];
            *vy = Pcurrent.Vel[1];
            *vz = Pcurrent.Vel[2];
            return 0;
        }
    }
    return -1;
}
int set_state_sph(int index, double mass, double x, double y, double z, double vx, double vy, double vz, double internal_energy){
    struct particle_data Pcurrent;
    struct sph_particle_data SphPcurrent;
    double a3;
    if(!(find_sph_particle(index, &SphPcurrent))){
        if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
        SphPcurrent.Entropy = GAMMA_MINUS1 * internal_energy / pow(SphPcurrent.Density / a3, GAMMA_MINUS1);
        if(!(find_particle(index, &Pcurrent))){
            Pcurrent.Mass = mass;
            Pcurrent.Pos[0] = x;
            Pcurrent.Pos[1] = y;
            Pcurrent.Pos[2] = z;
            Pcurrent.Vel[0] = vx;
            Pcurrent.Vel[1] = vy;
            Pcurrent.Vel[2] = vz;
            return 0;
        }
    }
    return -1;
}
int get_acceleration(int index, double * ax, double * ay, double * az){
    struct particle_data Pcurrent;
    struct sph_particle_data SphPcurrent;
    if(!(find_particle(index, &Pcurrent))){
        *ax = Pcurrent.GravAccel[0];
        *ay = Pcurrent.GravAccel[1];
        *az = Pcurrent.GravAccel[2];
        if(Pcurrent.Type == 0){
            if(!(find_sph_particle(index, &SphPcurrent))){
                *ax += SphPcurrent.HydroAccel[0];
                *ay += SphPcurrent.HydroAccel[1];
                *az += SphPcurrent.HydroAccel[2];
            }
        }
        return 0;
    }
    return -1;
}
int set_acceleration(int index, double ax, double ay, double az){
    return -1;
}
int get_internal_energy(int index, double *internal_energy){
    struct sph_particle_data Pcurrent;
    double a3;
    if(!(find_sph_particle(index, &Pcurrent))){
        if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
        *internal_energy = Pcurrent.Entropy * pow(Pcurrent.Density / a3, GAMMA_MINUS1) / GAMMA_MINUS1;
        return 0;
    }
    return -1;
}
int set_internal_energy(int index, double internal_energy){
    struct sph_particle_data Pcurrent;
    double a3;
    if(!(find_sph_particle(index, &Pcurrent))){
        if(All.ComovingIntegrationOn){a3 = All.Time * All.Time * All.Time;}else{a3 = 1;}
        Pcurrent.Entropy = GAMMA_MINUS1 * internal_energy / pow(Pcurrent.Density / a3, GAMMA_MINUS1);
        return 0;
    }
    return -1;
}
int get_smoothing_length(int index, double *smoothing_length){
    struct sph_particle_data Pcurrent;
    if (!(find_sph_particle(index, &Pcurrent))){
        if (!density_up_to_date){
            density();
            density_up_to_date = true;
        }
        *smoothing_length = Pcurrent.Hsml;
        return 0;
    }
    return -1;
}
int get_density(int index, double *density_out){
    struct sph_particle_data Pcurrent;
    if (!(find_sph_particle(index, &Pcurrent))){
        if (!density_up_to_date){
            density();
            density_up_to_date = true;
        }
        *density_out = Pcurrent.Density;
        return 0;
    }
    return -1;
}
int get_n_neighbours(int index, double *n_neighbours){
    struct sph_particle_data Pcurrent;
    if (!(find_sph_particle(index, &Pcurrent))){
        if (!density_up_to_date){
            density();
            density_up_to_date = true;
        }
        *n_neighbours = Pcurrent.NumNgb;
        return 0;
    }
    return -1;
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
    *time = All.Time;
    return 0;
}
int get_total_radius(double *radius){
    double r_squared;
    int i, j;
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    *radius = 0;
    for (i = 0; i < NumPart; i++){
        for (r_squared = 0, j = 0; j < 3; j++)
            r_squared += (SysState.CenterOfMass[j]-P[i].Pos[j])*(SysState.CenterOfMass[j]-P[i].Pos[j]);
        if (r_squared > *radius)
            *radius = r_squared;
    }
    *radius = sqrt(*radius);
    return 0;
}
int get_total_mass(double *mass){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    *mass = SysState.Mass;
    return 0;
}
int get_potential(double x, double y, double z, double *V){
    return -2;
}
int get_kinetic_energy(double *kinetic_energy){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    *kinetic_energy = SysState.EnergyKin;
    return 0;
}
int get_potential_energy(double *potential_energy){
    if (!(global_quantities_of_system_up_to_date && potential_energy_also_up_to_date))
        update_global_quantities(true);
    *potential_energy = SysState.EnergyPot;
    return 0;
}
int get_thermal_energy(double *thermal_energy){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    *thermal_energy = SysState.EnergyInt;
    return 0;
}
int get_number_of_particles(int *number_of_particles){
    *number_of_particles = All.TotNumPart;
    return 0; 
}
int get_indices_of_colliding_particles(int *index_of_particle1, int *index_of_particle2){
    return -1;
}
int get_center_of_mass_position(double *x, double *y, double *z){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    *x = SysState.CenterOfMass[0];
    *y = SysState.CenterOfMass[1];
    *z = SysState.CenterOfMass[2];
    return 0;
}
int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
    if (!global_quantities_of_system_up_to_date)
        update_global_quantities(false);
    *vx = SysState.Momentum[0]/SysState.Mass;
    *vy = SysState.Momentum[1]/SysState.Mass;
    *vz = SysState.Momentum[2]/SysState.Mass;
    return 0;
}
int get_gravity_at_point(double eps, double x, double y, double z, double *forcex, double *forcey, double *forcez){
    return -1;
}
int get_potential_at_point(double eps, double x, double y, double z, double * phi){
    return -1;
}
int get_hydro_state_at_point(double eps, double x, double y, double z, double * rho, double * rhovx, double * rhovy, double * rhovz, double * rhoe){
    double pos[3], vel[3];
    double h_out, ngb_out, dhsml_out, rho_out, rhov_out[3], rhov2_out, rhoe_out;
    int error;
    error = contruct_tree_if_needed();
    if (error)
        return error;
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    vel[0] = 0.0;
    vel[1] = 0.0;
    vel[2] = 0.0;
    hydro_state_at_point(pos, vel, &h_out, &ngb_out, &dhsml_out, &rho_out, rhov_out, &rhov2_out, &rhoe_out);
    if (debug) cout << ngb_out << endl;
    *rho   = rho_out;
    *rhovx = rhov_out[0];
    *rhovy = rhov_out[1];
    *rhovz = rhov_out[2];
#ifdef ISOTHERM_EQS
    *rhoe = rhoe_out + (rhov_out[0]*rhov_out[0] + rhov_out[1]*rhov_out[1] + rhov_out[2]*rhov_out[2]) / rho_out;
#else
    *rhoe = rhoe_out + rhov2_out;
#endif
    return 0;
}


