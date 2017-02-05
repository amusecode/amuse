/* SecularMultiple */
/* Adrian Hamers June 2016 */

#include "evolve.h"


#define FOUND_ROOT ((roots_found[i_root] == 1) || (roots_found[i_root] == -1))

void set_initial_ODE_variables(ParticlesMap *particlesMap, N_Vector &y, N_Vector &y_abs_tol, double abs_tol_spin_vec, double abs_tol_e_vec, double abs_tol_h_vec)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 0) // particle is a body
        {
            spin_vec[0] = P_p->spin_vec_x;
            spin_vec[1] = P_p->spin_vec_y;
            spin_vec[2] = P_p->spin_vec_z;

            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y,         k + k_component) = spin_vec[k_component];
                Ith(y_abs_tol, k + k_component) = abs_tol_spin_vec;
            }
            k=k+3;
        }
        if (P_p->is_binary == 1) // particle is a binary
        {
            e_vec[0] = P_p->e_vec_x;
            e_vec[1] = P_p->e_vec_y;
            e_vec[2] = P_p->e_vec_z;
            h_vec[0] = P_p->h_vec_x;
            h_vec[1] = P_p->h_vec_y;
            h_vec[2] = P_p->h_vec_z;
            double h = norm3(h_vec);
    
            for (k_component=0; k_component<3; k_component++)
            {
                Ith(y,k + k_component) = e_vec[k_component];
                Ith(y,k + k_component + 3) = h_vec[k_component];
                
                Ith(y_abs_tol,k + k_component) = abs_tol_e_vec;
                Ith(y_abs_tol,k + k_component + 3) = abs_tol_h_vec;
//                Ith(y_abs_tol,k + k_component + 3) = relative_tolerance*h;
            }
            k=k+6;
        }
    }
}

void extract_final_ODE_variables(ParticlesMap *particlesMap, N_Vector &y_out)
{
    ParticlesMapIterator it_p;
    int k=1;
    int k_component;

    double spin_vec[3],e_vec[3],h_vec[3];
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 0) // particle is a body
        {
            for (k_component=0; k_component<3; k_component++)
            {
                spin_vec[k_component] = Ith(y_out,k + k_component);
            }
            P_p->spin_vec_x = spin_vec[0];
            P_p->spin_vec_y = spin_vec[1];
            P_p->spin_vec_z = spin_vec[2];

            k=k+3;
        }
        if (P_p->is_binary == 1) // particle is a binary
        {
            for (k_component=0; k_component<3; k_component++)
            {
                e_vec[k_component] = Ith(y_out,k + k_component);
                h_vec[k_component] = Ith(y_out,k + k_component + 3);
            }

            P_p->e_vec_x = e_vec[0];
            P_p->e_vec_y = e_vec[1];
            P_p->e_vec_z = e_vec[2];
            P_p->h_vec_x = h_vec[0];
            P_p->h_vec_y = h_vec[1];
            P_p->h_vec_z = h_vec[2];
            k=k+6;
        }
    }
}


int evolve(ParticlesMap *particlesMap, double start_time, double time_step, double *output_time, double *hamiltonian, int *output_flag, int *error_code)
{
    int N_particles = particlesMap->size();
    int N_bodies, N_binaries;
    int N_root_finding;
    
    determine_binary_parents_levels_and_masses(particlesMap,&N_bodies,&N_binaries,&N_root_finding);

//    printf("N_bodies %d N_binaries %d N_particles %d N_root_finding %d\n",N_bodies,N_binaries,N_particles,N_root_finding);


    /*********************
     * setup of UserData *
     ********************/
     
	UserData data;
	data = NULL;
	data = (UserData) malloc(sizeof *data);
	data->particlesMap = particlesMap;
    data->N_root_finding = N_root_finding;


    /********************************
     * set ODE tolerances   *
     ********************************/
    if (relative_tolerance <= 0.0)
    {
        printf("relative tolerance cannot be zero; setting default value of 1e-16\n");
        relative_tolerance = 1.0e-16;
    }

    /* Note: some hardcoded parameters for ODE solver */
    double abs_tol_spin_vec = absolute_tolerance_spin_vectors;
    double abs_tol_e_vec = absolute_tolerance_eccentricity_vectors;
    double abs_tol_h_vec = absolute_tolerance_angular_momentum_vectors;
    double initial_ODE_timestep = 1.0e-6; /* one year */
    int maximum_number_of_internal_ODE_steps = 5e8;
    int maximum_number_of_convergence_failures = 100;    
    double maximum_ODE_integration_time = 13.8e10;


    /***************************
     * setup of ODE variables  *
     **************************/    
	N_Vector y, y_out, y_abs_tol;
	void *cvode_mem;
	int flag;

	y = y_out = y_abs_tol = NULL;
	cvode_mem = NULL;

    int number_of_ODE_variables = N_bodies*3 + N_binaries*6; // spin vectors for each body + e & h vectors for each binary
//    printf("N_ODE %d\n",number_of_ODE_variables);
    
    y = N_VNew_Serial(number_of_ODE_variables);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return 1;
    y_out = N_VNew_Serial(number_of_ODE_variables);
	if (check_flag((void *)y_out, "N_VNew_Serial", 0)) return 1;
    y_abs_tol = N_VNew_Serial(number_of_ODE_variables); 
	if (check_flag((void *)y_abs_tol, "N_VNew_Serial", 0)) return 1;         

    set_initial_ODE_variables(particlesMap, y, y_abs_tol,abs_tol_spin_vec,abs_tol_e_vec,abs_tol_h_vec);

    /***************************
     * setup of ODE integrator *
     **************************/    

    /* use Backward Differentiation Formulas (BDF)
        scheme in conjunction with Newton iteration --
        these choices are recommended for stiff ODEs
        in the CVODE manual                          
    */

    cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return 1;    
    
    /* essential initializations */
    flag = CVodeInit(cvode_mem, compute_y_dot, start_time, y);
	if (check_flag(&flag, "CVodeInit", 1)) return 1;    

	flag = CVodeSetUserData(cvode_mem, data);
	if (check_flag(&flag, "CVodeSetUsetData", 1)) return 1;

    flag = CVodeSVtolerances(cvode_mem, relative_tolerance, y_abs_tol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return 1;

	flag = CVDense(cvode_mem, number_of_ODE_variables);
	if (check_flag(&flag, "CVDense", 1)) return 1;

	flag = CVodeSetInitStep(cvode_mem, initial_ODE_timestep);
	if (check_flag(&flag, "CVodeSetInitStep", 1)) return 1;

    /* optional initializations */
//	flag = CVodeSetErrHandlerFn(cvode_mem, ehfun, eh_data); // error handling function
//	if (check_flag(&flag, "CVodeSetErrHandlerFn", 1)) return;
	  		
	flag = CVodeSetMaxNumSteps(cvode_mem, maximum_number_of_internal_ODE_steps);
	if (check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return 1;

//	flag = CVodeSetMinStep(cvode_mem, 0.1); // minimum step size
//	if (check_flag(&flag, "CVodeSetMinStep", 1)) return 1;

	flag = CVodeSetMaxHnilWarns(cvode_mem, 1);
	if (check_flag(&flag, "CVodeSetMaxHnilWarns", 1)) return 1;
			
//	flag = CVodeSetStopTime(cvode_mem, MAXTIME); // maximum time
//	if (check_flag(&flag, "CVodeSetStopTime", 1)) return 1;

	flag = CVodeSetMaxConvFails(cvode_mem, maximum_number_of_convergence_failures);
	if (check_flag(&flag, "CVodeSetMaxConvFails", 1)) return 1;

    /* initialization of root finding */
    int roots_found[N_root_finding];
	flag = CVodeRootInit(cvode_mem, N_root_finding, root_finding_functions);
	if (check_flag(&flag, "CVodeRootInit", 1)) return 1;	


    /***************************
     * ODE integration         *
     **************************/ 
    
	double user_end_time = start_time + time_step;
	double integrator_end_time;

	flag = CVode(cvode_mem, user_end_time, y_out, &integrator_end_time, CV_NORMAL);	

    if (check_for_initial_roots(particlesMap) > 0)
    {
        flag = CV_ROOT_RETURN;
    }

	if (flag == CV_SUCCESS)
	{
		*output_flag = CV_SUCCESS;
		*error_code = 0;
        *output_time = integrator_end_time;
	}
	else if (flag == CV_ROOT_RETURN) // a root was found during the integration
	{
		CVodeGetRootInfo(cvode_mem,roots_found);
        read_root_finding_data(particlesMap,roots_found);
        *output_flag = CV_ROOT_RETURN;
        *output_time = integrator_end_time;
    }
    else if (flag == CV_WARNING) // a warning has occurred during the integration
    {
		*output_flag = 99;
		*error_code = flag;
    }
	else // an error has occurred during the integration
    {
		*output_flag = flag;
		*error_code = flag;
    }

    /***************************
     * y_out -> particlesMap   *
     * ************************/

    extract_final_ODE_variables(particlesMap,y_out);

    *hamiltonian = data->hamiltonian;
    
    N_VDestroy_Serial(y);
    N_VDestroy_Serial(y_out);
    N_VDestroy_Serial(y_abs_tol);
    CVodeFree(&cvode_mem);

	return 0;
    
    
}

int read_root_finding_data(ParticlesMap *particlesMap, int *roots_found)
{
    ParticlesMapIterator it_p;
    
    int i_root = 0;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 1)
        {
            if (P_p->check_for_secular_breakdown == 1)
            {
                if FOUND_ROOT
                {
                    P_p->secular_breakdown_has_occurred = 1;
                }
                i_root++;

            }
            if (P_p->check_for_dynamical_instability == 1)
            {
                if FOUND_ROOT
                {
                    P_p->dynamical_instability_has_occurred = 1;
                }
                i_root++;                
            }
            if (P_p->check_for_physical_collision_or_orbit_crossing == 1)
            {
                if FOUND_ROOT
                {
                    P_p->physical_collision_or_orbit_crossing_has_occurred = 1;
                }
                i_root++;                
            }
            if (P_p->check_for_minimum_periapse_distance == 1)
            {
                if FOUND_ROOT
                {
                    P_p->minimum_periapse_distance_has_occurred = 1;
                }
                i_root++;                
            }
        }
        else /* P_p not a binary */
        {
            if (P_p->check_for_RLOF_at_pericentre == 1)
            {
                if FOUND_ROOT
                {
                    P_p->RLOF_at_pericentre_has_occurred = 1;
                }
                i_root++;
            }
        }
    }
    return 0;
}

int check_for_initial_roots(ParticlesMap *particlesMap)
{
    ParticlesMapIterator it_p;
    
    int N_root_found = 0;
    
    for (it_p = particlesMap->begin(); it_p != particlesMap->end(); it_p++)
    {
        Particle *P_p = (*it_p).second;
        if (P_p->is_binary == 1)
        {
            if (P_p->check_for_secular_breakdown == 1)
            {
                if (P_p->secular_breakdown_has_occurred == 1)
                {
                    N_root_found++;
                }

            }
            if (P_p->check_for_dynamical_instability == 1)
            {
                if (P_p->dynamical_instability_has_occurred == 1)
                {
                    N_root_found++;
                }
            }
            if (P_p->check_for_physical_collision_or_orbit_crossing == 1)
            {
                if (P_p->physical_collision_or_orbit_crossing_has_occurred == 1)
                {
                    N_root_found++;
                }
            }
            if (P_p->check_for_minimum_periapse_distance == 1)
            {
                if (P_p->minimum_periapse_distance_has_occurred == 1)
                {
                    N_root_found++;
                }
            }            
        }
        else /* P_p not a binary */
        {
            if (P_p->check_for_RLOF_at_pericentre == 1)
            {
                if (P_p->RLOF_at_pericentre_has_occurred == 1)
                {
                    N_root_found++;
                }
            }
        }
    }
    return N_root_found;
}

/* function to check ODE solver-related function return values */
static int check_flag(void *flagvalue, char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return 0;
}
