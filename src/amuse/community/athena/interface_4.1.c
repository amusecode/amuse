#define MAIN_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <math.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>

#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include <math.h>

//AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>
#include <time.h>

static MeshS mesh;
static VDFun_t Integrate;
static int is_dt_set_by_script=0;
static int is_restart = 0;
static int evolve_to_exact_time = 1;
static int has_external_gravitational_potential=0;
static int mpi_comm_size_interface = 0;
static Real last_dt_above_zero = 0.0;


#ifdef SELF_GRAVITY
  VDFun_t SelfGrav;      /* function pointer to self-gravity, set at runtime */
#endif


static inline void ijk_pos(
    const GridS *pG,
    const Real x1, const Real x2, const Real x3,
    Real * i, Real * j, Real * k
)
{
  if(x1 < pG->MinX[0] || x1 > pG->MaxX[0])
  {
      *i = *j = *k = 0.0;
      return;
  }
  if(x2 < pG->MinX[1] || x2 > pG->MaxX[1])
  {
      *i = *j = *k = 0.0;
      return;
  }
  if(x3 < pG->MinX[2] || x3 > pG->MaxX[2])
  {
      *i = *j = *k = 0.0;
      return;
  }

  if(pG->dx1 == 0) {
    *i = 0;
  } else {
    *i = ((x1 - pG->MinX[0])/pG->dx1) - 0.5 + pG->Disp[0];
    //fprintf(stderr, "index of pos (i) %d\n", i);
    modf(round(*i), i);
  }
  if(pG->dx2 == 0) {
     *j = 0;
  } else {
    *j = ((x2 - pG->MinX[1])/pG->dx2) - 0.5 + pG->Disp[1];
    modf(round(*j), j);

  }

  if(pG->dx3 == 0) {
    *k = 0.0;
  } else {
    *k = ((x3 - pG->MinX[2])/pG->dx3) - 0.5 + pG->Disp[2];
    modf(round(*k), k);
  }
}

static DomainS * get_domain_structure_with_index(int index_of_grid)
{
    DomainS * dom = 0;
    int ii = 0, level = 0;
    for(level = 0; level < mesh.NLevels; level++)
    {
        for(ii = 0; ii < mesh.DomainsPerLevel[level]; ii++)
        {
            dom = (DomainS*)&(mesh.Domain[level][ii]);
            //printf("get_position_of_index: %d, %d, %d\n", level, dom->InputBlock, index_of_grid);
            if(dom->InputBlock == index_of_grid)
            {
                return dom;
            }
        }
    }
    return dom;
}

static void boundary_left_x1_copy(GridS *pGrid)
{
    int is = pGrid->is;
    int js = pGrid->js, je = pGrid->je;
    int ks = pGrid->ks, ke = pGrid->ke;
    int i,j,k;
    #ifdef MHD
    int ju, ku; /* j-upper, k-upper */
    #endif

    if(pGrid->boundary->LeftX1 == 0) {fprintf(stderr, "LEFT X1 boundary is requested, but no data set!\n");return;}
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=1; i<=nghost; i++) {
                //fprintf(stderr, "i,j,k: %d, %d, %d (%d, %d, %d)\n", k, j, is -i, k-ks, j -js, nghost - i);
                pGrid->U[k][j][is-i] = pGrid->boundary->LeftX1[k-ks][j-js][nghost - i];
            }
        }
    }
  
}

static void boundary_right_x1_copy(GridS *pGrid)
{
    int ie = pGrid->ie;
    int js = pGrid->js, je = pGrid->je;
    int ks = pGrid->ks, ke = pGrid->ke;
    int i,j,k;
    #ifdef MHD
    int ju, ku; /* j-upper, k-upper */
    #endif

    if(pGrid->boundary->RightX1 == 0) {fprintf(stderr, "RIGHT X1 boundary is requested, but no data set!\n");return;}
    for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
            for (i=1; i<=nghost; i++) {
                //fprintf(stderr, "i,j,k: %d, %d, %d (%d, %d, %d)\n", k, j, is -i, k-ks, j -js, nghost - i);
                pGrid->U[k][j][i+ie] = pGrid->boundary->RightX1[k-ks][j-js][i-1];
            }
        }
    }
}

static void boundary_left_x2_copy(GridS *pGrid)
{
    int isl = pGrid->is - nghost, iel = pGrid->ie + nghost;
    int js = pGrid->js;
    int ks = pGrid->ks, ke = pGrid->ke;
    int i,j,k;

    if(pGrid->boundary->LeftX2 == 0) {fprintf(stderr, "LEFT X2 boundary is requested, but no data set!\n");return;}
    
    for (k=ks; k<=ke; k++) {
        for (j=1; j<=nghost; j++) {
            for (i=isl; i<=iel; i++) {
                //fprintf(stderr, "i,j,k: %d, %d, %d (%d, %d, %d) %d, %d\n", k, js-j, i, k-ks, nghost - j, i - isl, isl, isl);
                pGrid->U[k][js-j][i] = pGrid->boundary->LeftX2[k-ks][nghost-j][i-isl];
            }
        }
    }
}

static void boundary_right_x2_copy(GridS *pGrid)
{
    int isl = pGrid->is - nghost, iel = pGrid->ie + nghost;
    int je = pGrid->je;
    int ks = pGrid->ks, ke = pGrid->ke;
    int i,j,k;

    if(pGrid->boundary->RightX2 == 0) {fprintf(stderr, "RIGHT X2 boundary is requested, but no data set!\n");return;}
   
    for (k=ks; k<=ke; k++) {
        for (j=1; j<=nghost; j++) {
            for (i=isl; i<=iel; i++) {
                pGrid->U[k][je+j][i] = pGrid->boundary->RightX2[k-ks][j-1][i-isl];
            }
        }
    }
}

static void boundary_left_x3_copy(GridS *pGrid)
{
    int isl = pGrid->is - nghost, iel = pGrid->ie + nghost;
    int jsl = pGrid->js - nghost, jel = pGrid->je + nghost;
    int ks = pGrid->ks, ke = pGrid->ke;
    int i,j,k;

    if(pGrid->boundary->LeftX3 == 0) {fprintf(stderr, "LEFT X3 boundary is requested, but no data set!\n");return;}
    

    for (k=1; k<=nghost; k++) {
        for (j=jsl; j<=jel; j++) {
            for (i=isl; i<=iel; i++) {
                pGrid->U[ks-k][j][i] = pGrid->boundary->LeftX3[nghost - k][j-jsl][i-isl];
            }
        }
    }
}

static void boundary_right_x3_copy(GridS *pGrid)
{
    int isl = pGrid->is - nghost, iel = pGrid->ie + nghost;
    int jsl = pGrid->js - nghost, jel = pGrid->je + nghost;
    int ks = pGrid->ks, ke = pGrid->ke;
    int i,j,k;

    if(pGrid->boundary->RightX3 == 0) {fprintf(stderr, "RIGHT X3 boundary is requested, but no data set!\n");return;}
    
    for (k=1; k<=nghost; k++) {
        for (j=jsl; j<=jel; j++) {
            for (i=isl; i<=iel; i++) {
                pGrid->U[ke+k][j][i] = pGrid->boundary->RightX3[k - 1][j-jsl][i-isl];
            }
        }
    }
}


/*
 * determine i,j,k for x1, x2, x3
 * x1 will be between i and i+1
 * x2 will be between j and j+1
 * x3 will be between k and k+1
 */


static inline void ijk_pos_dom(
    const Real x1, const Real x2, const Real x3,
    Real * i, Real * j, Real * k,
    Real * dx, Real * dy, Real * dz
)
{
  DomainS * domain = &mesh.Domain[0][0];

  if(domain->dx[0] == 0) {
    *i = 0;
    *dx = 0;
  } else {
    *i = ((x1 - domain->MinX[0])/domain->dx[0]) - 0.5 + nghost;
    *dx = modf(*i, i);
  }
  if(domain->Nx[1] == 1) {
     *j = 0;
     *dy = 0.0;
  } else {
    *j = ((x2 - domain->MinX[1])/domain->dx[1]) - 0.5 + nghost;
    *dy = modf(*j, j);
  }
  
  if(domain->Nx[2] == 1) {
    *k = 0.0;
    *dz = 0;
  } else {
    *k = ((x3 - domain->MinX[2])/domain->dx[2]) - 0.5 + nghost;
    *dz = modf(*k, k);
  }
}

static Real ***Potentials=NULL;

static Real grav_pot(const Real x1, const Real x2, const Real x3)
{
    Real i, j, k;
    Real dx, dy, dz;
    ijk_pos_dom(x1,x2,x3, &i, &j, &k, &dx, &dy, &dz);
    int ii, jj, kk;
    ii = i;
    jj = j;
    kk = k;
    
    Real potential000 = Potentials[kk][jj][ii];
    
    Real potential001 = potential000;
    Real potential100 = potential000;
    Real potential101 = potential000;
    Real potential010 = potential000;
    Real potential011 = potential000;
    Real potential110 = potential000;
    Real potential111 = potential000;
    if(dx > 0) potential001 = Potentials[kk][jj][ii+1];
    if(dz > 0) potential100 = Potentials[kk+1][jj][ii];
    if(dz > 0 && dx > 0) potential101 = Potentials[kk+1][jj][ii+1];
    if(dy > 0) {
        potential010 = Potentials[kk][jj+1][ii];
        if(dx > 0) potential011 = Potentials[kk][jj+1][ii+1];
        if(dz > 0) potential110 = Potentials[kk+1][jj+1][ii];
        if(dz > 0 && dx > 0) potential111 = Potentials[kk+1][jj+1][ii+1];
    }
    
    Real potential = 
        (potential000 * (1 - dz) * (1 - dy) * (1 - dx)) +
        (potential100 * dz * (1 - dy) * (1 - dx)) +
        (potential010 * (1 - dz) * dy * (1 - dx)) +
        (potential001 * (1 - dz) * (1 - dy) * dx) +
        (potential101 * dz * (1 - dy) * dx) +
        (potential011 * (1 - dz) * dy * dx) +
        (potential110 * dz * dy * (1 - dx)) +
        (potential111 * dz * dy * dx );
    
    
    //printf("dx,dy,dx %f,%f,%f\n",dx,dy,dz);
    
    return potential;
}


int cleanup_code(){
  return 0;
}

int recommit_parameters(){
  return 0;
}

int set_has_external_gravitational_potential(int value) {
    has_external_gravitational_potential = value;
    return 0;
}

int get_has_external_gravitational_potential(int * value) {
    *value = has_external_gravitational_potential;
    return 0;
}



int get_number_of_grids(int * value) {
    DomainS * dom = 0;
    int ii = 0, level = 0;
    int result = -1;
    for(level = 0; level < mesh.NLevels; level++)
    {
        for(ii = 0; ii < mesh.DomainsPerLevel[level]; ii++)
        {
            dom = (DomainS*)&(mesh.Domain[level][ii]);
            if(dom->InputBlock > result)
            {
                result = dom->InputBlock;
            }
        }
    }
    *value = result;
    return 0;
}


int initialize_code(){

#ifdef SELF_GRAVITY
  grav_mean_rho = 1.0;
  four_pi_G = -1.0;
#endif
  last_dt_above_zero = 0.0;
#ifdef _WIN32
  par_open("nul"); 
#else
  par_open("/dev/null"); /* to trick athena into thinking it has opened a parameter file */
#endif

  par_sets("job","problem_id", "amuse", "all amuse runs");
  is_restart = 0;
  show_config_par();   /* Add the configure block to the parameter database */


  // AMUSE STOPPING CONDITIONS SUPPORT
  set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
  set_support_for_condition(TIMEOUT_DETECTION);

#ifdef MPI_PARALLEL
    /* Get my task id (rank in MPI) */
    if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD,&myID_Comm_world))
    {
        ath_error("Error on calling MPI_Comm_rank\n");
    }
    if(MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD,&(mpi_comm_size_interface)))
    {
        ath_error("Error on calling mpi_comm_size_interface\n");
    }
#else
    myID_Comm_world = 0;
    mpi_comm_size_interface = 0;
#endif
  return 0;
}

int get_timestep(double * value) {
    *value = mesh.dt;
    return 0;
}

int set_timestep(double value) {
    mesh.dt = value;
    is_dt_set_by_script = 1;
    return 0;
}

int get_nghost(int * value) {
    *value = nghost;
    return 0;
}

int get_evolve_to_exact_time(int * value) {
    *value = evolve_to_exact_time;
    return 0;
}

int set_evolve_to_exact_time(int value) {
    evolve_to_exact_time = value;
    return 0;
}



int commit_parameters(){
  int nl, nd;
  
  int out_level = par_geti_def("log","out_level",0);
  int err_level = par_geti_def("log","err_level",0);

  ath_log_set_level(out_level, err_level);
/*
  if(has_external_gravitational_potential) {
    StaticGravPot = grav_pot;
  }
*/
  CourNo = par_getd("time","cour_no");

#ifdef SELF_GRAVITY
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = par_getd("problem","grav_mean_rho");
#endif

#ifdef ISOTHERMAL
  Iso_csound = par_getd("problem","iso_csound");
  Iso_csound2 = Iso_csound*Iso_csound;
#else
  Gamma = par_getd("problem","gamma");
  Gamma_1 = Gamma - 1.0;
  Gamma_2 = Gamma - 2.0;
#endif

  if(has_external_gravitational_potential) {
    StaticGravPot = grav_pot;
  }

  init_mesh(&mesh);
  init_grid(&mesh);
  
  
    
  for (nl=0; nl<(mesh.NLevels); nl++){
    for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){
        if (mesh.Domain[nl][nd].Grid != NULL){
            GridS * currentGrid = mesh.Domain[nl][nd].Grid;
            currentGrid->boundary = (BoundaryCellS *) calloc(sizeof(BoundaryCellS), 1);
            
            if(par_geti( "domain1", "bc_ix1") == 10) {
                if (mesh.Domain[nl][nd].Disp[0] == 0) {  
                    if(currentGrid->Disp[0] == 0) {              
                        currentGrid->boundary->LeftX1 = (ConsS***) calloc_3d_array(
                            currentGrid->Nx[2],
                            currentGrid->Nx[1],
                            nghost,
                            sizeof(ConsS)
                        );
                    }
                    bvals_mhd_fun(&mesh.Domain[nl][nd], left_x1,  boundary_left_x1_copy);
                }     
            }     
            if( par_geti( "domain1", "bc_ox1") == 10) {
                if (mesh.Domain[nl][nd].MaxX[0] == mesh.Domain[nl][nd].RootMaxX[0])  {
                    if(currentGrid->Disp[0] + currentGrid->Nx[0] == mesh.Domain[nl][nd].Nx[0]) {             
                        currentGrid->boundary->RightX1 = (ConsS***) calloc_3d_array(
                            currentGrid->Nx[2],
                            currentGrid->Nx[1],
                            nghost,
                            sizeof(ConsS)
                        );
                    }
                    bvals_mhd_fun(&mesh.Domain[nl][nd], right_x1,  boundary_right_x1_copy);
                }
            }
            if(par_geti( "domain1", "bc_ix2") == 10) {
                if (mesh.Domain[nl][nd].Disp[1] == 0) {   
                    if(currentGrid->Disp[1] == 0) { 
                        currentGrid->boundary->LeftX2 = (ConsS***) calloc_3d_array(
                            currentGrid->Nx[2],
                            nghost,
                            currentGrid->Nx[0] + 2 * nghost,
                            sizeof(ConsS)
                        );
                    }
                    bvals_mhd_fun(&mesh.Domain[nl][nd], left_x2,  boundary_left_x2_copy);
                }     
            }       
            if( par_geti( "domain1", "bc_ox2") == 10) {
                if (mesh.Domain[nl][nd].MaxX[1] == mesh.Domain[nl][nd].RootMaxX[1])  {
                    
                    if(currentGrid->Disp[1] + currentGrid->Nx[1] == mesh.Domain[nl][nd].Nx[1]) {    
                        currentGrid->boundary->RightX2 = (ConsS***) calloc_3d_array(
                            currentGrid->Nx[2],
                            nghost,
                            currentGrid->Nx[0] + 2 * nghost,
                            sizeof(ConsS)
                        );
                    }
                    bvals_mhd_fun(&mesh.Domain[nl][nd], right_x2,  boundary_right_x2_copy);
                }
            }  
            int jfactor = currentGrid->Nx[1] > 1 ? 2 : 0;
            if(par_geti( "domain1", "bc_ix3") == 10) {
                if (mesh.Domain[nl][nd].Disp[2] == 0) {
                    if(currentGrid->Disp[2] == 0) {            
                        currentGrid->boundary->LeftX3 = (ConsS***) calloc_3d_array(
                            nghost,
                            currentGrid->Nx[1] + jfactor * nghost,
                            currentGrid->Nx[0] + 2 * nghost,
                            sizeof(ConsS)
                        );
                    }
                    bvals_mhd_fun(&mesh.Domain[nl][nd], left_x3,  boundary_left_x3_copy);
                }     
            }   
            if( par_geti( "domain1", "bc_ox3") == 10) {
                if (mesh.Domain[nl][nd].MaxX[2] == mesh.Domain[nl][nd].RootMaxX[2])  {
                    
                    if(currentGrid->Disp[2] + currentGrid->Nx[2] == mesh.Domain[nl][nd].Nx[2]) {   
                        currentGrid->boundary->RightX3 = (ConsS***) calloc_3d_array(
                            nghost,
                            currentGrid->Nx[1] + jfactor * nghost,
                            currentGrid->Nx[0] + 2 * nghost,
                            sizeof(ConsS)
                        );
                    }
                
                    bvals_mhd_fun(&mesh.Domain[nl][nd], right_x3,  boundary_right_x3_copy);
                }
            }  
        }
    }
  }
  
 
  
  if ((Potentials = (Real***)calloc_3d_array(
    mesh.Domain[0][0].Nx[2] + 2 + (2 * nghost),
    mesh.Domain[0][0].Nx[1] + 2 + (2 * nghost),
    mesh.Domain[0][0].Nx[0] + 2 + (2 * nghost),
    sizeof(Real))) == NULL)
  {
    return -1;
  }


  return 0;
}

int initialize_grid()
{

  int nl, nd;
/* restrict initial solution so grid hierarchy is consistent */
#ifdef STATIC_MESH_REFINEMENT
  SMR_init(&mesh);
  RestrictCorrect(&mesh);
#endif

  bvals_mhd_init(&mesh);

#ifdef SELF_GRAVITY
  bvals_grav_init(&mesh);
#endif
#if defined(SHEARING_BOX) || (defined(FARGO) && defined(CYLINDRICAL))
  bvals_shear_init(&mesh);
#endif
#ifdef PARTICLES
  bvals_particle_init(&mesh);
#endif


  for (nl=0; nl<(mesh.NLevels); nl++){
    for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){
      if (mesh.Domain[nl][nd].Grid != NULL){
        bvals_mhd(&(mesh.Domain[nl][nd]));
#ifdef PARTICLES
        bvals_particle(&(Mesh.Domain[nl][nd]));
#ifdef FEEDBACK
        exchange_feedback_init(&(Mesh.Domain[nl][nd]));
#endif
#endif
      }
    }
  }

#ifdef STATIC_MESH_REFINEMENT
  Prolongate(&mesh);
#endif

  if(!is_dt_set_by_script) {
    par_setd("time","tlim", "%f", 100.0, "-");
    new_dt(&mesh);
  }

  lr_states_init(&mesh);

  Integrate = integrate_init(&mesh);

#ifdef SELF_GRAVITY
  SelfGrav = selfg_init(&mesh);
  for (nl=0; nl<(mesh.NLevels); nl++){
    for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){
      if(mesh.Domain[nl][nd].Grid != NULL){
        (*SelfGrav)(&(mesh.Domain[nl][nd]));
        bvals_grav(&(mesh.Domain[nl][nd]));
      }
    }
  }
#endif
#if defined(RESISTIVITY) || defined(VISCOSITY) || defined(THERMAL_CONDUCTION)
  integrate_diff_init(&mesh);
#endif
  return 0;
}

static inline int is_on_grid(GridS * grid, int i0, int j0, int k0)
{
    if (grid->Nx[0] > 1 && (i0 < (grid->Disp[0])  || i0 >= (grid->Disp[0] + grid->Nx[0])))
    {
        return 0;
    }
    else if (grid->Nx[1] > 1 && (j0 < (grid->Disp[1])  || j0 >= (grid->Disp[1] + grid->Nx[1])))
    {
        return 0;
    }
    else if (grid->Nx[2] > 1 && (k0 < (grid->Disp[2])  || k0 >= (grid->Disp[2] + grid->Nx[2])))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}

static inline int is_on_extended_grid(DomainS * dom, GridS * grid, int i0, int j0, int k0)
{
    if (grid->Nx[0] > 1 && ((grid->Disp[0]>0 && i0 < grid->Disp[0]) || (grid->Disp[0] + grid->Nx[0]<dom->Nx[0] && i0 >= (grid->Disp[0] + grid->Nx[0]))) )
    {
        return 0;
    }
    else if (grid->Nx[1] > 1 && ((grid->Disp[1]>0 && j0 < grid->Disp[1])  || (grid->Disp[1] + grid->Nx[1]<dom->Nx[1] && j0 >= (grid->Disp[1] + grid->Nx[1]))))
    {
        return 0;
    }
    else if (grid->Nx[2] > 1 && ((grid->Disp[2]>0 && k0 < grid->Disp[2])  || (grid->Disp[2] + grid->Nx[2]<dom->Nx[2] && k0 >= (grid->Disp[2] + grid->Nx[2]))))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}


static inline int is_on_magnetic_grid(GridS * grid, int i0, int j0, int k0)
{
    if (grid->Nx[0] > 1 && (i0 < (grid->Disp[0])  || i0 >= (grid->Disp[0] + grid->Nx[0] + 1)))
    {
        return 0;
    }
    else if (grid->Nx[1] > 1 && (j0 < (grid->Disp[1])  || j0 >= (grid->Disp[1] + grid->Nx[1] + 1)))
    {
        return 0;
    }
    else if (grid->Nx[2] > 1 && (k0 < (grid->Disp[2])  || k0 >= (grid->Disp[2] + grid->Nx[2] + 1)))
    {
        return 0;
    }
    else
    {
        return 1;
    }
}
static inline int ijk_on_grid(GridS * grid, int * i0, int * j0, int * k0)
{
    *i0 -= grid->Disp[0] - grid->is;
    if(grid->Nx[1] > 1)
    {
        *j0 -= grid->Disp[1] - grid->js;
    }
    else
    {
        *j0 = 0;
    }
    if(grid->Nx[2] > 1)
    {
        *k0 -= grid->Disp[2] - grid->ks;
    }
    else
    {
        *k0 = 0;
    }
}


static inline int is_on_boundary_grid(DomainS * dom, GridS * grid, int boundary, int i0, int j0, int k0, int isset)
{
    //printf("DISP: %d, %d, %d \n", grid->Disp[0] , grid->Disp[1] , grid->Disp[2] );
    /*
    printf("I: %d, %d, %d, %d \n", grid->is, grid->ie,  grid->ie -  grid->is + 1, grid->Nx[0] + grid->Disp[0]);
    printf("J: %d, %d, %d, %d \n", grid->js, grid->je,  grid->je -  grid->js + 1, grid->Nx[1] + grid->Disp[1]);
    printf("K: %d, %d, %d, %d \n", grid->ks, grid->ke,  grid->ke -  grid->ks + 1, grid->Nx[2] + grid->Disp[2]);
    */
    int il = isset || grid->Disp[0] == 0 ? grid->Disp[0] : grid->Disp[0] + nghost;
    int ir = grid->Disp[0] + grid->Nx[0];
    ir = isset || ir  == dom->Nx[0] ? ir + (nghost*2) : ir + nghost;

    int jl = isset || grid->Disp[1] == 0 ? grid->Disp[1] : grid->Disp[1] + nghost;
    int jr = grid->Disp[1] + grid->Nx[1];
    jr = isset || jr  == dom->Nx[1] ? jr + (nghost*2) : jr + nghost;
    
    switch(boundary)
    {
        case 1:
        case 2:
            if (i0 < 0 || i0 >= nghost)
            {
                return 0;
            }
            else if (grid->Nx[1] > 1 &&  (j0 < (grid->Disp[1])  || j0 >= (grid->Disp[1] + grid->Nx[1])))
            {
                return 0;
            }
            else if (grid->Nx[2] > 1 && (k0 < (grid->Disp[2])  || k0 >= (grid->Disp[2] + grid->Nx[2])))
            {
                return 0;
            }
            else
            {
                return 1;
            }
            return 0;
        case 3:
        case 4:
            if (i0 < il || i0 >= ir)
            {
                return 0;
            }
            else if (j0 < 0 || j0 >= nghost)
            {
                return 0;
            }
            else if (grid->Nx[2] > 1 && (k0 < (grid->Disp[2])  || k0 >= (grid->Disp[2] + grid->Nx[2])))
            {
                return 0;
            }
            else
            {
                return 1;
            }
        case 5:
        case 6:
            if (i0 < il || i0 >= ir)
            {
                return 0;
            }
            else if (j0 < jl || j0 >= jr)
            {
                return 0;
            }
            else if (k0 < 0 || k0 >= nghost)
            {
                return 0;
            }
            else
            {
                return 1;
            }
        default:
            return 0;
    }
    
}


static inline int ijk_on_boundary_grid(GridS * grid, int boundary, int * i0, int * j0, int * k0)
{
    int il = grid->Disp[0] == 0 ? 0 : grid->Disp[0];
    int jl = grid->Disp[1] == 0 ? 0 : grid->Disp[1];
    int kl = grid->Disp[2] == 0 ? 0 : grid->Disp[2];
    
    switch(boundary)
    {
        case 1:
        case 2:
            if(grid->Nx[1] > 1)
            {
                *j0 -= grid->Disp[1];
            }
            else
            {
                *j0 = 0;
            }
            if(grid->Nx[2] > 1)
            {
                *k0 -= grid->Disp[2];
            }
            else
            {
                *k0 = 0;
            }
            return -1;
        case 4:
        case 3:
            *i0 -= il;
            
            if(grid->Nx[2] > 1)
            {
                *k0 -= grid->Disp[2];
            }
            else
            {
                *k0 = 0;
            }
            return 0;
        case 5:
        case 6:
            *i0 -= il;
            *j0 -= jl;
        default:
            return -1;
    }
}



int get_position_of_index(int *i, int *j, int *k, int *index_of_grid, double * x, double * y,
  double * z, int number_of_points){


    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    
    DomainS * dom = 0;
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        x[l] = 0.0;
        y[l] = 0.0;
        z[l] = 0.0;
        
        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_extended_grid(dom,grid, i0, j0, k0))
            {
                ijk_on_grid(grid, &i0, &j0, &k0);
                cc_pos(
                  grid,
                  i0,
                  j0,
                  k0,
                  &x[l] ,
                  &y[l] ,
                  &z[l]
                );
            }

        }
    }

#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(x, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(y, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(z, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, x, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, y, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, z, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}

int get_interpolated_gravitational_potential(double x, double y, double z,   double * potential) {
    if (Potentials == NULL) {
        return -1;
    }



    *potential = grav_pot(x,y,z);

    return 0;
}

int get_index_of_position(double x, double y,
  double z,int index_of_grid, double *i , double * j, double * k)
  {

    double pos[3] = {0,0,0};

    if (mesh.NLevels == 0) {
        return -1;
    }

    
    DomainS * dom = get_domain_structure_with_index(index_of_grid);
   
    if(dom == 0){
        return -3;
    }
    
    if(dom->Grid == NULL)
    {
        pos[0] = 0.0;
        pos[1] = 0.0;
        pos[2] = 0.0;
    }
    else
    {
        GridS * grid = dom->Grid;
        ijk_pos(grid, x, y, z, &pos[0], &pos[1], &pos[2]);
    }



#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(pos, NULL, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, pos, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    *i = pos[0];
    *j = pos[1];
    *k = pos[2];

    return 0;
}

int test() {
    return 1;
}

int get_grid_state(
    int * i, int * j, int * k,
    int * index_of_grid,
    double * rho,
    double * rhovx, double * rhovy, double * rhovz,
    double * en,
    int number_of_points)
{
    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    int ii = 0;
    DomainS * dom = 0;
    if (mesh.NLevels == 0) {
        return -1;
    }
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        
        rho[l] = 0;
        rhovx[l] = 0;
        rhovy[l] = 0;
        rhovz[l] = 0;
        en[l] = 0;
                
        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            if (is_on_grid(grid, i0, j0, k0))
            {
                ijk_on_grid(grid, &i0, &j0, &k0);


                rho[l] = grid->U[k0][j0][i0].d;
                rhovx[l] = grid->U[k0][j0][i0].M1;
                rhovy[l] = grid->U[k0][j0][i0].M2;
                rhovz[l] = grid->U[k0][j0][i0].M3;
                en[l] = grid->U[k0][j0][i0].E;
            }

        }
    }


#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(rho, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovx, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovy, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovz, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(en, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, rho, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovx, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovy, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovz, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, en, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}


int get_grid_gravitational_potential(
    int * i, int * j, int * k,
    int * index_of_grid,
    double * phi,
    int number_of_points)
{
    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    int ii = 0;
    DomainS * dom = 0;
    if (mesh.NLevels == 0) {
        return -1;
    }
#ifndef SELF_GRAVITY
    for(l=0; l < number_of_points; l++) {
        phi[l] = 0.0;
    }
    
#else

    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        current_index_of_grid = index_of_grid[l];
        
        phi[l] = 0;
        
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {
                
    //fprintf(stderr, "GRID i,j,k %d,%d,%d\n", i0, j0, k0);
                ijk_on_grid(grid, &i0, &j0, &k0);
                phi[l] = grid->Phi[k0][j0][i0];
            }

        }
    }



#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(phi, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     } else {
        MPI_Reduce(MPI_IN_PLACE, phi, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

#endif
    return 0;
}


int get_potential_at_point(double eps, double x, double y, double z, double *phi)
{
    if(myID_Comm_world)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    
    double ii, jj, kk;
    double dx, dy, dz;
    int i, j, k;
    int index_of_grid = 1; // supports only one grid!
    DomainS * dom = get_domain_structure_with_index(1);

    if(dom == 0)
    {
        return -1;
    }
    if(dom->Grid == NULL)
    {
        return -1;
    }

    ijk_pos(dom->Grid, x, y, z, &ii, &jj, &kk);//(x, y, z, &ii, &jj, &kk, &dx, &dy, &dz);
    i = ii; j = jj; k = kk;
    //fprintf(stderr, "XYZ i,j,k %f,%f,%f\n", ii, jj, kk);
    get_grid_gravitational_potential(&i, &j, &k, &index_of_grid,  phi, 1);
    

    return 0;
}

int get_gravity_at_point(double eps, double x,double y, double z,
                         double *fx, double *fy, double *fz)
{
    if(myID_Comm_world)     { // calculate only on the root mpi process, not on others
        return 0;
    }
    double ii, jj, kk;
    double dx, dy, dz;
    int i, j, k;
    int index_of_grid = 1; // supports only one grid!
    
    
    ijk_pos_dom(x, y, z, &ii, &jj, &kk, &dx, &dy, &dz);
    i = ii; j = jj; k = kk;
    get_grid_gravitational_acceleration(&i, &j, &k, &index_of_grid, fx, fy, fz, 1);

    return 0;
}


int get_grid_gravitational_acceleration(
    int * i, int * j, int * k,
    int * index_of_grid,
    double * fx, double * fy, double * fz,
    int number_of_points)
{
    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    int ii = 0;
    DomainS * dom = 0;
    if (mesh.NLevels == 0) {
        return -1;
    }
#ifndef SELF_GRAVITY
    for(l=0; l < number_of_points; l++) {
        fx[l] = 0;
        fy[l] = 0;
        fz[l] = 0;
    }
#else
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        current_index_of_grid = index_of_grid[l];
        
        fx[l] = 0;
        fy[l] = 0;
        fz[l] = 0;
        
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {

                
                ijk_on_grid(grid, &i0, &j0, &k0);
                fx[l] = (grid->Phi[k0][j0][i0-1] - grid->Phi[k0][j0][i0+1]) / grid->dx1;
                fy[l] = (grid->Phi[k0][j0-1][i0] - grid->Phi[k0][j0+1][i0]) / grid->dx2;
                fz[l] = (grid->Phi[k0-1][j0][i0] - grid->Phi[k0+1][j0][i0]) / grid->dx3;
            }

        }
    }



#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(fx, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(fy, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(fz, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
     } else {
        MPI_Reduce(MPI_IN_PLACE, fx, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, fy, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, fz, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

#endif
    return 0;
}

int get_grid_density(
    int * i, int * j, int * k,
    int * index_of_grid,
    double * rho,
    int number_of_points)
{
    int l=0;
    int i0,j0,k0 = 0;
    int ii = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    DomainS * dom = 0;
    if (mesh.NLevels == 0) {
        return -1;
    }
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        rho[l] = 0.0;
        
        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {

                
                ijk_on_grid(grid, &i0, &j0, &k0);
                rho[l] = grid->U[k0][j0][i0].d;
            }

        }
    }



#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(rho, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, rho, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}

int get_grid_scalar(
    int * i, int * j, int * k,
    int * index_of_grid,
    double * scalar,
    int number_of_points)
{
    int l=0;
    int i0,j0,k0 = 0;
    int ii = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    DomainS * dom = 0;
    if (mesh.NLevels == 0) {
        return -1;
    }
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        scalar[l] = 0.0;
        
        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {

                ijk_on_grid(grid, &i0, &j0, &k0);
#if (NSCALARS > 0)
                scalar[l] = grid->U[k0][j0][i0].s[0];
#endif
            }

        }
    }



#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(scalar, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, scalar, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}


int get_grid_momentum_density(
    int * i, int * j, int * k,
    int * index_of_grid,
    double * rhovx, double * rhovy, double * rhovz,
    int number_of_points)
{
    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    DomainS * dom = 0;
    if (mesh.NLevels == 0) {
        return -1;
    }
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        rhovx[l] = rhovy[l]= rhovz[l] = 0.0;
        
        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
           
            if (is_on_grid(grid, i0, j0, k0))
            {

                
                ijk_on_grid(grid, &i0, &j0, &k0);

                rhovx[l] = grid->U[k0][j0][i0].M1;
                rhovy[l] = grid->U[k0][j0][i0].M2;
                rhovz[l] = grid->U[k0][j0][i0].M3;
            }

        }
    }



#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(rhovx, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovy, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovz, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, rhovx, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovy, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovz, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}

int get_grid_energy_density(
    int * i, int * j, int * k,
    int * index_of_grid,
    double * en,
    int number_of_points)
{
    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    DomainS * dom = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        en[l] = 0.0;
        
        
        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {

                
                ijk_on_grid(grid, &i0, &j0, &k0);


                en[l] = grid->U[k0][j0][i0].E;
            }

        }
    }



#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(en, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, en, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}




int set_grid_state(
    int * i,
    int * j,
    int * k,
    double * rho,
    double * rhovx, double * rhovy, double * rhovz,
    double * en,
    int * index_of_grid,
    int number_of_points)
{

    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    DomainS * dom = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];

        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {
                ijk_on_grid(grid, &i0, &j0, &k0);
                
                grid->U[k0][j0][i0].d = rho[l];
                grid->U[k0][j0][i0].M1 = rhovx[l];
                grid->U[k0][j0][i0].M2 = rhovy[l];
                grid->U[k0][j0][i0].M3 = rhovz[l];
                grid->U[k0][j0][i0].E = en[l];
            }

        }
    }

    return 0;
}


int set_grid_density(
    int * i,
    int * j,
    int * k,
    double * rho,
    int * index_of_grid,
    int number_of_points)
{

    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    DomainS * dom = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];

        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {
                ijk_on_grid(grid, &i0, &j0, &k0);
                
                grid->U[k0][j0][i0].d = rho[l];
            }

        }
    }

    return 0;
}



int set_grid_scalar(
    int * i,
    int * j,
    int * k,
    double * scalar,
    int * index_of_grid,
    int number_of_points)
{

    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    DomainS * dom = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];

        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {
                ijk_on_grid(grid, &i0, &j0, &k0);

#if (NSCALARS > 0)
                grid->U[k0][j0][i0].s[0] = scalar[l];
#endif
            }

        }
    }

    return 0;
}


int set_grid_energy_density(
    int * i,
    int * j,
    int * k,
    double * en,
    int * index_of_grid,
    int number_of_points)
{

    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    DomainS * dom = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];

        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {
                ijk_on_grid(grid, &i0, &j0, &k0);
                
                grid->U[k0][j0][i0].E = en[l];
            }

        }
    }

    return 0;
}




int set_grid_momentum_density(
    int * i,
    int * j,
    int * k,
    double * rhovx, double * rhovy, double * rhovz,
    int * index_of_grid,
    int number_of_points)
{

    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    DomainS * dom = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];

        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_grid(grid, i0, j0, k0))
            {
                ijk_on_grid(grid, &i0, &j0, &k0);
                
                grid->U[k0][j0][i0].M1 = rhovx[l];
                grid->U[k0][j0][i0].M2 = rhovy[l];
                grid->U[k0][j0][i0].M3 = rhovz[l];
            }

        }
    }

    return 0;
}

int get_potential(
    int * i, int * j, int * k,
    double * potential,
    int number_of_points)
{
    if (mesh.NLevels == 0) {
        return -1;
    }
    if (Potentials == NULL) {
        return -2;
    }

    int imin = 0;
    int imax = mesh.Nx[0]+nghost;
    int jmin, jmax;
    if (mesh.Nx[1] > 1)
    {
        jmin = 0;
        jmax = mesh.Nx[1]+nghost;
    }
    else
    {
        jmin = jmax = 0;
    }
    int kmin, kmax;
    if (mesh.Nx[2] > 1)
    {
        kmin = 0;
        kmax = mesh.Nx[2]+nghost;
    }
    else
    {
        kmin = kmax = 0;
    }
    int l=0;
    int i0,j0,k0 = 0;

    for(l=0; l < number_of_points; l++)
    {
        i0 = i[l] + nghost;
        j0 = j[l];
        k0 = k[l];
        potential[l] = 0.0;

        if(mesh.Nx[1] > 1) {j0 += nghost;}
        if(mesh.Nx[2] > 1) {k0 += nghost;}

        if (
            (i0 >= imin && i0 <= imax) &&
            (j0 >= jmin && j0 <= jmax) &&
            (k0 >= kmin && k0 <= kmax)
        )
        {
            potential[l] = Potentials[k0][j0][i0];
        }
    }


#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(potential, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, potential, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}


int set_potential(
    int * i, int * j, int * k,
    double * potential,
    int number_of_points)
{

    if (mesh.NLevels == 0) {
        return -1;
    }
    if (Potentials == NULL) {
        return -2;
    }

    int imin = 0;
    int imax = mesh.Nx[0] + nghost;
    int jmin, jmax;
    if (mesh.Nx[1] > 1)
    {
        jmin = 0;
        jmax = mesh.Nx[1] + nghost;
    }
    else
    {
        jmin = jmax = 0;
    }
    int kmin, kmax;
    if (mesh.Nx[2] > 1)
    {
        kmin = 0;
        kmax = mesh.Nx[2] + nghost;
    }
    else
    {
        kmin = kmax = 0;
    }
    int l=0;
    int i0,j0,k0 = 0;

    for(l=0; l < number_of_points; l++)
    {
        i0 = i[l] + nghost;
        j0 = j[l];
        k0 = k[l];

        if(mesh.Nx[1] > 1) {j0 += nghost;}
        if(mesh.Nx[2] > 1) {k0 += nghost;}
        if (
            (i0 >= imin && i0 <= imax) &&
            (j0 >= jmin && j0 <= jmax) &&
            (k0 >= kmin && k0 <= kmax)
        )
        {
            Potentials[k0][j0][i0] = potential[l];
        }
    }

    return 0;
}

int esys_roe_adb_hydro(int * index, double * u, double * v, double * w,
  double * h, double * ev, double * rem0, double * rem1, double * rem2,
  double * rem3, double * rem4, double * lem0, double * lem1,
  double * lem2, double * lem3, double * lem4){
  return -2;
}


//static Gas ***Soln=NULL;

#include <math.h>

void _fill_grid_linearwave_1d(GridS *pGrid, DomainS *pDomain, int wave_flag, double amp, double vflow, int  wave_dir){
  ConsS ***Soln;
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,n,m,nx1,nx2,nx3;
  Real d0,p0,u0,v0,w0,h0;
  Real x1,x2,x3,r,ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
#ifdef MHD
  Real bx0,by0,bz0,xfact,yfact;
#endif /* MHD */
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;


/* allocate memory for solution on this level.  If this is root level
 * also allocate memory for RootSoln */

  if ((Soln = (ConsS***)calloc_3d_array(nx3,nx2,nx1,sizeof(ConsS)))==NULL)
    ath_error("[problem]: Error allocating memory for Soln\n");


/* Get eigenmatrix, where the quantities u0 and bx0 are parallel to the
 * wavevector, and v0,w0,by0,bz0 are perpendicular.  Background state chosen
 * carefully so wave speeds are well-separated: Va=1, Vf=2, Vs=0.5 */

  d0 = 1.0;
#ifndef ISOTHERMAL
  p0 = 1.0/Gamma;
  u0 = vflow*sqrt(Gamma*p0/d0);
#else
  u0 = vflow*Iso_csound;
#endif
  v0 = 0.0;
  w0 = 0.0;
#ifdef MHD
  bx0 = 1.0;
  by0 = sqrt(2.0);
  bz0 = 0.5;
  xfact = 0.0;
  yfact = 1.0;
#endif

  for (n=0; n<NWAVE; n++) {
    for (m=0; m<NWAVE; m++) {
      rem[n][m] = 0.0;
      lem[n][m] = 0.0;
    }
  }

#ifdef HYDRO
#ifdef ISOTHERMAL
  esys_roe_iso_hyd(u0,v0,w0,   ev,rem,lem);
#else
  h0 = ((p0/Gamma_1 + 0.5*d0*(u0*u0+v0*v0+w0*w0)) + p0)/d0;
  esys_roe_adb_hyd(u0,v0,w0,h0,ev,rem,lem);
  ath_pout(0,"Ux - Cs = %e, %e\n",ev[0],rem[0][wave_flag]);
  ath_pout(0,"Ux      = %e, %e\n",ev[1],rem[1][wave_flag]);
  ath_pout(0,"Ux + Cs = %e, %e\n",ev[4],rem[4][wave_flag]);
#endif /* ISOTHERMAL */
#endif /* HYDRO */

#ifdef MHD
#if defined(ISOTHERMAL)
  esys_roe_iso_mhd(d0,u0,v0,w0,bx0,by0,bz0,xfact,yfact,ev,rem,lem);
  ath_pout(0,"Ux - Cf = %e, %e\n",ev[0],rem[0][wave_flag]);
  ath_pout(0,"Ux - Ca = %e, %e\n",ev[1],rem[1][wave_flag]);
  ath_pout(0,"Ux - Cs = %e, %e\n",ev[2],rem[2][wave_flag]);
  ath_pout(0,"Ux + Cs = %e, %e\n",ev[3],rem[3][wave_flag]);
  ath_pout(0,"Ux + Ca = %e, %e\n",ev[4],rem[4][wave_flag]);
  ath_pout(0,"Ux + Cf = %e, %e\n",ev[5],rem[5][wave_flag]);
#else
  h0 = ((p0/Gamma_1+0.5*(bx0*bx0+by0*by0+bz0*bz0)+0.5*d0*(u0*u0+v0*v0+w0*w0))
               + (p0+0.5*(bx0*bx0+by0*by0+bz0*bz0)))/d0;
  esys_roe_adb_mhd(d0,u0,v0,w0,h0,bx0,by0,bz0,xfact,yfact,ev,rem,lem);
  ath_pout(0,"Ux - Cf = %e, %e\n",ev[0],rem[0][wave_flag]);
  ath_pout(0,"Ux - Ca = %e, %e\n",ev[1],rem[1][wave_flag]);
  ath_pout(0,"Ux - Cs = %e, %e\n",ev[2],rem[2][wave_flag]);
  ath_pout(0,"Ux      = %e, %e\n",ev[3],rem[3][wave_flag]);
  ath_pout(0,"Ux + Cs = %e, %e\n",ev[4],rem[4][wave_flag]);
  ath_pout(0,"Ux + Ca = %e, %e\n",ev[5],rem[5][wave_flag]);
  ath_pout(0,"Ux + Cf = %e, %e\n",ev[6],rem[6][wave_flag]);
#endif /* ISOTHERMAL */
#endif /* MHD */

/* Now initialize 1D solution vector  */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {

/* Set background state */

    cc_pos(pGrid,i,j,k,&x1,&x2,&x3);
    Soln[k][j][i].d = d0;
#ifndef ISOTHERMAL
    Soln[k][j][i].E = p0/Gamma_1 + 0.5*d0*u0*u0;
#ifdef MHD
    Soln[k][j][i].E += 0.5*(bx0*bx0 + by0*by0 + bz0*bz0);
#endif /* MHD */
#endif /* ISOTHERMAL */

/* Select appropriate solution based on direction of wavevector */

    switch(wave_dir){

/* wave in x1-direction */

    case 1:
      Soln[k][j][i].d += amp*sin(2.0*PI*x1)*rem[0][wave_flag];
#ifndef ISOTHERMAL
      Soln[k][j][i].E += amp*sin(2.0*PI*x1)*rem[4][wave_flag];
#endif /* ISOTHERMAL */
      Soln[k][j][i].M1 = d0*vflow;
      Soln[k][j][i].M2 = 0.0;
      Soln[k][j][i].M3 = 0.0;
      Soln[k][j][i].M1 += amp*sin(2.0*PI*x1)*rem[1][wave_flag];
      Soln[k][j][i].M2 += amp*sin(2.0*PI*x1)*rem[2][wave_flag];
      Soln[k][j][i].M3 += amp*sin(2.0*PI*x1)*rem[3][wave_flag];
#ifdef MHD
      Soln[k][j][i].B1c = bx0;
      Soln[k][j][i].B2c = by0 + amp*sin(2.0*PI*x1)*rem[NWAVE-2][wave_flag];
      Soln[k][j][i].B3c = bz0 + amp*sin(2.0*PI*x1)*rem[NWAVE-1][wave_flag];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Soln[k][j][i].s[n] = amp*(1.0 + sin(2.0*PI*x1));
#endif
      break;

/* Wave in x2-direction */

    case 2:
      Soln[k][j][i].d += amp*sin(2.0*PI*x2)*rem[0][wave_flag];
#ifndef ISOTHERMAL
      Soln[k][j][i].E += amp*sin(2.0*PI*x2)*rem[4][wave_flag];
#endif /* ISOTHERMAL */
      Soln[k][j][i].M1 = 0.0;
      Soln[k][j][i].M2 = d0*vflow;
      Soln[k][j][i].M3 = 0.0;
      Soln[k][j][i].M1 += amp*sin(2.0*PI*x2)*rem[3][wave_flag];
      Soln[k][j][i].M2 += amp*sin(2.0*PI*x2)*rem[1][wave_flag];
      Soln[k][j][i].M3 += amp*sin(2.0*PI*x2)*rem[2][wave_flag];
#ifdef MHD
      Soln[k][j][i].B1c = bz0 + amp*sin(2.0*PI*x2)*rem[NWAVE-1][wave_flag];
      Soln[k][j][i].B2c = bx0;
      Soln[k][j][i].B3c = by0 + amp*sin(2.0*PI*x2)*rem[NWAVE-2][wave_flag];
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Soln[k][j][i].s[n] = amp*(1.0 + sin(2.0*PI*x2));
#endif
      break;

/* Wave in x3-direction */

    case 3:
      Soln[k][j][i].d += amp*sin(2.0*PI*x3)*rem[0][wave_flag];
#ifndef ISOTHERMAL
      Soln[k][j][i].E += amp*sin(2.0*PI*x3)*rem[4][wave_flag];
#endif /* ISOTHERMAL */
      Soln[k][j][i].M1 = 0.0;
      Soln[k][j][i].M2 = 0.0;
      Soln[k][j][i].M3 = d0*vflow;
      Soln[k][j][i].M1 += amp*sin(2.0*PI*x3)*rem[2][wave_flag];
      Soln[k][j][i].M2 += amp*sin(2.0*PI*x3)*rem[3][wave_flag];
      Soln[k][j][i].M3 += amp*sin(2.0*PI*x3)*rem[1][wave_flag];
#ifdef MHD
      Soln[k][j][i].B1c = by0 + amp*sin(2.0*PI*x3)*rem[NWAVE-2][wave_flag];
      Soln[k][j][i].B2c = bz0 + amp*sin(2.0*PI*x3)*rem[NWAVE-1][wave_flag];
      Soln[k][j][i].B3c = bx0;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        Soln[k][j][i].s[n] = amp*(1.0 + sin(2.0*PI*x3));
#endif
      break;
    default:
      ath_error("[linear_wave1d]: wave direction %d not allowed\n",wave_dir);
    }
  }}}

/* Now set initial conditions to wave solution */

  for (k=ks; k<=ke; k++) {
  for (j=js; j<=je; j++) {
  for (i=is; i<=ie; i++) {
    pGrid->U[k][j][i].d = Soln[k][j][i].d;
#ifndef ISOTHERMAL
    pGrid->U[k][j][i].E = Soln[k][j][i].E;
#endif /* ISOTHERMAL */
    pGrid->U[k][j][i].M1 = Soln[k][j][i].M1;
    pGrid->U[k][j][i].M2 = Soln[k][j][i].M2;
    pGrid->U[k][j][i].M3 = Soln[k][j][i].M3;
#ifdef MHD
    pGrid->U[k][j][i].B1c = Soln[k][j][i].B1c;
    pGrid->U[k][j][i].B2c = Soln[k][j][i].B2c;
    pGrid->U[k][j][i].B3c = Soln[k][j][i].B3c;
    pGrid->B1i[k][j][i] = Soln[k][j][i].B1c;
    pGrid->B2i[k][j][i] = Soln[k][j][i].B2c;
    pGrid->B3i[k][j][i] = Soln[k][j][i].B3c;
#endif /* MHD */
#if (NSCALARS > 0)
      for (n=0; n<NSCALARS; n++)
        pGrid->U[k][j][i].s[n] = Soln[k][j][i].s[n];
#endif
  }}}
#ifdef MHD
  for (k=ks; k<=ke; k++) {
    for (j=js; j<=je; j++) {
      pGrid->B1i[k][j][ie+1] = pGrid->B1i[k][j][ie];
    }
  }
  if (pGrid->Nx[1] > 1) {
    for (k=ks; k<=ke; k++) {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][je+1][i] = pGrid->B2i[k][je][i];
      }
    }
  }
  if (pGrid->Nx[2] > 1) {
    for (j=js; j<=je; j++) {
      for (i=is; i<=ie; i++) {
        pGrid->B3i[ke+1][j][i] = pGrid->B3i[ke][j][i];
      }
    }
  }
#endif /* MHD */

/* For self-gravitating problems, read 4\piG and compute mean density */

#ifdef SELF_GRAVITY
  four_pi_G = par_getd("problem","four_pi_G");
  grav_mean_rho = d0;
#endif /* SELF_GRAVITY */

/* With viscosity and/or resistivity, read eta_R and nu_V */

#ifdef OHMIC
  eta_Ohm = par_getd("problem","eta");
#endif
#ifdef HALL_MHD
  eta_Ohm = par_getd("problem","eta_O");
  eta_Hall = par_getd("problem","eta_H");
#endif
#ifdef NAVIER_STOKES
  nu_V = par_getd("problem","nu");
#endif


  free_3d_array(Soln);
  return;
}


int fill_grid_linearwave_1d(int wave_flag, double amp, double vflow, int wave_dir)
{
    int level = 0;
    int domain = 0;
    for(level = 0 ; level < mesh.NLevels; level++)
    {
        for(domain = 0; domain < mesh.DomainsPerLevel[level]; domain++)
        {
            if(!(mesh.Domain[level][domain].Grid == NULL))
            {
                _fill_grid_linearwave_1d( mesh.Domain[level][domain].Grid,&mesh.Domain[level][domain], wave_flag, amp, vflow, wave_dir);
            }
        }
    }
    return 0;

}

int get_time(double * value){
    *value = mesh.time;
    return 0;
}

int evolve_model(double tlim) {
    int nl, nd;
    //AMUSE STOPPING CONDITIONS
    int is_number_of_steps_detection_enabled;
    int is_timeout_detection_enabled;
    int number_of_steps_innerloop = 0;
    int max_number_of_steps;
    double timeout;
    time_t clock_current, clock_init;
    int error;
    
    if(evolve_to_exact_time) {
        par_setd("time","tlim", "%.15e", tlim, "-");
        tlim=par_getd("time","tlim"); /* this fixes accuracy problem by using the par_setd stuff */
    } else {
        par_setd("time","tlim", "%.15e", 10.0 * tlim, "-"); /* far away */
    }

    //AMUSE STOPPING CONDITIONS SUPPORT
    error = is_stopping_condition_enabled(NUMBER_OF_STEPS_DETECTION,
                                          &is_number_of_steps_detection_enabled);
    error = is_stopping_condition_enabled(TIMEOUT_DETECTION, 
					  &is_timeout_detection_enabled);
    get_stopping_condition_number_of_steps_parameter(&max_number_of_steps);
    get_stopping_condition_timeout_parameter(&timeout);    
    time(&clock_init);
    if(evolve_to_exact_time && mesh.time + mesh.dt > tlim) {
        if(mesh.time < tlim)
        {
            mesh.dt = (tlim - mesh.time) / 2.0;
            fprintf(stderr, "set mesh.time %f, %f, %f\n", mesh.time, mesh.dt, tlim);
        }
        else 
        {
            return 0;
        }
    }
    if(mesh.dt == 0.0)
    {
        mesh.dt = last_dt_above_zero;
        new_dt(&mesh); 
        fprintf(stderr, "new mesh.time %f, %f (%f), %f\n", mesh.time, mesh.dt, last_dt_above_zero, tlim);
 
    }
    
    for (nl=0; nl<(mesh.NLevels); nl++){
      for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){
        if (mesh.Domain[nl][nd].Grid != NULL){
          bvals_mhd(&(mesh.Domain[nl][nd]));
#ifdef PARTICLES
          bvals_particle(&level0_Grid, &level0_Domain);
#endif
        }
      }
    }
    
#ifdef STATIC_MESH_REFINEMENT
    Prolongate(&mesh);
#endif

    while (mesh.time < tlim) {
        //fprintf(stderr, "mesh.time %g, %g, %g\n", mesh.time, mesh.dt, tlim);
        if(mesh.dt == 0.0)
        {
            break;
        }
    /*--- Step 9b. ---------------------------------------------------------------*/
    /* operator-split explicit diffusion: thermal conduction, viscosity, resistivity
     * Done first since CFL constraint is applied which may change dt  */

#if defined(RESISTIVITY) || defined(VISCOSITY) || defined(THERMAL_CONDUCTION)
            integrate_diff(&mesh);
            for (nl=0; nl<(mesh.NLevels); nl++){
              for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){
                if (mesh.Domain[nl][nd].Grid != NULL){
                  bvals_mhd(&(mesh.Domain[nl][nd]));
                }
              }
            }
#endif

    /*--- Step 9c. ---------------------------------------------------------------*/
    /* Loop over all Domains and call Integrator */

        for (nl=0; nl<(mesh.NLevels); nl++){
          for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){
            if (mesh.Domain[nl][nd].Grid != NULL){
              (*Integrate)(&(mesh.Domain[nl][nd]));

#ifdef FARGO
              Fargo(&(mesh.Domain[nl][nd]));
#ifdef PARTICLES
              advect_particles(&level0_Grid, &level0_Domain);
#endif
#endif /* FARGO */
            }
          }
        }

    /*--- Step 9d. ---------------------------------------------------------------*/
    /* With SMR, restrict solution from Child --> Parent grids  */

#ifdef STATIC_MESH_REFINEMENT
        RestrictCorrect(&mesh);
#endif

    /*--- Step 9e. ---------------------------------------------------------------*/
    /* User work (defined in problem()) */

        Userwork_in_loop(&mesh);

    /*--- Step 9f. ---------------------------------------------------------------*/
    /* Compute gravitational potential using new density, and add second-order
     * correction to fluxes for accelerations due to self-gravity. */

#ifdef SELF_GRAVITY
        for (nl=0; nl<(mesh.NLevels); nl++){
          for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){
            if (mesh.Domain[nl][nd].Grid != NULL){
              (*SelfGrav)(&(mesh.Domain[nl][nd]));
              bvals_grav(&(mesh.Domain[nl][nd]));
              selfg_fc(&(mesh.Domain[nl][nd]));
            }
          }
        }
#endif

    /*--- Step 9g. ---------------------------------------------------------------*/
    /* Update mesh time, and time in all Grid's.  Compute new dt */

        mesh.nstep++;
        mesh.time += mesh.dt;
        for (nl=0; nl<(mesh.NLevels); nl++){
          for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){
            if (mesh.Domain[nl][nd].Grid != NULL){
              mesh.Domain[nl][nd].Grid->time = mesh.time;
            }
          }
        }

        last_dt_above_zero = mesh.dt;
        if(mesh.time < tlim)
        {
            new_dt(&mesh);
        }
    /*--- Step 9h. ---------------------------------------------------------------*/
    /* Boundary values must be set after time is updated for t-dependent BCs.
     * With SMR, ghost zones at internal fine/coarse boundaries set by Prolongate */

        for (nl=0; nl<(mesh.NLevels); nl++){
          for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){
            if (mesh.Domain[nl][nd].Grid != NULL){
              bvals_mhd(&(mesh.Domain[nl][nd]));
#ifdef PARTICLES
              bvals_particle(&level0_Grid, &level0_Domain);
#endif
            }
          }
        }

#ifdef STATIC_MESH_REFINEMENT
        Prolongate(&mesh);
#endif
        //AMUSE STOPPING CONDITIONS SUPPORT
        if (is_number_of_steps_detection_enabled) {
            number_of_steps_innerloop++;
            if (number_of_steps_innerloop >= max_number_of_steps) {
                int stopping_index  = next_index_for_stopping_condition();
                set_stopping_condition_info(stopping_index, NUMBER_OF_STEPS_DETECTION);
            }
        }
        if (is_timeout_detection_enabled) {
	    time(&clock_current);
	    if ((clock_current - clock_init) > timeout) {
                int stopping_index  = next_index_for_stopping_condition();
                set_stopping_condition_info(stopping_index, TIMEOUT_DETECTION);
            }
        }

	if(set_conditions & enabled_conditions) {
	    break;
	}
    }

    
    return 0;
}


ConsS *** get_boundary_with_index(GridS * grid, int index_of_boundary)
{
    BoundaryCellS * boundary = grid->boundary;
    switch(index_of_boundary) {
        case 1:
            return boundary->LeftX1;
        case 2:
            return boundary->RightX1;
        case 3:
            return boundary->LeftX2;
        case 4:
            return boundary->RightX2;
        case 5:
            return boundary->LeftX3;
        case 6:
            return boundary->RightX3;
        default:
            fprintf(stderr, "Error, incorrect boundary index");
            return NULL;
    }
}

int set_boundary_state(
    int * i,
    int * j,
    int * k,
    double * rho,
    double * rhovx, double * rhovy, double * rhovz,
    double * en,
    int * index_of_boundary,
    int * index_of_grid,
    int number_of_points)
{

    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    DomainS * dom = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];

        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            ConsS *** U = get_boundary_with_index(grid, index_of_boundary[l]);
            if(U == 0) {    
                continue;
            }
            
            if (is_on_boundary_grid(dom, grid, index_of_boundary[l], i0, j0, k0, 1)) 
            {
                ijk_on_boundary_grid(grid, index_of_boundary[l], &i0, &j0, &k0);
                //fprintf(stderr, "i,j,k %d, %d, %d - %p\n", i0, j0, k0, U);
                U[k0][j0][i0].d = rho[l];
                U[k0][j0][i0].M1 = rhovx[l];
                U[k0][j0][i0].M2 = rhovy[l];
                U[k0][j0][i0].M3 = rhovz[l];
                U[k0][j0][i0].E = en[l];
            }

        }
    }

    return 0;
}


int get_boundary_state(
    int * i,
    int * j,
    int * k,
    int * index_of_boundary,
    int * index_of_grid,
    double * rho,
    double * rhovx, double * rhovy, double * rhovz,
    double * en,
    int number_of_points)
{

    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    DomainS * dom = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];

        rho[l] =  0.0;
        rhovx[l] = 0.0;
        rhovy[l] = 0.0;
        rhovz[l] = 0.0;
        en[l] = 0.0;
        
        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            ConsS *** U = get_boundary_with_index(grid, index_of_boundary[l]);
            //fprintf(stderr, "bndry %d - %p\n", index_of_boundary[l], U);
            if(U == 0) {
                continue;
            }
            
            if (is_on_boundary_grid(dom, grid, index_of_boundary[l], i0, j0, k0, 0)) 
            {
                ijk_on_boundary_grid(grid, index_of_boundary[l], &i0, &j0, &k0);
                //fprintf(stderr, "i,j,k %d, %d, %d - %p\n", i0, j0, k0, U);
                rho[l] = U[k0][j0][i0].d;
                rhovx[l] = U[k0][j0][i0].M1;
                rhovy[l] = U[k0][j0][i0].M2;
                rhovz[l] = U[k0][j0][i0].M3;
                en[l] = U[k0][j0][i0].E;
            }

        }
    }

#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(rho, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovx, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovy, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovz, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(en, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, rho, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovx, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovy, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovz, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, en, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif
    return 0;
}


int get_grid_magnetic_field(
    int * i, int * j, int * k,
    int * index_of_grid,
    double * B1i, double * B2i, double * B3i,
    int number_of_points)
{
    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    DomainS * dom = 0;
    if (mesh.NLevels == 0) {
        return -1;
    }
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        B1i[l] = B2i[l]= B3i[l] = 0.0;
        
        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
           
            if (is_on_magnetic_grid(grid, i0, j0, k0))
            {
                ijk_on_grid(grid, &i0, &j0, &k0);

#ifdef MHD
                B1i[l] = grid->B1i[k0][j0][i0];
                B2i[l] = grid->B2i[k0][j0][i0];
                B3i[l] = grid->B3i[k0][j0][i0];
#endif
            }

        }
    }



#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(B1i, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(B2i, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(B3i, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, B1i, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, B2i, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, B3i, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}

int set_grid_magnetic_field(
    int * i,
    int * j,
    int * k,
    double * B1i,
    double * B2i,
    double * B3i,
    int * index_of_grid,
    int number_of_points)
{

    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    DomainS * dom = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];

        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            
            if (is_on_magnetic_grid(grid, i0, j0, k0))
            {
                ijk_on_grid(grid, &i0, &j0, &k0);

#ifdef MHD								
                grid->B1i[k0][j0][i0] = B1i[l];
                grid->B2i[k0][j0][i0] = B2i[l];
                grid->B3i[k0][j0][i0] = B3i[l];
#endif
            }

        }
    }

    return 0;
}  

static char * boundary_names[] = {"bc_ix1","bc_ox1","bc_ix2","bc_ox2","bc_ix3","bc_ox3"};

int get_boundary_index_range_inclusive(
    int index_of_boundary,
    int index_of_grid,
    int * minx,
    int * maxx,
    int * miny,
    int * maxy,
    int * minz,
    int * maxz
)
{
    if(index_of_boundary < 1) {return -1;}
    if(index_of_boundary > 6) {return -1;}
    const char * boundary_name = boundary_names[index_of_boundary-1];
    int boundary_type = par_geti( "domain1", (char *) boundary_name);
    
    if(index_of_grid < 1) {return -1;}
    DomainS * dom = get_domain_structure_with_index(index_of_grid);
    
    if(dom == 0)
    {
        return -2;
    }
        
    //printf("boundary name: %s , type: %d\n", boundary_name, boundary_type);
    if(boundary_type == 10) {
        *minx = 0;
        *miny = 0;
        *minz = 0;
        switch(index_of_boundary) {
            case 1:
            case 2:
                *maxx = nghost - 1;
                *maxy = dom->Nx[1] - 1;
                *maxz = dom->Nx[2] - 1;
                break;
            case 3:
            case 4:
                *maxx = dom->Nx[0] - 1 + (2 * nghost); 
                *maxy = nghost - 1;
                *maxz = dom->Nx[2] - 1;
                break;
            case 5:
            case 6:
                *maxx = dom->Nx[0] - 1 + (2 * nghost);
                *maxy = dom->Nx[1] - 1 + (dom->Nx[1]  > 1 ? 2 * nghost : 0);
                *maxz = nghost - 1;
                break;
            default:
                *maxx = 0;
                *maxy = 0;
                *maxz = 0;
            
        }
        
    } else {
        *minx = 0;
        *maxx = 0;
        *miny = 0;
        *maxy = 0;
        *minz = 0;
        *maxz = 0;
    }
    return 0;
}



int get_boundary_position_of_index(
    int *i, int *j, int *k, 
    int *index_of_boundary, 
    int *index_of_grid, 
    double * x, double * y, double * z, 
    int number_of_points){


    int l=0;
    int i0,j0,k0 = 0;
    int previous_index_of_grid = -1, current_index_of_grid = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    DomainS * dom = 0;
    if (mesh.NLevels == 0) {
        return -1;
    }
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        x[l] = 0.0;
        y[l] = 0.0;
        z[l] = 0.0;
        
        current_index_of_grid = index_of_grid[l];
        if (current_index_of_grid != previous_index_of_grid)
        {
            dom = get_domain_structure_with_index(current_index_of_grid);
        }
        if(dom == 0)
        {
            continue;
        }
        

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            if (is_on_boundary_grid(dom,grid, index_of_boundary[l], i0, j0, k0, 0))
            {
                switch(index_of_boundary[l]) {
                    case 1:
                        i0 = i0 - nghost;
                        break;
                    case 2:
                        i0 = i0 + (grid->ie-grid->is) + 1 ;
                        break;
                    case 3:
                        i0 = i0 - nghost;
                        j0 = j0 - nghost;
                        break;
                    case 4:
                        i0 = i0 - nghost;
                        j0 = j0 + (grid->je-grid->js) + 1 ;
                        break;
                    case 5:
                        i0 = i0 - nghost;
                        j0 = j0 - nghost;
                        k0 = k0 - nghost;
                        break;
                    case 6:
                        i0 = i0 - nghost;
                        j0 = j0 - nghost;
                        k0 = k0 + (grid->ke-grid->ks) + 1 ;
                        break;
                    default:
                        i0 = 1;
                        j0 = 1;
                        k0 = 1;
                }
                ijk_on_grid(grid, &i0, &j0, &k0);
                cc_pos(
                    grid,
                    i0,
                    j0,
                    k0,
                    &x[l],
                    &y[l],
                    &z[l]
                );
            }

        }
    }

#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(x, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(y, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(z, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, x, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, y, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, z, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}


/* k,j,i offsets for weighing */
static int weighing_offsets[8][3] = {
    {0,0,0},
    {0,0,1},
    {0,1,0},
    {0,1,1},
    {1,0,0},
    {1,0,1},
    {1,1,0},
    {1,1,1},
};

int get_hydro_state_at_point(
        double * x, double * y, double *z,
        double * vx, double * vy, double * vz,
        double * rho,
        double * rhovx, double * rhovy, double * rhovz,
        double * rhoen,
        int number_of_points)
{
    DomainS * dom = 0;
    int ii = 0, level = 0, i = 0;
    int igrid, jgrid, kgrid;
    int nweigh = 0;
    for(i=0; i < number_of_points; i++) {
    for(level = mesh.NLevels - 1; level >= 0; level--)
    {
        rho[i]   = 0.0;
        rhovx[i] = 0.0;
        rhovy[i] = 0.0;
        rhovz[i] = 0.0;
        rhoen[i] = 0.0;
        for(ii = 0; ii < mesh.DomainsPerLevel[level]; ii++)
        {   
            dom = (DomainS*)&(mesh.Domain[level][ii]);
            if(dom->Grid == NULL)
            {
                continue;
            }
            else
            {
                GridS * grid = dom->Grid;
                if(
                    (grid->MinX[0] <= x[i] && x[i] < grid->MaxX[0]) &&
                    ((dom->Nx[1] <= 1) || (grid->MinX[1] <= y[i] && y[i] < grid->MaxX[1])) &&
                    ((dom->Nx[2] <= 1) || (grid->MinX[2] <= z[i] && z[i] < grid->MaxX[2])) 
                )
                {
                    double yval = 0.0;
                    double zval = 0.0;
                    double xval = (x[i] - grid->MinX[0] - (0.5 * grid->dx1))/grid->dx1;
                    igrid = floor(xval) + grid->is;
                    nweigh = 2;
                    if(dom->Nx[1] > 1) {
                        yval = (y[i] - grid->MinX[1] - (0.5 * grid->dx2))/grid->dx2;
                        jgrid = floor(yval) + grid->js;
                        nweigh = 4;
                    } else {
                        jgrid = 0;
                    }
                    if(dom->Nx[2] > 1) {
                        zval = (z[i] - grid->MinX[2] - (0.5 * grid->dx3))/grid->dx3;
                        kgrid = floor(zval) + grid->ks;
                        nweigh = 8;
                    } else {
                        kgrid = 0;
                    }
                    double x1dx = xval - floor(xval);
                    double x0dx = 1.0 - x1dx;
                    double y1dy = yval - floor(yval);
                    double y0dy = 1.0 - y1dy;
                    double z1dz = zval - floor(zval);
                    double z0dz = 1.0 - z1dz;
                    double weighingfactors[8];
                    
                    weighingfactors[0] = z0dz * y0dy * x0dx;
                    weighingfactors[1] = z0dz * y0dy * x1dx;
                    weighingfactors[2] = z0dz * y1dy * x0dx;
                    weighingfactors[3] = z0dz * y1dy * x1dx;
                    weighingfactors[4] = z1dz * y0dy * x0dx; 
                    weighingfactors[5] = z1dz * y0dy * x1dx; 
                    weighingfactors[6] = z1dz * y1dy * x0dx; 
                    weighingfactors[7] = z1dz * y1dy * x1dx; 
                    int k = 0;
                    for(k = 0; k < nweigh; k ++) {
                        int ioffset = weighing_offsets[k][2];
                        int joffset = weighing_offsets[k][1];
                        int koffset = weighing_offsets[k][0];
                        rho[i]   += weighingfactors[k] * grid->U[kgrid+koffset][jgrid+joffset][igrid+ioffset].d;
                        rhovx[i] += weighingfactors[k] * grid->U[kgrid+koffset][jgrid+joffset][igrid+ioffset].M1;
                        rhovy[i] += weighingfactors[k] * grid->U[kgrid+koffset][jgrid+joffset][igrid+ioffset].M2;
                        rhovz[i] += weighingfactors[k] * grid->U[kgrid+koffset][jgrid+joffset][igrid+ioffset].M3;
                        rhoen[i] += weighingfactors[k] * grid->U[kgrid+koffset][jgrid+joffset][igrid+ioffset].E;
                    }
                }
            }
        }
    }
    }
#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(rho  , NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovx, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovy, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhovz, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(rhoen, NULL, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, rho  , number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovx, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovy, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhovz, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(MPI_IN_PLACE, rhoen, number_of_points, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif

    return 0;
}
