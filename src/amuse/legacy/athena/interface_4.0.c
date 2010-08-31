#define MAIN_C

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>

#include "defs.h"
#include "athena.h"
#include "globals.h"
#include "prototypes.h"

#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

#include <math.h>
#include <stopcond.h>

static MeshS mesh;
static VDFun_t Integrate;
static int is_dt_set_by_script=0;
static int is_restart = 0;
static int has_external_gravitational_potential=0;
static int mpi_comm_size = 0;

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
  }
  if(pG->dx2 == 0) {
     *j = 0;
  } else {
    *j = ((x2 - pG->MinX[1])/pG->dx2) - 0.5 + pG->Disp[1];
  }
  
  if(pG->dx3 == 0) {
    *k = 0.0;
  } else {
    *k = ((x3 - pG->MinX[2])/pG->dx3) - 0.5 + pG->Disp[2];
  }
}


static inline void ijk_pos_dom(
    const Real x1, const Real x2, const Real x3,
    Real * i, Real * j, Real * k
)
{
  DomainS * domain = &mesh.Domain[0][0];
  
  if(domain->dx[0] == 0) {
    *i = 0;
  } else {
    *i = ((x1 - domain->MinX[0])/domain->dx[0]) - 0.5 + nghost;
  }
  if(domain->dx[1] == 0) {
     *j = 0;
  } else {
    *j = ((x2 - domain->MinX[1])/domain->dx[1]) - 0.5 + nghost;
  }
  
  if(domain->dx[2] == 0) {
    *k = 0.0;
  } else {
    *k = ((x3 - domain->MinX[2])/domain->dx[2]) - 0.5 + nghost;
  }
}

static Real ***Potentials=NULL;

static Real grav_pot(const Real x1, const Real x2, const Real x3)
{
    Real i, j, k;
    ijk_pos_dom(x1,x2,x3, &i, &j, &k);
    int ii, jj, kk;
    ii = i;
    jj = j;
    kk = k;
    
    Real potential0 = Potentials[kk][jj][ii];
    Real potential1 = potential0;
    if(ii < i) {
        potential1 = Potentials[kk][jj][ii + 1];
    } else if(jj < j) {
        potential1 = Potentials[kk][jj + 1][ii];
    } else if(kk < k) {
        potential1 = Potentials[kk + 1][jj][ii];
    } else {
        potential1 = potential0;
    }

    return (potential0 + potential1) / 2.0;
    
    return 0;
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

int initialize_code(){

  par_open("/dev/null"); /* to trick athena into thinking it has opened a parameter file, will not work on windows */
  
  par_sets("job","problem_id", "amuse", "all amuse runs");
  is_restart = 0;
  show_config_par();   /* Add the configure block to the parameter database */
  

  // AMUSE STOPPING CONDITIONS SUPPORT
  set_support_for_condition(NUMBER_OF_STEPS_DETECTION);
  
#ifdef MPI_PARALLEL
    /* Get my task id (rank in MPI) */
    if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD,&myID_Comm_world))
    {
        ath_error("Error on calling MPI_Comm_rank\n");
    }
    if(MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD,&(mpi_comm_size)))
    {
        ath_error("Error on calling MPI_Comm_size\n");
    }
#else
    myID_Comm_world = 0;
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

int commit_parameters(){
  
  int out_level = par_geti_def("log","out_level",0);
  int err_level = par_geti_def("log","err_level",0);
  
  ath_log_set_level(out_level, err_level);
  
  /*
   * if(has_external_gravitational_potential) {
    StaticGravPot = grav_pot;
  }
  *
  */
  
  CourNo = par_getd("time","cour_no");

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
  SelfGrav = selfg_init(&Mesh);
  for (nl=0; nl<(mesh.NLevels); nl++){ 
    for (nd=0; nd<(mesh.DomainsPerLevel[nl]); nd++){  
      ifmesh.Domain[nl][nd].Grid != NULL){
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

int get_position_of_index(int i, int j, int k, int level, int domain, double * x, double * y, 
  double * z){
 
    double pos[4] = {0,0,0,0};
    
    if (mesh.NLevels == 0) {
        return -1;
    }

    if (level >= mesh.NLevels){
        return -2;
    }

    if (domain >= mesh.DomainsPerLevel[level]){
        return -3;
    }
    DomainS * dom = (DomainS*)&(mesh.Domain[level][domain]);

    if(dom->Grid == NULL)
    {
        pos[0] = 0.0;
        pos[1] = 0.0;
        pos[2] = 0.0;
        pos[3] = 1;
    }
    else
    {
        GridS * grid = dom->Grid;
        //fprintf(stderr, "%d, %d, %d, %d\n", grid->is, grid->ie, grid->Nx[0], grid->Disp[0]);
        
        if (grid->Nx[0] > 1 && ( (grid->Disp[0] > 0 && i < (grid->Disp[0]))  || (((grid->Disp[0] + grid->Nx[0]) < dom->Nx[0]) && i >= (grid->Disp[0] + grid->Nx[0]))))
        {
            pos[0] = 0.0;
            pos[1] = 0.0;
            pos[2] = 0.0;
            pos[3] = 1;
        }
        else if (grid->Nx[1] > 1 && ((grid->Disp[1] > 0 && j < (grid->Disp[1]))  || (((grid->Disp[1] + grid->Nx[1]) < dom->Nx[1]) &&j >= (grid->Disp[1] + grid->Nx[1]))))
        {
            pos[0] = 0.0;
            pos[1] = 0.0;
            pos[2] = 0.0;
            pos[3] = 1;
        }
        else if (grid->Nx[2] > 1 && ((grid->Disp[2] > 0 && k < (grid->Disp[2]))  || (((grid->Disp[2] + grid->Nx[2]) < dom->Nx[2]) &&k >= (grid->Disp[2] + grid->Nx[2]))))
        {
            pos[0] = 0.0;
            pos[1] = 0.0;
            pos[2] = 0.0;
            pos[3] = 1;
        }
        else
        {
            cc_pos(
                grid,
                i + grid->is - grid->Disp[0],
                j + grid->js - grid->Disp[1],
                k + grid->ks - grid->Disp[2],
                &pos[0] ,
                &pos[1] ,
                &pos[2] 
            );
            pos[3] = 0;
        }
    }
    
    
#ifdef MPI_PARALLEL
    if(myID_Comm_world) {
        MPI_Reduce(pos, NULL, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    } else {
        MPI_Reduce(MPI_IN_PLACE, pos, 4, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }
#endif
    if(pos[3] == mpi_comm_size) {
        return -4;
    }
    *x = pos[0];
    *y = pos[1];
    *z = pos[2];
    
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
  double z,int level, int domain, double *i , double * j, double * k)
  {
      
    double pos[3] = {0,0,0};
    
    if (mesh.NLevels == 0) {
        return -1;
    }

    if (level >= mesh.NLevels){
        return -2;
    }

    if (domain >= mesh.DomainsPerLevel[level]){
        return -3;
    }
    DomainS * dom = (DomainS*)&(mesh.Domain[level][domain]);

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

int fill_grid_state(
  int i, int j, int k, 
  double rho, double rhovx, double rhovy, double rhovz, 
  double en){
  return -1;
}

int get_grid_state_mpi(
    int * i, int * j, int * k, 
    int * levels, int * domains,
    double * rho, 
    double * rhovx, double * rhovy, double * rhovz, 
    double * en, 
    int number_of_points)
{
    int l=0;
    int i0,j0,k0 = 0;
    int level = 0;
    int domain = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        level = levels[l];
        domain = domains[l];
        
        rho[l] = rhovx[l] = rhovy[l]= rhovz[l] = en[l] = 0.0;
            
        if (level >= mesh.NLevels){
            continue;
        }

        if (domain >= mesh.DomainsPerLevel[level]){
            continue;
        }
        
        DomainS * dom = (DomainS*)&(mesh.Domain[level][domain]);

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            if (grid->Nx[0] > 1 && (i0 < (grid->Disp[0])  || i0 >= (grid->Disp[0] + grid->Nx[0])))
            {
                
            }
            else if (grid->Nx[1] > 1 && (j0 < (grid->Disp[1])  || j0 >= (grid->Disp[1] + grid->Nx[1])))
            {
                
            }
            else if (grid->Nx[2] > 1 && (k0 < (grid->Disp[2])  || k0 >= (grid->Disp[2] + grid->Nx[2])))
            {
            }
            else
            {
                
                i0 -= grid->Disp[0] - grid->is;
                if(grid->Nx[1] > 1)
                {
                    j0 -= grid->Disp[1] - grid->js;
                }
                else
                {
                    j0 = 0;
                }
                if(grid->Nx[2] > 1)
                {
                    k0 -= grid->Disp[2] - grid->ks;
                }
                else
                {
                    k0 = 0;
                }
                
                
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


int fill_grid_state_mpi(
    int * i, 
    int * j, 
    int * k,
    double * rho, 
    double * rhovx, double * rhovy, double * rhovz, 
    double * en, 
    int * levels,
    int * domains, 
    int number_of_points)
{
    
    int l=0;
    int i0,j0,k0 = 0;
    int level = 0;
    int domain = 0;
    
    if (mesh.NLevels == 0) {
        return -1;
    }
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        level = levels[l];
        domain = domains[l];
        
        if (level >= mesh.NLevels){
            continue;
        }

        if (domain >= mesh.DomainsPerLevel[level]){
            continue;
        }
        
        DomainS * dom = (DomainS*)&(mesh.Domain[level][domain]);

        if(dom->Grid == NULL)
        {
            continue;
        }
        else
        {
            GridS * grid = dom->Grid;
            if (grid->Nx[0] > 1 && (i0 < (grid->Disp[0])  || i0 >= (grid->Disp[0] + grid->Nx[0])))
            {
                
            }
            else if (grid->Nx[1] > 1 && (j0 < (grid->Disp[1])  || j0 >= (grid->Disp[1] + grid->Nx[1])))
            {
                
            }
            else if (grid->Nx[2] > 1 && (k0 < (grid->Disp[2])  || k0 >= (grid->Disp[2] + grid->Nx[2])))
            {
            }
            else
            {   
                if(grid->Nx[0] > 1)
                {
                    i0 -= grid->Disp[0] - grid->is;
                }
                else
                {
                    i0 = 0;
                }
                if(grid->Nx[1] > 1)
                {
                    j0 -= grid->Disp[1] - grid->js;
                }
                else
                {
                    j0 = 0;
                }
                if(grid->Nx[2] > 1)
                {
                    k0 -= grid->Disp[2] - grid->ks;
                }
                else
                {
                    k0 = 0;
                }
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


int get_grid_state(
  int i, int j, int k, 
  double * rho, double * rhovx, double * rhovy, double * rhovz, 
  double * en){
      
  return -2;
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

int evolve(double tlim) {
    int nl, nd;
    Real dt_done;
    //AMUSE STOPPING CONDITIONS
    int is_number_of_steps_detection_enabled;
    int number_of_steps_innerloop = 0;
    int max_number_of_steps;
    int error;

    par_setd("time","tlim", "%.15e", tlim, "-");
    
    error = is_stopping_condition_enabled(NUMBER_OF_STEPS_DETECTION, 
					  &is_number_of_steps_detection_enabled);
    get_stopping_condition_number_of_steps_parameter(&max_number_of_steps);    

    while (mesh.time < tlim) {
        
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

        dt_done = mesh.dt;
        new_dt(&mesh);

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

	//SC
	if (is_number_of_steps_detection_enabled) {
	    number_of_steps_innerloop++;
	    if (number_of_steps_innerloop > max_number_of_steps) {
		int stopping_index  = next_index_for_stopping_condition();
		set_stopping_condition_info(stopping_index, NUMBER_OF_STEPS_DETECTION);
	    }
	}
	
    }

    
    return 0;
}

