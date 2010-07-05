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

static Grid level0_Grid;      /* only level0 Grid and Domain in this version */
static Domain level0_Domain;
static VGFun_t Integrate;
static int is_dt_set_by_script=0;
static int is_restart = 0;
static int has_external_gravitational_potential=0;

static inline void ijk_pos( const Grid *pG,
    const Real x1, const Real x2, const Real x3,
    Real * i, Real * j,Real * k
)
{
  *i = ((x1 - pG->x1_0)/pG->dx1) - 0.5 - pG->idisp;
  *j = ((x2 - pG->x2_0)/pG->dx2) - 0.5 - pG->jdisp;
  if(pG->dx3 == 0) {
    *k = 0.0;
  } else {
    *k = ((x3 - pG->x3_0)/pG->dx3) - 0.5 - pG->kdisp;
  }
}

static Real ***Potentials=NULL;

static Real grav_pot(const Real x1, const Real x2, const Real x3)
{
    Real i, j, k;
    ijk_pos(&level0_Grid, x1, x2, x3, &i, &j, &k);
    
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
#ifdef MPI_PARALLEL
/* Get my task id (rank in MPI) */
  if(MPI_SUCCESS != MPI_Comm_rank(MPI_COMM_WORLD,&(level0_Grid.my_id)))
    ath_error("Error on calling MPI_Comm_rank\n");

/* Get the number of processes */
  if(MPI_SUCCESS != MPI_Comm_size(MPI_COMM_WORLD,&(level0_Grid.nproc)))
    ath_error("Error on calling MPI_Comm_size\n");
#else
  level0_Grid.my_id = 0;
  level0_Grid.nproc = 1;
#endif

  par_open("/dev/null"); /* to trick athena into thinking it has opened a parameter file, will not work on windows */
  is_restart = 0;
  show_config_par();   /* Add the configure block to the parameter database */
  
  
  return 0;
}

int get_timestep(double * value) {
    *value = level0_Grid.dt;
    return 0;
}

int set_timestep(double value) {
    level0_Grid.dt = value;
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
  
  if(has_external_gravitational_potential) {
    StaticGravPot = grav_pot;
  }
  
  CourNo = par_getd("time","cour_no");

#ifdef ISOTHERMAL
  Iso_csound = par_getd("problem","iso_csound");
  Iso_csound2 = Iso_csound*Iso_csound;
#else
  Gamma = par_getd("problem","gamma");
  Gamma_1 = Gamma - 1.0;
  Gamma_2 = Gamma - 2.0;
#endif
  
  
  init_domain(&level0_Grid, &level0_Domain);
  
  init_grid(&level0_Grid, &level0_Domain);
  
  
  if ((Potentials = (Real***)calloc_3d_array(
    level0_Grid.Nx3 + 2 * nghost, 
    level0_Grid.Nx2 + 2 * nghost,
    level0_Grid.Nx1 + 2 * nghost,  
    sizeof(Real))) == NULL)
  {
    return -1;
  }

  
  return 0;
}

int initialize_grid()
{
  
  set_bvals_init(&level0_Grid, &level0_Domain);
  
  set_bvals(&level0_Grid, 0);
  
  if(!is_dt_set_by_script) {
    new_dt(&level0_Grid);
  }
  lr_states_init(level0_Grid.Nx1,level0_Grid.Nx2,level0_Grid.Nx3);
  
  Integrate = integrate_init(level0_Grid.Nx1,level0_Grid.Nx2,level0_Grid.Nx3);

#ifdef SELF_GRAVITY
  SelfGrav = selfg_init(&level0_Grid, &level0_Domain);
  if(is_restart == 0) {
      (*SelfGrav)(&level0_Grid, &level0_Domain);
  }
  /* Only bvals for Phi set when last argument of set_bvals = 1  */
  set_bvals(&level0_Grid, 1);
#endif
    
  return 0;
}

int get_position_of_index(int i, int j, int k, double * x, double * y, 
  double * z){

  if (level0_Domain.grid_block == NULL) {
    return -1;
  }
  
  if (i < (level0_Domain.ixs - nghost)  || i > (level0_Domain.ixe + nghost)) {
    return -1;
  }
  if (j < (level0_Domain.jxs - nghost)  || j > (level0_Domain.jxe + nghost)) {
    return -1;
  }
  if (k < (level0_Domain.kxs - nghost) || k > (level0_Domain.kxe + nghost)) {
    return -1;
  }
  
  
  cc_pos(
    &level0_Grid,
    i + level0_Grid.is,
    j + level0_Grid.js,
    k + level0_Grid.ks,
    x,y,z
  );
  
  return 0;
}

int get_interpolated_gravitational_potential(double x, double y, double z, double * potential) {
    
    if (level0_Domain.grid_block == NULL) {
        return -1;
    }

    *potential = grav_pot(x,y,z);
    
    return 0;
}

int get_index_of_position(double x, double y, 
  double z, double *i , double * j, double * k){

  if (level0_Domain.grid_block == NULL) {
    return -1;
  }
  
  ijk_pos(&level0_Grid,x,y,z,i,j,k);
  *i -= level0_Grid.is;
  *j -= level0_Grid.js;
  *k -= level0_Grid.ks;
  return 0;
}

int test() {
    return 1;
}

int fill_grid_state(
  int i, int j, int k, 
  double rho, double rhovx, double rhovy, double rhovz, 
  double en){
      
  if (level0_Domain.grid_block == NULL) {
    return -1;
  }
  if (i < level0_Domain.ixs  || i > level0_Domain.ixe) {
    return -1;
  }
  if (j < level0_Domain.jxs  || j > level0_Domain.jxe) {
    return -1;
  }
  if (k < level0_Domain.kxs  || k > level0_Domain.kxe) {
    return -1;
  }
  
  
  
  level0_Grid.U[k][j][i].d = rho;
  level0_Grid.U[k][j][i].M1 = rhovx;
  level0_Grid.U[k][j][i].M2 = rhovy;
  level0_Grid.U[k][j][i].M3 = rhovz;
  level0_Grid.U[k][j][i].E = en;
  
  return 0;
}

int get_grid_state_mpi(
    int * i, int * j, int * k, 
    double * rho, 
    double * rhovx, double * rhovy, double * rhovz, 
    double * en, 
    int number_of_points)
{
    int imin = level0_Grid.is + level0_Grid.idisp;
    int imax = level0_Grid.ie + level0_Grid.idisp;
    int jmin = level0_Grid.js + level0_Grid.jdisp;
    int jmax = level0_Grid.je + level0_Grid.jdisp;
    int kmin = level0_Grid.ks + level0_Grid.kdisp;
    int kmax = level0_Grid.ke + level0_Grid.kdisp;
    int l=0;
    int i0,j0,k0 = 0;
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        if ( 
            (i0 >= imin && i0 <= imax) &&
            (j0 >= jmin && j0 <= jmax) &&
            (k0 >= kmin && k0 <= kmax)
        ) {
            i0 -= level0_Grid.idisp;
            j0 -= level0_Grid.jdisp;
            k0 -= level0_Grid.kdisp;
            rho[l] = level0_Grid.U[k0][j0][i0].d;
            rhovx[l] = level0_Grid.U[k0][j0][i0].M1;
            rhovy[l] = level0_Grid.U[k0][j0][i0].M2;
            rhovz[l] = level0_Grid.U[k0][j0][i0].M3;
            en[l] = level0_Grid.U[k0][j0][i0].E;
        } else {
            rho[l] = rhovx[l] = rhovy[l]= rhovz[l] = en[l] = 0.0;
        } 
    }
    
#ifdef MPI_PARALLEL
    if(level0_Grid.my_id) {
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
    int * i, int * j, int * k, 
    double * rho, 
    double * rhovx, double * rhovy, double * rhovz, 
    double * en, 
    int number_of_points)
{
    int imin = level0_Grid.is + level0_Grid.idisp;
    int imax = level0_Grid.ie + level0_Grid.idisp;
    int jmin = level0_Grid.js + level0_Grid.jdisp;
    int jmax = level0_Grid.je + level0_Grid.jdisp;
    int kmin = level0_Grid.ks + level0_Grid.kdisp;
    int kmax = level0_Grid.ke + level0_Grid.kdisp;
    int l=0;
    int i0,j0,k0 = 0;
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        if ( 
            (i0 >= imin && i0 <= imax) &&
            (j0 >= jmin && j0 <= jmax) &&
            (k0 >= kmin && k0 <= kmax)
        ) {
            i0 -= level0_Grid.idisp;
            j0 -= level0_Grid.jdisp;
            k0 -= level0_Grid.kdisp;
            level0_Grid.U[k0][j0][i0].d = rho[l];
            level0_Grid.U[k0][j0][i0].M1 = rhovx[l];
            level0_Grid.U[k0][j0][i0].M2 = rhovy[l];
            level0_Grid.U[k0][j0][i0].M3 = rhovz[l];
            level0_Grid.U[k0][j0][i0].E = en[l];
        }
    }
    
    return 0;
}




int get_potential(
    int * i, int * j, int * k, 
    double * potential,
    int number_of_points)
{
      
    if (level0_Domain.grid_block == NULL) {
        return -1;
    }
    if (Potentials == NULL) {
        return -1;
    }

    int imin = level0_Grid.is + level0_Grid.idisp - 1;
    int imax = level0_Grid.ie + level0_Grid.idisp + 1;
    int jmin = level0_Grid.js + level0_Grid.jdisp;
    int jmax = level0_Grid.je + level0_Grid.jdisp;
    if(jmin != jmax) {
        jmin--;
        jmax++;
    }
    int kmin = level0_Grid.ks + level0_Grid.kdisp;
    int kmax = level0_Grid.ke + level0_Grid.kdisp;
    if(kmin != kmax) {
        kmin--;
        kmax++;
    }
    int l=0;
    int i0,j0,k0 = 0;
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        if ( 
            (i0 >= imin && i0 <= imax) &&
            (j0 >= jmin && j0 <= jmax) &&
            (k0 >= kmin && k0 <= kmax)
        ) {
            i0 -= level0_Grid.idisp;
            j0 -= level0_Grid.jdisp;
            k0 -= level0_Grid.kdisp;
            
            potential[l] = Potentials[k0][j0][i0];
        } else {
            potential[l] = 0.0;
        } 
    }
    
#ifdef MPI_PARALLEL
    if(level0_Grid.my_id) {
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
    if (level0_Domain.grid_block == NULL) {
        return -1;
    }
    if (Potentials == NULL) {
        return -1;
    }

    int imin = level0_Grid.is + level0_Grid.idisp - 1;
    int imax = level0_Grid.ie + level0_Grid.idisp + 1;
    int jmin = level0_Grid.js + level0_Grid.jdisp;
    int jmax = level0_Grid.je + level0_Grid.jdisp;
    if(jmin != jmax) {
        jmin--;
        jmax++;
    }
    int kmin = level0_Grid.ks + level0_Grid.kdisp - 1;
    int kmax = level0_Grid.ke + level0_Grid.kdisp + 1;
    if(kmin != kmax) {
        kmin--;
        kmax++;
    }
    int l=0;
    int i0,j0,k0 = 0;
    
    for(l=0; l < number_of_points; l++) {
        i0 = i[l];
        j0 = j[l];
        k0 = k[l];
        
        if ( 
            (i0 >= imin && i0 <= imax) &&
            (j0 >= jmin && j0 <= jmax) &&
            (k0 >= kmin && k0 <= kmax)
        ) {
            i0 -= level0_Grid.idisp;
            j0 -= level0_Grid.jdisp;
            k0 -= level0_Grid.kdisp;
            Potentials[k0][j0][i0] = potential[l];
        }
    }
    
    return 0;
}


int get_grid_state(
  int i, int j, int k, 
  double * rho, double * rhovx, double * rhovy, double * rhovz, 
  double * en){
      
  if (level0_Domain.grid_block == NULL) {
    return -1;
  }
  if (i < level0_Domain.ixs  || i > level0_Domain.ixe) {
    return -1;
  }
  if (j < level0_Domain.jxs  || j > level0_Domain.jxe) {
    return -1;
  }
  if (k < level0_Domain.kxs  || k > level0_Domain.kxe) {
    return -1;
  }
  
  *rho = level0_Grid.U[k][j][i].d;
  *rhovx = level0_Grid.U[k][j][i].M1;
  *rhovy = level0_Grid.U[k][j][i].M2;
  *rhovz = level0_Grid.U[k][j][i].M3;
  *en = level0_Grid.U[k][j][i].E;
  
  return 0;
}

int esys_roe_adb_hydro(int * index, double * u, double * v, double * w, 
  double * h, double * ev, double * rem0, double * rem1, double * rem2, 
  double * rem3, double * rem4, double * lem0, double * lem1, 
  double * lem2, double * lem3, double * lem4){
  return 0;
}


static Gas ***Soln=NULL;

#include <math.h>
    
void _fill_grid_linearwave_1d(Grid *pGrid, Domain *pDomain, int wave_flag, double amp, double vflow, int  wave_dir){
  int i=0,j=0,k=0;
  int is,ie,js,je,ks,ke,n,m,nx1,nx2,nx3;
  Real d0,p0,u0,v0,w0,h0;
  Real x1,x2,x3,r,ev[NWAVE],rem[NWAVE][NWAVE],lem[NWAVE][NWAVE];
  Real x1max,x2max,lambda;
#ifdef MHD
  Real bx0,by0,bz0,xfact,yfact;
#endif /* MHD */
  is = pGrid->is; ie = pGrid->ie;
  js = pGrid->js; je = pGrid->je;
  ks = pGrid->ks; ke = pGrid->ke;
  nx1 = (ie-is)+1 + 2*nghost;
  nx2 = (je-js)+1 + 2*nghost;
  nx3 = (ke-ks)+1 + 2*nghost;


  if ((Soln = (Gas***)calloc_3d_array(nx3,nx2,nx1,sizeof(Gas))) == NULL)
    ath_error("[linear_wave1d]: Error allocating memory\n");
    
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
  if (pGrid->Nx1 > 1) {
    for (k=ks; k<=ke; k++) {
      for (j=js; j<=je; j++) {
        pGrid->B1i[k][j][ie+1] = pGrid->B1i[k][j][ie];
      }
    }
  }
  if (pGrid->Nx2 > 1) {
    for (k=ks; k<=ke; k++) {
      for (i=is; i<=ie; i++) {
        pGrid->B2i[k][je+1][i] = pGrid->B2i[k][je][i];
      }
    }
  }
  if (pGrid->Nx3 > 1) {
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

  return;
}


int fill_grid_linearwave_1d(int wave_flag, double amp, double vflow, int wave_dir)
{
    
    if (level0_Domain.grid_block == NULL) {
        return -1;
    } 
    if (level0_Grid.U == NULL) {
        return -1;
    } 
    _fill_grid_linearwave_1d(&level0_Grid, &level0_Domain, wave_flag, amp, vflow, wave_dir);
    return 0;
}

int get_time(double * value){
    *value = level0_Grid.time;
}

int evolve(double tlim) {
    
    while (level0_Grid.time < tlim) {
        if ((tlim-level0_Grid.time) < level0_Grid.dt) {
            level0_Grid.dt = (tlim-level0_Grid.time);
        }
        
        (*Integrate)(&level0_Grid);
        
        Userwork_in_loop(&level0_Grid, &level0_Domain);
        
        
#ifdef SELF_GRAVITY
        (*SelfGrav)(&level0_Grid, &level0_Domain);
        /* Only bvals for Phi set when last argument of set_bvals = 1  */
        set_bvals(&level0_Grid, 1);
        selfg_flux_correction(&level0_Grid);
#endif
        
        
        level0_Grid.nstep++;
        level0_Grid.time += level0_Grid.dt;
        new_dt(&level0_Grid);
        set_bvals(&level0_Grid, 0);
    }
    
}

