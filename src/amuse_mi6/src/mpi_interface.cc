
#include"mpi_interface.h"

static const int NMAXPROC = 1024;

static int *int_list_send0;
static double *dbl_list_send0;
static Vector3 *vec_list_send0;
static Vector3 *vec_list_send1;
static Vector3 *vec_list_send2;
static int *int_list_recv0;
static double *dbl_list_recv0;
static Vector3 *vec_list_recv0;
static Vector3 *vec_list_recv1;
static Vector3 *vec_list_recv2;

struct dblint{
  double val;
  int rank;
};
static dblint *dblint_list_send0;
static dblint *dblint_list_recv0;

void allocate_mpi_buffer(const int &N){
  int_list_send0 = new int[N];
  dbl_list_send0 = new double[N];
  vec_list_send0 = new Vector3[N];
  vec_list_send1 = new Vector3[N];
  vec_list_send2 = new Vector3[N];
  dblint_list_send0 = new dblint[N];

  int_list_recv0 = new int[N];
  dbl_list_recv0 = new double[N];
  vec_list_recv0 = new Vector3[N];
  vec_list_recv1 = new Vector3[N];
  vec_list_recv2 = new Vector3[N];
  dblint_list_recv0 = new dblint[N];
}

void divide_proc(int Narray[], 
		 int Ndisp[], 
		 const int &Nprt,
		 const int &Nproc){
  for(int i=0; i<Nproc; i++){
    Narray[i] = Nprt / Nproc;
    if(i < Nprt%Nproc){Narray[i]++;}
  }
  for(int i=0; i<Nproc+1; i++){  Ndisp[i] = 0;  }
  for(int i=1; i<Nproc+1; i++){
    Ndisp[i] += Ndisp[i-1] + Narray[i-1];
  }
}


void divide_proc(const int &Nprt,
		 int &Njp_loc,
		 int &first_address){
  int Nproc = 0;
  int myrank = 0;
  first_address = 0;
#ifndef NOMPI
  MPI_Comm_size(MPI_COMM_WORLD, &Nproc);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
    Nproc = 1;
    myrank = 0;
#endif
  Njp_loc  = Nprt / Nproc;
  if(myrank < Nprt%Nproc){
    Njp_loc++;
  }
  first_address = (int)(Nprt/Nproc)*myrank;
  if(myrank < Nprt % Nproc){
    first_address += myrank;
  }

}

void reduce_force(Particle prt[],
		  int address[],
		  const int &Ni0,
		  const int &Ni,
		  const MPI_Comm &icomm){
#ifndef NOMPI
  int myrank;
  MPI_Comm_rank(icomm, &myrank);
  for(int i=Ni0; i<Ni0+Ni; i++){
    int addi = address[i];
    vec_list_send0[i-Ni0] = prt[addi].acc;
    vec_list_send1[i-Ni0] = prt[addi].jrk;
    vec_list_send2[i-Ni0] = prt[addi].acc2;
    dbl_list_send0[i-Ni0] = prt[addi].phi;
    dblint_list_send0[i-Ni0].rank = myrank;
    dblint_list_send0[i-Ni0].val =  prt[addi].r2min;
  }

  MPI_Allreduce(vec_list_send0, vec_list_recv0,
		3*Ni, MPI_DOUBLE, MPI_SUM, icomm);
  MPI_Allreduce(vec_list_send1, vec_list_recv1,
		3*Ni, MPI_DOUBLE, MPI_SUM, icomm);
  MPI_Allreduce(vec_list_send2, vec_list_recv2,
		3*Ni, MPI_DOUBLE, MPI_SUM, icomm);
  MPI_Allreduce(dbl_list_send0, dbl_list_recv0,
		Ni, MPI_DOUBLE, MPI_SUM, icomm);


  MPI_Allreduce(dblint_list_send0,  dblint_list_recv0, 
		Ni,  MPI_DOUBLE_INT,  MPI_MINLOC,  icomm);

  for(int i=Ni0; i<Ni0+Ni; i++){
    int addi = address[i];
    prt[addi].acc = vec_list_recv0[i-Ni0];
    prt[addi].jrk = vec_list_recv1[i-Ni0];
    prt[addi].acc2 = vec_list_recv2[i-Ni0];
    prt[addi].phi = dbl_list_recv0[i-Ni0];
    prt[addi].r2min = dblint_list_recv0[i-Ni0].val;
    MPI_Bcast(&(prt[addi].ngb_index),  1,  MPI_INT,  dblint_list_recv0[i-Ni0].rank,  icomm);
  }
#endif
}

void sum_double(double &x, 
		double &xsum,
		const int &nwords,
		const MPI_Comm &comm){
#ifndef NOMPI
  MPI_Allreduce(&x, &xsum, nwords, MPI_DOUBLE, MPI_SUM, comm); 
#else
    double * xsum_p = &xsum;
    double * x_p = &x;
    for(int i = 0; i < nwords; i++) {
        *(xsum_p + i) = *(x_p + i);
    }
#endif
}

void mpi_abort(){
  int error = 0;
#ifndef NOMPI
  MPI_Abort(MPI_COMM_WORLD, error);
#else
    exit(1);
#endif
}


