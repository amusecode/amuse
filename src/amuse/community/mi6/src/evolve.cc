#include"evolve.h"

extern double Tcal_grav0,  Tcal_grav1;
extern double Tcal_comm0,  Tcal_comm1;
extern double Tcal_fix0, Tcal_fix1;


extern int  EX_FLAG;

static int NSTEP = 0;

void set_NSTEP(const int &_NSTEP){
  NSTEP = _NSTEP;
} 

int get_NSTEP(){
  return NSTEP;
} 


static int NBH = 0;
static int NJP = 0;
void set_NBH_NJP(const int &_NBH, 
		 const int &_NJP){
  NBH = _NBH;
  NJP = _NJP;
}

void set_NJP(const int &_NJP){
  NJP = _NJP;
}

static double ETA_S = 0;
static double ETA_FS = 0;
static double ETA_SMBH = 0;
static double ETA_IMBH = 0;
void set_eta(const double &eta_s,
	     const double &eta_fs,
	     const double &eta_smbh,
	     const double &eta_imbh){
  ETA_S = eta_s;
  ETA_FS = eta_fs;
  ETA_SMBH = eta_smbh;
  ETA_IMBH = eta_imbh;
}



void setj_to_sapporo(Particle prt[],
		     int address[],
		     const int &Nip){
  static double eps2 = 0.0;
  for(int i=0; i<Nip; i++){
    int addi = address[i];
    Particle *prti = prt+addi;
    if(0 <= prti->address){
      set_j_particle(prti->address, (double*)&prti->pos,
		     (double*)&prti->vel,  (double*)&prti->acc,
		     (double*)&prti->jrk,  (double*)&prti->acc2,
		     (double*)&prti->acc3,  prti->mass,
		     prti->time,    prti->index,
		     eps2);
    }
  }
}

static const int NPRTMAX = 1048576+100;

static int ID_IN[NPRTMAX];
static double MASS_IN[NPRTMAX];
static double EPS2_IN[NPRTMAX];
static double POS_IN[NPRTMAX][3];
static double VEL_IN[NPRTMAX][3];
static double ACC_IN[NPRTMAX][3];

static double ACC_OUT[NPRTMAX][3];
static double JRK_OUT[NPRTMAX][3];
static double SNP_OUT[NPRTMAX][3];
static double CRK_OUT[NPRTMAX][3];
static double PHI_OUT[NPRTMAX];
static int NNB_OUT[NPRTMAX];
static double NNB_R2_OUT[NPRTMAX];

static double EPS2_FS_FS;
static double EPS2_BH_BH;
static double EPS2_BH_FS;
static double EPS2_FS_BH;

#ifdef DEBUG
static int ID_PRE[NPRTMAX];
static double MASS_PRE[NPRTMAX];
static double EPS2_PRE[NPRTMAX];
static double POS_PRE[NPRTMAX][3];
static double VEL_PRE[NPRTMAX][3];
static double ACC_PRE[NPRTMAX][3];
#endif //DEBUG

void set_eps2(const double &eps2_fs_fs,
	      const double &eps2_bh_fs,
	      const double &eps2_bh_bh){
  EPS2_FS_FS = eps2_fs_fs;
  EPS2_BH_FS = eps2_bh_fs;
  EPS2_BH_BH = eps2_bh_bh;
  EPS2_FS_BH = eps2_bh_fs;
}

void calc_force_using_sapporo(Particle prt[],
			      int address[],
			      const int &Nip,
			      const double &Tsys,
			      const int &mode){

  if(mode == 0){
    no_predict_all(Tsys, NJP);
    for(int i=0; i<Nip; i++){
      int addi = address[i];
      for(int k=0; k<3; k++){
	POS_IN[i][k] = prt[addi].pos[k];
	VEL_IN[i][k] = prt[addi].vel[k];
	ACC_IN[i][k] = prt[addi].acc[k];
      }
    }
  }
  else if(mode == 1){
    predict_all(Tsys, NJP);
    for(int i=0; i<Nip; i++){
      int addi = address[i];
      for(int k=0; k<3; k++){
	POS_IN[i][k] = prt[addi].pos_pre[k];
	VEL_IN[i][k] = prt[addi].vel_pre[k];
	ACC_IN[i][k] = prt[addi].acc_pre[k];
      }
    }
  }

  for(int i=0; i<Nip; i++){
    int addi = address[i];
    ID_IN[i]=prt[addi].index;
    MASS_IN[i]=prt[addi].mass;
    if(prt[addi].index < NBH){
      EPS2_IN[i] = EPS2_BH_FS;
    }
    else{
      EPS2_IN[i] = EPS2_FS_FS;
    }
    for(int k=0; k<3; k++){
      ACC_OUT[i][k] = 0.0;
      JRK_OUT[i][k] = 0.0;
      SNP_OUT[i][k] = 0.0;
    }
    PHI_OUT[i] = 0.0;
    NNB_OUT[i] = -1;
    NNB_R2_OUT[i] = 0.0;
  }


  calc_force_on_predictors(Nip, NJP,
			   ID_IN, POS_IN, 
			   VEL_IN, ACC_IN, 
			   MASS_IN, EPS2_IN,
			   ACC_OUT, JRK_OUT,
			   SNP_OUT, CRK_OUT,PHI_OUT,
			   NNB_OUT, NNB_R2_OUT);
  /*
  cout<<"Tsys="<<Tsys<<endl;
  for(int i=0; i<NJP; i++){
    pick_up_predictor_2(i, ID_PRE[i], MASS_PRE[i], EPS2_PRE[i], POS_PRE[i], VEL_PRE[i], ACC_PRE[i]);
    cout<<"ID_PRE[i]="<<ID_PRE[i]<<endl;
    cout<<"MASS_PRE[i]="<<MASS_PRE[i]<<endl;
    cout<<"POS_PRE[i]="<<POS_PRE[i][0]<<"   "<<POS_PRE[i][1]<<"   "<<POS_PRE[i][2]<<endl;
    cout<<"VEL_PRE[i]="<<VEL_PRE[i][0]<<"   "<<VEL_PRE[i][1]<<"   "<<VEL_PRE[i][2]<<endl;
    cout<<"ACC_PRE[i]="<<ACC_PRE[i][0]<<"   "<<ACC_PRE[i][1]<<"   "<<ACC_PRE[i][2]<<endl;
    cout<<endl;
  }
  */
  for(int i=0; i<Nip; i++){
    int addi = address[i];
    prt[addi].acc = ACC_OUT[i];
    prt[addi].jrk = JRK_OUT[i];
    prt[addi].acc2 = SNP_OUT[i];
#ifdef SAP
    prt[addi].phi = -PHI_OUT[i];
#else
    prt[addi].phi = PHI_OUT[i];
#endif
    prt[addi].ngb_index = NNB_OUT[i];
    if(prt[addi].index < NBH){
      prt[addi].r2min = NNB_R2_OUT[i] - EPS2_BH_FS;
    }
    else{
      prt[addi].r2min = NNB_R2_OUT[i] - EPS2_FS_FS;
    }
    /*
#ifdef DEBUG
    cout<<"prt[addi].index="<<prt[addi].index<<endl;
    cout<<"prt[addi].acc_new="<<prt[addi].acc<<endl;
    cout<<"prt[addi].jrk_new="<<prt[addi].jrk<<endl;
    cout<<"prt[addi].snp_new="<<prt[addi].acc2<<endl;
    cout<<"prt[addi].phi_new="<<prt[addi].phi<<endl;
    cout<<endl;
#endif //DEBUG
    */
  }
}

void evolve_initialize(Particle prt[],
		       int address[],
		       const int &Ntot,
		       const int &NBH,
		       const int &Njp,
		       const double &Tsys){
  static int first = 1;
  if(first){
    allocate_mpi_buffer(Ntot);
    initialize();
    first = 0;
  }
  set_NBH_NJP(NBH, Njp);
  setj_to_sapporo(prt, address, Ntot);
  for(int i=0; i<Ntot; i++){
    prt[address[i]].copyold();
    prt[address[i]].clear();
  }
  calc_force_using_sapporo(prt, address, Ntot, Tsys, 0);


  reduce_force(prt, address, 0, Ntot, MPI_COMM_WORLD);

  if(EX_FLAG == 1){
    calc_force_from_point_mass_to_array(prt, address, Ntot, 0);
  }
  for(int i=0; i<Ntot; i++){
    prt[i].accumulate_force();
  }
}

static double next_time[NPRTMAX];

// merge: return 1 
int evolve_onestep(Particle prt[],
		   int address[],
		   int &Nip,
		   const int &Ntot,
		   const int &NBH,
		   double &Tsys,
		   const double &Tmerge,
		   const double &dt_max,
		   const int &first_loop,
		   double &Egr,
		   const int &itr){

  int merge = 0;
  int ret_val = 0;
  double eta = 0.0;
  //cerr<<"first_loop="<<first_loop<<endl;
  for(int i=0; i<Nip; i++){
    int addi = address[i];
    //if(prt[addi].address != -1){
    if(prt[addi].mass != 0.0){
      if(first_loop){
	eta = ETA_S;
	//prt[addi].set_dt_2nd(dt_max, eta, Tmerge);
	prt[addi].set_dt_2nd(dt_max*pow(0.5, itr), eta, Tmerge);
      }
      else{
	if(NBH <= prt[addi].index){
	  eta = ETA_FS;
	}
	else{
	  eta = ETA_SMBH;
	}
	//prt[addi].set_dt_4th(dt_max, eta, Tmerge);
	//prt[addi].set_dt_6th(dt_max, eta, Tmerge);
	prt[addi].set_dt_6th(dt_max*pow(0.5, itr), eta, Tmerge);
      }
    }
    next_time[addi] = prt[addi].time + prt[addi].dtime;
  }

  /*
  if(Tsys == 0.0){
    for(int i=0; i<Ntot; i++){
      cout<<"address[i]="<<address[i]<<endl;
      prt[i].dump();
    }
  }
  */

  /*
  if(first_loop){
    for(int i=0; i<Ntot; i++){
      next_time[i] = prt[i].time + prt[i].dtime;
    }
    Qsort_index(next_time, address, 0, Ntot-1);
  }
  */

  /*
  for(int i=0; i<Ntot; i++){
    next_time[i] = prt[i].time + prt[i].dtime;
  }
  Qsort_index(next_time, address, 0, Ntot-1);
  */


  sort_select_iparticle(prt,  Nip,  Ntot,  address,  next_time);
  Tsys = next_time[address[0]];


  /*
  cout<<"Tsys="<<Tsys<<endl;
  cout<<"Ntot="<<Ntot<<endl;
  cout<<"Nip="<<Nip<<endl;
  for(int i=0; i<Nip; i++){
    cout<<"next_time[address[i]]="<<next_time[address[i]]<<endl;
    prt[address[i]].dump();
  }
  cout<<endl;
  cout<<endl;
  */

  NSTEP += Nip;
  for(int i=0; i<Nip; i++){
    prt[address[i]].predict(Tsys);
  }
  for(int i=0; i<NBH; i++){
    prt[i].predict(Tsys);
  }


  merge = merge_check(prt,  address,  Nip,  NBH);
  if( merge ){
    // if particles will be merge
    // all particles are predicted
    //Tmerge = Tsys;
    Nip = 0;
    for(int i=0; i<Ntot; i++){
      if(prt[i].mass != 0.0){
	prt[i].dtime = Tsys - prt[i].time;
	next_time[i] = Tsys;
	Nip++;
      }
    }
    for(int i=0; i<Nip; i++){
      prt[address[i]].predict(Tsys);
    }
  }


  for(int i=0; i<Nip; i++){
    prt[address[i]].copyold();
    prt[address[i]].clear();
  }
  Tcal_grav0 = MPI_Wtime();
  calc_force_using_sapporo(prt, address, Nip, Tsys, 1);
  Tcal_grav1 += MPI_Wtime() - Tcal_grav0;

  Tcal_comm0 = MPI_Wtime();
  reduce_force(prt,  address,  0,  Nip,  MPI_COMM_WORLD);
  Tcal_comm1 += MPI_Wtime() - Tcal_comm0;

  Tcal_fix0 = MPI_Wtime();
  if(EX_FLAG == 1){
    calc_force_from_point_mass_to_array(prt, address, Nip, 1);
  }
  Tcal_fix1 += MPI_Wtime() - Tcal_fix0;

#ifndef PECEC
  for(int i=0; i<Nip; i++){
    int addi = address[i];
    prt[addi].accumulate_force();
    prt[addi].correct_mix();
    prt[addi].interpolate();
    //prt[addi].correct();
    prt[addi].time = Tsys; 
    Egr += prt[addi].calc_Egr_6th();
  }
#else
  for(int i=0; i<Nip; i++){
    int addi = address[i];
    prt[addi].accumulate_force();
    prt[addi].correct_mix();

    prt[addi].pos_pre = prt[addi].pos;
    prt[addi].vel_pre = prt[addi].vel;
    prt[addi].acc_pre = prt[addi].acc;
    prt[addi].jrk_pre = prt[addi].jrk;

    prt[addi].acc -= (prt[addi].acc_ex + prt[addi].acc_PN_ex);
    prt[addi].jrk -= (prt[addi].jrk_ex + prt[addi].jrk_PN_ex);
    prt[addi].acc2 -= (prt[addi].acc2_ex + prt[addi].acc2_PN_ex);
    prt[addi].acc3 -= prt[addi].acc3_ex;
    prt[addi].phi -= prt[addi].phi_ex;

    prt[addi].acc_ex = prt[addi].acc_PN_ex = 0.0;
    prt[addi].jrk_ex = prt[addi].jrk_PN_ex = 0.0;
    prt[addi].acc2_ex = prt[addi].acc2_PN_ex = 0.0;
    prt[addi].acc3_ex = 0.0;
    prt[addi].phi_ex = 0.0;
  }
  Tcal_fix0 = MPI_Wtime();
  if(EX_FLAG == 1){
    calc_force_from_point_mass_to_array(prt, address, Nip, 1);
  }
  Tcal_fix1 += MPI_Wtime() - Tcal_fix0;
  for(int i=0; i<Nip; i++){
    int addi = address[i];
    prt[addi].accumulate_force();
    prt[addi].correct_mix();
    prt[addi].interpolate();
    prt[addi].time = Tsys; 
    Egr += prt[addi].calc_Egr_6th();
  }
#endif //PECEC


  setj_to_sapporo(prt, address, Nip);

  if(merge){
    ret_val = 1;
  }

  return ret_val;
}


void sort_time_all(Particle prt[], 
		   int address[], 
		   int Ntot){
  for(int i=0; i<Ntot; i++){
    next_time[i] = prt[i].time + prt[i].dtime;
  }
  Qsort_index(next_time, address, 0, Ntot-1);
}
