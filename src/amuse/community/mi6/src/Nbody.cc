/*
J-paralel
*/
#include<iostream>
#include<fstream>
#include<cmath>
#include<unistd.h>
#include"6thorder.h"
#include"IO.h"
#include"mpi_interface.h"
#include"evolve.h"
#include"energy.h"
#include"merge.h"
#include"external_field.h"

using namespace std;

double Tcal_grav0 = 0.0;
double Tcal_grav1 = 0.0;
double Tcal_comm0 = 0.0;
double Tcal_comm1 = 0.0;
double Tcal_fix0 = 0.0;
double Tcal_fix1 = 0.0;
double Tcal_tot0 = 0.0;
double Tcal_tot1 = 0.0;
double Tcal_all = 0.0;

int EX_FLAG = 1;
//int EX_FLAG = 0;


void iteration(Particle prt[],
	       Particle prt_old[],
	       int address[],
	       int address_old[],
	       const int &Ntot,

	       const int &Ndead_old,
	       const int &Nmerge_old,
	       const int &Naccrete_old,

	       const int &Nip_tot_old,
	       const int &Njp_old,

	       const double &E1_old,
	       const double &Ep1_old,
	       const double &Ek1_old,
	       const double &Egr_old,
	       const double &Emerge_old,

	       const double &Tsys_old,
	       const double &Tmerge_old,

	       const int &Nloop_old,
	       const int &Nstep_old,
	       const int &first_loop_old,
	       const double &eta_s,
	       const double &eta_fs,
	       const double &eta_smbh,
	       const double &eta_imbh,

	       int &Ndead,
	       int &Nmerge,
	       int &Naccrete,

	       int &Nip_tot,


	       int &Njp,

	       double &E1,
	       double &Ep1,
	       double &Ek1,
	       double &Egr,
	       double &Emerge,

	       double &Tsys,
	       double &Tmerge,

	       int &Nloop,
	       int &Nstep,
	       int &first_loop,
	       int &itr,
	       int &Nip){
  itr++;
  cerr<<"energy error is too large: iteration="<<itr<<endl;
  for(int i=0; i<Ntot; i++){
    prt[i] = prt_old[i];
    address[i] = address_old[i];
  }

  Ndead = Ndead_old;
  Nip_tot = Nip_tot_old;
  Njp = Njp_old;

  set_NJP(Njp);
  setj_to_sapporo(prt, address, Ntot);
  //setj_to_sapporo(prt, address, Njp);

  E1 = E1_old;
  Ek1 = Ek1_old;
  Ep1 = Ep1_old;
  Egr = Egr_old;
  Emerge = Emerge_old;

  Tsys = Tsys_old;
  Tmerge = Tmerge_old;

  Nloop = Nloop_old;
  Nstep = Nstep_old;
  first_loop = first_loop_old;
  //set_eta(eta_s*pow(0.5, itr), eta_fs*pow(0.5, itr), eta_smbh*pow(0.5, itr), eta_imbh*pow(0.5, itr));
  Tcal_grav1 = Tcal_comm1 = Tcal_fix1 = 0.0;
  Tcal_tot0 = MPI_Wtime();
  set_NSTEP(Nstep);

  Nmerge = Nmerge_old;
  Naccrete = Naccrete_old;
  set_Ndead(Nmerge, Naccrete);
  Nip = Nip_tot;

  copy_SMBH_OLD_TO_NEW();


}

int main(int argc, char *argv[]){

  set_speed_of_light(2285.604);
  set_SMBH(1.0, (Vector3)0.0,  (Vector3)0.0); 

  
  //const double dEcrit = 1e-10;
  const double dEcrit = 5e-5;
  //const double dEcrit = 1e30;

  int Nstep = 0;
  int Nstep_old = 0;
  int Nloop = 0;
  int Nloop_old = 0;

  int myrank = 0;
  int Nproc = 0;
  int Ntot = 0;
  int NFS = 0;
  int NBH = 0;
  int NSMBH = 0;
  int NIMBH = 0;
  int Nmerge = 0;
  int Nmerge_loop = 0;
  int Naccrete = 0;
  int Naccrete_loop = 0;
  double eps2_fs_fs = 0.0;
  double eps2_bh_bh = 0.0;
  double eps2_fs_smbh = 0.0;
  double eps2_fs_imbh = 0.0;
  double eta_s = 0.0;
  double eta_fs = 0.0;
  double eta_smbh = 0.0;
  double eta_imbh = 0.0;
  double Tsys = 0.0;
  double Tmerge = 0.0;
  double Tend = 0.0;
  double Egr = 0.0;
  double Emerge = 0.0;




  int Ndead = 0; // Nmerge + Naccrete (every node have same number)
  int Nip_tot = 0; // Ntot - Ndead     (every node have same number)
  int Nip = 0; // Nip  (every node have same number)
  int Njp_org = 0; // Njp + Ndead 
  int Njp = 0; // Njp (0 < Njp< Njp_org)


  int first_address = 0;

  cout<<setprecision(15);
  cerr<<setprecision(15);

  int READ_FLAG = 0; // 0:read nemo_ascii, 1:read snap shot

  MPI_Init(&argc,&argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Comm_rank (comm, &myrank); 
  MPI_Comm_size (comm, &Nproc); 

  int snpid = 0;
  char sinput[1024];
  char dirname[1024];
  char smerge[1024]; 
  char slog[1024];
  char sBH[1024];
  ofstream fout_merge; 
  ofstream fout_log;
  ofstream fout_BH;
  Particle *prt = NULL;
  Particle *prt_old = NULL;

  int c;
  if(myrank==0){
    while((c=getopt(argc,argv,"i:I:o:s:h")) != -1){
      switch(c){
      case 'i':
	READ_FLAG=0;
	sprintf(sinput,optarg);
	break;
      case 'I':
	READ_FLAG=1;
	sprintf(sinput,optarg);
	break;
      case 'o':
	sprintf(dirname,optarg);
	break;
      case 's':
	snpid=atoi(optarg);
	break;
      case 'h':
	cerr<<"i: input file name (nemo ascii)"<<endl;
	cerr<<"I: input file name (snapshoe binary)"<<endl;
	cerr<<"o: output directory name"<<endl;
	cerr<<"s: snap file id"<<endl;
	return 0;
      }
    }
    sprintf(smerge,"%s/merge.dat",dirname);
    sprintf(slog,"%s/log.dat",dirname);
    sprintf(sBH,"%s/BH.dat",dirname);
    fout_merge.open(smerge);
    fout_log.open(slog);
    fout_BH.open(sBH);

    fout_merge<<setprecision(15);
    fout_log<<setprecision(15);
    fout_BH<<setprecision(15);

    if(READ_FLAG==0){
      readpara0(prt, Ntot, NSMBH, NIMBH, NBH, NFS,
		eps2_fs_fs, eps2_bh_bh, eps2_fs_smbh, eps2_fs_imbh, 
		eta_s, eta_fs, eta_smbh, eta_imbh, 
		Tsys, Tmerge, Tend, sinput);
    }
    if(READ_FLAG==1){

      readpara1(prt, Ntot, NSMBH, NIMBH, NBH, NFS, Ndead,
		eps2_fs_fs, eps2_bh_bh, eps2_fs_smbh, eps2_fs_imbh, 
		eta_s, eta_fs, eta_smbh, eta_imbh,
		Egr,
		Tsys, Tmerge, Tend, sinput, EX_FLAG);

    }
  }

  bcast(&Ntot, 1, comm);
  bcast(&NSMBH, 1, comm);
  bcast(&NIMBH, 1, comm);
  bcast(&NBH, 1, comm);
  bcast(&NFS, 1, comm);
  bcast(&Ndead, 1, comm);
  bcast(&eps2_fs_fs, 1, comm);
  bcast(&eps2_bh_bh, 1, comm);
  bcast(&eps2_fs_smbh, 1, comm);
  bcast(&eps2_fs_imbh, 1, comm);
  bcast(&eta_s, 1, comm);
  bcast(&eta_fs, 1, comm);
  bcast(&eta_smbh, 1, comm);
  bcast(&eta_imbh, 1, comm);
  bcast(&Egr, 1, comm);
  bcast(&Tsys, 1, comm);
  bcast(&Tmerge, 1, comm);
  bcast(&Tend, 1, comm);
  if(myrank != 0){
    prt = new Particle[Ntot];
  }
  bcast(prt, Ntot, comm);

  // use SSE
  for(int i=0; i<Ntot; i++){
    if(prt[i].index < NBH){
      prt[i].radius = 0.0;
      prt[i].type = SMBH;
    }
    else{
      prt[i].calc_radius();
    }
  }

  Nip_tot = Ntot - Ndead;

  divide_proc(Ntot, Njp_org, first_address);


  int *address = new int[Ntot];
  int *address_old = new int[Ntot];
  for(int i=0; i<Ntot; i++){
    address[i] = i;
  }
  set_eps2(eps2_fs_fs, eps2_fs_smbh, eps2_bh_bh); 

  Njp = 0;
  for(int j=0; j<Ntot; j++){
    if(first_address <= j && 
       j < first_address+Njp_org && 
       prt[j].mass != 0.0){
      prt[j].address = Njp;
      Njp++;
    }
    else{
      prt[j].address = -1;
    }
  }

  Nip = Ntot; 
  set_eta(eta_s, eta_fs, eta_smbh, eta_imbh);
  evolve_initialize(prt, address, Ntot, NBH, Njp, Tsys);

  prt_old = new Particle[Ntot];
  for(int i=0; i<Ntot; i++){
    prt_old[i] = prt[i];
    address_old[i] = address[i];
  }

  double E0, Ek0, Ep0;
  calc_energy(prt, address, Nip_tot, E0, Ek0, Ep0, 0);
  double E1 = E0;
  double Ek1 = Ek0;
  double Ep1 = Ep0;
  double E1_old = E0;
  double Ek1_old = Ek0;
  double Ep1_old = Ep0;
  double Tmerge_old = Tmerge;
  int Nip_tot_old = Nip_tot;
  int Njp_old = Njp;

  int Ndead_old = Ndead;
  int Nmerge_old = Nmerge;
  int Naccrete_old = Naccrete;

  double Tsys_old = Tsys;
  double Egr_old = Egr;
  double Emerge_old = Emerge;

  copy_SMBH_NEW_TO_OLD();

  if(myrank == 0){
    cerr<<"E0="<<E0<<endl;
    cerr<<"Ek0="<<Ek0<<endl;
    cerr<<"Ep0="<<Ep0<<endl;
    write0(prt, Ntot, NBH, Ndead, Tsys, Tmerge, Egr, dirname, snpid, EX_FLAG);
  }

  int state = 0;
  //double dt_max = 1.0/8.0;
  double dt_max = 1.0/1024.0;
  //double dt_max = 1.0/16384.0;
  double dt_snp = 1.0/32.0;
  double Tsnp = Tsys + dt_snp;
  int first_loop = 1;
  int first_loop_old = first_loop;
  Particle *(prt_merged[10000]);
  Particle *(prt_accreted[10000]);
  Tcal_tot0 = MPI_Wtime();
  int itr = 0;

  while(Tsys < Tend){
    if(Nloop % 1000 == 0){
      cerr<<"Tsys="<<Tsys<<",  Nloop="<<Nloop<<",  Nstep="<<get_NSTEP()<<endl;
    }
    Nloop++;
    if(state == 1){
      state = 0;
      first_loop = 1;
    }

    set_eta(eta_s*pow(0.5, itr), eta_fs*pow(0.5, itr), eta_smbh*pow(0.5, itr), eta_imbh*pow(0.5, itr));
    state = evolve_onestep(prt, address, Nip, Ntot, NBH, Tsys, Tmerge, dt_max, first_loop, Egr, itr);

    first_loop = 0;
    if(state == 1){
      double E1_tmp = 0.0; 
      double Ek1_tmp = 0.0;
      double Ep1_tmp = 0.0;
      calc_energy(prt, Ntot, E1_tmp, Ek1_tmp, Ep1_tmp, 0);

      merge_prt();
      Tmerge = Tsys;
      Njp = 0;

      for(int j=first_address; j<first_address+Njp_org; j++){
	if(prt[j].address == -1){continue;}
	prt[j].address = Njp;
	Njp++;
      }

      // do somthing to evolve merged stars using SSE
      get_merged_prt(prt_merged, Nmerge_loop);
      Nmerge += Nmerge_loop;
      get_accreted_prt(prt_accreted, Naccrete_loop);
      Naccrete += Naccrete_loop;
      Ndead += Nmerge_loop + Naccrete_loop;
      Nip_tot = Ntot - Ndead;

      sort_time_all(prt, address, Ntot);


      if(myrank == 0){
	fout_merge<<"--------------------------"<<endl;
	fout_merge<<"Tsys="<<Tsys<<endl;
	fout_merge<<"Ndead="<<Ndead<<endl;
	fout_merge<<"merged particles"<<endl;
	fout_merge<<"Nmerge_loop="<<Nmerge_loop<<endl;
	fout_merge<<"Nmerge="<<Nmerge<<endl;
	for(int i=0; i<Nmerge_loop; i++){
	  fout_merge<<"index="<<prt_merged[i]->index<<endl;
	  fout_merge<<"mass="<<prt_merged[i]->mass<<endl;
	  fout_merge<<"pos="<<prt_merged[i]->pos<<endl;
	  fout_merge<<"vel="<<prt_merged[i]->vel<<endl;
	  fout_merge<<"phi="<<prt_merged[i]->phi<<endl;
	  fout_merge<<endl;
	}
	fout_merge<<endl;
	fout_merge<<"accreted particles"<<endl;
	fout_merge<<"Naccrete_loop="<<Naccrete_loop<<endl;
	fout_merge<<"Naccrete="<<Naccrete<<endl;
	for(int i=0; i<Naccrete_loop; i++){
	  fout_merge<<"index="<<prt_accreted[i]->index<<endl;
	  fout_merge<<"mass="<<prt_accreted[i]->mass<<endl;
	  fout_merge<<"pos="<<prt_accreted[i]->pos<<endl;
	  fout_merge<<"vel="<<prt_accreted[i]->vel<<endl;
	  fout_merge<<"phi="<<prt_accreted[i]->phi<<endl;
	  fout_merge<<endl;
	}
	fout_merge<<endl;
      }

      for(int i=0; i<Ntot; i++){
	prt[i].clear();
	prt[i].acc4 = 0.0;
	prt[i].acc5 = 0.0;
      }
      evolve_initialize(prt,  address,  Ntot,  NBH,  Njp,  Tsys);
      calc_energy(prt, Ntot, E1, Ek1, Ep1, 0);

      // accumulate dissipation energy through merger
      Emerge += E1 - E1_tmp;
      Nip = Nip_tot;
      Tmerge = Tsys;
    }

    if(Tsnp <= Tsys){
      write0(prt, Ntot, NBH, Ndead, Tsys, Tmerge, Egr, dirname, snpid, EX_FLAG);
      Tsnp += dt_snp;
    }

    //if( fmod(Tsys-Tmerge, dt_max) == 0.0 || state == 1){
    if( fmod(Tsys-Tmerge, dt_max) == 0.0 && state == 0){
      Tcal_tot1 = MPI_Wtime() - Tcal_tot0;
      //calc_energy(prt, address, Nip_tot, E1, Ek1, Ep1, 0);
      calc_energy(prt, Ntot, E1, Ek1, Ep1, 0);
      if(myrank == 0){
	cerr<<endl;
	cerr<<"Tsys="<<Tsys<<endl;
	cerr<<"E1="<<E1<<endl;
	cerr<<"E1_old="<<E1_old<<endl;
	cerr<<"Ek1="<<Ek1<<endl;
	cerr<<"Ek1_old="<<Ek1_old<<endl;
	cerr<<"Ep1="<<Ep1<<endl;
	cerr<<"Ep1_old="<<Ep1_old<<endl;
	cerr<<"Emerge="<<Emerge<<endl;
	cerr<<"Emerge_old="<<Emerge_old<<endl;
	cerr<<"Egr="<<Egr<<endl;
	cerr<<"Egr_old="<<Egr_old<<endl;
	cerr<<"(E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old/dt_max="<<(E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old/dt_max<<endl;
	cerr<<"(E1-E0-Egr-Emerge)/E0="<<(E1-E0-Egr-Emerge)/E0<<endl;
      }
      if( dEcrit < fabs( (E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old/dt_max ) ){
	iteration(prt,  prt_old,  address,  address_old,  
		  Ntot,  Ndead_old,  Nmerge_old,  Naccrete_old,
		  Nip_tot_old,  Njp_old,
		  E1_old,  Ep1_old,  Ek1_old,  Egr_old,  Emerge_old,
		  Tsys_old,  Tmerge_old,  
		  Nloop_old,  Nstep_old,  first_loop_old,
		  eta_s, eta_fs, eta_smbh, eta_imbh, 
		  Ndead,   Nmerge,  Naccrete,
		  Nip_tot,  Njp,
		  E1,  Ep1,  Ek1,  Egr,  Emerge,
		  Tsys,  Tmerge,
		  Nloop,  Nstep,  first_loop,
		  itr, 
		  Nip);

	goto ITERATION;
      }
      else{
	if(myrank == 0){
	  cerr<<"no iteration"<<endl;
	  cerr<<endl;
	}
      }

      if(myrank == 0){
	Nstep = get_NSTEP();
	cerr<<"Tsys="<<Tsys<<" ( = "<<Tsys*7.45e3<<"[yr] = "<<Tsys*7.45e3*365<<"[day]"<<endl;
	cerr<<"Tcal_tot="<<Tcal_tot1<<"[sec]"<<endl;
	cerr<<"Tcal_grav="<<Tcal_grav1<<"[sec]"<<endl;
	cerr<<"Tcal_fix="<<Tcal_fix1<<"[sec]"<<endl;
	cerr<<"Tcal_comm="<<Tcal_comm1<<"[sec]"<<endl;
	cerr<<"Nip_tot="<<Nip_tot<<endl;
	cerr<<"Ndead="<<Ndead<<endl;
	cerr<<"Nstep="<<Nstep<<endl;
	cerr<<"Nloop="<<Nloop<<endl;
	cerr<<((double)(Nstep-Nstep_old)*(double)(Nip_tot-Ndead)*97.0)*1e-9/Tcal_tot1<<"[Gflops]"<<endl;
	cerr<<"Nave="<<(Nstep-Nstep_old)/(Nloop-Nloop_old)<<endl;
	cerr<<"Egr="<<Egr<<endl;
	cerr<<"Emerge="<<Emerge<<endl;
	cerr<<"fabs((E1-E0-Egr-Emerge)/E0)="<<fabs((E1-E0-Egr-Emerge)/E0)<<endl;
	cerr<<"fabs((E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old)="<<fabs((E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old)<<endl;
	cerr<<endl;

	fout_log<<"Tsys="<<Tsys<<endl;
	fout_log<<"Tcal_tot1="<<Tcal_tot1<<endl;
	fout_log<<((double)(Nstep-Nstep_old)*(double)(Nip_tot-Ndead)*97.0)*1e-9/Tcal_tot1<<"[Gflops]"<<endl;
	fout_log<<"Nave="<<(Nstep-Nstep_old)/(Nloop-Nloop_old)<<endl;
	fout_log<<"Nloop="<<Nloop<<endl;
	fout_log<<"Egr="<<Egr<<endl;
	fout_log<<"Emerge="<<Emerge<<endl;
	fout_log<<"fabs((E1-E0-Egr-Emerge)/E0)="<<fabs((E1-E0-Egr-Emerge)/E0)<<endl;
	fout_log<<"fabs((E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E0)="<<fabs((E1-E1_old-(Egr-Egr_old)-(Emerge-Emerge_old))/E1_old)<<endl;
	fout_log<<"get_NSTEP()="<<get_NSTEP()<<endl;
	fout_log<<endl;
	fout_BH<<Tsys;
	for(int i=0; i<NBH; i++){
	  fout_BH<<"   "<<prt[i].pos<<"   "<<prt[i].vel<<"   "<<prt[i].phi;
	}
	fout_BH<<endl;
      }

      if(state == 0){
	/*
	cerr<<endl;
	cerr<<endl;
	cerr<<endl;
	cerr<<"Tsys="<<Tsys<<endl;
	cerr<<endl;
	cerr<<endl;
	cerr<<endl;
	*/
	for(int i=0; i<Ntot; i++){
	  prt_old[i] = prt[i];
	  address_old[i] = address[i];
	}

	Ndead_old = Ndead;
	Nmerge_old = Nmerge;
	Naccrete_old = Naccrete;

	Nip_tot_old = Nip_tot;
	Njp_old = Njp;

	E1_old = E1;
	Ek1_old = Ek1;
	Ep1_old = Ep1;

	Egr_old = Egr;
	Emerge_old = Emerge;

	Tcal_tot0 = MPI_Wtime();
	Nstep_old = Nstep;


	Nloop_old = Nloop;
	first_loop_old = first_loop;
	Tcal_grav1 = Tcal_comm1 = Tcal_fix1 = 0.0;
	itr = 0;

	Tsys_old = Tsys;
	Tmerge_old = Tmerge;

	copy_SMBH_NEW_TO_OLD();

      }
    }
  ITERATION: ;
  }
  //for(int j=0; j<Ntot; j++){prt[j].copyold();}

  MPI_Finalize();

  return 0;
}
