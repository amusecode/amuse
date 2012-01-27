#include"merge.h"

struct p2p{
  Particle *prt0;
  Particle *prt1;
};

static p2p merge_candidate[1000];
static Particle *accrete_candidate[1000];

static int Nmerge = 0;
static int Nmerge_loop = 0;
static int Naccrete = 0;
static int Naccrete_loop = 0;
//static int Ndestroy = 0;
//static int Ndestroy_loop = 0;

int merge_check(Particle prt[],
		int address[],
		const int &Nip,
		const int &NBH){
  Nmerge_loop = 0;
  Naccrete_loop = 0;
  int merge_flag = 0;
  double Tsys = prt[address[0]].time + prt[address[0]].dtime;
  for(int i=0; i<Nip; i++){
    int addi = address[i];
    Particle *prti = prt + addi;
    int addj = prti->ngb_index;
    Particle *prtj = prt + addj;
    if(prti->mass == 0.0 || prtj->mass == 0.0){continue;}
    prtj->predict(Tsys);
    Vector3 rij = (prti->pos_pre - prtj->pos_pre);
    double rij2 = rij*rij;
    double rij6 = rij2*rij2*rij2;
    double rcol2 = 0.0;
    double rcol6 = 0.0;
    if(prti->index < NBH && NBH <= prtj->index){
      rcol2 = prtj->radius*prtj->radius;
      rcol6 = rcol2*rcol2*rcol2;
      rcol6 *= (prti->mass*prti->mass)/(prtj->mass*prtj->mass);
    }
    else if(prtj->index < NBH && NBH <= prti->index){
      rcol2 = prti->radius*prti->radius;
      rcol6 = rcol2*rcol2*rcol2;
      rcol6 *= (prtj->mass*prtj->mass)/(prti->mass*prti->mass);
    }
    else{
      rcol2 = 2.0*(prti->radius + prtj->radius);
      rcol2 *= rcol2; 
      rcol6 = rcol2*rcol2*rcol2;
    }
    if(rij2 < rcol2){
    //if(0){
      if(prti-> flag != 0 || prtj->flag !=0){continue;}
      merge_candidate[Nmerge_loop].prt0 = prti;
      merge_candidate[Nmerge_loop].prt1 = prtj;
      prti->flag = 1;
      prtj->flag = 1;
      Nmerge_loop++; 
      merge_flag = 1;
    }
    else{
      double r2 = prti->pos * prti->pos;
      double r6 = r2*r2*r2;
      double rtide6 = calc_rtide_cu(prti->mass, prti->radius);
      rtide6 *= rtide6 ;
      if(r6 < rtide6){
      //if(r2 < 5e-5*5e-5){
      //if(0){
	accrete_candidate[Naccrete_loop] = prti;
	prti->flag = 1;
	merge_flag = 1;
	Naccrete_loop++; 
      }
    }
  }
  return merge_flag;
}

void merge_prt(){

  cerr<<"Naccrete_loop="<<Naccrete_loop<<endl;
  cerr<<"Nmerge_loop="<<Nmerge_loop<<endl;

  for(int i=0; i<Naccrete_loop; i++){
    accrete(accrete_candidate[i]);
  }

  Particle *prt0, *prt1;
  int Nmerge_loop_max = Nmerge_loop;
  Nmerge_loop = 0;
  for(int i=0; i<Nmerge_loop_max; i++){
    prt0 = merge_candidate[i].prt0;
    prt1 = merge_candidate[i].prt1;

    // decision making ... it is not good
    Vector3 rij = prt0->pos - prt1->pos;
    Vector3 vij = prt0->vel - prt1->vel;
    double E2body = 0.5*vij*vij - (prt0->mass+prt1->mass)/sqrt(rij*rij);
    if(get_mpi_rank() == 0){
      cerr<<"E2body="<<E2body<<endl;
    }
    //if(E2body < 0.0){
    if(1){
      merge(prt0, prt1);
      Nmerge_loop++;
    }
    else{
      /*
      if(get_mpi_rank() == 0){
	cerr<<"destroy..."<<endl;
      }
      destroy(prt0, prt1);
      Ndestroy++;
      */
    }
  }
  Naccrete += Naccrete_loop;
  Nmerge += Nmerge_loop;
}


void merge(Particle *prt0,
	   Particle *prt1){
  Particle *prt_light = prt0;
  Particle *prt_heavy = prt1;
  if(prt1->mass < prt0->mass){
    prt_light = prt1;
    prt_heavy = prt0;
  }
  double newmass = prt0->mass + prt1->mass; 
  prt_heavy->mass = newmass;
  prt_heavy->pos = ( (prt0->mass * prt0->pos) + (prt1->mass * prt1->pos) ) / newmass;
  prt_heavy->vel = ( (prt0->mass * prt0->vel) + (prt1->mass * prt1->vel) ) / newmass;
  prt_heavy->calc_radius();
  prt_heavy->flag = 0;
  prt_light->dead();
  if(get_mpi_rank() == 0){
    cerr<<"merge particles..."<<endl;
    cerr<<"prt_heavy->index="<<prt_heavy->index<<endl;
    cerr<<"prt_light->index="<<prt_light->index<<endl;
  }
}

void destroy(Particle *prt0,
	     Particle *prt1){
  prt0->dead();
  prt1->dead();
}

void accrete(Particle *prt0){
  if(get_mpi_rank() == 0){
    cerr<<"accrete onto SMBH"<<endl;
    cerr<<"index="<<prt0->index<<endl;
  }
  accrete_mass(prt0->mass);
  prt0->dead();
}

void get_merged_prt(Particle *(prt_merged[]), int &_Nmerge_loop){
  Particle *prt0, *prt1;
  for(int i=0; i<Nmerge_loop; i++){
    prt0 = merge_candidate[i].prt0;
    prt1 = merge_candidate[i].prt1;
    if(prt0->mass == 0.0){
      prt_merged[i] = prt1;
    }
    else{
      prt_merged[i] = prt0;
    }
  }
  _Nmerge_loop = Nmerge_loop;
}

void get_accreted_prt(Particle *(prt_accreted[]), int &_Naccrete_loop){
  for(int i=0; i<Naccrete_loop; i++){
    prt_accreted[i] = accrete_candidate[i];
  }
  _Naccrete_loop = Naccrete_loop;
}


void set_Ndead(const int &_Nmerge,
	       const int &_Naccrete){
  Nmerge = _Nmerge;
  Naccrete = _Naccrete;
}
