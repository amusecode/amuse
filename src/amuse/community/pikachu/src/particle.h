/*
for PN hybrid scheme 
for 4th order Hermite
PN is evaluated between BH and BH
 */
#ifndef PARTICLE_H
#define PARTICLE_H

#include<cmath>
#include<iostream>
#include"Vector3.h"
//#include"/home/masaki/code/include/Vector/Vector3.h"

using namespace std;

#define dump_cout(x) cout<<#x<<" = "<<x<<endl;
#define dump_cerr(x) cerr<<#x<<" = "<<x<<endl;

const double over3 = 1.0/3.0;
const double over4 = 1.0/4.0;
const double over5 = 1.0/5.0;
const double over6 = 1.0/6.0;
const double over12 = 1.0/12.0;


enum prt_type{star, blackhole, dead, pseudo};

class Neighbour_List;

class Particle
{
private:
    //Particle(Particle &obj);
public:
    prt_type type;
    int have_ngh;

    int index;
    int Nj;
    int idx_ngh_FS;
    int idx_ngh_BH;
    int node_org;

    double mass;
    double radius;
    double pot;
    double pot_short;
    double r2_ngh_FS;
    double r2_ngh_BH;

    // acc, jrk ...etc include PN force also
    // e.g. acc = accFS + accPN
    Vector3 pos;
    Vector3 vel;
    Vector3 acc;
    Vector3 acc_long;
    Vector3 acc_short;

    //for tree
    Particle *prt_next;

    Particle(){
	have_ngh = 0;
	index = 0;
	mass = 0.0;
	radius = 0.0;
	pot = 0.0;
	pot_short = 0.0;
	pos = 0.0;
	vel = 0.0;
	acc = 0.0;
	acc_long = 0.0;
	acc_short = 0.0;
	type = star;
	idx_ngh_FS = -1;
	r2_ngh_FS = -0.0;
	idx_ngh_BH = -1;
	r2_ngh_BH = -0.0;
	Nj = 0;
	node_org = -1;

        prt_next = NULL;
    }

    ~Particle(){}
    
    Particle& operator = (const Particle& p){
	have_ngh = p.have_ngh;
	index = p.index;
	mass = p.mass;
	radius = p.radius;
	pot = p.pot;
	pot_short = p.pot_short;
	pos = p.pos;
	vel = p.vel;
	acc = p.acc;
	acc_long = p.acc_long;
	acc_short = p.acc_short;
	type = p.type;
	Nj = p.Nj;

	idx_ngh_FS = p.idx_ngh_FS;
	r2_ngh_FS = p.r2_ngh_FS;
	idx_ngh_BH = p.idx_ngh_BH;
	r2_ngh_BH = p.r2_ngh_BH;
	node_org = p.node_org;

	return *this;
    }

    void clear(){
	acc = 0.0;
	acc_long = 0.0;
	acc_short = 0.0;
	pot = 0.0;
	pot_short = 0.0;
	have_ngh = 0;
	idx_ngh_FS = -1;
	r2_ngh_FS = 0.0;
	idx_ngh_BH = -1;
	r2_ngh_BH = 0.0;
    }

    void kick(double dt){vel +=  acc*dt;}
    void kick_for_tree(double dt){vel +=  acc_long*dt;}
    void drift(double dt){pos +=  vel*dt;}

    void dump(int flag = 0, ostream& fout=cout){
	fout<<setprecision(15);
	if(flag == 0){
	    fout<<endl;
	    fout<<"index = "<<index<<endl;
	    fout<<"mass = "<<mass<<endl;
	    fout<<"pot = "<<pot<<endl;
	    fout<<"pos = "<<pos<<endl;
	    fout<<"vel = "<<vel<<endl;
	    fout<<"type = "<<type<<endl;
	    fout<<"node_org = "<<node_org<<endl;
	    fout<<"have_ngh = "<<have_ngh<<endl;
	    fout<<endl;
	}
	else if(flag == 1){
	    fout<<endl;
	    fout<<"index = "<<index<<endl;
	    fout<<"mass = "<<mass<<endl;
	    fout<<"pot = "<<pot<<endl;
	    fout<<"pot_short = "<<pot_short<<endl;
	    fout<<"pos = "<<pos<<endl;
	    fout<<"vel = "<<vel<<endl;
	    fout<<"acc = "<<acc<<endl;
	    fout<<"acc_short = "<<acc_short<<endl;
	    fout<<"type = "<<type<<endl;
	    fout<<"node_org = "<<node_org<<endl;
	    fout<<"Nj = "<<Nj<<endl;
	    fout<<"have_ngh = "<<have_ngh<<endl;
	    fout<<"idx_ngh_FS = "<<idx_ngh_FS<<endl;
	    fout<<"r2_ngh_FS = "<<r2_ngh_FS<<endl;
	    fout<<"idx_ngh_BH = "<<idx_ngh_BH<<endl;
	    fout<<"r2_ngh_BH = "<<r2_ngh_BH<<endl;
	    fout<<endl;
	}
	else if(flag == 2){
	    fout<<"index = "<<index<<endl;
	    fout<<"mass = "<<mass<<endl;
	    fout<<"pot = "<<pot<<endl;
	    fout<<"pot_short = "<<pot_short<<endl;
	    fout<<"pos = "<<pos<<endl;
	    fout<<"vel = "<<vel<<endl;
	    fout<<"acc = "<<acc<<endl;
	    fout<<"acc_long = "<<acc_long<<endl;
	    fout<<"acc_short = "<<acc_short<<endl;
	    //fout<<"accPN_short = "<<accPN_short<<endl;
	    fout<<"type = "<<type<<endl;
	    fout<<"idx_ngh_FS = "<<idx_ngh_FS<<endl;
	    fout<<"r2_ngh_FS = "<<r2_ngh_FS<<endl;
	    fout<<"idx_ngh_BH = "<<idx_ngh_BH<<endl;
	    fout<<"r2_ngh_BH = "<<r2_ngh_BH<<endl;
	    fout<<"Nj = "<<Nj<<endl;
	    fout<<"node_org = "<<node_org<<endl;
	    fout<<"have_ngh = "<<have_ngh<<endl;
	    fout<<endl;
	}
    }
};

class Particle_Short
{
public:
    prt_type type;
    Neighbour_List *ngh_list_first;

  int index;
  int Nj;
  int idx_ngh_FS;
  int idx_ngh_BH;

  double time;
  double time_pre;
  double delta_time;
  double mass;
  double pot;
  double pot_short;
  double r2_ngh_FS;
  double r2_ngh_BH;

  // acc, jrk ...etc include PN force also
  // e.g. acc = accFS + accPN
  Vector3 pos;
  Vector3 pos_pre;
  Vector3 pos_old;
  Vector3 vel;
  Vector3 vel_pre;
  Vector3 vel_old;
  Vector3 acc;
  Vector3 acc_short;
  Vector3 acc_short_old;
  Vector3 jrk_short;
  Vector3 jrk_short_old;
  Vector3 acc2_short;
  Vector3 acc3_short;


#ifdef POST_NEWTON
  Vector3 acc_short_pre;
  Vector3 accPN;
  Vector3 accPN_short;
  Vector3 accPN_short_old;
  Vector3 jrkPN;
  Vector3 jrkPN_short;
  Vector3 jrkPN_short_old;
  Vector3 acc2PN_short;
  Vector3 acc3PN_short;
#endif //POST_NEWTON

  double radius;

    // to measur the performance
    int step;

  Particle_Short(){
    index = 0;

    time = 0.0;
    time_pre = 0.0;
    delta_time = 0.0;

    mass = 0.0;
    pot = 0.0;
    pot_short = 0.0;
    pos = 0.0;
    pos_pre = 0.0;
    pos_old = 0.0;
    vel = 0.0;
    vel_pre = 0.0;
    vel_old = 0.0;
    acc = 0.0;
    acc_short = 0.0;
    acc_short_old = 0.0;
    jrk_short = 0.0;
    jrk_short_old = 0.0;
    acc2_short = 0.0;
    acc3_short = 0.0;
    type = star;
    idx_ngh_FS = -1;
    r2_ngh_FS = -0.0;
    idx_ngh_BH = -1;
    r2_ngh_BH = -0.0;
    Nj = 0;

#ifdef POST_NEWTON
    acc_short_pre = 0.0;
    accPN = 0.0;
    accPN_short = 0.0;
    accPN_short_old = 0.0;
    jrkPN = 0.0;
    jrkPN_short = 0.0;
    jrkPN_short_old = 0.0;
    acc2PN_short = 0.0;
    acc3PN_short = 0.0;
#endif

    radius = 0.0;

    step = 0;

  }

  ~Particle_Short(){}

  Particle_Short& operator = (const Particle_Short& p){

    index = p.index;
    time = p.time;
    time_pre = p.time_pre;
    delta_time = p.delta_time;
    mass = p.mass;
    pot = p.pot;
    pot_short = p.pot_short;
    pos = p.pos;
    pos_pre = p.pos_pre;
    pos_old = p.pos_old;
    vel = p.vel;
    vel_pre = p.vel_pre;
    vel_old = p.vel_old;
    acc = p.acc;
    acc_short = p.acc_short;
    acc_short_old = p.acc_short_old;
    jrk_short = p.jrk_short;
    jrk_short_old = p.jrk_short_old;
    acc2_short = p.acc2_short;
    acc3_short = p.acc3_short;
    type = p.type;
    Nj = p.Nj;

    idx_ngh_FS = p.idx_ngh_FS;
    r2_ngh_FS = p.r2_ngh_FS;
    idx_ngh_BH = p.idx_ngh_BH;
    r2_ngh_BH = p.r2_ngh_BH;

    //node_org = p.node_org;


#ifdef POST_NEWTON
    acc_short_pre = p.acc_short_pre;
    accPN = p.accPN;
    accPN_short = p.accPN_short;
    accPN_short_old = p.accPN_short_old;
    jrkPN = p.jrkPN;
    jrkPN_short = p.jrkPN_short;
    jrkPN_short_old = p.jrkPN_short_old;
    acc2PN_short = p.acc2PN_short;
    acc3PN_short = p.acc3PN_short;
#endif


    radius = p.radius;

    step = p.step;

    return *this;
  }

    void clear(){
	acc = 0.0;
	acc_short = 0.0;
	jrk_short = 0.0;
	pot = 0.0;
	pot_short = 0.0;
	idx_ngh_FS = -1;
	r2_ngh_FS = 0.0;
	idx_ngh_BH = -1;
	r2_ngh_BH = 0.0;

#ifdef POST_NEWTON
	acc_short_pre = 0.0;
	accPN = 0.0;
	accPN_short = 0.0;
	accPN_short_old = 0.0;
	jrkPN = 0.0;
	jrkPN_short = 0.0;
	jrkPN_short_old = 0.0;
	acc2PN_short = 0.0;
	acc3PN_short = 0.0;
#endif

    }

  void kick(double dt){vel +=  acc*dt;}
  void drift(double dt){pos +=  vel*dt;}

  void predict_short_h4(const double &dt){
    pos_pre = ( ( jrk_short*dt*over3 + acc_short ) *dt*0.5 + vel ) *dt + pos;  // 8ops
    vel_pre =   ( jrk_short*dt*0.5   + acc_short ) *dt + vel;  // 5ops

    //pos_pre = ( ( ( acc2_short*dt*over4 + jrk_short ) *dt*over3 + acc_short ) *dt*0.5 + vel ) *dt + pos; // 11ops
    //vel_pre =   ( ( acc2_short*dt*over3 + jrk_short ) *dt*0.5   + acc_short ) *dt + vel;  // 8ops

    //pos_pre = ( ( ( ( acc3_short*dt*over5 + acc2_short ) *dt*over4 + jrk_short ) *dt*over3 + acc_short ) *dt*0.5 + vel ) *dt + pos; // 14 ops
    //vel_pre =   ( ( ( acc3_short*dt*over4 + acc2_short ) *dt*over3 + jrk_short ) *dt*0.5   + acc_short ) *dt + vel;  // 11 ops

  }

    void correct_h2(){
        vel +=  ( acc_short_old + acc_short ) * 0.5 * delta_time;
        pos +=  ( vel_old + vel ) * 0.5 * delta_time;
    }

    void correct_h4(){
        vel = vel_old + ( ( jrk_short_old - jrk_short ) *delta_time*over6 + ( acc_short_old + acc_short ) ) *0.5*delta_time; // 7ops
        pos = pos_old + ( ( acc_short_old - acc_short ) *delta_time*over6 + ( vel_old + vel ) ) *0.5*delta_time;  // 7ops
        //pos = pos_old + ( ( ( ( jrk_short_old + jrk_short ) *delta_time*over12 + ( acc_short_old - acc_short ) ) *delta_time*over5 ) + (vel_old + vel)) *delta_time*0.5;    // 11ops
  }

  void interpolate_h4(){
    double DT = 1.0 / delta_time;
    double DT2 = DT * DT;
    double DT3 = DT2 * DT;
    Vector3 acc2_short_old = -2.0*DT2* ( 3.0* ( acc_short_old - acc_short ) + ( 2.0*jrk_short_old + jrk_short ) * delta_time);
    Vector3 acc3_short_old =  6.0*DT3* ( 2.0* ( acc_short_old - acc_short ) + ( jrk_short_old + jrk_short ) * delta_time);
    acc2_short = acc3_short_old*delta_time + acc2_short_old;
    acc3_short = acc3_short_old;
  }


    void set_dt_2nd(const double &maxdt, 
		    const double &eta, 
		    const double &Toffset,
		    const double &Tsync,
		    const double acc_offset_sq = 0.0){
	static double log2 = log(2.0);
	double a0_2 = acc_short * acc_short;
	double a1_2 = jrk_short * jrk_short;
	double dt_tmp = 0.0;
	if(a0_2 == 0.0 || a1_2 == 0.0){
	    dt_tmp = maxdt;
	}
	else{
	    a0_2 += acc_offset_sq;
	    double a1_2 = jrk_short * jrk_short;
	    dt_tmp = eta * sqrt( a0_2/a1_2 );
	}
	
	if(time == Toffset){
	    //if(maxdt < dt_tmp){
	    if(maxdt <= dt_tmp){
		delta_time = maxdt;
	    }
	    else{
		int power = log(dt_tmp)/log2;
		delta_time = pow(2.0, (double)(power-1));
	    }
	}
	else{
	    if( dt_tmp < delta_time ){
		int power = log(dt_tmp)/log2;
		delta_time = pow(2.0, (double)(power-1));
	    }
	    /*
	    else if( (dt_tmp >= 2.0*delta_time)  && 
		     (dt_tmp <= maxdt)  && 
		     (fmod((time - Toffset), 2.0*delta_time) == 0.0) ){
	    */

	    else if( (dt_tmp >= 2.0*delta_time)  && 
		     (2.0*delta_time <= maxdt)  && 
		     (fmod((time - Toffset), 2.0*delta_time) == 0.0) ){

		delta_time *= 2.0;
	    }
	}
	if(Tsync < delta_time + time){
	    delta_time = Tsync -time;
	}
    }


    void set_dt_4th(const double &maxdt, 
		    const double &eta, 
		    const double &Toffset,
		    const double &Tsync,
		    const double acc_offset_sq = 0.0){
	static double log2 = log(2.0);
	double a0_2 = acc_short * acc_short;
	double a1_2 = jrk_short * jrk_short;
	double dt_tmp = 0.0;
	if(a0_2 == 0.0 || a1_2 == 0.0){
	    dt_tmp = maxdt;
	}
	else{
	    a0_2 += acc_offset_sq;
	    double a0 = sqrt(a0_2);
	    double a1 = sqrt(a1_2);
	    double a2_2 = acc2_short * acc2_short;
	    double a2 = sqrt(a2_2);
	    double a3_2 = acc3_short * acc3_short;
	    double a3 = sqrt(a3_2);
	    double A0 = a0 * a2 + a1_2;
	    double A1 = a1 * a3 + a2_2;
	    dt_tmp = eta*sqrt(A0/A1);
	}

	if(time == Toffset){
	    //if(maxdt < dt_tmp){
	    if(maxdt <= dt_tmp){
		delta_time = maxdt;
	    }
	    else{
		int power = log(dt_tmp)/log2;
		delta_time = pow(2.0, (double)(power-1));
	    }
	}
	else{
	    if( dt_tmp < delta_time ){
		int power = log(dt_tmp)/log2;
		delta_time = pow(2.0, (double)(power-1));
	    }
	    /*
	    else if( (dt_tmp > 2.0*delta_time)  && 
		     (dt_tmp <= maxdt)  && 
		     (fmod((time - Toffset), 2.0*delta_time) == 0.0) ){
	    */
	    else if( (dt_tmp >= 2.0*delta_time)  && 
		     (2.0*delta_time <= maxdt)  && 
		     (fmod((time - Toffset), 2.0*delta_time) == 0.0) ){
		delta_time *= 2.0;
	    }
	}
	if(Tsync < delta_time + time){
	    delta_time = Tsync -time;
	}
    }


    void dump(int flag = 0, ostream &fout = cout){
    fout<<setprecision(15);
    if(flag == 0){
      fout<<endl;
      fout<<"index = "<<index<<endl;
      fout<<"time = "<<time<<endl;
      fout<<"delta_time = "<<delta_time<<endl;
      fout<<"mass = "<<mass<<endl;
      fout<<"pot = "<<pot<<endl;
      fout<<"pos = "<<pos<<endl;
      fout<<"vel = "<<vel<<endl;
      fout<<"type = "<<type<<endl;
      //fout<<"node_org = "<<node_org<<endl;
      fout<<endl;
    }

    if(flag == 1){
      fout<<endl;
      fout<<"index = "<<index<<endl;
      fout<<"time = "<<time<<endl;
      fout<<"mass = "<<mass<<endl;
      fout<<"pot = "<<pot<<endl;
      fout<<"pos = "<<pos<<endl;
      fout<<"vel = "<<vel<<endl;
      fout<<"acc_short = "<<acc_short<<endl;
      fout<<"jrk_short = "<<jrk_short<<endl;
      fout<<"type = "<<type<<endl;
      //fout<<"node_org = "<<node_org<<endl;
      fout<<"Nj = "<<Nj<<endl;
      fout<<"idx_ngh_FS = "<<idx_ngh_FS<<endl;
      //fout<<"adr_ngh_FS="<<adr_ngh_FS<<endl;
      fout<<"r2_ngh_FS = "<<r2_ngh_FS<<endl;
      fout<<"idx_ngh_BH = "<<idx_ngh_BH<<endl;
      //fout<<"adr_ngh_BH="<<adr_ngh_BH<<endl;
      fout<<"r2_ngh_BH = "<<r2_ngh_BH<<endl;
      fout<<endl;
    }

    if(flag == 2){
    fout<<"index = "<<index<<endl;
    fout<<"time = "<<time<<endl;
    fout<<"delta_time = "<<delta_time<<endl;
    fout<<"time_pre = "<<time_pre<<endl;
    fout<<"mass = "<<mass<<endl;
    fout<<"pot = "<<pot<<endl;
    fout<<"pot_short = "<<pot_short<<endl;
    fout<<"pos = "<<pos<<endl;
    fout<<"pos_pre = "<<pos_pre<<endl;
    fout<<"pos_old = "<<pos_old<<endl;
    fout<<"vel = "<<vel<<endl;
    fout<<"vel_pre = "<<vel_pre<<endl;
    fout<<"vel_old = "<<vel_old<<endl;
    fout<<"acc = "<<acc<<endl;
    fout<<"acc_short = "<<acc_short<<endl;
    fout<<"acc_short_old = "<<acc_short_old<<endl;
    fout<<"jrk_short = "<<jrk_short<<endl;
    fout<<"jrk_short_old = "<<jrk_short_old<<endl;
    fout<<"acc2_short = "<<acc2_short<<endl;
    fout<<"acc3_short = "<<acc3_short<<endl;

#ifdef POST_NEWTON
    fout<<"acc_short_pre = "<<acc_short_pre<<endl;
    fout<<"accPN_short = "<<accPN_short<<endl;
    fout<<"accPN_short_old = "<<accPN_short_old<<endl;
    fout<<"jrkPN_short = "<<jrkPN_short<<endl;
    fout<<"jrkPN_short_old = "<<jrkPN_short_old<<endl;
    fout<<"acc2PN_short = "<<acc2PN_short<<endl;
    fout<<"acc3PN_short = "<<acc3PN_short<<endl;
#endif

    fout<<"type = "<<type<<endl;

    fout<<"idx_ngh_FS = "<<idx_ngh_FS<<endl;
    fout<<"r2_ngh_FS = "<<r2_ngh_FS<<endl;

    fout<<"idx_ngh_BH = "<<idx_ngh_BH<<endl;
    fout<<"r2_ngh_BH = "<<r2_ngh_BH<<endl;

    fout<<"Nj = "<<Nj<<endl;
    fout<<endl;
    }
  }


#ifdef POST_NEWTON
    double calc_Egr_2nd(){
        return 0.5*mass*delta_time*(accPN_short*vel 
                                    + accPN_short_old*vel_old);
    }

    double calc_Egr_4th(){
        return 0.5*mass*delta_time*( (
                                         accPN_short*vel 
                                         + accPN_short_old*vel_old
                                         )
                                     + delta_time*over6*( 
                                         jrkPN_short_old*vel_old 
                                         + accPN_short_old*acc_short_old 
                                         - jrkPN_short*vel 
                                         - accPN_short*acc_short
                                         ) 	 );
    }

    Vector3 calc_Lgr_2nd(){
        return 0.5*mass*delta_time*( (pos^accPN_short) + (pos_old^accPN_short_old) );
    }
 
    Vector3 calc_Lgr_4th(){
        return 0.5*mass*delta_time*( (
                                         (pos^accPN_short) + (pos_old^accPN_short_old) 
                                         )
                                     + delta_time*over6*( 
                                         (vel_old^accPN_short_old)
                                         +(pos_old^jrkPN_short_old)
                                         -(vel^accPN_short)
                                         -(pos^jrkPN_short)
                                         )   );
    }
#endif

  void calc_radius(){
    //radius = invMunit * get_radius_zams0(mass*Munit);
    //radius = (Rsolsi/Pcsi) * get_radius_zams0(mass*Munit);
    radius = 0.0;
  }

  double calc_rtide(const double &SMBHmass){
    //rtide = radius * pow(SMBHmass/mass, over3);
    //return radius * pow(SMBHmass/mass, over3);
    return 0.0;
  } 

};







class Neighbour_List{
public:
    Particle_Short *prti;
    Particle_Short *prtj;
    Neighbour_List(){
	prti = NULL;
	prtj = NULL;
    }
    void init(){
	prti = NULL;
	prtj = NULL;
    }
    void dump(ostream& fout = cout){
	int idxi = prti->index;
	int idxj = prtj->index;
	fout<<idxi<<"   "<<idxj<<endl;
    }
};


class Neighbour_List_Index{
public:
    int idx_i;
    int idx_j;
    void set(const Neighbour_List &_ngh_list){
	idx_i = _ngh_list.prti->index;
	idx_j = _ngh_list.prtj->index;
    }
    void dump(ostream& fout = cout){
	fout<<idx_i<<"   "<<idx_j<<endl;
    }
};

#endif //PARTICLE_H
