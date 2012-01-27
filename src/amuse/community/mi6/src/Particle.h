#ifndef PARTICLE_H
#define PARTICLE_H

#include<cmath>
#include<iostream>
#include"Vector3.h"
#include"unit.h"
#include"stellar_evolution.h"

using namespace std;

enum prt_type{MS, RG, WR, WD, NS, BH, SMBH, IMBH, DEAD};

//static double over3 = 1.0/3.0;

class Particle
{
 private:
  //Particle(Particle &obj);
 public:
  int index;

  double radius;

  double mass;
  double time;
  double dtime;

  double phi;
  double phi_old;

  Vector3 pos;
  Vector3 pos_old;
  Vector3 pos_pre;
  Vector3 vel;
  Vector3 vel_old;
  Vector3 vel_pre;
  Vector3 acc;
  Vector3 acc_old;
  Vector3 acc_pre;
  Vector3 jrk;
  Vector3 jrk_old;
  Vector3 jrk_pre;
  Vector3 acc2;
  Vector3 acc2_old;
  Vector3 acc3;
  Vector3 acc3_old;
  Vector3 acc4;
  Vector3 acc4_old;
  Vector3 acc5;
  Vector3 acc5_old;
  /*
  Vector3 acc6;
  Vector3 acc6_old;
  Vector3 acc7;
  Vector3 acc7_old;
  */
  Vector3 acc_ex;
  Vector3 acc_ex_old;
  Vector3 jrk_ex;
  Vector3 jrk_ex_old;
  Vector3 acc2_ex;
  Vector3 acc2_ex_old;
  Vector3 acc3_ex;
  Vector3 acc3_ex_old;
  double phi_ex;
  double phi_ex_old;

  Vector3 acc_PN_ex;
  Vector3 acc_PN_ex_old;
  Vector3 jrk_PN_ex;
  Vector3 jrk_PN_ex_old;
  Vector3 acc2_PN_ex;
  Vector3 acc2_PN_ex_old;


  prt_type type;

  int address;

  int flag;
  double rtide;

  int ngb_index;
  double r2min;

  Particle(){
    index = 0;

    radius = 0.0;

    mass = 0.0;
    time = 0.0;
    dtime = 0.0;

    phi = 0.0;
    phi_old = 0.0;

    pos = 0.0;
    pos_old = 0.0;
    pos_pre = 0.0;
    vel = 0.0;
    vel_old = 0.0;
    vel_pre = 0.0;
    acc = 0.0;
    acc_old = 0.0;
    acc_pre = 0.0;
    jrk = 0.0;
    jrk_old = 0.0;
    jrk_pre = 0.0;
    acc2 = 0.0;
    acc3 = 0.0;
    acc4 = 0.0;
    acc5 = 0.0;
    acc4_old = 0.0;
    acc5_old = 0.0;
    /*
    acc6 = 0.0;
    acc7 = 0.0;
    acc6_old = 0.0;
    acc7_old = 0.0;
    */
    acc_ex = 0.0;
    acc_ex_old = 0.0;
    jrk_ex = 0.0;
    jrk_ex_old = 0.0;
    acc2_ex = 0.0;
    acc2_ex_old = 0.0;
    acc3_ex = 0.0;
    acc3_ex_old = 0.0;
    phi_ex = 0.0;
    phi_ex_old = 0.0;

    acc_PN_ex = 0.0;
    acc_PN_ex_old = 0.0;
    jrk_PN_ex = 0.0;
    jrk_PN_ex_old = 0.0;
    acc2_PN_ex = 0.0;
    acc2_PN_ex_old = 0.0;

    type = MS;

    address = 0;
    flag = 0;

    rtide = 0.0;

    ngb_index = 0;
    r2min = 0.0;
  }

  /*
  Particle& operator = (const Particle& p){
    index = p.index;

    radius = p.radius;

    mass = p.mass;
    time = p.time;
    dtime = p.dtime;

    phi = p.phi;
    phi_old = p.phi_old;

    pos = p.pos;
    pos_old = p.pos_old;
    pos_pre = p.pos_pre;
    vel = p.vel;
    vel_old = p.vel_old;
    vel_pre = p.vel_pre;
    acc = p.acc;
    acc_old = p.acc_old;
    acc_pre = p.acc_pre;
    jrk = p.jrk;
    jrk_old = p.jrk_old;
    jrk_pre = p.jrk_pre;
    acc2 = p.acc2;
    acc2_old = p.acc2_old;
    acc3 = p.acc3;
    acc3_old = p.acc3_old;
    acc4 = p.acc4;
    acc4_old = p.acc4_old;
    acc5 = p.acc5;
    acc5_old = p.acc5_old;

    acc_ex = p.acc_ex;
    acc_ex_old = p.acc_ex_old;
    jrk_ex = p.jrk_ex;
    jrk_ex_old = p.jrk_ex_old;
    acc2_ex = p.acc2_ex;
    acc2_ex_old = p.acc2_ex_old;
    acc3_ex = p.acc3_ex;
    acc3_ex_old = p.acc3_ex_old;
    phi_ex = p.phi_ex;
    phi_ex_old = p.phi_ex_old;

    acc_PN_ex = p.acc_PN_ex;
    acc_PN_ex_old = p.acc_PN_ex_old;
    jrk_PN_ex = p.jrk_PN_ex;
    jrk_PN_ex_old = p.jrk_PN_ex_old;
    acc2_PN_ex = p.acc2_PN_ex;
    acc2_PN_ex_old = p.acc2_PN_ex_old;

    type = p.type;
    address = p.address;

    flag = p.flag;

    rtide = p.rtide;

    ngb_index = p.ngb_index;
    r2min = r2min;

    return *this;
  }
  */


  void copyold(){ 
    pos_old = pos;  
    vel_old = vel;
    acc_old = acc;  
    jrk_old = jrk;  
    acc2_old = acc2;  
    acc3_old = acc3; 
    phi_old = phi;

    acc_ex_old = acc_ex; 
    jrk_ex_old = jrk_ex; 
    acc2_ex_old = acc2_ex; 
    acc3_ex_old = acc3_ex; 
    phi_ex_old = phi_ex; 

    acc_PN_ex_old = acc_PN_ex; 
    jrk_PN_ex_old = jrk_PN_ex; 
    acc2_PN_ex_old = acc2_PN_ex;
  }

  void clear(){ 
    acc = 0.0; 
    jrk = 0.0; 
    acc2 = 0.0; 
    acc3 = 0.0; 
    phi = 0.0;

    acc_ex = 0.0; 
    jrk_ex = 0.0; 
    acc2_ex = 0.0; 
    acc3_ex = 0.0; 
    phi_ex = 0.0;

    acc_PN_ex = 0.0; 
    jrk_PN_ex = 0.0; 
    acc2_PN_ex = 0.0;
  }

  void set_dt_2nd(const double &maxdt, 
		  const double &eta, 
		  double time_offset = 0.0){
    static double log2 = log(2.0);
    double a0_2 = acc * acc;
    double a1_2 = jrk * jrk;
    double dt_tmp = eta * sqrt( a0_2/a1_2 );
    if(time == time_offset){
      if(maxdt < dt_tmp){
	dtime = maxdt;
      }
      else{
	int power = log(dt_tmp)/log2;
	dtime = pow(2.0, (double)(power-1));
      }
    }
    else{
      if( dt_tmp < dtime ){
	int power = log(dt_tmp)/log2;
	dtime = pow(2.0, (double)(power-1));
      }
      else if( (dt_tmp >= 2.0*dtime)  && 
	       (dt_tmp <= maxdt)  && 
	       (fmod((time - time_offset), 2.0*dtime) < 1e-18) && 
	       (time - time_offset) - dtime != 0.0){
	dtime = dtime*2.0;
      }
    }
  }

  void set_dt_4th(const double &maxdt, 
		  const double &eta, 
		  double time_offset=0.0){
    double a0_2 = acc * acc;
    double a0 = sqrt(a0_2);
    double a1_2 = jrk * jrk;
    double a1 = sqrt(a1_2);
    double a2_2 = acc2 * acc2;
    double a2 = sqrt(a2_2);
    double a3_2 = acc3 * acc3;
    double a3 = sqrt(a3_2);

    double A0 = a0 * a2 + a1_2;
    double A1 = a1 * a3 + a2_2;
    double dt_tmp = eta*sqrt(A0/A1);

    static double log2 = log(2.0);
    if(time == time_offset){
      if(maxdt < dt_tmp){
	dtime = maxdt;
      }
      else{
	int power = log(dt_tmp)/log2;
	dtime = pow(2.0, (double)(power-1));
      }
    }
    else{
      if( dt_tmp < dtime ){
	int power = log(dt_tmp)/log2;
	dtime = pow(2.0, (double)(power-1));
      }
      else if( (dt_tmp > 2.0*dtime)  && 
	       (dt_tmp <= maxdt)  && 
	       (fmod((time - time_offset), 2.0*dtime) < 1e-18) && 
	       (time - time_offset) - dtime != 0.0){
	dtime = dtime*2.0;
      }
    }
  }


  void set_dt_6th(const double &maxdt, 
		  const double &eta, 
		  double time_offset=0.0){
    double a0_2 = acc * acc;
    double a0 = sqrt(a0_2);
    double a1_2 = jrk * jrk;
    double a2_2 = acc2 * acc2;
    double a2 = sqrt(a2_2);
    double a3_2 = acc3 * acc3;
    double a3 = sqrt(a3_2);
    double a4_2 = acc4 * acc4;
    double a5_2 = acc5 * acc5;
    double a5 = sqrt(a5_2);
    double A0 = a0 * a2 + a1_2;
    double A1 = a3 * a5 + a4_2;
    double dt_tmp = eta*pow(A0/A1, 0.1666666666666666666);
    static double log2 = log(2.0);
    if(time == time_offset){
      if(maxdt < dt_tmp){
	dtime = maxdt;
      }
      else{
	int power = log(dt_tmp)/log2;
	dtime = pow(2.0, (double)(power-1));
      }
    }
    else{
      if( dt_tmp < dtime ){
	int power = log(dt_tmp)/log2;
	dtime = pow(2.0, (double)(power-1));
      }
      else if( (dt_tmp >= 2.0*dtime)  && 
	       (dt_tmp <= maxdt)  && 
	       (fmod((time - time_offset), 2.0*dtime) < 1e-18) && 
	       (time - time_offset) - dtime != 0.0){
	dtime = dtime*2.0;
      }
    }
  }



  void predict(double system_time){
    double dt_tmp = system_time - time;
    /*
    pos_pre = ((((((((acc7*dt_tmp*n9 + acc6)*dt_tmp*n8 + acc5)*dt_tmp*n7  + acc4)*dt_tmp*n6  + acc3)*dt_tmp*n5  + acc2)*dt_tmp*n4  + jrk)*dt_tmp*n3   + acc)*dt_tmp*0.5 + vel)*dt_tmp + pos;
    vel_pre =  (((((((acc7*dt_tmp*n8 + acc6)*dt_tmp*n7 + acc5)*dt_tmp*n6  + acc4)*dt_tmp*n5  + acc3)*dt_tmp*n4  + acc2)*dt_tmp*n3  + jrk)*dt_tmp*0.5  + acc)*dt_tmp + vel;
    acc_pre =   ((((((acc7*dt_tmp*n7 + acc6)*dt_tmp*n6 + acc5)*dt_tmp*n5  + acc4)*dt_tmp*n4  + acc3)*dt_tmp*n3  + acc2)*dt_tmp*0.5  + jrk)*dt_tmp  + acc;
    jrk_pre =    (((((acc7*dt_tmp*n6 + acc6)*dt_tmp*n5 + acc5)*dt_tmp*n4  + acc4)*dt_tmp*n3  + acc3)*dt_tmp*0.5 + acc2)*dt_tmp  + jrk;
    */

    pos_pre = ((((((acc5*dt_tmp*n7  + acc4)*dt_tmp*n6  + acc3)*dt_tmp*n5  + acc2)*dt_tmp*n4  + jrk)*dt_tmp*n3  + acc)*dt_tmp*0.5 + vel)*dt_tmp + pos;
    vel_pre =  (((((acc5*dt_tmp*n6  + acc4)*dt_tmp*n5  + acc3)*dt_tmp*n4  + acc2)*dt_tmp*n3  + jrk)*dt_tmp*0.5 + acc)*dt_tmp     + vel;
    acc_pre =   ((((acc5*dt_tmp*n5  + acc4)*dt_tmp*n4  + acc3)*dt_tmp*n3  + acc2)*dt_tmp*0.5 + jrk)*dt_tmp     + acc;
    jrk_pre =    (((acc5*dt_tmp*n4  + acc4)*dt_tmp*n3  + acc3)*dt_tmp*0.5 + acc2)*dt_tmp     + jrk;

    /*
    pos_pre = ((((acc3*dt_tmp*n5  + acc2)*dt_tmp*n4  + jrk)*dt_tmp*n3  + acc)*dt_tmp*0.5 + vel)*dt_tmp + pos;
    vel_pre =  (((acc3*dt_tmp*n4  + acc2)*dt_tmp*n3  + jrk)*dt_tmp*0.5 + acc)*dt_tmp     + vel;
    acc_pre =   ((acc3*dt_tmp*n3  + acc2)*dt_tmp*0.5 + jrk)*dt_tmp     + acc;
    */
    /*
    pos_pre = ((jrk*dt_tmp*n3  + acc)*dt_tmp*0.5 + vel)*dt_tmp + pos;
    vel_pre =  (jrk*dt_tmp*0.5 + acc)*dt_tmp     + vel;
    */

  }

  void correct(){
    /*
    pos = ((((((((acc7_old*dtime*n9 + acc6_old)*dtime*n8 + acc5_old)*dtime*n7  + acc4_old)*dtime*n6  + acc3_old)*dtime*n5  + acc2_old)*dtime*n4  + jrk_old)*dtime*n3  + acc_old)*dtime*0.5 + vel_old)*dtime + pos_old;
    vel =  (((((((acc7_old*dtime*n8 + acc6_old)*dtime*n7 + acc5_old)*dtime*n6  + acc4_old)*dtime*n5  + acc3_old)*dtime*n4  + acc2_old)*dtime*n3  + jrk_old)*dtime*0.5  + acc_old)*dtime + vel_old;
    */
    /*
    pos = ((((((acc5_old*dtime*n7  + acc4_old)*dtime*n6  + acc3_old)*dtime*n5  + acc2_old)*dtime*n4  + jrk_old)*dtime*n3  + acc_old)*dtime*0.5 + vel_old)*dtime + pos_old;
    vel =  (((((acc5_old*dtime*n6  + acc4_old)*dtime*n5  + acc3_old)*dtime*n4  + acc2_old)*dtime*n3  + jrk_old)*dtime*0.5 + acc_old)*dtime     + vel_old;
    */

    /*
    vel += ((( (acc3_old - acc3)*n1680*dtime + (acc2_old + acc2)*n84 )*dtime + (jrk_old - jrk)*n3_28 )*dtime + (acc_old + acc)*0.5)*dtime;
    pos += ((( (acc2_old - acc2)*n1680*dtime + (jrk_old + jrk)*n84   )*dtime + (acc_old - acc)*n3_28 )*dtime + (vel_old + vel)*0.5)*dtime;
    */
    /*
    vel += ((  (acc2_old + acc2)*n120*dtime + (jrk_old - jrk)*n10 )*dtime + (acc_old + acc)*0.5 )*dtime;
    pos += ((( (acc2_old - acc2)*n1680*dtime + (jrk_old + jrk)*n84 )*dtime + (acc_old - acc)*n3_28 )*dtime + (vel_old + vel)*0.5)*dtime;
    */

    vel +=  ( (jrk_old - jrk)*n12*dtime  + (acc_old + acc)*0.5 )*dtime;
    pos +=  ( (acc_old - acc)*n12*dtime  + (vel_old + vel)*0.5 )*dtime;

  }


  void correct_mix(){
    vel = vel_old + ((( (acc3_ex_old - acc3_ex)*n1680*dtime + (acc2_ex_old + acc2_ex)*n84 )*dtime + (jrk_ex_old - jrk_ex)*n3_28 )*dtime + (acc_ex_old + acc_ex)*0.5)*dtime
      + ((  ( (acc2_old-acc2_ex_old) + (acc2-acc2_ex))*n120*dtime  + ((jrk_old-jrk_ex_old) - (jrk-jrk_ex))*n10 )*dtime + ((acc_old-acc_ex_old) + (acc-acc_ex))*0.5 )*dtime;
    pos = pos_old + ((( (acc2_old - acc2)*n1680*dtime + (jrk_old + jrk)*n84 )*dtime + (acc_old - acc)*n3_28 )*dtime + (vel_old + vel)*0.5)*dtime;
    /*
    vel += ((( (acc3_ex_old - acc3_ex)*n1680*dtime + (acc2_ex_old + acc2_ex)*n84 )*dtime + (jrk_ex_old - jrk_ex)*n3_28 )*dtime + (acc_ex_old + acc_ex)*0.5)*dtime;
    vel += ((  ( (acc2_old-acc2_ex_old) + (acc2-acc2_ex))*n120*dtime  + ((jrk_old-jrk_ex_old) - (jrk-jrk_ex))*n10 )*dtime + ((acc_old-acc_ex_old) + (acc-acc_ex))*0.5 )*dtime;
    pos += ((( (acc2_old - acc2)*n1680*dtime + (jrk_old + jrk)*n84 )*dtime + (acc_old - acc)*n3_28 )*dtime + (vel_old + vel)*0.5)*dtime;
    */
  }



  void interpolate(){
    /*
    double DT = 1.0/dtime; double DT2 = DT*DT; double DT3 = DT2*DT;  double DT4 = DT2*DT2;
    acc4_old = -4.0*(((210.0*(acc_old - acc)*DT + (120.0*jrk_old + 90.0*jrk))*DT +  (30.0*acc2_old - 15.0*acc2))*DT + (4.0*acc3_old + acc3))*DT;
    acc5_old = 60.0*(((168.0*(acc_old - acc)*DT + (90.0*jrk_old + 78.0*jrk))*DT +   (20.0*acc2_old - 14.0*acc2))*DT + (2.0*acc3_old + acc3))*DT2;
    acc6_old = -120.0*(((420.0*(acc_old - acc)*DT + (216.0*jrk_old + 204.0*jrk))*DT + (45.0*acc2_old - 39.0*acc2))*DT + (4.0*acc3_old + 3.0*acc3))*DT3;
    acc7_old = 840.0*(((120.0*(acc_old - acc)*DT + 60.0*(jrk_old + jrk))*DT +   12.0*(acc2_old - acc2))*DT + (acc3_old + acc3))*DT4;
    acc4 = ((acc7_old*dtime*n3 + acc6_old)*dtime*0.5 + acc5_old)*dtime + acc4_old;
    acc5 = (acc7_old*dtime*0.5 + acc6_old)*dtime + acc5_old;
    acc6 = acc7_old*dtime + acc6_old;
    acc7 = acc7_old;
    */

    double DT=1.0/dtime; double DT2=DT*DT; double DT3=DT2*DT; 
    acc3_old = -3.0*( ((20.0*(acc_old - acc)*DT + (12.0*jrk_old + 8.0*jrk))*DT + (3.0*acc2_old - acc2))*DT);
    acc4_old = 12.0*( ((30.0*(acc_old - acc)*DT + (16.0*jrk_old + 14.0*jrk))*DT + (3.0*acc2_old - 2.0*acc2))*DT2);
    acc5_old = -60.0*(((12.0*(acc_old - acc)*DT + 6.0*(jrk_old + jrk))*DT + (acc2_old - acc2))*DT3);
    acc3 = (acc5_old*dtime*0.5 + acc4_old)*dtime + acc3_old; // O(dt^3)
    acc4 = acc5_old*dtime + acc4_old;  // O(dt^2)
    acc5 = acc5_old;  // O(dt^1)

    /*
    double DT=1.0/dtime; double DT2=DT*DT; double DT3=DT2*DT; 
    acc2_old = -2.0*DT2*(3.0*(acc_old - acc) + (2.0*jrk_old + jrk)*dtime);
    acc3_old = 6.0*DT3*(2.0*(acc_old - acc) + (jrk_old + jrk)*dtime);
    acc2 = acc3_old*dtime + acc2_old;
    acc3 = acc3_old;
    */
  }

  void accumulate_force(){
    acc += (acc_ex + acc_PN_ex);
    jrk += (jrk_ex + jrk_PN_ex);
    acc2 += (acc2_ex + acc2_PN_ex);
    acc3 += acc3_ex;
    phi += phi_ex;
  }


  void dump(){
    cout<<setprecision(15);
    cout<<"index = "<<index<<endl;
    cout<<"address = "<<address<<endl;
    cout<<"radius = "<<radius<<endl;
    cout<<"mass = "<<mass<<endl;
    cout<<"time = "<<time<<endl;
    cout<<"dtime = "<<dtime<<endl;
    cout<<"phi = "<<phi<<endl;
    cout<<"phi_old = "<<phi_old<<endl;
    cout<<"pos = "<<pos<<endl;
    cout<<"pos_pre = "<<pos_pre<<endl;
    cout<<"pos_old = "<<pos_old<<endl;
    cout<<"vel = "<<vel<<endl;
    cout<<"vel_pre = "<<vel_pre<<endl;
    cout<<"vel_old = "<<vel_old<<endl;
    cout<<"acc = "<<acc<<endl;
    cout<<"acc_pre = "<<acc_pre<<endl;
    cout<<"acc_old = "<<acc_old<<endl;
    cout<<"jrk = "<<jrk<<endl;
    cout<<"jrk_pre = "<<jrk_pre<<endl;
    cout<<"jrk_old = "<<jrk_old<<endl;
    cout<<"acc2 = "<<acc2<<endl;
    cout<<"acc2_old = "<<acc2_old<<endl;
    cout<<"acc3 = "<<acc3<<endl;
    cout<<"acc3_old = "<<acc3_old<<endl;
    cout<<"acc4 = "<<acc4<<endl;
    cout<<"acc4_old = "<<acc4_old<<endl;
    cout<<"acc5 = "<<acc5<<endl;
    cout<<"acc5_old = "<<acc5_old<<endl;



    cout<<"phi_ex = "<<phi_ex<<endl;
    cout<<"acc_ex = "<<acc_ex<<endl;
    cout<<"acc_ex_old = "<<acc_ex_old<<endl;
    cout<<"jrk_ex = "<<jrk_ex<<endl;
    cout<<"jrk_ex_old = "<<jrk_ex_old<<endl;
    cout<<"acc2_ex = "<<acc2_ex<<endl;
    cout<<"acc2_ex_old = "<<acc2_ex_old<<endl;
    cout<<"acc3_ex = "<<acc3_ex<<endl;
    cout<<"acc3_ex_old = "<<acc3_ex_old<<endl;


    cout<<"ngb_index="<<ngb_index<<endl;
    cout<<"r2min="<<r2min<<endl;

    cout<<endl;
  }

  void dead(){
    mass = 0.0;
    dtime = 999999999999999999999999.9;
    radius = 0.0;
    rtide = 0.0;
    address = -1;
    type = DEAD;
    dtime = LARGE_DOUBLE;
  }


  double calc_Egr_2nd(){
    return 0.5*mass*dtime*(acc_PN_ex*vel 
			   + acc_PN_ex_old*vel_old);
  }
  double calc_Egr_4th(){
    return 0.5*mass*dtime*((acc_PN_ex*vel 
			    + acc_PN_ex_old*vel_old)
			   +dtime*n6*( jrk_PN_ex_old*vel_old 
				       + acc_PN_ex_old*acc_old 
				       - jrk_PN_ex*vel 
				       - acc_PN_ex*acc) );
  }
  double calc_Egr_6th(){
    return 0.5*mass*dtime*((acc_PN_ex*vel 
			    + acc_PN_ex_old*vel_old)
			   + dtime*n5*( (jrk_PN_ex_old*vel_old 
					 + acc_PN_ex_old*acc_old 
					 - jrk_PN_ex*vel 
					 - acc_PN_ex*acc) 
					+ dtime*n12*(acc2_PN_ex_old*vel_old 
						     + 2.0*jrk_PN_ex_old*acc_old 
						     + acc_PN_ex_old*jrk_old 
						     + acc2_PN_ex*vel
						     + 2.0*jrk_PN_ex*acc
						     + acc_PN_ex*jrk)));
  }

  void calc_radius(){
    radius = (Rsolsi/Pcsi)*get_radius_zams0(mass*4e6);
    /*
    cerr<<"mass="<<mass<<endl;
    cerr<<"mass*4e6="<<mass*4e6<<endl;
    cerr<<"radius="<<radius<<endl;
    cerr<<"radius*Pcsi/Rsolsi="<<radius*Pcsi/Rsolsi<<endl;
    cerr<<endl;
    */
    //radius = 2.254e-8; // Rsun[pc]
  }

};

#endif //PARTICLE_H


