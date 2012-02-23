#include"external_field.h"

static double SMBH_MASS = 1.0; 
static Vector3 SMBH_POS = 0.0;
static Vector3 SMBH_VEL = 0.0; 
static Vector3 SMBH_ACC = 0.0; 
static Vector3 SMBH_JRK = 0.0; 

static double SMBH_MASS_OLD = 1.0; 
static Vector3 SMBH_POS_OLD = 0.0;
static Vector3 SMBH_VEL_OLD = 0.0; 

//static double speed_of_light = 999999999999999999999999999999999.9; 
static double speed_of_light = 2285.604;
//static double speed_of_light = 7.94719414;
static double inv_c2 = 1.0/(speed_of_light*speed_of_light);
static double inv_c5 = inv_c2*inv_c2/speed_of_light;

static double EPS_SQ = 0.0;

static int calculate_postnewtonian = 1;
/////////////////////////////
///// forces from SMBH //////
/////////////////////////////

void calc_force_from_point_mass(Particle &prti,
				const int &mode){

  Vector3 posi, veli, acci, jrki;
  if(mode == 0){
    posi = prti.pos; 
    veli = prti.vel;
    acci = prti.acc;
    jrki = prti.jrk;
  }
  else if(mode == 1){
    posi = prti.pos_pre; 
    veli = prti.vel_pre;
    acci = prti.acc_pre;
    jrki = prti.jrk_pre;
    //cerr<<"posi="<<posi<<endl;
  }
  double eps2 = 0.0;
  double r2 = 0.0;

  /*
  calc_acc_jrk(posi,  veli, 
	       SMBH_MASS,  SMBH_POS,  SMBH_VEL,
	       eps2,
	       prti.acc_ex,  prti.jrk_ex,  prti.phi_ex, r2);
  */

  /*
  calc_acc_jrk_acc2(posi, veli, acci, 
		    SMBH_MASS, SMBH_POS, SMBH_VEL, SMBH_ACC, 
		    eps2,
		    prti.acc_ex, prti.jrk_ex, prti.acc2_ex, prti.phi_ex, r2);
  */

  if (!calculate_postnewtonian){
      calc_acc_jrk_acc2_acc3(posi, veli, acci, jrki, 
			 SMBH_MASS, SMBH_POS, SMBH_VEL, SMBH_ACC, SMBH_JRK,
			 eps2,
			 prti.acc_ex, prti.jrk_ex, prti.acc2_ex, prti.acc3_ex,  prti.phi_ex,
			 r2);



  } else {
      calc_acc0_jrk_acc2_acc3_PN0_PN1_PN25(prti.mass,
				       posi, veli, acci, jrki,
				       SMBH_MASS, SMBH_POS, SMBH_VEL, SMBH_ACC, SMBH_JRK,
				       prti.phi_ex,
				       prti.acc_ex, prti.jrk_ex, prti.acc2_ex, prti.acc3_ex, 
				       prti.acc_PN_ex, prti.jrk_PN_ex, prti.acc2_PN_ex,
				       prti.acc_PN_ex, prti.jrk_PN_ex, prti.acc2_PN_ex,
				       inv_c2, inv_c5, EPS_SQ);
  }
}

void calc_force_from_point_mass_to_array(Particle prt[],
					 int address[],
					 const int &Nip,
					 const int &mode){
  for(int i=0; i<Nip; i++){
    calc_force_from_point_mass(prt[address[i]], mode);
  }
}




void accrete_mass(const double &mass){
  SMBH_MASS += mass;
}

double calc_rtide_cu(const double &mass,
		     const double &radius){
		       
  return radius*radius*radius*(SMBH_MASS/mass);
}


void get_SMBH(double &mass, 
	      Vector3 &pos,
	      Vector3 &vel){
  mass = SMBH_MASS;
  pos = SMBH_POS;
  vel = SMBH_VEL;
}

void set_SMBH(const double &mass, 
	      const Vector3 &pos,
	      const Vector3 &vel){
  SMBH_MASS = mass;
  SMBH_POS = pos;
  SMBH_VEL = vel;
}

void copy_SMBH_OLD_TO_NEW(){
  SMBH_MASS = SMBH_MASS_OLD;
  SMBH_POS = SMBH_POS_OLD;
  SMBH_VEL = SMBH_VEL_OLD;
}

void copy_SMBH_NEW_TO_OLD(){
  SMBH_MASS_OLD = SMBH_MASS;
  SMBH_POS_OLD = SMBH_POS;
  SMBH_VEL_OLD = SMBH_VEL;
}

void set_speed_of_light(const double &_c){
  speed_of_light = _c;

}

int get_calculate_postnewtonian(int *value){
    *value = calculate_postnewtonian;
    return 0;
}
int set_calculate_postnewtonian(int value){
    calculate_postnewtonian = value;
    return 0;
}
