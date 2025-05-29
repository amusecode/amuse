#ifndef FORCE_H
#define FORCE_H

#include<cmath>
#include<iostream>
#include"Vector3.h"

using namespace std;

//////////////////////////
///// pure Newton   //////
//////////////////////////


inline void calc_acc_jrk(const Vector3 &posi, 
			 const Vector3 &veli, 
			 const double &massj, 
			 const Vector3 &posj, 
			 const Vector3 &velj, 
			 const double &eps2,
			 Vector3 &acc_out,
			 Vector3 &jrk_out,
			 double &phi_out,
			 double &r2_out){

  Vector3 rij = posi - posj;
  r2_out = rij*rij;
  double R2 = 1.0/(r2_out + eps2); 
  double R = sqrt(R2); 
  double R3 = R*R2; 

  phi_out -= massj*R;
  double mjR3 = massj*R3;
  Vector3 F0ij = -mjR3*rij; 
  acc_out += F0ij;

  Vector3 vij = veli - velj;
  double A1 = rij*vij*R2;
  Vector3 F1ij = -mjR3*vij - 3.0*A1*F0ij;  
  jrk_out += F1ij;

}

inline void calc_acc_jrk_acc2(const Vector3 &posi, 
			      const Vector3 &veli, 
			      const Vector3 &acci, 
			      const double &massj, 
			      const Vector3 &posj, 
			      const Vector3 &velj, 
			      const Vector3 &accj,
			      const double &eps2,
			      Vector3 &acc_out,
			      Vector3 &jrk_out,
			      Vector3 &acc2_out,
			      double &phi_out,
			      double &r2_out){
  Vector3 rij = posi - posj;
  r2_out = rij*rij;
  double R2 = 1.0/(r2_out + eps2); 
  double R = sqrt(R2); 
  double R3 = R*R2; 

  phi_out -= massj*R;
  double mjR3 = massj*R3;
  Vector3 F0ij = -mjR3*rij; 
  acc_out += F0ij;

  Vector3 vij = veli - velj;
  double A1 = rij*vij*R2;
  Vector3 F1ij = -mjR3*vij - 3.0*A1*F0ij;  
  jrk_out += F1ij;

  Vector3 aij = acci - accj;
  double A1_sq = A1*A1;
  double A2 = (vij*vij + rij*aij)*R2 + A1_sq;
  Vector3 F2ij = -mjR3*aij - 6.0*A1*F1ij - 3.0*A2*F0ij;
  acc2_out += F2ij;

}


inline void calc_acc_jrk_acc2_acc3(const Vector3 &posi, 
				   const Vector3 &veli, 
				   const Vector3 &acci, 
				   const Vector3 &jrki, 
				   const double &massj, 
				   const Vector3 &posj, 
				   const Vector3 &velj, 
				   const Vector3 &accj,
				   const Vector3 &jrkj,
				   const double &eps2,
				   Vector3 &acc_out,
				   Vector3 &jrk_out,
				   Vector3 &acc2_out,
				   Vector3 &acc3_out,
				   double &phi_out,
				   double &r2_out){
  
  Vector3 rij = posi - posj;
  r2_out = rij*rij;
  double R2 = 1.0/(r2_out + eps2); 
  double R = sqrt(R2); 
  double R3 = R*R2; 

  phi_out -= massj*R;
  double mjR3 = massj*R3;
  Vector3 F0ij = -mjR3*rij; 
  acc_out += F0ij;

  Vector3 vij = veli - velj;
  double A1 = rij*vij*R2;
  Vector3 F1ij = -mjR3*vij - 3.0*A1*F0ij;  
  jrk_out += F1ij;

  Vector3 aij = acci - accj;
  double A1_sq = A1*A1;
  double A2 = (vij*vij + rij*aij)*R2 + A1_sq;
  Vector3 F2ij = -mjR3*aij - 6.0*A1*F1ij - 3.0*A2*F0ij;
  acc2_out += F2ij;

  Vector3 jij = jrki - jrkj;
  double A3 = (3.0*vij*aij + rij*jij)*R2 + A1*(3.0*A2 - 4.0*A1_sq);
  Vector3 F3ij = -mjR3*jij - 9.0*A1*F2ij - 9.0*A2*F1ij - 3.0*A3*F0ij;
  acc3_out += F3ij;

}



/////////////////////
///// calc PN ///////
/////////////////////

// PN to acc2
// Newton to acc3
inline void calc_acc0_jrk_acc2_acc3_PN0_PN1_PN25(const double &massi,
						 const Vector3 &posi,
						 const Vector3 &veli,
						 const Vector3 &acci,
						 const Vector3 &jrki,
						 const double &massj,
						 const Vector3 &posj,
						 const Vector3 &velj,
						 const Vector3 &accj,
						 const Vector3 &jrkj,
						 double &phi,
						 Vector3 &accPN0,
						 Vector3 &jrkPN0,
						 Vector3 &acc2PN0,
						 Vector3 &acc3PN0,
						 Vector3 &accPN1,
						 Vector3 &jrkPN1,
						 Vector3 &acc2PN1,
						 Vector3 &accPN25,
						 Vector3 &jrkPN25,
						 Vector3 &acc2PN25,
						 const double c2_inv,
						 const double c5_inv,
						 const double &eps2 = 0.0){

  Vector3 rij = posi - posj;
  Vector3 vij = veli - velj;
  Vector3 aij = acci - accj;
  Vector3 jij = jrki - jrkj;
  double R2 = 1.0/((rij*rij) + eps2); 
  double R = sqrt(R2); 
  double R3= R*R2;
  double A0 = (rij*vij)*R2;
  double B0 = (vij*vij +  rij*aij)*R2 + A0*A0;
  double C0 = (3.0*vij*aij + rij*jij)*R2 + A0*(3.0*B0 - 4.0*A0*A0);
  double mjR3 = massj*R3;
  Vector3 Aij = -mjR3*rij;
  Vector3 Jij = -mjR3*vij - 3.0*A0*Aij;
  Vector3 Sij = -mjR3*aij - 6.0*A0*Jij - 3.0*B0*Aij;
  Vector3 Cij = -mjR3*jij - 9.0*A0*Sij - 9.0*B0*Jij - 3.0*C0*Aij;
  accPN0 += Aij;
  jrkPN0 += Jij;
  acc2PN0 += Sij;
  acc3PN0 += Cij;
  phi -= massj*R;


  // PN 1st
  double A1 = c2_inv*(5.0*massi*massj + 4.0*massj*massj);
  double B1 = c2_inv*1.5*massj;
  double C1 = c2_inv*massj;
  double rvj = rij*velj;
  double rvj_sq = rvj*rvj;
  double fb1 = rvj_sq;
  double vi_sq = veli*veli;
  double vj_sq = velj*velj;
  double vivj = veli*velj;
  double fc1 = (-vi_sq + 4.0*vivj - 2.0*vj_sq);
  double R4 = R2*R2;
  double R5 = R3*R2;
  double F1 = 
    A1*R4 + 
    B1*R5*fb1 + 
    C1*R3*fc1;

  double D1 = c2_inv*4.0*massj;
  double E1 = c2_inv*-3.0*massj;
  double rvi = rij*veli;
  double G1 = D1*R3*rvi + E1*R3*rvj;

  accPN1 += F1*rij + G1*vij;



  double vvj = vij*velj;
  double raj = rij*accj;
  double fb1dot = 2.0*(vvj + raj);
  double viai = veli*acci;
  double vjaj = velj*accj;
  double aivj = acci*velj;
  double viaj = veli*accj;
  double fc1dot = -2.0*viai + 4.0*(aivj + viaj) - 4.0*vjaj;
  double R6 = R4*R2;
  double R7 = R4*R3;
  double rv = rij*vij;
  double F1dot = 
    A1*( -4.0*R6*rv ) + 
    B1*( -5.0*R7*rv*fb1 + R5*fb1dot )+
    C1*( -3.0*R5*rv*fc1 + R3*fc1dot );

  double vvi = vij*veli;
  double rai = rij*acci;
  double G1dot = 
    D1*( -3.0*R5*rv*rvi + R3*(vvi+rai) ) +
    E1*( -3.0*R5*rv*rvj + R3*(vvj+raj) );

  jrkPN1 += F1dot*rij +F1*vij + G1dot*vij + G1*aij;

  double avj = aij*velj;
  double vaj = vij*accj;
  double rjj = rij*jrkj;
  double fb1dot2 = 2.0*(avj + 2.0*vaj + rjj);
  double ai_sq = acci*acci;
  double viji = veli*jrki;
  double jivj = jrki*velj;
  double aiaj = acci*accj;
  double vijj = veli*jrkj;
  double aj_sq = accj*accj;
  double vjjj = velj*jrkj;
  double fc1dot2 = -2.0*(ai_sq + viji) + 4.0*(jivj + 2.0*aiaj + vijj) - 4.0*(aj_sq + vjjj);
  double ra = rij*aij;
  double v_sq = vij*vij;
  double rvdot = ra + v_sq;
  double rv_sq = rv*rv;
  double R8 = R6*R2;
  double R9 = R6*R3;
  double F1dot2 = 
    A1*( 24.0*R8*rv_sq - 4.0*R6*rvdot ) +
    B1*( 35.0*R9*rv_sq*fb1 - 5.0*R7*(rvdot*fb1 + rv*fb1dot) +
	-5.0*R7*rv*fb1dot + R5*fb1dot2) +
    C1*( 15.0*R7*rv_sq*fc1 - 3.0*R5*(rvdot*fc1 + rv*fc1dot) +
	-3.0*R5*rv*fc1dot + R3*fc1dot2);

  double avi = aij*veli;
  double vai = vij*acci;
  double rji = rij*jrki;
  double G1dot2 = 
    D1*( (15.0*R7*rv_sq*rvi - 3.0*R5*(rvdot*rvi + rv*(vvi+rai))) +
	(-3.0*R5*rv*(vvi+rai) + R3*(avi + vai + vai + rji) ) ) +
    E1*( (15.0*R7*rv_sq*rvj - 3.0*R5*(rvdot*rvj + rv*(vvj+raj))) +
	(-3.0*R5*rv*(vvj+raj) + R3*(avj + vaj + vaj + rjj) ) );

  acc2PN1 += 
    F1dot2*rij + 2.0*F1dot*vij + F1*aij + 
    G1dot2*vij + 2.0*G1dot*aij + G1*jij;



  // PN 2.5th
  static double N208_15 = 208.0/15.0;
  static double N24_5 = 24.0/5.0;
  static double N12_5 = 12.0/5.0;
  static double N8_5 = 8.0/5.0;
  static double N32_5 = 32.0/5.0;
  static double N4_5 = 4.0/5.0;
  double A25 = c5_inv*massi*massj*(N208_15*massj - N24_5*massi);
  double B25 = c5_inv*massi*massj*N12_5;
  double C25 = c5_inv*massi*massj*(N8_5*massi - N32_5*massj);
  double D25 = c5_inv*-N4_5*massi*massj;

  double fb25 = rv*v_sq;
  double F25 = A25*R6*rv + B25*R5*fb25;

  double G25 = C25*R4 + D25*R3*v_sq;

  accPN25 += F25*rij + G25*vij;


  double va = vij*aij;
  double fb25dot = rvdot*v_sq + rv*2.0*va;
  double F25dot = 
    A25*(-6.0*R8*rv_sq + R6*rvdot) +
    B25*(-5.0*R7*rv*fb25 + R5*fb25dot);

  double G25dot = 
    C25*( -4.0*R6*rv ) +
    D25*( -3.0*R5*rv*v_sq + R3*2.0*va );

  jrkPN25 += F25dot*rij + F25*vij + G25dot*vij + G25*aij;


  double vj = vij*jij;
  double a_sq = aij*aij;
  double rj = rij*jij;
  double rvdot2 = rj + 3.0*va;
  double fb25dot2 = rvdot2*v_sq + rvdot*2.0*va + + 2.0*(rvdot*va+rv*(a_sq+vj));
  double R10 = R8*R2;
  double F25dot2 = 
    A25*( 48.0*R10*rv*rv_sq - 6.0*R8*(v_sq + ra) - 6.0*R8*rv*rvdot + R6*rvdot2 ) +
    B25*( 35.0*R9*rv_sq*fb25 - 5.0*R7*(rvdot*fb25 + rv*fb25dot) - 5.0*R7*rv*fb25dot + R5*fb25dot2);

  double G25dot2 = 
    C25*( 24.0*R8*rv_sq - 4.0*R6*rvdot ) + 
    D25*( 15.0*R7*rv_sq*v_sq - 3.0*R5*( rvdot*v_sq + rv*2.0*va) +
	  -6.0*R5*rv*va + 2.0*R3*(a_sq+vj) );
	   
  acc2PN25 += 
    F25dot2*rij + 2.0*F25dot*vij + F25*aij + 
    G25dot2*vij + 2.0*G25dot*aij + G25*jij;



}


/*
inline void calc_PN0_PN1(const int &idi,
			 const double &massi,
			 const Vector3 &posi,
			 const Vector3 &veli,
			 const Vector3 &acci,
			 const int &idj,
			 const double &massj,
			 const Vector3 &posj,
			 const Vector3 &velj,
			 const Vector3 &accj,
			 Particle prt[],
			 const int &Nj,
			 double &phi,
			 Vector3 &accPN0,
			 Vector3 &accPN1,
			 const double c2_inv){

  Vector3 rij = posi - posj;
  Vector3 vij = veli - velj;
  double R2 = 1.0/(rij*rij);
  double R = sqrt(R2); 
  double R3= R*R2;
  double A0 = (rij*vij)*R2;
  double mjR3 = massj*R3;
  Vector3 Aij = -mjR3*rij;
  Vector3 Jij = -mjR3*vij - 3.0*A0*Aij;
  accPN0 += Aij;
  jrkPN0 += Jij;
  phi -= massj*R;


  double C0;
  for(int j=0; j<Nj; j++){
    if(idi == prt[j].index){continue;}
    Vector3 r_tmp = posi - prt[j].pos;
    double R_tmp = 1.0/sqrt(r_tmp*r_tmp);
    C0 += prt[j].mass*R_tmp;
  }

  double C1;
  for(int j=0; j<Nj; j++){
    if(idj == prt[j].index){continue;}
    Vector3 r_tmp = posj - prt[j].pos;
    double R_tmp = 1.0/sqrt(r_tmp*r_tmp);
    C1 += prt[j].mass*R_tmp;
  }
  Vector3 H1 = c2_inv*massj*R3(4.0*C0 + C1)*rij;

  Vector3 I1 = c2_inv*massj*R3*0.5*(rij*accj)*rij;

  Vector3 J1 = c2_inv*massj*3.5*R*accj


  // PN 1st
  double A1 = c2_inv*(5.0*massi*massj + 4.0*massj*massj);
  double B1 = c2_inv*1.5*massj;
  double C1 = c2_inv*massj;
  double rvj = rij*velj;
  double rvj_sq = rvj*rvj;
  double fb1 = rvj_sq;
  double vi_sq = veli*veli;
  double vj_sq = velj*velj;
  double vivj = veli*velj;
  double fc1 = (-vi_sq + 4.0*vivj - 2.0*vj_sq);
  double R4 = R2*R2;
  double R5 = R3*R2;
  double F1 = 
    A1*R4 + 
    B1*R5*fb1 + 
    C1*R3*fc1;

  double D1 = c2_inv*4.0*massj;
  double E1 = c2_inv*-3.0*massj;
  double rvi = rij*veli;
  double G1 = D1*R3*rvi + E1*R3*rvj;

  accPN1 += F1*rij + G1*vij;



  double vvj = vij*velj;
  double raj = rij*accj;
  double fb1dot = 2.0*(vvj + raj);
  double viai = veli*acci;
  double vjaj = velj*accj;
  double aivj = acci*velj;
  double viaj = veli*accj;
  double fc1dot = -2.0*viai + 4.0*(aivj + viaj) - 4.0*vjaj;
  double R6 = R4*R2;
  double R7 = R4*R3;
  double rv = rij*vij;
  double F1dot = 
    A1*( -4.0*R6*rv ) + 
    B1*( -5.0*R7*rv*fb1 + R5*fb1dot )+
    C1*( -3.0*R5*rv*fc1 + R3*fc1dot );

  double vvi = vij*veli;
  double rai = rij*acci;
  double G1dot = 
    D1*( -3.0*R5*rv*rvi + R3*(vvi+rai) ) +
    E1*( -3.0*R5*rv*rvj + R3*(vvj+raj) );

  jrkPN1 += F1dot*rij +F1*vij + G1dot*vij + G1*aij;

  double avj = aij*velj;
  double vaj = vij*accj;
  double rjj = rij*jrkj;
  double fb1dot2 = 2.0*(avj + 2.0*vaj + rjj);
  double ai_sq = acci*acci;
  double viji = veli*jrki;
  double jivj = jrki*velj;
  double aiaj = acci*accj;
  double vijj = veli*jrkj;
  double aj_sq = accj*accj;
  double vjjj = velj*jrkj;
  double fc1dot2 = -2.0*(ai_sq + viji) + 4.0*(jivj + 2.0*aiaj + vijj) - 4.0*(aj_sq + vjjj);
  double ra = rij*aij;
  double v_sq = vij*vij;
  double rvdot = ra + v_sq;
  double rv_sq = rv*rv;
  double R8 = R6*R2;
  double R9 = R6*R3;
  double F1dot2 = 
    A1*( 24.0*R8*rv_sq - 4.0*R6*rvdot ) +
    B1*( 35.0*R9*rv_sq*fb1 - 5.0*R7*(rvdot*fb1 + rv*fb1dot) +
	-5.0*R7*rv*fb1dot + R5*fb1dot2) +
    C1*( 15.0*R7*rv_sq*fc1 - 3.0*R5*(rvdot*fc1 + rv*fc1dot) +
	-3.0*R5*rv*fc1dot + R3*fc1dot2);

  double avi = aij*veli;
  double vai = vij*acci;
  double rji = rij*jrki;
  double G1dot2 = 
    D1*( (15.0*R7*rv_sq*rvi - 3.0*R5*(rvdot*rvi + rv*(vvi+rai))) +
	(-3.0*R5*rv*(vvi+rai) + R3*(avi + vai + vai + rji) ) ) +
    E1*( (15.0*R7*rv_sq*rvj - 3.0*R5*(rvdot*rvj + rv*(vvj+raj))) +
	(-3.0*R5*rv*(vvj+raj) + R3*(avj + vaj + vaj + rjj) ) );

  acc2PN1 += 
    F1dot2*rij + 2.0*F1dot*vij + F1*aij + 
    G1dot2*vij + 2.0*G1dot*aij + G1*jij;

}
*/
#endif //FORCE_H
