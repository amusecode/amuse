#ifndef FORCE_H
#define FORCE_H

#include"Matrix3.h"

inline double duncan_function_4th(const double &rij, 
				  const double &rout, 
				  const double &rin){
  //double K = 0.0;

    double x = (rij - rin)/(rout - rin);
    double K = 0.0;
    if(x <= 0.0){
	K = 0.0;
    }
    else if(1.0 <= x){
	K = 1.0;
    }
    else{
	double x2 = x*x;
	double x4 = x2*x2;
	double x5 = x4*x;
	double x6 = x4*x2;
	double x7 = x6*x;
	K = -20.0*x7 + 70.0*x6 - 84.0*x5 + 35.0*x4;
    }
    
    return K;
}

inline double duncan_function_4th_dot(const double &rij, 
				      const double &rijvij,
				      const double &rout, 
				      const double &rin){
    //double Kdot = 0.0;
    double x = (rij - rin)/(rout - rin);
    double xdot = rijvij/(rij*(rout - rin));
    double Kdot = 0.0;
    if(x <= 0.0){
	Kdot = 0.0;
    }
    else if(1.0 <= x){
	Kdot = 0.0;
    }
    else{
	double x2 = x*x;
	double x3 = x2*x;
	double x4 = x2*x2;
	double x5 = x4*x;
	double x6 = x4*x2;
	Kdot = (-140.0*x6 + 420.0*x5 - 420.0*x4 + 140.0*x3) * xdot;
    }
    return Kdot;
}



inline void pairwise_acc(const double massi,
			 const Vector3 &posi,
			 Vector3 &acci,
			 double &poti,
			 const double &massj,
			 const Vector3 &posj,
			 Vector3 &accj,
			 double &potj,
			 const double &eps2,
			 double &r2){
    Vector3 rij = posi - posj;
    r2 = rij * rij;
    double R2 = 1.0 / (r2 + eps2); 
    double R = sqrt(R2); 
    double R3 = R * R2;
    acci -= massj * R3 * rij;
    accj += massi * R3 * rij;
    poti -= massj * R;
    potj -= massi * R;
}


inline void calc_acc_duncan4(const Vector3 &posi,
			     Vector3 &acci,
			     double &poti,
			     const Vector3 &posj, 
			     const double &massj, 
			     const double &eps2, 
			     const double &rcut_in,  
			     const double &rcut_out,
			     double &r2){

    Vector3 rij = posi - posj;
    r2 = rij * rij;
    if(r2 <= rcut_out*rcut_out){
	double r = sqrt(r2);
	double r2_eps = r2 + eps2;
	double R = 1.0/sqrt(r2_eps);
	double R2 = R*R;
	double R3 = R2*R;
	double K = duncan_function_4th(r,  rcut_out,  rcut_in);
	Vector3 F0 = -massj * R3 * rij * (1.0-K);
	acci += F0;
	poti -= massj * R * (1.0-K);
    }
}


inline void calc_acc_jrk_duncan4(const Vector3 &posi,
				 const Vector3 &veli,
				 Vector3 &acci,
				 Vector3 &jrki,
				 double &poti,
				 const Vector3 &posj, 
				 const Vector3 &velj, 
				 const double &massj, 
				 const double &eps2, 
				 const double &rcut_in,  
				 const double &rcut_out,
				 double &r2){

    Vector3 rij = posi - posj;
    r2 = rij * rij;
    if(r2 <= rcut_out*rcut_out){
	Vector3 vij = veli - velj;
	double rijvij = rij * vij;
	double r = sqrt(r2);
	double r2_eps = r2 + eps2;
	double R = 1.0/sqrt(r2_eps);
	double R2 = R*R;
	double R3 = R2*R;
	double A = (rijvij)*R2;
	double K = duncan_function_4th(r,  rcut_out,  rcut_in);
	double Kdot = duncan_function_4th_dot(r,  rijvij, rcut_out,  rcut_in);
	Vector3 F0 = -massj * R3 * rij * (1.0-K);
	Vector3 F1 = -massj * R3 * vij * (1.0-K) - 3.0 * A * F0 + massj * R3 * rij * Kdot;
	acci += F0;
	jrki += F1;
	poti -= massj * R * (1.0-K);
    }
}


#endif //FORCE_H
