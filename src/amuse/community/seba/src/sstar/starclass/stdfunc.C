//#INCLUDE "constants.h"
//#include "stdinc.h"

#include "stdfunc.h"

#ifndef TOOLBOX

// Standard functions for various implementations.

// Tidal energy dissipation during close encounter.
// Portegies Zwart, SF & Meinen AT, 1992, AA 280, 174.
// for polytrope n=1.5
real tf2_energy_diss(const real eta, const stellar_type type) {

//		Polytropic index 1.5 is default
      real coeff2[] = {-0.397, 1.678, 1.277, -12.42, 9.446, -5.550};

/*
      switch (type) {
         case Brown_Dwarf:		//polytropic index = 3
         case Horizontal_Branch:	//fully radiative stars.
         case Helium_Star:
         case Wolf_Rayet:
         case White_Dwarf:
         case Neutron_Star:
         case Black_Hole:
              coeff2[0] =  -1.124;
              coeff2[1] =   0.877;
              coeff2[2] = -13.37;
              coeff2[3] =  21.55;
              coeff2[4] = -16.48;
              coeff2[5] =   4.124;
              break;
         case Main_Sequence:		//polytropic index = 2
         case Hertzsprung_Gap:		
              coeff2[0] =  -0.517;
              coeff2[1] =  -0.906;
              coeff2[2] =  23.88;
              coeff2[3] = -93.49;
              coeff2[4] = 112.3;
              coeff2[5] = -44.15;
              break;
         case Sub_Giant:		//polytropic index = 1.5
         case Super_Giant:              //fully convective stars.
              coeff2[0] =  -0.397;
              coeff2[1] =   1.678;
              coeff2[2] =   1.277;
              coeff2[3] = -12.42;
              coeff2[4] =   9.446;
              coeff2[5] =  -5.550;
              break;
         default:		//SAFETY polytropic index = 1.5
              coeff2[0] =  -0.397;
              coeff2[1] =   1.678;
              coeff2[2] =   1.277;
              coeff2[3] = -12.42;
              coeff2[4] =   9.446;
              coeff2[5] =  -5.550;
      }
*/

      real y = log10(eta);
      y = min(1.0, max(y, 0.0));
      real logT = ((((coeff2[5]*y + coeff2[4])*y + coeff2[3])*y
                   + coeff2[2])*y + coeff2[1])*y + coeff2[0];
//      logT = min(0, max(logT, -10));

      return pow(10., logT);
   }

// Tidal energy dissipation during close encounter.
// Portegies Zwart, SF & Meinen AT, 1992, AA 280, 174.
// for polytrope n=1.5
real tf3_energy_diss(const real eta, const stellar_type type) {

//		Initialization with polytrope index 1.5
      real coeff3[] = {-0.909, 1.574, 12.37, -57.40, 80.10, -46.43};
/*
      switch (type) {
         case Brown_Dwarf:              //polytropic index = 3
         case Helium_Star:		//fully convecitve stars.
         case Wolf_Rayet:
         case White_Dwarf:
         case Neutron_Star:
         case Black_Hole:
              coeff3[0] =  -1.703;
              coeff3[1] =   2.653;
              coeff3[2] = -14.34;
              coeff3[3] =  12.85;
              coeff3[4] =  -0.492;
              coeff3[5] =  -3.600;
//              Do nothing just polytrope index 1.5.
              break;
         case Main_Sequence:		//polytropic index = 2
         case Hertzsprung_Gap:         
         case Horizontal_Branch:
              coeff3[0] =   -1.040;
              coeff3[1] =   -1.354;
              coeff3[2] =   37.64;
              coeff3[3] = -139.9;
              coeff3[4] =  168.2;
              coeff3[5] =  -66.53;
              break;
         case Sub_Giant:                //polytropic index = 1.5
         case Super_Giant:		//Convective stars.
              coeff3[0] =  -0.909;
              coeff3[1] =   1.574;
              coeff3[2] =  12.37;
              coeff3[3] = -57.40;
              coeff3[4] =  80.10;
              coeff3[5] = -46.43;
              break;
         default:               //SAFETY polytropic index = 1.5
              coeff3[0] =  -0.909;
              coeff3[1] =   1.574;
              coeff3[2] =  12.37;
              coeff3[3] = -57.40;
              coeff3[4] =  80.10;
              coeff3[5] = -46.43;
      }
*/

      real y = log10(eta);
      y = min(1.0, max(y, 0.0));
      real logT = ((((coeff3[5]*y + coeff3[4])*y + coeff3[3])*y
                   + coeff3[2])*y + coeff3[1])*y + coeff3[0];
//      logT = min(0, max(logT, -10));

      return pow(10., logT);
   }

//              Standard lineair interpolation routine.
real lineair_interpolation(const real x,
			   const real x1, const real x2,
			   const real y1, const real y2) {

        real a = (y2-y1)/(x2-x1);
	real b = y1 - a*x1;

	real y = a*x + b;
	return y;
}

//		Super nova utilities.

//	Check!!
real post_sn_cm_velocity(const real a_init,     const real e_init,
                         const real separation,
                         const real m1_0, const real m2_0,
                         const real m1, const real m2,
                         const real v_kick,
                         const real theta,      const real phi) {

//	Check!!
     real v_orb = sqrt((cnsts.physics(G)*cnsts.parameters(solar_mass)/cnsts.parameters(solar_radius))*(m1_0+ m2_0)*(2/separation - 1/a_init));
//          v_orb /= cnsts.physics(km_per_s);
     real vr_k = cnsts.physics(km_per_s)*v_kick/v_orb;

     real mu = m1_0*m2_0/(m1_0+m2_0);
          mu *= cnsts.parameters(solar_mass);
     real dm = (m1_0+m2_0) - (m1+m2);
          dm *= cnsts.parameters(solar_mass);
     real v_cm = v_orb/(cnsts.parameters(solar_mass)*(m1 + m2));
     v_cm *= sqrt(pow(mu*dm/(cnsts.parameters(solar_mass)*m1_0), 2) - 2*(mu*dm*m1/m1_0)
           * vr_k*sin(theta)*cos(phi) + pow(cnsts.parameters(solar_mass)*m1*vr_k, 2));
     v_cm /= cnsts.physics(km_per_s);

//	Check!!
     return v_cm;
   }


real post_sn_semi_major_axis(const real a_init,     const real e_init,
                             const real separation,
                             const real m1_0, const real m2_0,
                             const real m1, const real m2,
                             const real v_kick,
                             const real theta,      const real phi) {

//              SPZ's method.
     real mu = (m1 + m2)/(m1_0 + m2_0);
     real v_orb = sqrt((cnsts.physics(G)*cnsts.parameters(solar_mass)/cnsts.parameters(solar_radius))*(m1_0+ m2_0)*(2/separation - 1/a_init));
          v_orb /= cnsts.physics(km_per_s);
     real vr_k = v_kick/v_orb;
     real vr2 = 1 + vr_k*vr_k + 2*vr_k*sin(theta)*cos(phi);
     real a_new = 1/(2/separation - (vr2/mu)*(2/separation - 1/a_init));
     real alpha = a_new/a_init;

     real epsilon = (1/(mu*alpha))*(vr2 - vr_k*vr_k*pow(cos(theta), 2));
     real e_new = sqrt(1 - (1-e_init*e_init)*epsilon);


//     cerr<<"final: "<< alpha*a_init<<" "<<e_new<<endl;

     return a_new;
}

real post_sn_eccentricity(const real a_init,     const real e_init,
                          const real separation,
                          const real m1_0, const real m2_0,
                          const real m1, const real m2,
                          const real v_kick,
                          const real theta,      const real phi) {

     real mu = (m1 + m2)/(m1_0 + m2_0);
     real v_orb = sqrt((cnsts.physics(G)*cnsts.parameters(solar_mass)/cnsts.parameters(solar_radius))*(m1_0+ m2_0)*(2/separation - 1/a_init));
          v_orb /= cnsts.physics(km_per_s);
     real vr_k = v_kick/v_orb;
     real vr2 = 1 + vr_k*vr_k + 2*vr_k*sin(theta)*cos(phi);
     real a_new = 1/(2/separation - (vr2/mu)*(2/separation - 1/a_init));
     real alpha = a_new/a_init;

     real epsilon = (1/(mu*alpha))*(vr2 - vr_k*vr_k*pow(cos(theta), 2));
     real e_new = sqrt(1 - (1-e_init*e_init)*epsilon);

     return e_new;
}

real kinetic_energy(const real mass, const real velocity) {
             real Mkm_s2 = cnsts.parameters(solar_mass)*cnsts.physics(km_per_s)*cnsts.physics(km_per_s);
             real k = 0.5*Mkm_s2*mass*velocity*velocity;
             return k;
        }
real potential_energy(const real sep, const real m1, const real m2) {

       real GM2_R = cnsts.physics(G)*cnsts.parameters(solar_mass)*cnsts.parameters(solar_mass)/cnsts.parameters(solar_radius);
       real u = GM2_R*m1*m2/sep;
       return -u;
   }

real velocity_at_infinity(const real vel, 
                          const real sep, 
                          const real m_prim, 
                          const real m_sec) {

      real e_kin = kinetic_energy(m_prim+m_sec, vel);
      real e_pot = potential_energy(sep, m_prim, m_sec);
      real m_tot = cnsts.parameters(solar_mass)*(m_prim + m_sec);
      real velocity = 0;
      if (e_kin>abs(e_pot)) 
         velocity = sqrt((e_kin + e_pot)/m_tot);

      return velocity/cnsts.physics(km_per_s);
   }

real random_angle(const real min, const real max) {

     extern real randinter(real, real);

     real rnd = randinter(min, max);

//     cerr << "rdn: " << rnd<< endl;
     return rnd;
}

real random_eccentric_anomaly(const real ecc) {

//      Solve Kepler equation by iterating: M = E - e sin E
//      Lit.: Sterne, T.E., 1960, An introdiction to Celestial
//            Mechanics, p. 13-14

/*
     real mean_anomaly = random_angle(0, cnsts.math.two_pi);
//cerr<<" mean_anomaly="<<mean_anomaly<<endl;
     real ecc_anomaly = mean_anomaly;     // first guess.

//              Interpolate for better result.
         ecc_anomaly = mean_anomaly
                           - (ecc_anomaly - ecc*sin(ecc_anomaly)
                              - mean_anomaly)
                           / (1 - ecc*cos(ecc_anomaly));

*/

      real mean_anomaly = random_angle(0, cnsts.mathematics(two_pi));
//      cerr<<" mean_anomaly="<<mean_anomaly<<endl;

//		First guess.
     real ecc_anomaly = mean_anomaly
                            + cnsts.mathematics(two_pi)*(ecc*sin(mean_anomaly)
                            + 0.5*ecc*ecc*sin(2*mean_anomaly));
//     cerr<<"ecc_anomaly=" <<ecc_anomaly<<endl;
     real m0 = ecc_anomaly - ecc*sin(ecc_anomaly);
     real de0 = (mean_anomaly-m0)
              / (1 - ecc*cos(ecc_anomaly));

     ecc_anomaly += de0;
     real m1, de1;
     do {
        m1 = ecc_anomaly - ecc*sin(ecc_anomaly);

        de1 = (mean_anomaly-m1)
            / (1 - ecc*cos(ecc_anomaly));
        
        ecc_anomaly += de1;

     }
     while(de1>=0.001);

//cerr<<" => ecc_anomaly="<<ecc_anomaly<<endl;
     return ecc_anomaly;
}

real eccentric_anomaly(const real ecc, const real mean_anomaly) {

//      Solve Kepler equation by iterating: M = E - e sin E
//      Lit.: Sterne, T.E., 1960, An introdiction to Celestial
//            Mechanics, p. 13-14

//      cerr<<" mean_anomaly="<<mean_anomaly<<endl;

//              First guess.
     real ecc_anomaly = mean_anomaly
                            + cnsts.mathematics(two_pi)*(ecc*sin(mean_anomaly)
                            + 0.5*ecc*ecc*sin(2*mean_anomaly));
//     cerr<<"ecc_anomaly=" <<ecc_anomaly<<endl;
     real m0 = ecc_anomaly - ecc*sin(ecc_anomaly);
     real de0 = (mean_anomaly-m0)
              / (1 - ecc*cos(ecc_anomaly));

     ecc_anomaly += de0;
     real m1, de1;
     do {
        m1 = ecc_anomaly - ecc*sin(ecc_anomaly);

        de1 = (mean_anomaly-m1)
            / (1 - ecc*cos(ecc_anomaly));

        ecc_anomaly += de1;

     }
     while(de1>=0.001);

//cerr<<" => ecc_anomaly="<<ecc_anomaly<<endl;
     return ecc_anomaly;
}


real random_separation(const real semi, const real ecc) {
//cerr<<"random_separation(a="<<semi<<", e="<<ecc<<")"<<endl;

     real ecc_anomaly = random_eccentric_anomaly(ecc);

//cerr<<"ecc_anomaly=" << ecc_anomaly;
     real separation = semi*(1-ecc*cos(ecc_anomaly));
//cerr<<" ==> separation="<<separation<<endl;

     return separation;
}

/*
real maxwellian_distribution() {

        static int iset=0;
        static real gset;
        real fac,rsq,v1,v2;

        if  (iset == 0) {
                do {
                        v1=2.0*randinter(0,1)-1.0;
                        v2=2.0*randinter(0,1)-1.0;
                        rsq=v1*v1+v2*v2;
                } while (rsq >= 1.0 || rsq == 0.0);
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
   }
*/

real maxwellian(const real velocity, const real v_disp) {


      real Mxwmax = 4/(sqrt(cnsts.mathematics(pi))*exp(1));
      real u2 = pow(velocity/v_disp, 2);
      real prob = (4./sqrt(cnsts.mathematics(pi))) * u2 / exp(u2);
      prob /= Mxwmax;

      return prob;
   }

real random_maxwellian_velocity(const real v_disp) {

     real rnd_gauss, v[3];
     real final_vel=0;
     real mean_vel = 0.;
     real v1d_disp = v_disp/sqrt(3.);

     for (int i=0; i<3; i++) {
         rnd_gauss = gauss();
         v[i] = mean_vel + v1d_disp*rnd_gauss;  // km/s
         final_vel += v[i]*v[i];
     }

     return sqrt(final_vel);
  }
  
  


//   Hobbs, Lorimer, Lyne & Kramer, 2005, 360, 974 - maxwellian
real random_hobbs_velocity() {
//   cerr<<"random_hobbs_velocity()"<<endl;
   return random_maxwellian_velocity(265*sqrt(3.));  
}


//    Arzoumanian ea 2002, 568, 289 - combination of two maxwellians
real random_arzoumanian_velocity() {
//   cerr<<"random_arzoumanian_velocity()"<<endl;
   real prob = randinter(0., 1.);
   if (prob < 0.4){
       return random_maxwellian_velocity(90*sqrt(3.));
   }else {
       return random_maxwellian_velocity(500*sqrt(3.));
   }
}


//    Verbunt, Igoshev & Cator, 2017, 608, 57 - combination of two maxwellians
real random_verbunt_velocity() {
//   cerr<<"random_verbunt_velocity()"<<endl;
   real prob = randinter(0., 1.);
   if (prob < 0.42){
       return random_maxwellian_velocity(75*sqrt(3.));
   }else {
       return random_maxwellian_velocity(316*sqrt(3.));
   }
}

  

real paczynski_distribution(const real velocity, const real v_disp) {

      real prob = (4./cnsts.mathematics(pi))/pow(1+pow(velocity/v_disp, 2.), 2.);

//		Normaization to unity for its maximum.
      prob /= (4./cnsts.mathematics(pi));

      return prob;
   }


real random_paczynski_velocity(const real v_disp) {
//	Velocity distribution used by Paczynski, B., 1990, ApJ 348, 485.
//	with a dispersion of 270 km/s.
//	Phinney uses the same distribution with a dispersion of 600 km/s.
//	The distribution:
//	P(u)du = \frac{4}{\pi} \cdot \frac{du}{(1+u^2)^2},
//	u=v/v_d,
// 	v_d is the dispersion velocity = 270 or 600 km/s respectively.

   cerr<<"random_paczynski_velocity()"<<endl;

//	The safe and clumsy method
	
      real prob, velocity, dist_value;
      real max_distr = 4./cnsts.mathematics(pi);
      real v_max = 4;
      do {
         prob = max_distr*randinter(0., 1.);
         velocity = v_max*randinter(0., 1.);
         dist_value = max_distr/pow(1+pow(velocity, 2.), 2.);
     }
     while(prob>dist_value);

     return velocity*v_disp;
   }
 
real cross_section(const real m_prim, const real m_sec,
                   const real r_min,  const real v_rel) {

      real cgs_const = cnsts.physics(G)*cnsts.parameters(solar_mass)/(cnsts.parameters(solar_radius)*pow(cnsts.physics(km_per_s), 2));
      real gr_foc = cgs_const*(m_prim+m_sec)/(r_min*pow(v_rel, 2));
      real sigma = cnsts.mathematics(pi)*pow(r_min*cnsts.parameters(solar_radius), 2) * (1 + gr_foc);
      sigma /= pow(cnsts.parameters(solar_radius), 2);

      return sigma;     // in R_SUN^2
   }

real eddington_limit(const real radius,
		     const real dt,
		     const real mu) {
  return cnsts.parameters(solar_radius)*radius*dt*mu*1.5e-08;
}

real gf_velocity_distribution_extremum(const real m_prim, const real m_sec,
                                     const real r_min, const real v_disp) {

      real abc_A = 3*cnsts.parameters(solar_radius)*r_min/pow(cnsts.physics(km_per_s)*v_disp, 2); 
      real abc_B = 6*cnsts.physics(G)*cnsts.parameters(solar_mass)*(m_prim+m_sec)/pow(cnsts.physics(km_per_s)*v_disp, 2) - 4*cnsts.parameters(solar_radius)*r_min; 
      real abc_C = -4*cnsts.physics(G)*cnsts.parameters(solar_mass)*(m_prim+m_sec);
      real discriminant = pow(abc_B, 2) - 4*abc_A*abc_C;
      real kappa1 = (-abc_B + sqrt(discriminant))/(2*abc_A);
      real kappa2 = (-abc_B - sqrt(discriminant))/(2*abc_A);
      real extremum1 = sqrt(abs(kappa1))/cnsts.physics(km_per_s);
      // real extremum2 = sqrt(abs(kappa2))/cnsts.physics(km_per_s);
     
//      cerr << " maxima: "<< extremum1 << " " << extremum2 << endl;
//      return max(v_disp, min(extremum1, extremum2));
      return extremum1;
   }

//The relative velocity between the two encountering
//objects must be chosen from the distribution $P(v)$ given by:
//        P(v) dv = \sigma v {\cal F}(v)dv.
//here {\cal F}(v)dv is the Maxwellian.
real gravitational_focussed_velocity(const real m_prim, const real m_sec,
                                     const real r_min,  const real vel,
				     const real v_disp) {

      real konst = 4*sqrt(cnsts.mathematics(pi))*pow(3./(2*pow(v_disp*cnsts.physics(km_per_s), 2)), 3./2.);
      real GMm = cnsts.physics(G)*cnsts.parameters(solar_mass)*(m_prim+m_sec);
      real e_x = 1./exp((3./2.)*pow(vel/v_disp, 2));
      real velocity = vel*cnsts.physics(km_per_s);
      real s1 = pow(velocity, 3)*pow(cnsts.parameters(solar_radius)*r_min, 2);
      real s2 = 2*GMm*velocity*cnsts.parameters(solar_radius)*r_min;

      real pvs = konst*e_x*(s1 + s2);

      return pvs;

   }

real random_focussed_maxwellian_velocity(const real m_prim, const real m_sec,
                                     const real r_min,  const real v_disp,
                                     const real v_esc) {

      real v_extremum = gf_velocity_distribution_extremum(m_prim, m_sec,
                                                        r_min, v_disp);
      real max_prob = gravitational_focussed_velocity(m_prim, m_sec,
                                                   r_min,  v_extremum, v_disp);

      real pprob, velocity, prob;
      do {
         pprob    = randinter(0, 1.25);
         velocity = randinter(0, v_esc);
         prob     = gravitational_focussed_velocity(m_prim, m_sec,
                                                    r_min,  velocity, v_disp)
                  / max_prob;
      }
      while (prob<pprob);

      return velocity;

   }

/*
real random_focussed_maxwellian_velocity(const real v_disp, const real v_esc) {

      real v_disp2  = pow(v_disp, 2);
      real vprob, pprob, prob;
      do {
         pprob = randinter(0, 1.);
         vprob = randinter(0, v_esc);
         prob = 1./exp(pow(vprob, 2)/(sqrt(2.)*v_disp2));
      }
      while (prob>pprob);

      return vprob;
   }
*/

real gauss(const real velocity, const real v_disp) {

      real u2 = pow(velocity/v_disp, 2);
      real prob = 1./exp(u2);

      return prob;
  }

real gauss() {

    extern real randinter(real, real);

    static int iset = 0;

    /* System generated locals */
    real ret_val, r_1, r_2;

    /* Builtin functions */
//    double log(), sqrt();

    /* Local variables */
    static real gset, r, v1, v2, fac;

//    srand(idum);

    if (iset == 0) {
        do {
            v1 = 2*randinter(0, 1) - 1.;
            v2 = 2*randinter(0, 1) - 1.;
/* Computing 2nd power */
            r_1 = v1;
/* Computing 2nd power */
            r_2 = v2;
            r = r_1 * r_1 + r_2 * r_2;
        }
        while(r >= 1.);
        fac = sqrt(-2*log(r)/r);
        gset = v1 * fac;
        ret_val = v2 * fac;
        iset = 1;
    } else {
        ret_val = gset;
        iset = 0;
    }
//cerr << "gauss: " << ret_val << endl;
    return ret_val;
} 

real turn_off_mass(const real time, const real z) {

     real mass = 0.5*cnsts.parameters(maximum_main_sequence);
     real mdot = mass;
     real to_time = main_sequence_time(mass, z);

     while(abs(time-to_time)>1.e-4) {
        mdot *= 0.5;
        if(to_time>time)
           mass += mdot;
        else
           mass -= mdot;
        to_time = main_sequence_time(mass, z);

        if(mass>=cnsts.parameters(maximum_main_sequence)
	        -cnsts.safety(minimum_mass_step) ||
           mass<=cnsts.parameters(minimum_main_sequence)
	        -cnsts.safety(minimum_mass_step))
           to_time = time;
     }

     return mass;
}

real main_sequence_time(const real mass, const real z) {

    // Eq.4
    real pow_mass_7 = pow(mass, 7);
    real teller = smc.a(1, z) +smc.a(2, z)*pow(mass, 4) 
                      +smc.a(3, z)*pow(mass, 5.5)  
                      +       pow_mass_7;
    real noemer =  smc.a(4, z)*pow(mass, 2) +smc.a(5, z)*pow_mass_7; 
    real t_bgb = teller/noemer;


    real zeta = log10(z/cnsts.parameters(solar_metalicity));	
    // Eq.6 identified as 'x' by Hurley
    real stars_without_main_sequence_hook = max(0.95, min(0.99,
                         0.95 - 0.03*(zeta + 0.30103)));
			 
    // Eq.7 identified as 'mu' by Hurley
    real stars_with_main_sequence_hook = max(0.5, 
		1.0 - 0.01*max(smc.a(6, z)/pow(mass, smc.a(7, z)),
			       smc.a(8, z) + smc.a(9, z)
			       /pow(mass, smc.a(10, z))));		 

    //Eq. 5
    real t_ms = t_bgb * max(stars_with_main_sequence_hook,
                  stars_without_main_sequence_hook); 
 
    return t_ms;
}


real zero_age_main_sequnece_radius(const real mass, const real z){

    real mx = pow(mass, 0.5);
    real teller = (smc.c(8,z)*pow(mass,2) + smc.c(9,z)*pow(mass,6))*mx + smc.c(10,z)*pow(mass,11) +(smc.c(11,z) +
                                   smc.c(12,z)*mx)*pow(mass,19);
    real noemer = smc.c(13,z) + smc.c(14,z)*pow(mass,2) + (smc.c(15,z)*pow(mass,8) + pow(mass,18) + 
                                                smc.c(16,z)*pow(mass,19))*mx;
    
    return teller/noemer;

}

real roche_radius(const real a, const real m1, const real m2) {

    real q = m1/m2;
    real q1_3 = pow(q, cnsts.mathematics(one_third));
    real q2_3 = pow(q1_3, 2);   //pow(mr, TWO_THIRD);
  
    real Rl =  a*0.49*q2_3/(0.6*q2_3 + log(1 + q1_3));

    return Rl;
}

/*
#define EPS 1.0e-4
#define JMAX 20
#define JMAXP (JMAX+1)
#define K 5

real integrate(const real a,
               const real b,
               real (*func)(const real)) {

        void polint(real xa[], real ya[], int n, real x, real *y, real *dy);
        real trapzd(real (*func)(const real), 
                    real a, real b, int n);
        real ss,dss;
        real s[JMAXP+1],h[JMAXP+1];
        int j;

        h[1]=1.0;
        for (j=1;j<=JMAX;j++) {
                s[j]=trapzd(func,a,b,j);
                if (j >= K) {
                        polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
                        if (abs(dss) < EPS*abs(ss)) return ss;
                }
                s[j+1]=s[j];
                h[j+1]=0.25*h[j];
        }
        cerr<<"Too many steps in routine integrate";
        return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

real trapzd(real (*func)(const real), 
            real a, real b, int n)
{
        real x,tnm,sum,del;
        static real s;
        int it,j;

        if (n == 1) {
                return (s=0.5*(b-a)*((*func)(a)+(*func)(b)));
        } else {
                for (it=1,j=1;j<n-1;j++) it <<= 1;
                tnm=it;
                del=(b-a)/tnm;
                x=a+0.5*del;
                for (sum=0.0,j=1;j<=it;j++,x+=del) sum += (*func)(x);
                s=0.5*(s+(b-a)*sum/tnm);
                return s;
        }
}

#define NRANSI

*/
void polint(real xa[], real ya[], int n, real x, real *y, real *dy)
{
        int i,m,ns=1;
        real den,dif,dift,ho,hp,w;
        real *c,*d;
        void free_vector(real *v, long nl, long nh);
        real *vector(long nl, long nh);

        dif=fabs(x-xa[1]);
        c=vector(1,n);
        d=vector(1,n);
        for (i=1;i<=n;i++) {
                if ( (dift=fabs(x-xa[i])) < dif) {
                        ns=i;
                        dif=dift;
                }
                c[i]=ya[i];
                d[i]=ya[i];
        }
        *y=ya[ns--];
        for (m=1;m<n;m++) {
                for (i=1;i<=n-m;i++) {
                        ho=xa[i]-x;
                        hp=xa[i+m]-x;
                        w=c[i+1]-d[i];
                        if ( (den=ho-hp) == 0.0) 
                           cerr<<"Error in routine polin t";
                        den=w/den;
                        d[i]=hp*den;
                        c[i]=ho*den;
                }
                *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
        }
        free_vector(d,1,n);
        free_vector(c,1,n);
}
#undef NRANSI
#define NR_END 1
#define FREE_ARG char*

void free_vector(real *v, long nl, long nh)
/* free a real vector allocated with vector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}

real *vector(long nl, long nh)
/* allocate a real vector with subscript range v[nl..nh] */
{
        real *v;

        v=(real *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(real)));
        if (!v) cerr<<"allocation failure in vector()";
        return v-nl+NR_END;
}

void fool(char * message) {

      cerr << "Fool: You canot " << message << "!" << endl;
   }

// stellar evolution timescales timescales.
// should the effective_radius be used??
real kelvin_helmholds_timescale(const real mass,
				const real radius,
				const real luminosity) {
    return 31.56*pow(mass,2.)/(radius*luminosity);
  }

real nucleair_evolution_timescale(const real mass,
				  const real luminosity) {
  // t_nuc = 10^10 [years] Msun/Lsun.
  // Assumed that 0.1 Msun is thermalized.

  real fused_mass = 0.1*mass;

  return cnsts.parameters(energy_to_mass_in_internal_units)
       * fused_mass/luminosity;
  //return cnsts.parameters(nucleair_efficiency)*1.4e+7*mass/luminosity;
}

real dynamic_timescale(const real mass,
		       const real radius) {

    return 5.08e-11*sqrt(pow(radius, 3.)/mass); 
}

#else

void main() {

  real v_disp = 8; // km/s
//  real v_disp = 4; // km/s
  for(int i=0; i< 1000; i++) {
//    real velocity = i/(4.*v_disp);
//    real velocity = random_maxwellian_velocity(v_disp);
//    cerr << " v= " << velocity << endl;
//         << " p= " << maxwellian(velocity, v_disp) << endl;
    real velocity = random_focussed_maxwellian_velocity(39.3, 15.3, 
                                     1.,  v_disp, 4*v_disp);
    cerr << " v= " << velocity << endl;

  }

}

#endif
