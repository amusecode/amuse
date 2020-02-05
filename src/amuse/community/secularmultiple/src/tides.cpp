/*
*/

#include "types.h"
#include "tides.h"
#include "newtonian.h" /* for orbital_period */
#include <stdio.h>

double const MINIMUM_MASS_FOR_RADIATIVE_DAMPING_MSUN = 1.2; // in future: make user-adjustable
int const MAIN_SEQUENCE = 1;
int const CHeB = 4;
int const HeMS = 7;
int const HeWD = 10;

/* start addition 13-09-2016 */
int const PreMS = 17;
/* end addition 13-09-2016 */

bool check_for_radiative_damping(int stellar_type, double mass, double convective_envelope_mass, double convective_envelope_radius)
{
    if ((stellar_type == MAIN_SEQUENCE) && (mass/CONST_MSUN >= MINIMUM_MASS_FOR_RADIATIVE_DAMPING_MSUN))
    {
        return true;
    }
    else if ((stellar_type == CHeB) || (stellar_type == HeMS))
    {
        return true;
    }
    else
    {
        return false;
    }
}   
bool check_for_convective_damping(int stellar_type)
{
    /* start addition 13-09-2016 */
    if (stellar_type == PreMS)
    {
        return true;
    }
    /* end addition 13-09-2016 */
    
    if (stellar_type < HeWD)
    {
        return true;
    }
    else
    {
        return false;
    }
}

double from_k_AM_div_T_to_t_V(double k_AM_div_T, double apsidal_motion_constant)
{
    return c_3div2*(2.0*apsidal_motion_constant + 1.0)/k_AM_div_T;
}

double compute_t_V(Particle *star, Particle *companion, double semimajor_axis)
{
    int tides_viscous_time_scale_prescription = star->tides_viscous_time_scale_prescription;
    double t_V = 1.0e10; /* large value by default, i.e. weak tides if tides_viscous_time_scale_prescription is not given the correct value */
    
    if (tides_viscous_time_scale_prescription == 0)
    {
        t_V = star->tides_viscous_time_scale;
    }
    else if (tides_viscous_time_scale_prescription == 1)
    {
        t_V = compute_t_V_hurley
        (
            star->stellar_type,
            star->mass,
            star->convective_envelope_mass,
            companion->mass,
            semimajor_axis,
            star->radius,
            star->convective_envelope_radius,
            star->luminosity,
            star->spin_vec_norm,
            star->tides_gyration_radius,
            star->tides_apsidal_motion_constant
        ); 
    }
    
    return t_V;
}

double compute_t_V_hurley
(
    int stellar_type,
    double mass,
    double convective_envelope_mass,
    double companion_mass,
    double semimajor_axis,
    double radius,
    double convective_envelope_radius,
    double luminosity,
    double spin_angular_frequency,
    double gyration_radius,
    double apsidal_motion_constant
)
{
    bool USE_RADIATIVE_DAMPING = check_for_radiative_damping(stellar_type,mass,convective_envelope_mass,convective_envelope_radius);
    bool USE_CONVECTIVE_DAMPING = check_for_convective_damping(stellar_type);
    double k_AM_div_T,t_V;
    
    if (USE_CONVECTIVE_DAMPING == true && ((convective_envelope_mass <= 0.0) || (convective_envelope_radius <= 0.0)))
    {
        //printf("to rad \n");
        USE_RADIATIVE_DAMPING = true;
    }
    if (radius <= 0.0)
    {
        return 1.0e100;
    }
    
    //printf("stellar_type %d \n",stellar_type);
    //printf("USE_RADIATIVE_DAMPING %d \n",USE_RADIATIVE_DAMPING);
    //printf("USE_CONVECTIVE_DAMPING %d \n",USE_CONVECTIVE_DAMPING);
    
    if (USE_RADIATIVE_DAMPING == true) // radiative damping
    {
        double E2 = 1.592e-09*pow(mass/CONST_MSUN,2.84); // Hurley prescription; Zahn, 1977, A&A, 57, 383 and 1975, A&A, 41, 329
        
#ifdef IGNORE
        // Izzard's prescription not yet implemented
		else if(stardata->preferences->E2_prescription==E2_IZZARD)
		{
			if(stardata->star[star_number].stellar_type<HERTZSPRUNG_GAP)
			{
				double fburn=1.0;
				double x=stardata->star[star_number].aj/stardata->star[star_number].tms;
				if(x<fburn)
				{
					/* log mass and Z */
					double logm=log10(stardata->star[star_number].mass);
					double logz=log10(stardata->common.metallicity);
 
					/* fits for am and E20 */
					double am = 0.15*sin(3.2*logm)+0.31*logm;
					double E20 = -1.23250e+01+1.04550e+01*logm-4.15420e-01*logz-7.18650e+00*logm*logm+1.97940e+00*logm*logm*logm;
			  		E20=pow(10.0,E20);

					/* calc log(E2/E20) */
					E2 = -pow(x+am,4.0)*pow(MAX(1,x/0.95),30.0);
			  
					/* hence E2 */
					E2 = E20 * pow(10.0,E2);
			  
					/* verbosity */
					/*
					if(star->starnum==1)
					{
						printf("E2 kw=%d I=%g (fburn=%g x=%g E20=%g) H=%g\n",
						star->stellar_type,
						E2,
						fburn,x,
						E20,E2_Hurley);
					}
					*/
				}
				else
				{
					/* no conv core */
					E2=0.0;
				}
			}
			else
			{
				E2=0.0;
			}
		}
#endif
        
        k_AM_div_T = E2*pow(1.0 + companion_mass/mass,5.0/6.0)*radius*sqrt(CONST_G*mass/(pow(semimajor_axis,5.0)));
        t_V = from_k_AM_div_T_to_t_V(k_AM_div_T,apsidal_motion_constant);
        return t_V;
        
    }
    else if (USE_CONVECTIVE_DAMPING == true) // convective damping
    {
        double P_orb = 2.0*M_PI*sqrt((semimajor_axis*semimajor_axis*semimajor_axis)/(CONST_G*(mass + companion_mass)));
//        printf("a %g\n",semimajor_axis);
        double P_spin,P_tid;
        
        if (spin_angular_frequency == 0.0)
        {
            P_tid = P_orb;
        }
        else
        {
            P_spin = 2.0*M_PI/spin_angular_frequency;
            P_tid = 1.0/( 1e-10 + fabs( 1.0/P_orb - 1.0/P_spin) );
        }

        double tau_convective = pow( (convective_envelope_mass*convective_envelope_radius*(radius - (1.0/2.0)*convective_envelope_radius))/(3.0*luminosity), 1.0/3.0);
        //double tau_convective = pow( (convective_envelope_mass*radius*radius)/(3.0*luminosity), 1.0/3.0);
//	print 'tau',envelope_mass,envelope_mass*envelope_radius*(radius - (1.0/2.0)*envelope_radius)/(3.0*luminosity)

        double f_convective = pow(P_tid/(2.0*tau_convective),2.0);
        f_convective = min(1.0,f_convective);

        k_AM_div_T = (2.0/21.0)*(f_convective/tau_convective)*(convective_envelope_mass/mass);
        t_V = from_k_AM_div_T_to_t_V(k_AM_div_T,apsidal_motion_constant);

        //printf("test %g %g %g %g %g %g %g  \n",P_spin,P_tid,P_orb,tau_convective,f_convective,k_AM_div_T,t_V);

        //if ((convective_envelope_mass <= 0.0) || (convective_envelope_radius <= 0.0))
        // {
        //     t_V = 1.0e100;
        // }
        
//        printf("test par conv %g %g %g %g %g \n",mass,radius,convective_envelope_mass,convective_envelope_radius,spin_angular_frequency);
//        printf("test conv %g %g %g %g %g \n",P_orb,tau_convective,P_tid,P_spin,f_convective);
        return t_V;

    }
    else // degenerate damping -- 1984MNRAS.207..433C
    {
        double seconds_in_year = 365.25*24.0*3600.0;
        double tau_degenerate = 1.3e7*seconds_in_year;
        k_AM_div_T = (1.0/(3.0*tau_degenerate))*gyration_radius*gyration_radius*pow((luminosity/CONST_L_SUN)/(mass/CONST_MSUN),5.0/7.0);
        t_V = from_k_AM_div_T_to_t_V(k_AM_div_T,apsidal_motion_constant);
        
        return t_V;
    }
}

double compute_EOM_equilibrium_tide_BO_full(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, int tides_method)
/* Barker & Ogilvie (2009; http://adsabs.harvard.edu/abs/2009MNRAS.395.2268B) */

/* NOTE: in SecularMultiple, the h-vector is defined as the orbital angular momentum vector,
 * NOT the SPECIFIC orbital angular momentum vector. Compared to the notation used by Eggleton,
 * h_vec_SecularMultiple = mu*h_vec_Eggleton where mu = m*M/(m+M) is the reduced mass.
 * In particular, note the line `star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I;' */
{
//    printf("tides BO full\n");
//    printf("TIDES %d %d %d\n",binary_index,star_index,companion_index);
    Particle *binary = (*particlesMap)[binary_index];
    Particle *star = (*particlesMap)[star_index];
    Particle *companion = (*particlesMap)[companion_index];
    
    /* orbit quantities */
    double e = binary->e;
    double e_p2 = binary->e_p2;
    double a = binary->a;
    double h = binary->h;
    double *e_vec = binary->e_vec;
    double *h_vec = binary->h_vec;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p3 = binary->j_p3;
    double j_p4 = binary->j_p4;
    double j_p8 = j_p4*j_p4;
    double j_p10 = j_p2*j_p8;  
    double j_p13 = j_p3*j_p10; 
    double j_p4_inv = 1.0/j_p4; 
    double j_p10_inv = 1.0/j_p10; 
    double j_p13_inv = 1.0/j_p13;
    double P_orb = compute_orbital_period(binary);
    double n = 2.0*M_PI/P_orb; /* mean motion */

    /* stellar properties */
    double *spin_vec = star->spin_vec;
    double M = star->mass;
    double m = companion->mass;
    double R = star->radius;
    double k_AM = star->tides_apsidal_motion_constant;
    double rg = star->tides_gyration_radius;

    double t_V = compute_t_V(star,companion,a);
    star->tides_viscous_time_scale = t_V;

    if (t_V!=t_V)
    {
        printf("ERRORRRR\n");
        printf("t_V %g \n",t_V);
        printf("st %d\n",star->stellar_type);
        printf("pr %d\n",star->tides_viscous_time_scale_prescription);
        printf("M %g\n",M);
        printf("star->convective_envelope_mass %g\n",star->convective_envelope_mass);
        printf("m %g\n",m);
        printf("a %g\n",a);
        printf("R %g\n",R);
        printf("star->convective_envelope_radius %g\n",star->convective_envelope_radius);
        printf("star->luminosity %g\n",star->luminosity);
    }
    if (1==0)
    {
        printf("t_V %g \n",t_V);
        printf("st %d\n",star->stellar_type);
        printf("pr %d\n",star->tides_viscous_time_scale_prescription);
        printf("M %g\n",M);
        printf("star->convective_envelope_mass %d\n",star->convective_envelope_mass);
        printf("m %g\n",m);
        printf("a %g\n",a);
        printf("R %g\n",R);
        printf("star->convective_envelope_radius %g\n",star->convective_envelope_radius);
        printf("star->luminosity %g\n",star->luminosity);
    }

    double tau = 3.0*(1.0 + 1.0/(2.0*k_AM))*R*R*R/(CONST_G*M*t_V);

    double I = rg*M*R*R; // moment of intertia

    double R_div_a = R/a;
	double R_div_a_p5 = pow(R_div_a,5.0);
    double t_f_inv = 3.0*k_AM*tau*n*n*(m/M)*R_div_a_p5;

    double f_tides1 = f_tides1_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides2 = f_tides2_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides3 = f_tides3_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides4 = f_tides4_function_BO(e_p2,j_p10_inv,j_p13_inv);
    double f_tides5 = f_tides5_function_BO(e_p2,j_p10_inv,j_p13_inv);

    double spin_vec_dot_e_vec = dot3(spin_vec,e_vec);
    double spin_vec_dot_h_vec = dot3(spin_vec,h_vec);

    double q_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);

    double spin_vec_dot_e_vec_unit = dot3(spin_vec,e_vec_unit);
    double spin_vec_dot_h_vec_unit = dot3(spin_vec,h_vec_unit);
    double spin_vec_dot_q_vec_unit = dot3(spin_vec,q_vec_unit);
    double mu = m*M/(m+M);
    double C = m*k_AM*R_div_a_p5/(mu*n);
    
    double C_rot;
    double X_rot = 0.0;
    double Y_rot = 0.0;
    double Z_rot = 0.0;
    double Z_TB = 0.0;
    if (include_rotation_precession_terms == 1)
    {
        if (tides_method == 1)
        {
            C_rot = C*j_p4_inv*spin_vec_dot_h_vec_unit;
            X_rot = -C_rot*spin_vec_dot_e_vec_unit;
            Y_rot = -C_rot*spin_vec_dot_q_vec_unit;
        }
//        printf("e %g q %g \n",spin_vec_dot_e_vec_unit/norm3(spin_vec),spin_vec_dot_q_vec_unit/norm3(spin_vec));

        Z_rot = C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit);
//        printf("1 %g 2 %g 3 %g\n",2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit ,- spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit ,- spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit);
//        printf("X %g Y %g Z %g\n",X_rot,Y_rot,Z_rot);
    }    
    if (include_tidal_bulges_precession_terms == 1)
    {
        Z_TB = C*15.0*n*n*(mu/M)*f_tides2;
    }
    
    double X = X_rot;
    double Y = Y_rot;
    double Z = Z_rot + Z_TB;
    
    //printf("M %g m %g R %g k_AM %g a %g e %g Z %g Z_rot %g Z_TB %g\n",M,R,k_AM,a,e,Z,Z_rot,Z_TB);


    double dh_vec_dt_star_i;
   
    for (int i=0; i<3; i++)
    {

        if (include_tidal_friction_terms == 1)
        {
            dh_vec_dt_star_i = -t_f_inv*( (h/(2.0*n))*(spin_vec_dot_e_vec*f_tides5*e_vec[i] - spin_vec[i]*f_tides3) \
                + h_vec[i]*(f_tides4 - spin_vec_dot_h_vec*f_tides2/(2.0*n*h)) );
            binary->dh_vec_dt[i] += dh_vec_dt_star_i;

            binary->de_vec_dt[i] += -(t_f_inv/h)*( spin_vec_dot_e_vec*f_tides2*h_vec[i]/(2.0*n) \
                + 9.0*e_vec[i]*(f_tides1*h - c_11div18*spin_vec_dot_h_vec*f_tides2/n) );
            star->dspin_vec_dt[i] += -dh_vec_dt_star_i/I;
            //printf("test %g %g\n",spin_vec_dot_e_vec*f_tides5*e_vec[i] - spin_vec[i]*f_tides3);
        }
        if (include_rotation_precession_terms == 1 || include_tidal_bulges_precession_terms == 1)
        {
            if (e >= minimum_eccentricity_for_tidal_precession)
            {
                binary->de_vec_dt[i] += e*(Z*q_vec_unit[i] - Y*h_vec_unit[i]);
                
                dh_vec_dt_star_i = h*(-X*q_vec_unit[i] + Y*e_vec_unit[i]);
                binary->dh_vec_dt[i] += dh_vec_dt_star_i;
                star->dspin_vec_dt[i] += -dh_vec_dt_star_i/I;
                
//                printf("ok %d %d\n",include_rotation_precession_terms,include_tidal_bulges_precession_terms);
            }
        }
    }

//    printf("I %g %g %g %g\n",rg,M,R,Q_prime);
//    printf("f dh_vec_dt %g %g %g\n",binary->dh_vec_dt[0],binary->dh_vec_dt[1],binary->dh_vec_dt[2]);
//    printf("f dspin_vec_dt %g %g %g\n",star->dspin_vec_dt[0],star->dspin_vec_dt[1],star->dspin_vec_dt[2]);
    return 0;
}


double f_tides1_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div4 + e_p2*(c_15div8 + e_p2*c_5div64)));
}
double f_tides2_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_3div2 + e_p2*c_1div8));
}
double f_tides3_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_9div2 + e_p2*c_5div8));
}
double f_tides4_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div2 + e_p2*(c_45div8 + e_p2*c_5div16)));
}
double f_tides5_function_BO(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(3.0 + c_1div2*e_p2);
}  


double compute_EOM_equilibrium_tide(ParticlesMap *particlesMap, int binary_index, int star_index, int companion_index, int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession)

/* Equilibrium tide in vector form adopted from Eggleton & Kisseleva 1998 */
 
/* NOTE: in SecularMultiple, the h-vector is defined as the orbital angular momentum vector,
 * NOT the SPECIFIC orbital angular momentum vector. Compared to the notation used by Eggleton,
 * h_vec_SecularMultiple = mu*h_vec_Eggleton where mu = m*M/(m+M) is the reduced mass.
 * In particular, note the line `star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I;' */

{
//    printf("tides EK \n");
//    printf("TIDES %d %d %d\n",binary_index,star_index,companion_index);
    Particle *binary = (*particlesMap)[binary_index];
    Particle *star = (*particlesMap)[star_index];
    Particle *companion = (*particlesMap)[companion_index];
    
    /* orbit quantities */
    double e = binary->e;
    double e_p2 = binary->e_p2;
    double a = binary->a;
    double h = binary->h;
    double *e_vec = binary->e_vec;
    double *h_vec = binary->h_vec;
    double *e_vec_unit = binary->e_vec_unit;
    double *h_vec_unit = binary->h_vec_unit;
    double j = binary->j;
    double j_p2 = binary->j_p2;
    double j_p3 = binary->j_p3;
    double j_p4 = binary->j_p4;
    double j_p8 = j_p4*j_p4;
    double j_p10 = j_p2*j_p8;  
    double j_p13 = j_p3*j_p10; 
    double j_p4_inv = 1.0/j_p4;
    double j_p10_inv = 1.0/j_p10; 
    double j_p13_inv = 1.0/j_p13;
    double P_orb = compute_orbital_period(binary);
    double n = 2.0*M_PI/P_orb; /* mean motion */
    
    /* stellar properties */
    double *spin_vec = star->spin_vec;
    double M = star->mass;
    double m = companion->mass;
    double mu = m*M/(m+M);
    double R = star->radius;
    double k_AM = star->tides_apsidal_motion_constant;
    double t_V = compute_t_V(star,companion,a);
    star->tides_viscous_time_scale = t_V;
    double rg = star->tides_gyration_radius;
    double I = rg*M*R*R; // moment of intertia
    
    double R_div_a = R/a;
	double R_div_a_p5 = pow(R_div_a,5.0);
    double R_div_a_p8 = pow(R_div_a,8.0);
    double t_f_inv = (9.0/t_V)*R_div_a_p8*((M+m)*m/(M*M))*(1.0 + 2.0*k_AM)*(1.0 + 2.0*k_AM);
    
    double q_vec_unit[3];
    cross3(h_vec_unit,e_vec_unit,q_vec_unit);

    double spin_vec_dot_e_vec_unit = dot3(spin_vec,e_vec_unit);
    double spin_vec_dot_h_vec_unit = dot3(spin_vec,h_vec_unit);
    double spin_vec_dot_q_vec_unit = dot3(spin_vec,q_vec_unit);
        
    double V,W,X,Y,Z;
    VWXYZ_tides_function(
        include_tidal_friction_terms, include_tidal_bulges_precession_terms, include_rotation_precession_terms, minimum_eccentricity_for_tidal_precession, \
        t_f_inv,k_AM, \
        m, M, mu, n, R_div_a_p5, \
        e, e_p2, j_p4_inv, j_p10_inv, j_p13_inv, \
        spin_vec_dot_h_vec_unit, spin_vec_dot_e_vec_unit, spin_vec_dot_q_vec_unit, \
        &V, &W, &X, &Y, &Z);

    //printf("t_f_inv V %g W %g X %g Y %g Z %g\n",t_f_inv,V,W,X,Y,Z);

//    printf("js %g %g %g\n",j_p10_inv,j_p13_inv);
//    printf("fs %g %g %g %g %g \n",f_tides1,f_tides2,f_tides3,f_tides4,f_tides5);

    double dh_vec_dt_star[3];
    
    for (int i=0; i<3; i++)
    {
        dh_vec_dt_star[i] = h*( Y*e_vec_unit[i] - W*h_vec_unit[i] - X*q_vec_unit[i] );
        binary->dh_vec_dt[i] += dh_vec_dt_star[i];

        binary->de_vec_dt[i] += e*( -V*e_vec_unit[i] - Y*h_vec_unit[i] + Z*q_vec_unit[i] );
        
        star->dspin_vec_dt[i] += -dh_vec_dt_star[i]/I; /* conservation of total angular momentum (orbit+spin) */




//        printf("test %g %g\n",I*spin_vec[2], h_vec[2]);
//        printf("test2 %g %g %g\n",I*spin_vec[0] + h_vec[0],I*spin_vec[1] + h_vec[1],I*spin_vec[2] + h_vec[2]);        
//        printf("test2 %g %g %g\n",I*spin_vec[0] + mu*h_vec[0],I*spin_vec[1] + mu*h_vec[1],I*spin_vec[2] + mu*h_vec[2]);        
    }

//    double omega = norm3(spin_vec);
//    double h2 = sqrt(CONST_G*(M+m)*a*(1.0-e*e));

//    printf("R %g\n",h/(mu*h2));
//    printf("H%g H_a %g H_b %g\n",mu*h2 + I*omega, mu*h2, I*omega);
//    printf("H%g H_a %g H_b %g\n",mu*h + I*omega, mu*h, I*omega);

//    printf("I %g %g %g %g\n",rg,M,R,Q_prime);
//    printf("f dh_vec_dt %g %g %g\n",binary->dh_vec_dt[0],binary->dh_vec_dt[1],binary->dh_vec_dt[2]);
//    printf("f dspin_vec_dt %g %g %g\n",star->dspin_vec_dt[0],star->dspin_vec_dt[1],star->dspin_vec_dt[2]);
    return 0;
}

double f_tides1_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div4 + e_p2*(c_15div8 + e_p2*c_5div64)));
}
double f_tides2_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_3div2 + e_p2*c_1div8));
}
double f_tides3_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p13_inv*(1.0 + e_p2*(c_15div2 + e_p2*(c_45div8 + e_p2*c_5div16)));
}
double f_tides4_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(3.0 + e_p2*c_3div8));
}
double f_tides5_function(double e_p2, double j_p10_inv, double j_p13_inv)
{
    return j_p10_inv*(1.0 + e_p2*(c_9div2 + e_p2*c_5div8));
}


int VWXYZ_tides_function
(
    int include_tidal_friction_terms, int include_tidal_bulges_precession_terms, int include_rotation_precession_terms, double minimum_eccentricity_for_tidal_precession, \
    double t_f_inv, double k_AM, \
    double m, double M, double mu, double n, double R_div_a_p5, \
    double e, double e_p2, double j_p4_inv, double j_p10_inv,double j_p13_inv, \
    double spin_vec_dot_h_vec_unit, double spin_vec_dot_e_vec_unit,double spin_vec_dot_q_vec_unit, \
    double* V, double* W, double* X, double* Y, double* Z
)
{
   
    *V = 0.0;
    *W = 0.0;
    *X = 0.0;
    *Y = 0.0;
    *Z = 0.0;

    double f2 = f_tides2_function(e_p2,j_p10_inv,j_p13_inv); /* needed for both pure tidal dissipation and pure tidal bulges terms */

//    if (e < minimum_eccentricity_for_tidal_precession)
//    {
//        return 0;
//    }


    if (include_tidal_friction_terms == 1)
    {
        double f1 = f_tides1_function(e_p2,j_p10_inv,j_p13_inv);
        double f3 = f_tides3_function(e_p2,j_p10_inv,j_p13_inv);
        double f4 = f_tides4_function(e_p2,j_p10_inv,j_p13_inv);
        double f5 = f_tides5_function(e_p2,j_p10_inv,j_p13_inv);
    
        *V += 9.0*t_f_inv*(f1    - c_11div18*(spin_vec_dot_h_vec_unit/n)*f2);
        *W += t_f_inv*(f3        - (spin_vec_dot_h_vec_unit/n)*f4);
        *X +=  -t_f_inv*spin_vec_dot_q_vec_unit*f5/(2.0*n);
        *Y += t_f_inv*spin_vec_dot_e_vec_unit*f2/(2.0*n);
//	printf("TF X %g Y %g V %g W  %g\n",*X,*Y,*Z,*V,*W);

    }

    if (e < minimum_eccentricity_for_tidal_precession)
    {
        return 0;
    }


    if ((include_tidal_bulges_precession_terms == 1) || (include_rotation_precession_terms) == 1)
    {
        double C = m*k_AM*R_div_a_p5/(mu*n);
    
        if (include_tidal_bulges_precession_terms == 1)
        {
            *Z += C*15.0*n*n*(mu/M)*f2;
//            printf("include_tidal_bulges_precession_terms Z %g\n",*Z);
        }    
        if (include_rotation_precession_terms == 1)
        {
            double C_XY = -C*spin_vec_dot_h_vec_unit*j_p4_inv;
            *X += C_XY*spin_vec_dot_e_vec_unit;
            *Y += C_XY*spin_vec_dot_q_vec_unit;            


//		printf("ROT X %g Y %g Z %g\n",C_XY*spin_vec_dot_e_vec_unit,C_XY*spin_vec_dot_q_vec_unit,C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit));
//        
//            printf("Z1 %g\n",*Z);
            *Z += C*c_1div2*j_p4_inv*(2.0*spin_vec_dot_h_vec_unit*spin_vec_dot_h_vec_unit - spin_vec_dot_q_vec_unit*spin_vec_dot_q_vec_unit - spin_vec_dot_e_vec_unit*spin_vec_dot_e_vec_unit);
//            printf("Z2 %g\n",*Z);
//            printf("O %g %g %g\n",spin_vec_dot_h_vec_unit,spin_vec_dot_e_vec_unit,spin_vec_dot_q_vec_unit);
//            printf("include_rotation_precession_terms Z %g\n",*Z);
        }
    }
    return 0;
}

