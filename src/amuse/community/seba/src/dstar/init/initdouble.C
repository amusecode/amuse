//// initdouble:  provides routines for the initialization of primordial
////              binaries.
////
////              Initialized parameters include:
////              - mass of the most massice component (primary star)
////              - mass of its binary companion (secondary star)
////              - semi-major axis of the binary system
////              - orbital eccentricity
////
////              routines incuded can be found in double_star.h.
////              The mass function routines are adoped from mkmass.C
////              and are defined in starbase.h
//// 
////              externally visible routines are:
////              -get_random_mass_ratio
////              -get_random_semi_major_axis
////              -get_random_eccentricity
////              The two utilities for the various parameter are:
////              -extract_...._distribution_type_string(....)
////              and
////              -type_string(char*)
////
////              The executable takes initial conditions (see Options)
////              and returns randomized binary parameters.                 
////                 
//// Options:   -M    upper primary mass limit [100[Msun]]
////            -m    lower limit to primary mass [0.1[Msun]]
////            -x    mass function exponent in case of power law [-2.35]
////            -F/f  mass function option: 0) Equal mass
////                                        1) Power-law [default]
////                                        2) Miller & Scalo
////                                        3) Scalo
////                                        4) Kroupa
////            Option -F requires one of the following strings:
////                      (mf_Power_Law, Miller_Scalo, Scalo, Kroupa)
////                   -f requires the appropriate interger (see mkmass.C)
////             -A   maximum semi-major axis limit [1000000[Rsun]]   
////             -a   minimum semi-major axis limit [0] 
////             -y   exponent for a power-law distribution  [0] (flat in log)
////             -G/g Semi major axis option: 0) Equal_sma
////                                          1) Power Law [default]
////                                          2) Duquennoy & Mayor (1987)
////                                          3) Eggleton (1999)
////            Option -G requires one of the following strings:
////                      (Equal_sma, sma_Power_Law, Duquennoy_Mayor, Eggleton)
////                   -g requires appropriate interger (see double_star.h)
////             -E   maximum eccentricity [1] 
////             -e   minimum eccentricity [0] 
////             -v   exponent for a power-law distribution 
////             -U/u eccentricity option: 0) Equal eccentricity
////                                       1) Power Law 
////                                       2) Thermal distribution [default]
////            Option -U requires one of the following strings:
////                      (Equal_ecc, ecc_Power_Law, Thermal_Distribution)
////                   -u requires appropriate interger (see double_star.h)
////             -Q   maximum mass ratio [1]
////             -q   minimum mass ratio [0.1 / selected primary mass]
////             -w   exponent for a power-law distribution  
////             -P/p eccentricity option: 0) constant mass ratio
////                                       1) Flat distribution
////                                       1) Power Law 
////                                       2) Hogeveen (1992)
////            Option -P requires one of the following strings:
////                      (Equal_q, Flat_q, qf_Power_Law, Hogeveen)
////                   -p requires appropriate interger (see double_star.h)
////
////             -n   number of output steps. [1]
////
//   Note:  libnode.a is referenced for the routines which produce the 
//          mass function
//
//	version 1.0	Simon Portegies Zwart, Utrecht, 1994
//	version 3.3	Simon Portegies Zwart, Cambridge, March 1999
//

#include "node.h"
#include "double_star.h"
#include "main_sequence.h"

#ifndef TOOLBOX

#define REPORT_ADD_DOUBLE        false
#define REPORT_DEL_DOUBLE        false
#define REPORT_EVOLVE_DOUBLE     false

#define SEED_STRING_LENGTH 255

void add_secondary(node* original, real mass_ratio) {

    node* primary = new node;
    node* secondary = new node;

    // Add new links.

    original->set_oldest_daughter(primary);

    primary->set_parent(original);
    secondary->set_parent(original);

    primary->set_younger_sister(secondary);
    secondary->set_elder_sister(primary);

    // Set new masses.

    primary->set_mass(original->get_mass());
    secondary->set_mass(mass_ratio*original->get_mass());
    original->inc_mass(secondary->get_mass());

    // Naming convention:

    if (original->get_name() == NULL)
        if (original->get_index() >= 0) {
            char tmp[64];
            sprintf(tmp, "%d", original->get_index());
            original->set_name(tmp);
        }

    primary->set_name(original->get_name());
    secondary->set_name(original->get_name());
    strcat(primary->get_name(), "a");
    strcat(secondary->get_name(), "b");

   }

void mksecondary(node* b, real binary_fraction, real lower_limit) {

    // For now, use a flat distribution in secondary mass ratio.
    // Assume that the probability of a star being the primary of
    // a binary is independent of mass.

    real sum = 0;
    b->set_mass(0);

    for_all_daughters(node, b, bi) {
        sum += binary_fraction;
        if (sum >= 1) {
            sum -= 1;

            real mass_ratio = randinter(lower_limit, 1);        // Quick fix...
            add_secondary(bi, mass_ratio);

        }
        b->inc_mass(bi->get_mass());
    }
}

#if 0
real random_exponential_mass(const real m_min,
				   const real m_max,
				   const real m_alpha) {

  real random    = randinter(0., 1.);
  real m_const   = pow(m_min, 1-m_alpha)
    - pow(m_max, 1-m_alpha);
  real m_prim    = pow(random*m_const
		       + pow(m_max, 1-m_alpha), 1/(1-m_alpha));

  return m_prim;
}
#endif


char* type_string(mass_ratio_distribution qf) {
    
    local char  qf_name[SEED_STRING_LENGTH];
    switch(qf) {
       case Equal_q:
            sprintf(qf_name, "Equal_q"); 
	    break;	    
       case Flat_q:
            sprintf(qf_name, "Flat_q"); 
	    break;	    
       case qf_Power_Law:
            sprintf(qf_name, "Power_Law"); 
	    break;	    
       case Thermal_Distribution:
            sprintf(qf_name, "Hogeveen"); 
	    break;	    
       default:
            sprintf(qf_name, "Unknown_qf"); 
	    break;
    }
    return qf_name;
}

mass_ratio_distribution extract_mass_ratio_distribution_type_string(char* type_string) {

     mass_ratio_distribution type = Unknown_qf;

     if (!strcmp(type_string, "Equal_q"))
        type = Equal_q;
     else if (!strcmp(type_string, "Unknown_qf")) 
        type = Unknown_qf;
     else if (!strcmp(type_string, "Flat_q")) 
        type = Flat_q;
     else if (!strcmp(type_string, "Power_Law")) 
        type = qf_Power_Law;
     else if (!strcmp(type_string, "Hogeveen")) 
        type = Hogeveen;
     else {
	 cerr << " in extract_eccentricity_type_string." << endl;
	 err_exit("No proper mass-ratio distribution indicated");
//	 exit(-1);
     }

     return type;
 }

//    Hogeveen S. PhD Thesis Amsterdam 1991.
local real qf_Hogeveen(real q_lower, real q_upper) { 

    real q;
    do {
	q = 1./pow(randinter(0.125, 1.), 1/3.) - 1;
    }
    while(q<q_lower || q>=q_upper);

    return q;
}

real get_random_mass_ratio(real q_lower, real q_upper, 
			   mass_ratio_distribution qf, 
			   real exponent) {
    real q;
    switch(qf) {
       case Equal_q:
	  if (q_lower==0) {
	     cerr << "get_random_mass_ratio:"<<endl;
	     cerr << "unambiguous choise of Equal_q."<<endl;
	     cerr << "Use -q option to set fixed mass ratio."<<endl;
	     exit(1);
	  }
	  q =  q_lower;
	       break;
       case Flat_q:
	  q =  randinter(q_lower, q_upper);
	       break;
       case qf_Power_Law:
	  q =  general_power_law(q_lower, q_upper, exponent);
	       break;
       case Hogeveen:
	  q = qf_Hogeveen(q_lower, q_upper);
	       break;
       default:
          cerr << "WARNING: \n"
	       << "        real get_random_mass_ratio:\n"
	       << "        parameters not properly defined.\n";
	  exit(1);
    }

    return q;
}

char* type_string(sma_distribution smaf) {
    
    local char  smaf_name[SEED_STRING_LENGTH];	
    switch(smaf) {
       case Equal_sma:
            sprintf(smaf_name, "Equal_sma"); 
	    break;	    
       case sma_Power_Law:
            sprintf(smaf_name, "Power_Law"); 
	    break;	    
       case Duquennoy_Mayor:
            sprintf(smaf_name, "Duquennoy_Mayor"); 
	    break;	    
       case Eggleton:
            sprintf(smaf_name, "Eggleton"); 
	    break;	    
       default:
            sprintf(smaf_name, "Unknown_smaf"); 
	    break;
    }
    return smaf_name;
}

sma_distribution 
    extract_semimajor_axis_distribution_type_string(char* type_string) {

     sma_distribution type = Unknown_smaf;

     if (!strcmp(type_string, "Equal_sma"))
        type = Equal_sma;
     else if (!strcmp(type_string, "Duquennoy_Mayor"))
        type = Duquennoy_Mayor;
     else if (!strcmp(type_string, "Eggleton"))
        type = Eggleton;
     else if (!strcmp(type_string, "Unknown_smaf")) 
        type = Unknown_smaf;
     else if (!strcmp(type_string, "Power_Law")) 
        type = sma_Power_Law;
     else {
	 cerr << "No proper semimajor axis distribution indicated"
	      << " in extract_semimajor_axis_type_string." << endl;
	 exit(1);
     }

     return type;
 }

// Duquennoy, A., \& Mayor, M., 1991,  AAp 248, 485.
local real smaf_Duquenoy_Mayor(real a_lower, real a_upper, 
			       real total_mass) {

    real lp_mean = 4.8;
    real lp_sigma = 2.3;
    real lp, sma;
    do {
	lp = lp_mean + lp_sigma*gauss();
	real a_tmp = pow(pow(10., lp)*cnsts.physics(seconds_per_day), 2)
	    * (cnsts.physics(G)*cnsts.parameters(solar_mass)*
	       total_mass)/4*pow(cnsts.mathematics(pi), 2);
	sma = pow(a_tmp, 1./3.)/cnsts.parameters(solar_radius);
    }
    while(sma<a_lower || sma>=a_upper);

    return sma;
}

// See Eggleton 1999 Equation 1.6.3 (in private copy of hist book).
local real smaf_Eggleton(real a_lower, real a_upper, 
			 real m_prim, real m_sec) {

    real p_lower = semi_to_period(a_lower, m_prim, m_sec);
    real p_upper = semi_to_period(a_upper, m_prim, m_sec);

//    PRC(p_lower);PRL(p_upper);

    real mpf     = pow(m_prim, 2.5)/5.e+4;
    real rnd_min = pow(p_upper * mpf, 1./3.3)
                 / (1 + pow(p_upper * mpf, 1./3.3));
    real rnd_max = pow(p_lower * mpf, 1./3.3) 
	         / (1 + pow(p_lower * mpf, 1./3.3)); 

    real rnd = randinter(rnd_min, rnd_max);
    real p   = pow(rnd/(1-rnd), 3.3)/mpf;

//    PRC(mpf);PRC(rnd_min);PRC(rnd_max);PRL(p);

    real sma = period_to_semi(p, m_prim, m_sec);

    return sma;
}

real get_random_semimajor_axis(real a_lower, real a_upper, 
			       sma_distribution smaf, 
			       real exponent, real m_prim, real m_sec) {
    real a;
    switch(smaf) {
       case Equal_sma:
	  if (a_lower!=a_upper) {
	     cerr << "get_random_semimajor_axis:"<<endl;
	     cerr << "unambiguous choise of Equal_sma."<<endl;
	     exit(1);
	  }
	  a =  a_lower;
	       break;
       case sma_Power_Law:
	  a =  general_power_law(a_lower, a_upper, exponent);
	       break;
       case Duquennoy_Mayor: 
	  a = smaf_Duquenoy_Mayor(a_lower, a_upper, m_prim+m_sec);
	       break;
       case Eggleton:
	  a = smaf_Eggleton(a_lower, a_upper, m_prim, m_sec);
	       break;
       default:
          cerr << "WARNING: \n"
	       << "        real get_random_semi_major_axis:\n"
	       << "        parameters not properly defined.\n";
	  exit(1);
    }

    return a;
}


char* type_string(ecc_distribution eccf) {
    
    local char  eccf_name[SEED_STRING_LENGTH];	// permanent
    switch(eccf) {
       case Equal_ecc:
            sprintf(eccf_name, "Equal_ecc"); 
	    break;	    
       case ecc_Power_Law:
            sprintf(eccf_name, "Power_Law"); 
	    break;	    
       case Thermal_Distribution:
            sprintf(eccf_name, "Thermal_Distribution"); 
	    break;	    
       default:
            sprintf(eccf_name, "Unknown_eccf"); 
	    break;
    }
    return eccf_name;
}

ecc_distribution 
    extract_eccentricity_distribution_type_string(char* type_string) {

     ecc_distribution type = Unknown_eccf;

     if (!strcmp(type_string, "Equal_ecc"))
        type = Equal_ecc;
     else if (!strcmp(type_string, "Unknown_eccf")) 
        type = Unknown_eccf;
     else if (!strcmp(type_string, "Power_Law")) 
        type = ecc_Power_Law;
     else {
	 cerr << "No proper eccentricity distribution indicated"
	      << " in extract_eccentricity_type_string." << endl;
	 exit(1);
     }

     return type;
 }

local real eccf_Thermal_Distribution(real e_lower, real e_upper) { 

    real ecc = sqrt(randinter(e_lower, e_upper));
    return ecc;
}

real get_random_eccentricity(real e_lower, real e_upper, 
			     ecc_distribution eccf, 
			       real exponent) {
    real e;
    switch(eccf) {
       case Equal_sma:
	  if (e_lower<0) {
	     cerr << "get_random_eccentricity:"<<endl;
	     cerr << "unambiguous choise of Equal_ecc."<<endl;
	     cerr << "Use -e option to set fixed eccentricity."<<endl;
	     exit(1);
	  }
	  e =  e_lower;
	       break;
       case ecc_Power_Law:
	  if (e_lower<0) e_lower=0;
	  e =  general_power_law(e_lower, e_upper, exponent);
	       break;
       case Thermal_Distribution:
	  if (e_lower<0) e_lower=0;
	  e = eccf_Thermal_Distribution(e_lower, e_upper);
	       break;
       default:
          cerr << "WARNING: \n"
	       << "        real get_random_eccentricity:\n"
	       << "        parameters not properly defined.\n";
	  exit(1);
    }

    return e;
}

local void determine_semi_major_axis_limits(real m_prim, real m_sec, real ecc,
                                            real &a_min, real &a_max, real z) {

    // This is very wrong, but has no effect, happily (SPZ+MS 9 July 2003)
    //    real a_prim=0, a_sec=0;
    //    if (a_min>0) {
    //	a_prim = roche_radius(a_min, m_prim, m_sec);
    //	a_sec  = roche_radius(a_min, m_sec, m_prim);
    //    }
    //    PRC(a_prim);PRL(a_sec);
    //    a_min = Starlab::max(a_min, Starlab::max(a_prim, a_sec));

    if(ecc<0) err_exit("eccentricity <0");

    real r_prim = zero_age_main_sequnece_radius(m_prim, z);
    real r_sec  = zero_age_main_sequnece_radius(m_sec, z);

    real peri_min = Starlab::max(a_min, r_prim+r_sec);
    a_min = peri_min/(1-ecc);

    //    real apo_min = Starlab::max(a_max, r_prim+r_sec);
    //    a_max = Starlab::max(a_max, a_min);

    return;
}

void mkrandom_binary( real m_min,  real m_max,
		      mass_function mf,  real m_exp,
		      real q_min,  real q_max,
		      mass_ratio_distribution qf,  real q_exp,
		      real a_min,  real a_max,
		      sma_distribution af,  real a_exp,
		      real e_min,  real e_max,
		      ecc_distribution ef,  real e_exp,
		     real &m_prim, real &m_sec, real &semi,
		     real &ecc, real z) {

    m_prim = get_random_stellar_mass(m_min, m_max, mf, m_exp);
    //    PRL(m_prim);

  // Initial secondary mass (selected between m_min and m_prim
  // with equal probability per unit mass.
    if(q_min<=0){
      	q_min =  0.1/m_prim; //minimum mass secondary = 0.1
      	//q_min = m_min/m_prim;
    }
    real q = get_random_mass_ratio(q_min, q_max, qf, q_exp);

    m_sec = q*m_prim;
    //    PRL(m_sec);

    // Assume for now that
    //	stellar radius [Rsun] = mass [Msun].
    // This is of course not correct, but for the moment
    // good enough.  mkbinary does not know much about
    // stars.
    real r_prim = zero_age_main_sequnece_radius(m_prim, z);
    real r_sec = zero_age_main_sequnece_radius(m_sec, z);

    if(e_max>=1 && ef!=Equal_ecc) 
	e_max = Starlab::max(e_min, Starlab::min(1., 1 - (r_prim+r_sec)/a_max));

    // SPZ+MS 9 July 2003
    // Make sure that binary is initialized between a_min (minimum
    // pericenter distance) and a_max (maximium orbital
    // separation). Eccentricity is rechosen until a_min<a_max.
    real a_min_org = a_min;
    do {
      a_min = a_min_org;
      ecc = get_random_eccentricity(e_min, e_max, ef, m_prim+m_sec);
      //      PRL(ecc);

      // The Initial orbital separation is chosen flat in log a.
      determine_semi_major_axis_limits(m_prim, m_sec, ecc, a_min, a_max, z);
      //      PRC(a_min);PRL(a_max);
    }
    while(a_min>a_max);

    semi = get_random_semimajor_axis(a_min, a_max, af, a_exp, m_prim, m_sec);
    //    PRL(semi);

    //    cerr << m_prim <<" "
    //    	 << m_sec <<" "
    //    	 << semi <<" "
    //    	 << ecc << endl;
}


void print_initial_binary_distributions(real m_min,  real m_max,
					mass_function mf,  real m_exp,
					real q_min,  real q_max,
					mass_ratio_distribution qf,  
					real q_exp,
					real a_min,  real a_max,
					sma_distribution af,  real a_exp,
					real e_min,  real e_max,
					ecc_distribution ef,  real e_exp) {


    cout << "Use the following initial distribution functions:" << endl;
    cout << "    -Mass function is " << type_string(mf);
    if(mf==mf_Power_Law)
        cout << " with exponent " << m_exp;

    if(mf==Equal_Mass)
        cout << " with value " << m_min << endl;
    else
        cout << "\n                      between "
             << m_min << " and " << m_max << endl;
    cout << "    -mass-ratio distribution is " << type_string(qf);
    if(qf==qf_Power_Law)
        cout << " with exponent " << q_exp;

    if (qf==Equal_q)
        cout << " with value " << q_min << endl;
    else if (qf!=Flat_q)
        cout << endl;
    else
        cout << "\n                      between "
             << q_min << " and " << q_max << endl;

    cout << "    -Semi-major axis distribution is " << type_string(af);
    if(af==sma_Power_Law)
        cout << " with exponent " << a_exp;

    if(af==Equal_sma)
        cout << " with value " << a_min << endl;
    else
        cout << "\n                      between "
             << a_min << " and " << a_max << endl;
    cout << "    -eccenctriciy distribution is " << type_string(ef);
    real e_lower = 0;
    if (e_min>=0) e_lower=e_min;
    if(ef==ecc_Power_Law)
        cout << " with exponent " << e_exp;

    if(ef==Equal_ecc)
        cout << " with value " << e_lower << endl;
    else
        cout << "\n                      between "
             << e_lower << " and " << e_max << endl;
  }


void  adddouble(node * b, real dyn_time,
		binary_type type,		// Defaults set in
		bool random_initialization,	// double_star.h
		real a_min, real a_max,
		real e_min, real e_max)
{
  if (REPORT_ADD_DOUBLE) {
    int p = cerr.precision(HIGH_PRECISION);
    cerr<<"adddouble: "<<b<<" "<<dyn_time;
    cerr.precision(p);
    cerr <<" "<<a_min<<" "<<a_max<<" "<<e_min<<" "<<e_max << endl;
  }
    
    real stellar_time = b->get_starbase()->conv_t_dyn_to_star(dyn_time);
cerr << "addstar(node... called from adddouble(node ..." << endl;
    addstar(b, stellar_time, Main_Sequence);

    real ecc, sma;
    // binary_type local_type = Unknown_Binary_Type;

    int id;
    real a_const = log(a_max) - log(a_min);
    for_all_nodes(node, b, bi) {
	if (bi->is_parent() && !bi->is_root()) {
            story * s = b->get_starbase()->get_star_story();
//            extract_story_chapter(local_type, sma, ecc, *s);
            b->get_starbase()->set_star_story(NULL);
            delete s;

            id = bi->get_index();

            if (random_initialization) {
               sma = a_min*exp(randinter(0., a_const));
               do {
                  ecc = sqrt(randinter(0., 1.));
               }
               while (ecc<e_min || ecc>=e_max);
            }
	    else { //if (local_type==Unknown_Binary_Type)
//               type = local_type;
//            else {
               id = bi->get_index();
               sma = a_min*exp(randinter(0., a_const));
               do {
                   ecc = sqrt(randinter(0., 1.));
               }
               while (ecc<e_min || ecc>=e_max);
            }
	    if (REPORT_ADD_DOUBLE) 
	      cerr<<"adddouble: iae: "<<id<<" "<<sma<<" "<<ecc<<endl;
	    
            double_star* new_double = new_double_star(bi, sma, ecc, stellar_time, id,
                         type);

	    // added spz:8Aug2002
	    new_double->dump("SeBa.data", true);

	    if (REPORT_ADD_DOUBLE) {
	      cerr<<"double: "<<new_double<<endl;
	      new_double->dump(cerr);
	      put_state(make_state(new_double));
	    }
         }
      }
   }
#else

void main(int argc, char ** argv) {

    bool F_flag = false;
    bool P_flag = false;
    bool U_flag = false;
    bool G_flag = false;
    char *mfc = new char[64];
    mass_function mf = mf_Power_Law;
    real m_min = 0.1;
    real m_max = 100;
    real m_exp = -2.35;
    char *qfc = new char[64];
    mass_ratio_distribution qf = Flat_q;
    real q_min = 0;
    real q_max = 1;
    real q_exp = 0;
    char *afc = new char[64];
    sma_distribution af = sma_Power_Law;
    real a_min = 0;
    real a_max = 1.e+6; 
    real a_exp = -1;
    char *efc = new char[64];
    ecc_distribution ef = Thermal_Distribution;
    real e_min = 0;    // allow detection of constant eccentricity
    real e_max = 1;
    real e_exp;

    int n = 1;

    int random_seed = 0;
    char seedlog[64];

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "M:m:x:F:f:A:a:y:G:g:E:e:v:U:u:Q:q:w:P:p:n:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1)
	switch(c) {
            case 'M': m_max = atof(poptarg);
		      break;
            case 'm': m_min = atof(poptarg);
		      break;
            case 'x': m_exp = atof(poptarg);
		      break;
	    case 'F': F_flag = true;
		      strcpy(mfc, poptarg);
	              break;
	    case 'f': mf = (mass_function)atoi(poptarg);
	              break;
            case 'A': a_max = atof(poptarg);
		      break;
            case 'a': a_min = atof(poptarg);
		      break;
            case 'y': a_exp = atof(poptarg);
		      break;
	    case 'G': G_flag = true;
		      strcpy(afc, poptarg);
	              break;
	    case 'g': af = (sma_distribution)atoi(poptarg);
	              break;
            case 'E': e_max = atof(poptarg);
		      break;
            case 'e': e_min = atof(poptarg);
		      break;
            case 'v': e_exp = atof(poptarg);
		      break;
	    case 'U': U_flag = true;
		      strcpy(efc, poptarg);
	              break;
	    case 'u': ef = (ecc_distribution)atoi(poptarg);
	              break;
            case 'Q': q_max = atof(poptarg);
		      break;
            case 'q': q_min = atof(poptarg);
		      break;
            case 'w': q_exp = atof(poptarg);
		      break;
	    case 'P': P_flag = true;
		      strcpy(qfc, poptarg);
	              break;
	    case 'p': qf = (mass_ratio_distribution)atoi(poptarg);
	              break;
	    case 'n': n = atoi(poptarg);
	              break;
	    case 's': random_seed = atoi(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    int actual_seed = srandinter(random_seed);
    printf("random number generator seed = %d\n",actual_seed);

//    if (binary_fraction < 0 || binary_fraction > 1)
//	err_exit("mkbinary: Illegal binary fraction");

    if(F_flag)
	mf = extract_mass_function_type_string(mfc);
    delete mfc;
    if(G_flag)
	af = extract_semimajor_axis_distribution_type_string(afc);
    delete afc;
    if(U_flag)
	ef = extract_eccentricity_distribution_type_string(efc);
    delete efc;
    if(P_flag)
	qf = extract_mass_ratio_distribution_type_string(qfc);
    delete qfc;

    real m_prim, m_sec, sma, ecc;
    cerr << "\tM\tm\ta\te"<<endl;
    for (int i=0; i<n; i++) {
	mkrandom_binary(m_min, m_max, mf, m_exp, 
			q_min, q_max, qf, q_exp, 
			a_min, a_max, af, a_exp, 
			e_min, e_max, ef, e_exp, 
			m_prim, m_sec, sma, ecc);
	cout << "\t" << m_prim 
	     << "\t" << m_sec 
	     << "\t" << sma
	     << "\t" << ecc << endl;
    }

}
#endif

