
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add a mass spectrum to an input snapshot.  Existing node masses
//// are overwritten.
////
//// If only one mass limit is set, the other is automatically forced
//// to the same value.  In the TwoComponent case, the lower and upper
//// limits specify the two component masses.
////
//// Usage: makemass [OPTIONS] < input > output
////
//// Options:
////         -C        output in 'col' format (dyn version only) [no]
////         -e/E/x/X  exponent [-2.35 (Salpeter)]
////         -F/f      mass function option:
////                        1) Power-law [default]
////                        2) Miller & Scalo
////                        3) Scalo
////                        4) Kroupa
////                        5) GdeMarchi
////                        6) Kroupa, Tout & Gilmore 1991
////                        7) TwoComponent (uses -h, -l and -u)
////                   option -F requires one of the following strings:
////                        Power_Law, Miller_Scalo, Scalo, Kroupa,
////                        GdeMarchi, KTG91, TwoComponent
////                   option -f requires the appropriate integer.
////         -h/H      fraction of stars in high-mass group (TwoComponent) [0]
////         -i        (re)number stellar index from highest to lowest mass.
////         -l/L      lower mass limit [1]
////         -m/M      scale to specified total mass [don't scale]
////         -u/U      upper mass limit [1]
////         -s        random seed
////
//// Written by Steve McMillan and Simon Portegies Zwart.
////
//// Report bugs to starlab@sns.ias.edu.

//++ Note: The conversion factor for scaling between dynamical and stellar
//++       mass is properly set in the output snapshot.

//		Steve McMillan, July 1996
//		Simon Portegies Zwart, Tokyo, December 1997

// Allow double use -- dyn tool will have TOOLBOX and DYN set and will
// include this file.  The only real difference is the "-C" option.

#ifndef DYN
#   include "node.h"
#   include "util_math.h"
#else
#   include "dyn.h"
#endif

#define  MAXIMUM_ZAMS_MASS 100
#define  SEED_STRING_LENGTH  256
static char  tmp_string[SEED_STRING_LENGTH];

#define TOL 1.e-4

#ifndef TOOLBOX

// Note from Steve (5/2005: These sampling methods do very poorly when m_lower
// and m_upper are similar, as the target interval is very small.  This can
// happen by accident, e.g. if a MF is specified but no range is given, which
// defaults to [1,1] and an infinite loop.

// See Eggleton, Fitchet & Tout 1989, ApJ 347, 998.

local real mf_Miller_Scalo(real m_lower, real m_upper) {

    if (twiddles(m_lower, m_upper, TOL))
	return 0.5*(m_lower+m_upper);		// Steve, 5/05

    real m, rnd;
    do {
	rnd = randinter(0,1);
	m = 0.19*rnd
	    / (pow(1-rnd, 0.75) + 0.032*pow(1-rnd, 0.25));
    }
    while (m < m_lower || m > m_upper);
    return m;
}

#if 0

// Tapered power-law (from Guido de Marchi private communication April. 2002)
// dN/dLogm= m^(-1.35)*(1-exp(-m/0.25)^2.4)

local real tapered_power_law(real m, 
			     const real x1 = -1.35, 
			     const real x2 = 2.4, 
			     const real mbreak = 0.25) {

    if (twiddles(m_lower, m_upper, TOL))
	return 0.5*(m_lower+m_upper);		// Steve, 5/05

    real m, rnd;
    do {
	rnd = randinter(0,1);
	//	m = 0.08 + (0.19*pow(rnd, 1.55) + 0.05*pow(rnd, 0.6))
	//	  /        pow(1-rnd, 0.58);
	real X = 1-rnd;
	m = 0.15 * (1 /(pow(X, 0.75) + 0.04*pow(X, 0.25)) - pow(X, 2)/1.04);
    }
    while (m < m_lower || m > m_upper);
    return m;
}
#endif

#if 0

// Tapered power-law (from Guido de Marchi private communication April. 2002)
// dN/dLogm= m^(-1.35)*(1-exp(-m/0.25)^2.4)

local real tapered_power_law(real m, 
			     const real x1 = -1.35, 
			     const real x2 = 2.4, 
			     const real mbreak = 0.25) {

    if (twiddles(m_lower, m_upper, TOL))
	return 0.5*(m_lower+m_upper);		// Steve, 5/05

  real x = randinter(1, 91.43);
  real ms = 0;
  if(x < 1) ms = randinter(5.62, 100);
  else if(x < 2) ms = randinter(3.55, 5.62);
  else if(x < 5.16) ms = randinter(2.24, 3.55);
  else if(x < 9.14) ms = randinter(1.41, 2.24);
  else if(x < 15.17) ms = randinter(0.89, 1.41);
  else if(x < 21.93) ms = randinter(0.56, 0.89);
  else if(x < 37.07) ms = randinter(0.35, 0.56);
  else if(x < 57.03) ms = randinter(0.22, 0.35);
  else if(x < 73.63) ms = randinter(0.14, 0.22);
  else if(x < 85.93) ms = randinter(0.09, 0.14);
  else if(x < 91.43) ms = randinter(0.04, 0.09);
  else ms = randinter(0.0, 0.04);
  
  // 562   3.55   2.24   1.41   0.89   0.56   0.35   0.22   0.14   0.09  0.04
  //     |      |      |      |      |      |      |      |      |      |
  //    1.00   3.16   3.98   6.03   6.76  15.14  19.96  16.60  12.30   5.5
  // 1    2.00   5.16   9.14   15.17  21.93 37.07  57.03  73.63  85.93  91.43

  return ms;
}
#endif

local real mf_GdeMarchi(real m_lower, real m_upper) {

    if (twiddles(m_lower, m_upper, TOL))
	return 0.5*(m_lower+m_upper);		// Steve, 5/05

    real m, X;
    do {
	X = randinter(0,1);
	// better fit to Guido's results at low mass end, 
        // but lacks high mass component.
	//m = 0.15 * (1 /(pow(X, 0.9) + 0.03*pow(X, 0.2)) - pow(X, 8)/1.03);

	// most satisfactory
	m = 0.15 * (1 /(pow(X, 0.75) + 0.03*pow(X, 0.25)) - pow(X, 8)/1.03);
    }
    while (m < m_lower || m > m_upper);
    return m;

#if 0

//dN/dLogm= m^(-1.35)*(1-exp(-m/0.15)^2.4)
// for now we do it the poor way by Monte Carlo

  real x1, x2, mbreak;
    real mass = tapered_power_law(m_lower, x1, x2, mbreak);
    return mass;

  real x1 = -2.35;
  real x2 = 2.40;
  real mbreak = 0.15;
    real dNdm_min = tapered_power_law(m_lower, x1, x2, mbreak);
    real dNdm_max = tapered_power_law(m_upper, x1, x2, mbreak);
//    PRC(dNdm_min);PRL(dNdm_max);
    real m, rnd;
    real dNdm_try, dNdm;
    do {
	rnd = randinter(0,1);
	m = m_lower + rnd*(m_upper-m_lower);
	dNdm_try = randinter(dNdm_min,dNdm_max);
	dNdm = tapered_power_law(m, x1, x2, mbreak);
    }
    while(dNdm_try>dNdm);
    return m;
#endif

}

// See de la Fuente Marcos, Aarseth, Kiseleva, Eggleton
// 1997, in Docobo, Elipe, McAlister (eds.), Visual Double
// Stars: Formation, dynamics and Evolutionary Tracks, KAP: ASSL
// Series vol. 223,  165.

local real mf_Scalo(real m_lower, real m_upper) {

    if (twiddles(m_lower, m_upper, TOL))
	return 0.5*(m_lower+m_upper);		// Steve, 5/05

    real m, rnd;
    do {
	rnd = randinter(0,1);
	m = 0.3*rnd
	    / pow(1-rnd, 0.55);
    }
    while (m < m_lower || m > m_upper);
    return m;
}

local real Kroupa_Tout_Gilmore(real m_lower, real m_upper) {

    if (twiddles(m_lower, m_upper, TOL))
	return 0.5*(m_lower+m_upper);		// Steve, 5/05

    real m, rnd;
    do {
	rnd = randinter(0,1);
	m = 0.08 + (0.19*pow(rnd, 1.55) + 0.05*pow(rnd, 0.6))
	  /        pow(1-rnd, 0.58);
	//real X = 1-rnd;
	//m = 0.15 * (1 /(pow(X, 0.75) + 0.04*pow(X, 0.25)) - pow(X, 2)/1.04);
    }
    while (m < m_lower || m > m_upper);
    return m;
}

local real Kroupa_Tout_Gilmore_1991(real m_lower, real m_upper) {

  real m, rnd;
  do {
    rnd = randinter(0,1);
    real X = 1-rnd;
    m = 0.33 * (1 /(pow(X, 0.75) + 0.04*pow(X, 0.25)) - pow(X, 2)/1.04);
  }
  while(m_lower>m || m>m_upper);
  return m;
}

static int ntot = 0, n_heavy = 0;

local real mf_TwoComponent(real m_lower, real m_upper, real heavyfrac)
{
    // Random choice between low and high mass, except that we do't
    // let the actual fraction deviate too far from the target.

    if (randinter(0,1) < heavyfrac) {
	ntot++;
	n_heavy++;
	if (n_heavy > (int)(heavyfrac*ntot) + 2) {
	    n_heavy--;
	    return m_lower;
	} else
	    return m_upper;
    } else {
	ntot++;
	if (n_heavy < (int)(heavyfrac*ntot) - 2) {
	    n_heavy++;
	    return m_upper;
	} else
	    return m_lower;
    }
}


//----------------------------------------------------------------------
// Global functions.
//----------------------------------------------------------------------

real general_power_law(real lowerl, real upperl, real exponent) {

    real x, norm;

    // Normalizing factor.

    if (exponent == -1)
	norm = log(upperl/lowerl);
    else
	norm = pow(upperl/lowerl, 1+exponent) - 1;

    if (exponent == -1)
	x = lowerl*exp(norm*randinter(0,1));
    else
	x = lowerl*pow(norm*randinter(0,1) + 1, 1/(1+exponent));

    return x;
}

real get_random_stellar_mass(real m_lower, real m_upper,
			     mass_function mf, real exponent) {

    real m;
    switch(mf) {
       case Equal_Mass:
	   if (m_lower==0) {
	       cerr << "get_random_stellar_mass:"<<endl;
	       cerr << "unambiguous choise of Equal_Mass."<<endl;
	       cerr << "Use -m option to set fixed stellar mass."<<endl;
	       exit(1);
	   }
	   m =  m_lower;
	   break;
       case mf_Power_Law:
	   m =  general_power_law(m_lower, m_upper, exponent);
	   break;
       case Miller_Scalo:
	   m = mf_Miller_Scalo(m_lower, m_upper);
	   break;
       case Scalo:
	   m = mf_Scalo(m_lower, m_upper);
	   break;
       case Kroupa: // Kroupa, Tout & Gilmore 1993, MNRAS 262, 545
	   m = Kroupa_Tout_Gilmore(m_lower, m_upper);
	   break;
       case GdeMarchi:
	  m = mf_GdeMarchi(m_lower, m_upper);
	      break;
       case KTG91:
	  m = Kroupa_Tout_Gilmore_1991(m_lower, m_upper);
       case TwoComponent:
	   m = mf_TwoComponent(m_lower, m_upper, exponent);   // use exponent to
							      // store heavyfrac
	   break;
       default:
	   cerr << "WARNING: \n"
		<< "        real get_random_stellar_mass:\n"
		<< "        Mass function parameters not properly defined.\n";
	   exit(1);
    }

    return m;
}

char* type_string(mass_function mf) {

    local char  mf_name[SEED_STRING_LENGTH];	// permanent
    switch(mf) {

       case Equal_Mass:
            sprintf(mf_name, "Equal_Mass");
	    break;	
       case mf_Power_Law:
            sprintf(mf_name, "Power_Law");
	    break;	
       case Miller_Scalo:
            sprintf(mf_name, "Miller_Scalo");
	    break;	
       case Scalo:
            sprintf(mf_name, "Scalo");
	    break;	
       case Kroupa:
            sprintf(mf_name, "Kroupa");
	    break;
       case GdeMarchi:
            sprintf(mf_name, "GdeMarchi");
	    break;
       case KTG91:
            sprintf(mf_name, "KTG91");
	    break;
       case TwoComponent:
            sprintf(mf_name, "TwoComponent");
	    break;
       default:
            sprintf(mf_name, "Unknown");
	    break;
    }
    return mf_name;
}

mass_function extract_mass_function_type_string(char* type_string) {

     mass_function type = Unknown_MF;

     if (!strcmp(type_string, "Equal_Mass"))
        type = Equal_Mass;
     else if (!strcmp(type_string, "Power_Law"))
        type = mf_Power_Law;
     else if (!strcmp(type_string, "Miller_Scalo"))
        type = Miller_Scalo;
     else if (!strcmp(type_string, "Scalo"))
        type = Scalo;
     else if (!strcmp(type_string, "Kroupa"))
        type = Kroupa;
     else if (!strcmp(type_string, "GdeMarchi"))
        type = GdeMarchi;
     else if (!strcmp(type_string, "KTG91"))
        type = KTG91;
     else if (!strcmp(type_string, "TwoComponent"))
        type = TwoComponent;
     else if (!strcmp(type_string, "Unknown"))
        type = Unknown_MF;
     else {
	 cerr << "No proper mass function indicated in makemass."<<endl;
	 exit(-1);
     }

     return type;
 }


//-----------------------------------------------------------------------------
//  compare_mass  --  compare the masses of two particles
//-----------------------------------------------------------------------------

typedef  struct
{
    node* str;
    real  mass;
} nm_pair, *nm_pair_ptr;

local int compare_mass(const void * pi, const void * pj)
{
    if (((nm_pair_ptr) pi)->mass < ((nm_pair_ptr) pj)->mass)
        return(1);
    else if (((nm_pair_ptr)pi)->mass > ((nm_pair_ptr)pj)->mass)
        return(-1);
    else
        return(0);
}

void makemass(node* b, mass_function mf,
	      real m_lower, real m_upper,
	      real exponent, real total_mass, bool renumber_stars)
{
    int n_target = 0;
    if (mf == TwoComponent) {
	ntot = n_heavy = 0;
	for_all_daughters(node, b, bi) n_target++;
	real x_target = n_target*exponent;
	n_target = (int)x_target;
    }

    real m, m_sum = 0;
    int n = 0;
    for_all_daughters(node, b, bi) {
	m = get_random_stellar_mass(m_lower, m_upper, mf, exponent);
	n++;
	bi->set_mass(m);
	m_sum += bi->get_mass();
    }

    if (mf == TwoComponent) {
	// cerr << "makemass: "; PRC(n_heavy); PRL(n_target);
	while (n_heavy < n_target) {	// need to add more heavy stars
	    // cerr << "makemass: adding heavy star" << endl;
	    int ntop = 0;
	    for_all_daughters(node, b, bi)
		if (bi->get_mass() != m_upper) ntop++;
	    real f = 1.0/ntop;
	    for_all_daughters(node, b, bi)
		if (bi->get_mass() != m_upper && randinter(0,1) < f) {
		    bi->set_mass(m_upper);
		    n_heavy++;
		    break;
		}
	}
	// cerr << "makemass: "; PRC(n_heavy); PRL(n_target);
	while (n_heavy > n_target) {	// need to remove heavy stars
	    // cerr << "makemass: removing heavy star" << endl;
	    int ntop = 0;
	    for_all_daughters(node, b, bi)
		if (bi->get_mass() == m_upper) ntop++;
	    real f = 1.0/ntop;
	    for_all_daughters(node, b, bi)
		if (bi->get_mass() == m_upper && randinter(0,1) < f) {
		    bi->set_mass(m_lower);
		    n_heavy--;
		    break;
		}
	}
	cerr << "makemass: "; PRC(n_heavy); PRL(n_target);
    }

    b->set_mass(m_sum);

    // Renumber the stars in order of mass.
    // Highest mass gets smallest number (strange choice, but...).

    if (renumber_stars) {
	int istart = 1;
	bool M_flag = true;
	renumber(b, istart, M_flag);
    }

#if 0
    nm_pair_ptr nm_table = new nm_pair[n];
    if (nm_table == NULL) {
	cerr << "makemass: "
	     << "not enough memory left for nm_table\n";
	return;
    }

    int i=0;
    for_all_daughters(node, b, bi) {
	nm_table[i].str = bi;
	nm_table[i].mass = bi->get_mass();
	i++;
    }

    qsort((void *)nm_table, (size_t)n, sizeof(nm_pair), compare_mass);

    for (i=0; i<n; i++) {
	nm_table[i].str->set_index(i+1);  // Ok, lets number from 1 to n
    }

    delete []nm_table;
#endif

    if (total_mass > 0) {
	real m_factor = total_mass/m_sum;
	for_all_daughters(node, b, bi)
	    bi->set_mass(m_factor*bi->get_mass());
	b->set_mass(total_mass);
    }

    if(m_lower>=m_upper)
	sprintf(tmp_string,
		"         %s mass function, total mass = %8.2f",
		type_string(mf), m_sum);
//		type_string(Equal_Mass), m_sum);
    else
	sprintf(tmp_string,
		"       %s mass function, total mass = %8.2f Solar",
		type_string(mf), m_sum);
    b->log_comment(tmp_string);
    cerr << "makemass: mass function is " << type_string(mf) << ", ";
    PRC(m_lower); PRL(m_upper);
}

#else

int main(int argc, char ** argv)
{
#ifdef DYN
    bool C_flag = false;
#endif

    bool F_flag   = false;                      // input mf via string
    mass_function mf = mf_Power_Law;            // default = Power-law
    char *mfc = new char[64];
    real m_lower  = 1, m_upper = 1;		// default = equal masses
    bool x_flag   = false;
    real exponent = -2.35;			// Salpeter mass function;
                                                // if no exponent is given
                                                // Miller & Scalo is used
    real hfrac = 0;
    bool hfrac_set = false;

    real m_total = -1;				// don't rescale

    bool renumber_stars = false;                // renumber stellar index
                                                // from high to low mass

    int random_seed = 0;
    char seedlog[64];

    bool lower_set = false, upper_set = false;

    check_help();

    extern char *poptarg;
    int c;

#ifndef DYN
    const char *param_string = "E:e:iF:f:h:L:l:M:m:s:U:u:X:x:";
#else
    const char *param_string = "CE:e:iF:f:h:L:l:M:m:s:U:u:X:x:";
#endif

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.24 $", _SRC_)) != -1)
	switch(c) {

#ifdef DYN
	    case 'C': C_flag = true;
		      break;
#endif
	    case 'F': F_flag = true;
		      strncpy(mfc, poptarg, 63);
	              break;
	    case 'f': mf = (mass_function)atoi(poptarg);
	              break;
	    case 'i': renumber_stars = true;
	              break;
	    case 'E':
	    case 'e':
	    case 'X':
	    case 'x': x_flag = true;
	              exponent = atof(poptarg);
		      break;
	    case 'H':
	    case 'h': hfrac = atof(poptarg);
		      hfrac_set = true;
		      break;
	    case 'L':
	    case 'l': m_lower = atof(poptarg);
		      lower_set = true;
		      break;
	    case 'M':
	    case 'm': m_total = atof(poptarg);
		      break;
	    case 's': random_seed = atoi(poptarg);
		      break;
	    case 'U':
	    case 'u': m_upper = atof(poptarg);
		      upper_set = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	    	      get_help();
	    	      exit(1);
	}

    if (lower_set && !upper_set) m_upper = m_lower;
    if (!lower_set && upper_set) m_lower = m_upper;

    if (m_lower <= 0 ||
	m_upper <= 0 ||
	m_lower > m_upper)
	err_exit("makemass: Illegal mass limits");

    if (F_flag)
      mf = extract_mass_function_type_string(mfc);
    delete [] mfc;

    if (mf == TwoComponent) {
	if (!hfrac_set)
	    warning("makemass: heavy fraction not set; default = 0");
	exponent = hfrac;	// avoid redefining/overloading functions.
    }

#ifndef DYN
    node *b;
    b = get_node();
#else
    dyn *b;
    b = get_dyn();
#endif

    b->log_history(argc, argv);

    int actual_seed = srandinter(random_seed);

    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    makemass(b, mf, m_lower, m_upper, exponent, m_total, renumber_stars);
    real initial_mass = getrq(b->get_log_story(), "initial_mass");

    if (initial_mass > -VERY_LARGE_NUMBER)
        putrq(b->get_log_story(), "initial_mass", b->get_mass(),
	      HIGH_PRECISION);

    if (b->get_mass() != initial_mass) {

	// Erase any info on virial radius or total energy, as we
	// do not recompute these quantities here.

	if (find_qmatch(b->get_log_story(), "initial_total_energy"))
	    rmq(b->get_log_story(), "initial_total_energy");
	if (find_qmatch(b->get_log_story(), "initial_rvirial"))
	    rmq(b->get_log_story(), "initial_rvirial");
    }

    // Add stellar conversion data, so add_star will work.

    real m_sum = b->get_mass();
    real old_mtot = b->get_starbase()->conv_m_dyn_to_star(1);
    if(old_mtot<=0) {
	real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
	real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
	b->get_starbase()->set_stellar_evolution_scaling(m_sum,
							 old_r_vir,
							 old_t_vir);
    }

#ifndef DYN
    put_node(b);
#else
    b->set_col_output(false);
    if (C_flag) b->set_col_output(true);
    put_dyn(b);
#endif

    rmtree(b);
    return 0;
}
#endif

