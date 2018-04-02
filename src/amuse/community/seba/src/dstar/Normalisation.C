//// SeBa:        Binary evolution program SeBa.
////              computes the evolution of a binary given any 
////              initial conditions (M, m, a, e).
////             
//// Output:      in the form of the following files:
////              -init.dat           contains selected initial conditons.
////              -SeBa.data           contains binary evolution histories
////              -binev.data          contains remnant formation information
////
////              in the form of standard output (cerr).
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
////            -x    mass function exponent in case of power law [-2.35] or mean mass 
////                                                                for log-normal IMF
////            -d        sigma for log-normal IMF
////            -F/f  mass function option: 0) Equal mass
////                                        1) Power-law [default]
////                                        2) Miller & Scalo
////                                        3) Scalo
////                                        4) Kroupa/
////                                        5) GdeMarchi
////                                        6) Kroupa, Tout & Gilmore 1991
////                                        7) TwoComponent (uses -h, -l and -u)
////                                        8) LogNormal (uses -x (for mu) and -d)
////            Option -F requires one of the following strings:
////                      (mf_Power_Law, Miller_Scalo, Scalo, Kroupa)
////                   -f requires the appropriate interger (see mkmass.C)
////             -A   maximum semi-major axis limit [1000000[Rsun]]   
////             -a   minimum semi-major axis limit [0] 
////             -y   exponent for a power-law distribution [0] (flat in log)
////             -G/g Semi major axis option: 0) Equal_sma
////                                          1) Power Law [default]
////                                          2) Duquennoy & Mayor (1987)
////                                          3) Raghavan (2010)
////                                          4) Eggleton (1999)
////            Option -G requires one of the following strings:
////                      (Equal_sma, sma_Power_Law, Duquennoy_Mayor, Raghavan, Eggleton)
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
////             -q   minimum mass ratio [0]
////                    extra option: q_min<0 : q_min=0.1/selected primary mass
////             -w   exponent for a power-law distribution  
////             -P/p mass ratio option: 0) constant mass ratio
////                                       1) Flat_q
////                                       2) Power Law 
////                                       3) Hogeveen (1992)
////            Option -P requires one of the following strings:
////                      (Equal_q, Flat_q, qf_Power_Law, Hogeveen)
////                   -p requires appropriate interger (see double_star.h)
////
////            -I select input file for reading initial conditions.
////               -uses: double_star::dump as input format.  [no default]
////            -R select random initial conditions    [false]
////               with parameters as discribed above.   
////            -n number of binaries to be simulated.  [1]
////               Options: -I all binaries in input file are computed.
////                        -R the number of binaries indicated.
////                        oterwise one binary is simulated with
////                        -M, -m, -a, -e as initial conditions.
////            -T or -t  binary end time. [13500] Myr
////            -s Random seed
////            -z select metallicity of binaries to be simulated. [0.02] Solar
////
//   Note:  libnode.a is referenced for the routines which produce the 
//          mass function
//
//	version 1.0	Simon Portegies Zwart, Utrecht, 1992
//                      -First version with class structure
//	version 2.0	Simon Portegies Zwart, Utrecht, 1994
//                      -Coupling to starlab
//	version 3.0	Simon Portegies Zwart, Amsterdam, June 1997
//	version 3.3	Simon Portegies Zwart, Cambridge, March 1999
//

#include "dyn.h" 
#include "double_star.h"
#include "main_sequence.h"
//#include "dstar_to_dyn.h"
//#include "seba.h"


#ifdef TOOLBOX

local bool read_binary_params(ifstream& in, real &m_prim, 
			      real &m_sec, real &sma, real &ecc, real &z) {

    m_prim = 1000;
    m_sec = 1000;
    sma = 0;
    ecc = 0;
    z = 1;
    while (m_prim>100.0 || m_sec>100.0 || ecc<0 || ecc>1 || z < 0.0001 || z > 0.03){   
      if(in.eof())
        return false;
    
      // reading from input file
    //  int id, tpp, tps; 
    //  real time, Rlp, Rls;
    
      in >> sma >> ecc >> m_prim >> m_sec >> z;
      //  in >>id >> time >> sma >> ecc 
      //     >> tpp >> m_prim >> Rlp >> tps >> m_sec >> Rls;
    }
    
    PRC(m_prim);PRC(m_sec);PRC(sma);PRC(ecc);PRL(z);
  return true;
}

/*-----------------------------------------------------------------------------
 *  binev  --
 *-----------------------------------------------------------------------------
 */
local bool  evolve_binary(dyn * bi,
                          real start_time, real end_time,
			  bool stop_at_merger_or_disruption,
			  bool stop_at_remnant_formation,
			  char* SeBa_outfile) { 


  double_star* ds = dynamic_cast(double_star*, 
				 bi->get_starbase());
  
  //		Setup star from input data.
  real dt, time=start_time;
  ds->evolve_element(time);
  
  //char SeBa_outfile[] = "SeBa.data";
  //char SeBa_outfile[] = SeBa_Filename();
  ds->dump(SeBa_outfile, true);

  if (!bi->is_root() &&
      bi->get_parent()->is_root()) 

    do {

      //dt = ds->get_evolve_timestep() + cnsts.safety(minimum_timestep);
      dt =
      Starlab::max(ds->get_evolve_timestep(),cnsts.safety(minimum_timestep));

      time = Starlab::min(time+dt, end_time);

      ds->evolve_element(time);

      if (stop_at_merger_or_disruption &&
	  (ds->get_bin_type() == Merged || 
	   ds->get_bin_type() == Disrupted))
	return false;

      if (stop_at_remnant_formation &&
	 (ds->get_primary()->remnant() || ds->get_secondary()->remnant()))
	return false;

    }
    while (time<end_time);

  ds->dump(SeBa_outfile, true);
  ds->set_star_story(NULL);
    
  rmtree(bi, false);
  return true;

}

int main(int argc, char ** argv) {

    bool e_flag = false;
    bool R_flag = false;
    bool F_flag = false;
    bool I_flag = false;
    bool O_flag = false;//doesn't work fully as "SeBa.data" is written in the code multiple times
    bool P_flag = false;
    bool U_flag = false;
    bool G_flag = false;

    bool stop_at_merger_or_disruption = false;
    bool stop_at_remnant_formation = false;
    bool random_initialization = false;
    real binary_fraction = 1.0;
    real metal = cnsts.parameters(Zsun);


    int n_init = 0;
    int n = 1;

    real  m_prim;
    real  m_sec;
    real  sma;
    real  ecc;

    real M_p = 0;
    real M_s = 0;
    real M_tot = 0;
    real M_mean_single = 0;
    real M_mean_double = 0;
    int n_run = 0;

    real M_p_all = 0;
    real M_s_all = 0;
    real M_tot_all = 0;
    real M_mean_single_all = 0;
    real M_mean_double_all = 0;

    char *mfc = new char[64];
    mass_function mf = mf_Power_Law;
    real m_min = 0.1, m_min_std = 0.08;
    real m_max = 100, m_max_std = 100;
    real m_exp = -2.35;
    real sigma = 0;
    char *qfc = new char[64];
    mass_ratio_distribution qf = Flat_q;
    real q_min = 0, q_min_std = 0;
    real q_max = 1, q_max_std = 1;
    real q_exp = 0;
    char *afc = new char[64];
    sma_distribution af = sma_Power_Law;
    real a_min = 0, a_min_std = 0;
    real a_max = 1.e+6, a_max_std = 1.e+6; 
    real a_exp = -1;                           
    char *efc = new char[64];
    ecc_distribution ef = Thermal_Distribution;
    real e_min = -1, e_min_std = -1; // allow detection of const eccentricity
    real e_max = 1, e_max_std = 1;
    real e_exp;

    real start_time = 0;
    real end_time   = 1e4;

    char* input_filename;
    char* output_filename;
    output_filename = "SeBa.data";

    int input_seed=0; 
    char seedlog[64];
    char paramlog[90];

    //check_help();

//    seba_counters* new_seba_counters = new seba_counters;

    extern char *poptarg;
    int c;
    const char *param_string = "n:N:RDSM:m:x:d:F:f:A:a:y:G:g:E:e:v:U:u:Q:q:T:t:I:O:w:P:p:n:s:z:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c) {
            case 'R': random_initialization = true;
		      break;
            case 'D': stop_at_merger_or_disruption = true;
		      break;
            case 'S': stop_at_remnant_formation = true;
		      break;
            case 'M': m_max = atof(poptarg);
		      break;
            case 'm': m_min = atof(poptarg);
		      break;
            case 'x': m_exp = atof(poptarg);
		      break;
            case 'd': sigma = atof(poptarg);
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
            case 't': 
            case 'T': end_time = atof(poptarg);
		      break;
            case 'I': I_flag = true;
		      input_filename = poptarg;
		      break;
            case 'O': O_flag = true;
		      output_filename = poptarg;
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
	    case 'N': n_init = atoi(poptarg);
	              break;
	    case 's': input_seed = atoi(poptarg);
		      break;
	    case 'z': metal = atof(poptarg);
	    	      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    int actual_seed = srandinter(input_seed);
//    cerr << "random number generator seed = " << actual_seed << endl;
    sprintf(paramlog, 
   "         alpha  = %3.1f\n         lambda = %3.1f\n         beta   = %3.1f\n         gamma  = %4.2f",
	    cnsts.parameters(common_envelope_efficiency),
	    cnsts.parameters(envelope_binding_energy),
	    cnsts.parameters(specific_angular_momentum_loss),
	    cnsts.parameters(dynamic_mass_transfer_gamma));


    if (n <= 0) err_exit("mknodes: N > 0 required!");

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
    if (q_min < 0){
        q_min_std = q_min;
    }

    actual_seed = srandinter(input_seed);
    sprintf(seedlog, "         random number generator seed = %d",actual_seed);

    ifstream infile(input_filename, ios::in);
    if(I_flag) {
      if (!infile) cerr << "error: couldn't create file "
	                << input_filename <<endl;
	cerr << "Reading input from file "<< input_filename <<endl;
    }



    print_initial_binary_distributions(m_min, m_max, mf, m_exp, //sigma,
			       q_min, q_max, qf, q_exp,
			       a_min, a_max, af, a_exp,
			       e_min, e_max, ef, e_exp);


    for (int i=0; i<n; i++) {

      if(I_flag) {
	if(read_binary_params(infile, m_prim, m_sec, sma, ecc, metal)) 
	  n=i+2;
	else
	  break;

      }
      else if (random_initialization) 
	mkrandom_binary(m_min_std, m_max_std, mf, m_exp, //sigma,
			q_min_std, q_max_std, qf, q_exp,
			a_min_std, a_max_std, af, a_exp,
			e_min_std, e_max_std, ef, e_exp,
			m_prim, m_sec, sma, ecc, metal);
      else {
	m_prim = m_max;
	m_sec  = m_min;
	sma    = a_min;
	ecc    = e_min;
	n = 1;
      }

      M_p_all += m_prim;
      M_s_all += m_sec;
      M_tot_all += m_prim + m_sec;
      
      real q = m_sec/m_prim;

      if (m_prim <= m_max && m_prim >= m_min &&
	  q      <= q_max && q      >= q_min &&
	  sma    <= a_max && sma    >= a_min &&
	  ecc    <= e_max && ecc    >= e_min ) {
	
	n_run ++;
	M_p += m_prim;
	M_s += m_sec;
	M_tot += m_prim + m_sec;
      }
      

    }

    M_mean_single_all = M_p_all/n;
    M_mean_double_all = M_tot_all/n;

    M_mean_single = M_p/n_run;
    M_mean_double = M_tot/n_run;

    cout << "In total mass spectrum (0.08 -- 100 Msun): " <<endl;
    cout << "Total number of binaries  : " << n << endl;
    cout << "Total mass in primaries   : " << M_p_all << endl;
    cout << "Total mass in secondaries : " << M_s_all << endl;
    cout << "Total mass in binaries    : " << M_tot_all << endl;
    cout << "Mean mass of primaries    : " << M_mean_single_all << endl;
    cout << "Mean mass of binaries     : " << M_mean_double_all << endl;

    cout << "In run (m_min = "<< m_min << ", m_max = "<< m_max<< ")"<<endl; 
    cout << "Total number of binaries  : " << n_run << " fraction: " 
	                                   << 1.*n_run/n <<endl;
    cout << "Total mass in primaries   : " << M_p << endl;
    cout << "Total mass in secondaries : " << M_s << endl;
    cout << "Total mass in binaries    : " << M_tot << endl;
    cout << "Mean mass of primaries    : " << M_mean_single << endl;
    cout << "Mean mass of binaries     : " << M_mean_double << endl;

}


#endif

