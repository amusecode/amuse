//
////  semi_to_period.C: transformed semi-major axis into orbital period
////                    and vice versa.
////
////	       Simon Portegies Zwart  & Gijs Nelemans, Tokyo, Sept 1998
////
////
#include "double_support.h"
#include "starlab_constants.h"

#ifndef TOOLBOX

real period_to_semi(real period,
		    real m_prim,
		    real m_sec) {

  real semi = pow(period*cnsts.physics(seconds_per_day), 2)
	    *     cnsts.physics(G)*cnsts.parameters(solar_mass)
	    *     (m_prim+m_sec)
            /     (4*pow(cnsts.mathematics(pi), 2));
  semi = pow(semi, 1./3.)
       / cnsts.parameters(solar_radius);
  
  return semi;
}

real semi_to_period(real semi,
		    real m_prim,
		    real m_sec) {

  real period = sqrt((4*pow(cnsts.mathematics(pi), 2)
		       *pow(cnsts.parameters(solar_radius)*semi, 3))
		       /(cnsts.physics(G)*cnsts.parameters(solar_mass)
			 *(m_prim+m_sec))); 

  return period/cnsts.physics(seconds_per_day);
}

#else

int main(int argc, char ** argv) {

    bool P_flag = false;
    bool q_flag = false;

    real m_prim = 13.1;
    real m_sec  =  9.8;
    real mass_ratio = m_sec/m_prim;

    real semi    = 138;
    real period = semi_to_period(semi, m_prim, m_sec);
    
    extern char *poptarg;
    int c;
    const char *param_string = "a:P:M:m:q:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.5 $", _SRC_)) != -1)
	switch(c) {
            case 'a': semi = atof(poptarg);
                      break;
            case 'P': P_flag = true;
	              period = atof(poptarg);
                      break;
            case 'M': m_prim = atof(poptarg);
                      break;
            case 'm': m_sec = atof(poptarg);
                      break;
            case 'q': q_flag = true;
	              mass_ratio = atof(poptarg);
                      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    if (q_flag) 
      m_sec = mass_ratio*m_prim;

    if (P_flag) {
      semi = period_to_semi(period, m_prim, m_sec);
      cerr << "(P [Days]; M, m [Msun]) = ("
	   << period << "; " << m_prim << ", " << m_sec << ")";
      cerr << "  ===>  a = " << semi << " [Rsun]." << endl;
    }
    else {
      period = semi_to_period(semi, m_prim, m_sec);
      cerr << "(a [Rsun]; M, m [Msun]) = ("
	   << semi << "; " << m_prim << ", " << m_sec << ")";
      cerr << "  ===>  P = " << period << " [Days]." << endl;
    }
    return 0;
}

#endif
