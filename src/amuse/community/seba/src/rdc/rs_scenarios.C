//// rs_population: Reads SeBa.data standard output format and 
////                integrates over time to find current 
////                number of binaries in the Galaxy
////          
//// Options:   
////             -t     snapshot time at which the observer views [0] Myr
////             -x     characteristic time of the birthrate function in
////                    Myr. Zero for flat birthrate [0]
////                    positive for exponential decay.
////             -N     normalize output to 100% [false]
////             -c     normalization constant in number of systems per Myr 
////                    [1] 
////
//// 
////
//----------------------------------------------------------------------------
//   version 1:  Sept 1998   Simon Portegies Zwart   spz@grape.c.u-tokyo.ac.jp
//                                                   University of Tokyo
//   version 2:  July 2000  Simon Portegies Zwart   spz@mit.edu
//                          Gijs Nelemans           gijsn@astro.uva.nl
//............................................................................
//   non-local functions: 
//----------------------------------------------------------------------------


#include "SeBa_hist.h"

#define MAX_NSCENARIOS 5000

#ifdef TOOLBOX
local int extract_scenarios(SeBa_hist *hi, const real end_time) {

  int N_read = 0;

  SeBa_hist *scenarios[MAX_NSCENARIOS];
  int scenario_frequency[MAX_NSCENARIOS];

  for(int i=0; i< MAX_NSCENARIOS; i++) {
    scenarios[i] = NULL;
    scenario_frequency[i] = 0;
  }
  int N_scenarios = 0;

  SeBa_hist* next_hi = get_history(hi, cin);
  N_read++;
//    cerr << "Print next_hi"<<endl;
//    for_all_SeBa_hist(SeBa_hist, hi, ha) {
//      cerr << *ha;
//    }
//    cerr << "Printed"<<endl;

  scenarios[N_scenarios] = hi;
  scenario_frequency[N_scenarios]++;
  N_scenarios++;
//  delete hi;
  hi = next_hi;
  next_hi = get_history(hi, cin);
  N_read++;
  bool new_entry;
  do {

    new_entry = true;

    for(int i = 0; i<N_scenarios; i++) {
      
      if(scenarios_identical(hi, scenarios[i])) {
	new_entry = false;
	scenario_frequency[i]++;
	break;
      }
    }

    if(new_entry) {
      scenarios[N_scenarios] = hi;
      scenario_frequency[N_scenarios]++;
      
      N_scenarios++;
      if(N_scenarios >= MAX_NSCENARIOS)
	break;
    }    

    if(!new_entry)
      delete hi;

    hi = next_hi;
    next_hi = get_history(hi, cin);
    N_read++;
  }
  while (next_hi);

  int i=0;
#if 0
  cerr << "print all scenarios"<<endl;
  do {

    for_all_SeBa_hist(SeBa_hist, scenarios[i], ha)
      cerr << *ha;
    i++;
  }
  while(scenarios[i]!=NULL);
#endif

  cerr << "Scn. # \t freq." << endl;
  for(i=0; i<N_scenarios; i++) {
    cerr << i << "\t" << scenario_frequency[i] << ":  "; 
    put_state(scenarios[i], cerr);
  }

  cerr << "Binaries read = " << N_read-1 << endl;
  return N_scenarios;
} 

//----------------------------------------------------------------------------
//  main  --  driver to reduce SeBa short dump data
//----------------------------------------------------------------------------

main(int argc, char ** argv) {

    real end_time = VERY_LARGE_NUMBER;

    char  *comment;
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "T:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
            case 'T': end_time = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            


    SeBa_hist* hi = new SeBa_hist;
    if (!hi->read_SeBa_hist(cin))
      exit(-1);

    if(end_time<VERY_LARGE_NUMBER)
      cerr << "End time = " << end_time << " [Myr]" << endl << endl;

    int N_scenarios = extract_scenarios(hi, end_time); 


    cerr << "Total number of scenarios recognized: " << N_scenarios << endl;
}

#endif // endof: TOOLBOX



