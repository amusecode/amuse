//// rs_population: Reads SeBa.data standard output format and 
////                integrates over time to find current 
////                number of binaries in the Galaxy
////          
//// Options:       no options
////
//// 
////
//----------------------------------------------------------------------------
//   version 1:  July 2001  Simon Portegies Zwart   spz@mit.edu
//............................................................................
//   non-local functions: 
//----------------------------------------------------------------------------


#include "SeBa_hist.h"

#define MAX_NSCENARIOS 5000

#ifdef TOOLBOX

bool binary_contains(char a[], char b[]) {
  bool identical = false;
  if(!strcmp(a, b))
    identical = true;

  return identical;
}

bool SeBa_hist_contains(SeBa_hist *hi, SeBa_hist *ha) {

//  cerr << "SeBa_hist_contains(SeBa_hist *hi, SeBa_hist *ha)" << endl;
//  PRI(10);PRL(*hi);
//  PRI(10);PRL(*ha);

  bool same_primary = false;
  bool same_secondary = false;
  for_all_SeBa_hist(SeBa_hist, hi->get_first(), ho) {
    for_all_SeBa_hist(SeBa_hist, ha->get_first(), he) {

//      PRC(ho->get_label_prim());PRL(he->get_label_prim());

      if(binary_contains(ho->get_label_prim(), he->get_label_prim())) {
	same_primary = true;
      }
      if(binary_contains(ho->get_label_sec(), he->get_label_sec())) {
	same_secondary = true;
      }
    }
  }

//  PRC(same_primary);
//  PRL(same_secondary);
  return same_primary;
}


local int reorder_binaries(SeBa_hist *hi, const real end_time) {

  int N_lread = 0;
  int N_scenarios = 0;

  SeBa_hist *scenarios[MAX_NSCENARIOS];
  int scenario_frequency[MAX_NSCENARIOS];

  for(int i=0; i< MAX_NSCENARIOS; i++) {
    scenarios[i] = NULL;
    scenario_frequency[i] = 0;
  }

  scenarios[N_scenarios] = hi;
  scenario_frequency[N_scenarios]++;
  N_scenarios++;

//  PRL(*hi);

  bool new_entry;
  do {
//    for(int i = 0; i<N_scenarios; i++) {
//      PRC(i);PRC(scenarios[i]);PRL(scenarios[i]->get_last());
//    }

    hi = new SeBa_hist();
    hi->read_SeBa_hist(cin);

    N_lread++;

    new_entry = true;
    for(int i = 0; i<N_scenarios; i++) {
      
      if(SeBa_hist_contains(scenarios[i], hi)) {

	if(!new_entry) {
	  cerr << "Same star is member of more thanone binary"<<endl;
	}
	new_entry = false;
	SeBa_hist *last = scenarios[i]->get_last();
	scenarios[i]->set_last(hi);
	hi->set_past(last);
	hi->set_number(scenarios[i]->get_number());

//	put_state(scenarios[i], cerr);
	scenario_frequency[i]++;
//	break;
      }
    }

    if(new_entry) {
      if(hi->get_number()==-1)
	hi->set_number(N_scenarios);

      scenarios[N_scenarios] = hi;
      scenario_frequency[N_scenarios]++;
      
      N_scenarios++;
      if(N_scenarios >= MAX_NSCENARIOS)
	break;
    }    

  }
  while (!cin.eof());

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
//    cerr << i << "\t" << scenario_frequency[i] << ":  "; 
//    put_state(scenarios[i], cerr);
    for_all_SeBa_hist(SeBa_hist, scenarios[i], ha) {
      cout << *ha;
    }
  }

  cerr << "Binaries read = " << N_lread-1 << endl;
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

    int N_binaries = reorder_binaries(hi, end_time); 

    cerr << "Total number of binaries found: " << N_binaries << endl;
}

#endif // endof: TOOLBOX



