//// rs_population: Reads SeBa.data standard output format and 
////                computes the birthrates of binaries in the Galaxy
////          
//// Options:   
////             -t     snapshot time at which the observer views [0] Myr
////             -x     characteristic time of the birthrate function in
////                    Myr. Zero for flat birthrate [0]
////                    positive for exponential decay.
////             -N     normalize output to 100% [N]
////             -S     Count supernova type [false]
////             -c     normalization constant in number of systems per Myr 
////                    [1] 
////
//// 
////
//----------------------------------------------------------------------------
//   version 1:  Sept 1998   Simon Portegies Zwart   spz@grape.c.u-tokyo.ac.jp
//                                                   University of Tokyo
//   version 2:  August 2000  Simon Portegies Zwart  spz@mit.edu
//............................................................................
//   non-local functions: 
//----------------------------------------------------------------------------


#include "SeBa_hist.h"

#ifdef TOOLBOX

real **mkrarray(int rows, int cols) {
    real ** ret = new real *[rows];
    for (int i=0; i<cols; i++)
	ret[i] = new real[cols];

    for (int i=0; i<rows; i++)
	for (int j=0; j<cols; j++)
	ret[i][j] = 0;

    return ret;
}

void print_binary_matrix(real **population) {

    real T_row = 0; 
    cerr << "    ";
    for (int j=0; j<no_of_stellar_type-1; j++) {
      cerr << " " << type_short_string((stellar_type)j);
      if(j==(int)Main_Sequence || 
	 j==(int)Super_Giant   || 
	 j==(int)Thorn_Zytkow)
        cerr << "  "; 
    }
    cerr << endl;
      
    for (int k = 0; k < no_of_stellar_type-1; k++) {
      T_row = 0;
      for (int i = 0; i < no_of_stellar_type-1; i++)
	T_row += population[i][k];
      
      if(T_row>0) {
	cerr << "   " << type_short_string((stellar_type)k);
	  
	for (int i = 0; i < no_of_stellar_type-1; i++) {
	  cerr << " " << population[i][k];
	  if(i==(int)Main_Sequence || 
	     i==(int)Super_Giant   || 
	     i==(int)Thorn_Zytkow)
	    cerr << "  "; 
	}
	cerr << endl;
      }
    }
  }

void print_binary_array(real **population) {

    real N_row = 0;

    int n_all = no_of_stellar_type*no_of_stellar_type;
    real *pri_sec = new real[n_all];
    for (int j=0; j<n_all; j++)
      pri_sec[j] = 0;

    int js, ks, ts;
    for (int j=0; j<no_of_stellar_type-1; j++)
      for (int k = 0; k < no_of_stellar_type-1; k++) {
        js = (int)stellar_type((stellar_type)j);
        ks = (int)stellar_type((stellar_type)k);
        if(js>ks) {
          ts = js;
          js=ks;
          ks = ts;
        }
        pri_sec[js + no_of_stellar_type*ks] += population[j][k];
      }

    for (js=0; js<no_of_stellar_type-1; js++)
      for (ks = 0; ks < no_of_stellar_type-1; ks++) {
        if(pri_sec[js + no_of_stellar_type*ks]>0)
          cerr << " (" << type_short_string((stellar_type)js)
               << ", " << type_short_string((stellar_type)ks) <<")\t"
               << pri_sec[js + no_of_stellar_type*ks] << endl;
      }

    delete []pri_sec;
  }

local real calculate_birthrate(real t_birth, real tau, real cnst) {

  real birthrate;
  if (tau==0) 
    birthrate = cnst;
  else
    birthrate = cnst * exp(-t_birth/tau);

  return birthrate;

}

local real contribution_to_population(real current_time, real t_start,
				      real tau, real cnst) {

  real t_birth = current_time - t_start;

  return calculate_birthrate(t_birth, tau, cnst);

}

local int extract_population(SeBa_hist *hi, real snap_time,
			     real tau, real imf_cnst,
			     real normalize, bool count_sn) {

    int N_sn = 0;
    int binaries_read = 0;
    real total_number_of_binaries = 0;

    real** detached_population;
    real** contact_population;
    real* single_population;
    real N_detached=0, N_contact=0, N_single=0;
    binary_type tpe_bin;
    stellar_type prim_type, sec_type, ts;
    real prim_mass, sec_mass;

      detached_population = mkrarray((int)no_of_stellar_type, 
				    (int)no_of_stellar_type);
      contact_population = mkrarray((int)no_of_stellar_type, 
				   (int)no_of_stellar_type);
      single_population = new real[(int)no_of_stellar_type];
      for(int i=0; i<(int)no_of_stellar_type; i++) 
	single_population[i]=0;

    int mergers = 0;	
    SeBa_hist* next_hi = get_history(hi, cin);
    do {

      binaries_read++;

      if (hi->get_last()->get_binary_type() == Merged) 
	mergers++;

//	SeBa_hist *hj = hi->get_SeBa_hist_at_time(snap_time);

      real tstart, tend, birthrate;
//      real current_time = 23;
//      PRC(snap_time);PRL(current_time);
      if (hi->get_last()->get_time() < snap_time) {
	cerr << "ERROR: snap_time not reached" << endl;
	exit(-1);
      }

      for_all_SeBa_hist(SeBa_hist, hi, hj) {

	// (SPZ+GN: 31 Jul 2000)
	// Evolution is stepped until last output.
	// Stop at t_max: i.e. hj->get_future()->get_time() <= t_max
	// and as soon as hj->get_time()>t_max.
	// must stop at hj == hi->get_last().

	tstart = hj->get_time();
	if(tstart>snap_time) {
	  break;
	}
	else {
	  tend = hj->get_future()->get_time();
	  if (tend>snap_time)
	    tend = snap_time;
	}

	// make sure the initial binary is counted.
	if(hj == hj->get_first()) 
	  birthrate = contribution_to_population(snap_time, tstart, 
						 tau, imf_cnst);
	// count binaries only once
	else if(*hj != *hj->get_past()) 
	  birthrate = contribution_to_population(snap_time, tstart, 
						 tau, imf_cnst);
	else
	  birthrate = 0;

	if(hj) {
	  tpe_bin     = hj->get_binary_type();
	  prim_type   = hj->get_primary_type();
	  prim_mass   = hj->get_parameter(primary_mass);
	  sec_type    = hj->get_secondary_type();
	  sec_mass    = hj->get_parameter(secondary_mass);

	  if(prim_type>=0 && sec_type>=0) {
	    if(tpe_bin==Detached) {
	      N_detached += birthrate; 
	      detached_population[(int)prim_type][(int)sec_type] 
		  += birthrate;  
	      total_number_of_binaries += birthrate;

	    }
	    else if(tpe_bin==Semi_Detached ||
		    tpe_bin==Contact) {
	      N_contact += birthrate; 
	      contact_population[(int)prim_type][(int)sec_type] += birthrate;  
	      total_number_of_binaries += birthrate;
	    }
	    else {
	      N_single += birthrate; 
	      single_population[(int)prim_type] += birthrate;  
	      total_number_of_binaries += birthrate;
	      
	      if(tpe_bin==Disrupted) {
		N_single += birthrate; 
		single_population[(int)sec_type] += birthrate;  
		total_number_of_binaries += birthrate;
	      }
	    }
	  }
	}
      }

      delete hi;
      
      hi = next_hi;
      next_hi = get_history(hi, cin);
    }
    while (next_hi);

    // renormalize to 100%

    PRL(total_number_of_binaries);
    if(normalize>0) {
	int p = cerr.precision(4);
	cerr << "Normalized to 100%"<<endl;
	total_number_of_binaries /= normalize;

	for (int j = 0; j<no_of_stellar_type-1; j++)  {

	    single_population[j] /= total_number_of_binaries;

	    for (int i = 0; i < no_of_stellar_type-1; i++) {
		detached_population[i][j] /= total_number_of_binaries;
		contact_population[i][j] /= total_number_of_binaries;
	    }
	}
	cerr.precision(p);
    }

    if (N_detached>0) {
	cerr << "     Detached population"<<endl;
	print_binary_matrix(detached_population);
        cerr << "\n     Summary:" << endl;
        print_binary_array(detached_population);
    }
    else
	cerr << "     ---No detached binaries---" <<endl;

    cerr <<  endl;

      if (N_contact>0) {
	cerr << "     Semi-detached/contact population"<<endl;
	print_binary_matrix(contact_population);
        cerr << "\n     Summary:" << endl;
        print_binary_array(contact_population);
      }
      else
	cerr << "     ---No contact Binaries---" <<endl;

    cerr <<  endl;

      cerr << "     Single population"<<endl;
      cerr << "    ";
      for (int j=0; j<no_of_stellar_type-1; j++) {
	cerr << " " << type_short_string((stellar_type)j);
	if(j==(int)Main_Sequence || 
	   j==(int)Super_Giant   || 
	   j==(int)Thorn_Zytkow)
	  cerr << "  "; 
      }

    cerr <<  endl;

      cerr << endl << "     ";
      for (int j = 0; j<no_of_stellar_type-1; j++) {
	cerr << " " << single_population[j];
	if(j==(int)Main_Sequence || 
	   j==(int)Super_Giant   || 
	   j==(int)Thorn_Zytkow)
	  cerr << "  "; 
      }
    cerr << endl;

    return binaries_read;

  }
//----------------------------------------------------------------------------
//  main  --  driver to reduce SeBa short dump data
//----------------------------------------------------------------------------

main(int argc, char ** argv) {

    real normalize = -1;
    real snap_time = 0;
    real tau = 0;
    real cnst = 1;
    bool count_sn = false;

    char  *comment;
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "N:St:x:c:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
            case 't': snap_time = atof(poptarg);
		      break;
            case 'x': tau = atof(poptarg);
		      break;
            case 'N': normalize = atof(poptarg);
		      break;
            case 'S': count_sn = !count_sn;
		      break;
            case 'c': cnst = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            


    SeBa_hist* hi = new SeBa_hist;
    if (!hi->read_SeBa_hist(cin))
      exit(-1);

    cerr << "Time = " << snap_time << " [Myr]" << endl << endl;
    

    int binaries_read = extract_population(hi, snap_time, 
					 tau, cnst, normalize, count_sn);

    cerr << "Total number of binaries read: " << binaries_read << endl;
}

#endif // endof: TOOLBOX



