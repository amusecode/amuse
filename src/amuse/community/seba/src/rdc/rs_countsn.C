//// rs_snapshot: Reads SeBa.data standard output format and creates
////              snapshot at time T
////          
//// Options:   -t     snapshot time in Myr [0]
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

#ifdef TOOLBOX


// Note that only one star can go SN at a time.
supernova_type count_supernovae(SeBa_hist *hi) {

  supernova_type sn_type = NAT;

  if (post_supernova_star(hi->get_primary_type())) 
    sn_type = type_of_supernova(hi->get_past()->get_primary_type());
    
  if(sn_type == NAT) {
    if (post_supernova_star(hi->get_secondary_type())) 
      sn_type = type_of_supernova(hi->get_past()->get_secondary_type());

    return sn_type;
  }
}

local int count_supernova_types(SeBa_hist *hi, bool normalize) {

    int binaries_read = 0;
    supernova_type sn_type;
    real n_sn = 0;

    real supernovae[no_of_supernova_type];
    for(int i=0; i<(int)no_of_supernova_type; i++) {
      supernovae[i] = 0;
    }

    int mergers = 0;	
    SeBa_hist* next_hi = get_history(hi, cin);
    do {

      binaries_read++;
      for_all_SeBa_hist(SeBa_hist, hi, ha) {

	sn_type = count_supernovae(ha);
	if(sn_type!=NAT) {
	  supernovae[(int)sn_type]++;
	  n_sn++;
	}
      }

      delete hi;
      hi = next_hi;
      next_hi = get_history(hi, cin);
    }
    while (next_hi);

    if (n_sn<=0) {
      cerr << "     ---No supernovae---" <<endl;
      exit(-1);
    }

    PRL(n_sn);
    if(normalize) {
	int p = cerr.precision(LOW_PRECISION);
	cerr << "Normalized to 100%"<<endl;
	n_sn /= 100;
      
	for (int i = 0; i<no_of_supernova_type-1; i++)  {
	    supernovae[i] /= n_sn;
	}
	cerr.precision(p);
    }

    cerr << "     Supernova population " << n_sn << endl;
    for(int i=NAT; i<no_of_supernova_type; i++) 
      cerr << type_string((supernova_type)i) << "\t";
    cerr << endl;

    for(int i=NAT; i<no_of_supernova_type; i++) 
      cerr <<  supernovae[i] << "\t";
    cerr << endl;

    return binaries_read;
  }
//----------------------------------------------------------------------------
//  main  --  driver to reduce SeBa short dump data
//----------------------------------------------------------------------------

main(int argc, char ** argv) {

  bool normalize = false;
    real snap_time = 0;

    char  *comment;
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "Nt:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
            case 't': snap_time = atof(poptarg);
		      break;
            case 'N': normalize = !normalize;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            


    SeBa_hist* hi = new SeBa_hist;
    if (!hi->read_SeBa_hist(cin))
      exit(-1);

    int binaries_read = count_supernova_types(hi, normalize);

    cerr << "Total number of binaries read: " << binaries_read << endl;
}

#endif // endof: TOOLBOX



