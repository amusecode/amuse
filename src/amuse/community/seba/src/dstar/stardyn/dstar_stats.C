
//  dstar_stats.C: Print out various diagnostic binary statistics on a system.

#include "dstar_to_dyn.h"

#include "hdyn.h"		// hdyn.h needed here only because of the
				// hdyn::set_use_dstar() reference in main()

#ifndef TOOLBOX

struct pop {
  int stypes[no_of_stellar_type];

  pop() {
    for (int i=0; i<no_of_stellar_type; i++)
      stypes[i]=0;
  }
};


struct nm_bin {
          int   no_of_stars;
          real  mass_limit;
          real  total_mass;
  nm_bin() {
    no_of_stars=0;
    mass_limit=total_mass=0;
  }
};

local void print_binary_content(dyn* b, binary_type which) {

  pop *content = new pop[no_of_stellar_type-1];

    int N_total = 0;
    double_state bst;
    for_all_daughters(dyn, b, bi) 

      // Count single stars and binaries separately.

      if (bi->get_oldest_daughter() != NULL &&
	  has_dstar(bi->get_oldest_daughter())) {
	
	bst = make_state(bi);
	if(bst.type==which) {
	  N_total++;	
	  content[bst.primary.type].stypes[bst.secondary.type]++;
	}
      }

    int N_row = 0; 
    if(N_total>0) {

      //      cerr << "    ";
      //      for (int j=0; j<no_of_stellar_type-1; j++)
      //	cerr << " " << type_short_string((stellar_type)j);
      //      cerr << endl;

      cerr << "    ";
      for (int j = 0; j <= dynamic_cast(int, Super_Giant); j++) 
	cerr << " " << type_short_string(dynamic_cast(stellar_type, j));

      cerr << "     ";
      for (int j = dynamic_cast(int, Carbon_Star); 
	   j < dynamic_cast(int, no_of_stellar_type-1); j++) 
	cerr << " " << type_short_string(dynamic_cast(stellar_type, j));
      cerr << endl;
      
      for (int k = 0; k < no_of_stellar_type-1; k++) {
	N_row = 0;
	for (int i = 0; i < no_of_stellar_type-1; i++)
	  N_row += content[i].stypes[k];

	if(N_row>0) {
	  cerr << "   " << type_short_string((stellar_type)k);

	  //	  for (int i = 0; i < no_of_stellar_type-1; i++) 
	  //	    cerr << " " << content[i].stypes[k];
	  //	  cerr << endl;
	  for (int i=0; i<=dynamic_cast(int, Super_Giant); i++)
	    cerr << " " << content[i].stypes[k];
	  cerr << "     ";
	  for (int i=Carbon_Star; i<no_of_stellar_type-1; i++)
	    cerr << " " << content[i].stypes[k];
	  cerr << endl;
	}
      }

    }

    else
      cerr << "     ---No "<< type_string(which) << " binaries---" <<endl;
    
    delete []content;
}

void dstar_stats(dyn* b, bool mass_spectrum, vec center,
		 bool verbose) { 

  sstar_stats(b, mass_spectrum, center, verbose); 
  
  if (verbose)
    cerr << "\n  Detached binary content:       (Primary)\n";
  print_binary_content(b, Detached);

  if (verbose)
    cerr << "\n  Synchronized binary content:       (Primary)\n";
  print_binary_content(b, Synchronized);

  if (verbose)
    cerr << "\n  Semi-detached binary content:  (Primary)\n";
  print_binary_content(b, Semi_Detached);

  if (verbose)
    cerr << "\n  Contact binary content:        (Primary)\n";
  print_binary_content(b, Contact);
      
}

#else

main(int argc, char **argv)
{
    bool binaries = true, verbose = true, out = false,
         n_sq = true, calc_e = true;
    int which_lagr = 0;

    extern char *poptarg;
    int c;
    const char *param_string = "bneost";	// Note: "v" removed because only the
					// "verbose" option currently works.

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.6 $", _SRC_)) != -1)
	switch(c) {

	    case 'b': binaries = !binaries;
		      break;
	    case 'e': calc_e = !calc_e;
		      break;
	    case 'n': n_sq = !n_sq;
		      break;
	    case 'o': out = !out;
		      break;
	    case 's': which_lagr = 2;
		      break;
	    case 't': which_lagr = 1;
		      break;
	    case 'v': verbose = !verbose;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    // Loop over input until no more data remain.

    hdyn *b;
    int i = 0;
    vec zero = 0;
    bool mass_spectrum = true;
    while (b = get_hdyn()) {

      if(find_qmatch(b->get_oldest_daughter()
		 ->get_starbase()->get_star_story(), "Type")) {
	addstar(b,                             // Note that T_start and
		b->get_system_time(),          // Main_Sequence are
		Main_Sequence,                 // defaults. They are
		true);                         // ignored if a star

	b->set_use_sstar(true);

	adddouble(b, b->get_system_time());

	b->set_use_dstar(true);
      }
      
      if (i++ > 0) cerr << endl;
      dstar_stats(b, mass_spectrum, zero,  verbose);

      if (out) put_hdyn(b);
      rmtree(b);
    }
}

#endif
