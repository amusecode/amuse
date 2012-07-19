//// rs_findtype: Reads SeBa.data standard output format and finds
////              specified binary types with specified components.
////          
//// Options:   -A    maximum separation [inf]
////            -a    minimum separation [0]
////            -M    maximum primary mass [inf]
////            -m    minimum primary mass [0]
////            -E    maximum eccentricity [1]
////            -e    minimum eccentricity [0]
////            -Q    maximum mass ratio (m_sec/m_prim) [inf]
////                  may be larger than unity!
////            -q    minimum mass ratio (m_sec/m_prim) [0]
////            -f    output only at first occasion [true]
////            -P/p  primary type (string or summary or 'any') [bla]
////                  options: static_star, SPZDCH_star, nas, 
////		               proto_star, planet, brown_dwarf,
////			       main_sequence, hyper_giant, hertzsprung_gap,
////			       sub_giant, horizontal_branch, super_giant,
////			       carbon_star, helium_star, helium_giant, 
////		               carbon_dwarf, helium_dwarf, oxygen_dwarf,
////			       thorn_zytkow,
////			       xray_pulsar, radio_pulsar, neutron_star, 
////                           black_hole, disintegrated, 
////                           double, no_of_stellar_type
////                 or: ZAMS, early_giant, late_giant,
////			   helium_remnant, white_dwarf, 
////			   neutron_remnant, inert_remnant,
////			   unspecified, undefined
////                 or: SS, PZHs, nas, ps, pl, bd, ms, wr, hg, gs, hb, 
////                     sg, co, he, gh, hd, cd, od, TZO, xp, rp, ns, bh
////                     di, bin
////            -S/s  secondary type (see primary type) [bla]
////            -B    Binary type [detached]
////                  options: detached, semi_detached, 
////                           contact, merged, disrupted
////            -C    Mass transfer type [Unknown]
////                  options: Unknown, Nuclear, 
////                           AML_Driven, Thermal, Dynamic
////            -t    time interval between formation of primary (-P/p)
////                  and secondary (-S/s) [10000]
////            -R    output entire scenario in standard format [false]
////            -v    verbose [false]
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


local void extract_binary_on_stellar_types(SeBa_hist *ha,
					   char* prim_string,
					   char* sec_string,
					   binary_type bt,
					   mass_transfer_type mt,
					   real a_min, real a_max,
					   real m_min, real m_max,
					   real e_min, real e_max,
					   real q_min, real q_max,
					   real dt,
					   bool first_occasion,
					   bool R_flag) {

  if (!strcmp(prim_string, sec_string)) {
    for_all_SeBa_hist(SeBa_hist, ha, hi) {
      if (hi->binary_contains(prim_string, sec_string, bt, mt) &&
	   hi->binary_limits(semi_major_axis, a_min, a_max) &&
	   hi->binary_limits(primary_mass, m_min, m_max) &&
	   hi->binary_limits(mass_ratio, q_min, q_max) &&
	   hi->binary_limits(eccentricity, e_min, e_max)) {
	 
	 if (R_flag) {
	   for_all_SeBa_hist(SeBa_hist, ha, ho) cout << *ho;	   
	 } 
	 else if (!strcmp(prim_string, "any")) {
	   
	   cout << *hi->get_first();
	   cout << *hi;
	 } else {

	   hi->put_first_formed_left(prim_string, dt);
	 }

	 if (first_occasion) 
	   return;
       }
     }
  }
  else {
    for_all_SeBa_hist(SeBa_hist, ha, hi) {
      if ((hi->binary_limits(semi_major_axis, a_min, a_max) &&
	   hi->binary_limits(eccentricity, e_min, e_max))) {
	
	if ((hi->binary_contains(prim_string, sec_string, bt, mt) &&
	     hi->binary_limits(mass_ratio, q_min, q_max) &&
	     hi->binary_limits(primary_mass, m_min, m_max))) {
	  

	  if (R_flag) {
	    for_all_SeBa_hist(SeBa_hist, ha, ho) cout << *ho;	   
	  } 
	  else {
	    cout << *hi->get_first();
	    cout << *hi;
	  }

	  if (first_occasion) 
	    return;

	} else if ((hi->binary_contains(sec_string, prim_string, bt, mt) &&
		    hi->binary_limits(mass_ratio, 1/q_max, 1/q_min) &&
		    hi->binary_limits(secondary_mass, m_min, m_max))) {
	  
	  if (R_flag) {
	    for_all_SeBa_hist(SeBa_hist, ha, ho) cout << *ho;	   
	  } 
	  else {
	    hi->get_first()->put_single_reverse(cout);
	    hi->put_single_reverse(cout);
	  }

	  if (first_occasion) 
	    return;

	}	
      }
    }
  }

}

//-----------------------------------------------------------------------------
//  main  --  driver to reduce SeBa short dump data
//-----------------------------------------------------------------------------


main(int argc, char ** argv) {

    bool  c_flag = false;
    bool  v_flag = false;
    bool  R_flag = false;
    bool  first_occasion = true;
    
    char * primary_type_string = "bla";
    char * secondary_type_string = "bla";

    binary_type bt = Detached;
    mass_transfer_type mt = Unknown;

    real a_max = VERY_LARGE_NUMBER;
    real a_min = 0;
    real m_max = VERY_LARGE_NUMBER;
    real m_min = 0;
    real e_max = 1;
    real e_min = 0;
    real q_max = VERY_LARGE_NUMBER;
    real q_min = 0;
    real dt = 1.e4;

    real snap_time = -1;
    bool T_flag = false;

    int binaries_read = 0;

    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "P:p:S:s:B:C:A:a:M:m:E:e:Q:q:ft:vR";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'A': a_max = atof(poptarg);
		      break;
	    case 'a': a_min = atof(poptarg);
		      break;
	    case 'M': m_max = atof(poptarg);
		      break;
	    case 'm': m_min = atof(poptarg);
		      break;
	    case 'E': e_max = atof(poptarg);
		      break;
	    case 'e': e_min = atof(poptarg);
		      break;
	    case 'Q': q_max = atof(poptarg);
		      break;
	    case 'q': q_min = atof(poptarg);
		      break;
	    case 'f': first_occasion = !first_occasion;
		      break;
            case 'p':
            case 'P': primary_type_string = poptarg;
		      break;
            case 's':
            case 'S': secondary_type_string = poptarg;
		      break;
            case 'B': bt = extract_binary_type_string(poptarg);
		      break;
            case 'C': mt = extract_mass_transfer_type_string(poptarg);
		      break;
            case 't': dt = atof(poptarg);
		      break;
            case 'R': R_flag = true;
		      break;
	    case 'v': v_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            

    SeBa_hist* hi = new SeBa_hist;
    if (!hi->read_SeBa_hist(cin))
      exit(-1);


    // (GN+SPZ Dec 13 1999) disrupted has e = 1
    if (bt == Disrupted) e_max = 2.;

    int mergers = 0;	
    SeBa_hist* next_hi = get_history(hi, cin);
    do {

      binaries_read++;

      if (hi->get_last()->get_binary_type() == Merged) 
	mergers++;

	extract_binary_on_stellar_types(hi, primary_type_string,
					secondary_type_string,
					bt, mt,
					a_min, a_max, m_min, m_max,
					e_min, e_max, q_min, q_max,
					dt, first_occasion,
					R_flag);

      if (v_flag)
	hi->put_history(cout, v_flag);
      

      delete hi;
      
      hi = next_hi;
      next_hi = get_history(hi, cin);
    }
    while (next_hi);

     extract_binary_on_stellar_types(hi->get_last(), primary_type_string,
				     secondary_type_string,
				     bt, mt,
				     a_min, a_max, m_min, m_max,
				     e_min, e_max, q_min, q_max, 
				     dt, first_occasion,
				     R_flag);

     if (v_flag) {
       hi->put_history(cout, v_flag);
       cerr <<"Number of Mergers: " << mergers << endl;
     }

    cerr << "Total number of binaries read: " << binaries_read << endl;
}

#endif // endof: TOOLBOX





