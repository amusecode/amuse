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

real **mkrarray(int rows, int cols) {
    real ** ret = new real *[rows];
    for (int i=0; i<cols; i++)
	ret[i] = new real[cols];

    for (int i=0; i<rows; i++)
	for (int j=0; j<cols; j++)
	ret[i][j] = 0;

    return ret;
}

void print_binary_matric(real **population) {

    real N_row = 0; 
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
      N_row = 0;
      for (int i = 0; i < no_of_stellar_type-1; i++)
	N_row += population[i][k];
      
      if(N_row>0) {
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

    int n_all = no_of_star_type_summ*no_of_star_type_summ;
    real *pri_sec = new real[n_all];
    for (int j=0; j<n_all; j++)
      pri_sec[j] = 0;

    int js, ks, ts;
    for (int j=0; j<no_of_stellar_type-1; j++)
      for (int k = 0; k < no_of_stellar_type-1; k++) {
        js = (int)summarize_stellar_type((stellar_type)j);
        ks = (int)summarize_stellar_type((stellar_type)k);
        if(js>ks) {
          ts = js;
          js=ks;
          ks = ts;
        }
        pri_sec[js + no_of_star_type_summ*ks] += population[j][k];
      }

    for (js=0; js<no_of_star_type_summ-1; js++)
      for (ks = 0; ks < no_of_star_type_summ-1; ks++) {
        if(pri_sec[js + no_of_star_type_summ*ks]>0)
          cerr << " (" << type_short_string((stellar_type_summary)js)
               << ", " << type_short_string((stellar_type_summary)ks) <<")\t"
               << pri_sec[js + no_of_star_type_summ*ks] << endl;
      }

    delete []pri_sec;
  }



local real extract_snapshot(SeBa_hist *hi, real snap_time,
			   bool normalize) {

    real binaries_read = 0;

    real ** detached_population;
    real** contact_population;
    real* single_population;
    real* blue_straggler_population;
    bool p_bss, s_bss, tb;
    int N_detached=0, N_contact=0, N_single=0;
    binary_type tpe_bin;
    stellar_type prim_type, sec_type, ts;
    real prim_mass, sec_mass;
    real m_to;
      m_to = turn_off_mass(snap_time);

      detached_population = mkrarray((int)no_of_stellar_type, 
				    (int)no_of_stellar_type);
      contact_population = mkrarray((int)no_of_stellar_type, 
				   (int)no_of_stellar_type);
      single_population = new real[(int)no_of_stellar_type];
      blue_straggler_population = new real[(int)no_of_stellar_type];
      for(int i=0; i<(int)no_of_stellar_type; i++) {
	single_population[i]=0;
        blue_straggler_population[i] = 0;
      }


    int mergers = 0;	
    SeBa_hist* next_hi = get_history(hi, cin);
    do {

      binaries_read++;

      if (hi->get_last()->get_time() < snap_time) {
	cerr << "ERROR: snap_time not reached" << endl;
	exit(-1);
      }

      if (hi->get_last()->get_binary_type() == Merged) 
	mergers++;
//	cerr << hi->get_parameter(identity) << endl;

	SeBa_hist *hj = hi->get_SeBa_hist_at_time(snap_time);

//	tpe_bin     = hi->binary_type_at_time(snap_time);
//	prim_type   = hi->primary_type_at_time(snap_time, p_bss);
//	sec_type    = hi->secondary_type_at_time(snap_time, s_bss);
	if(hj) {
	  tpe_bin     = hj->get_binary_type();
	  prim_type   = hj->get_primary_type();
	  prim_mass   = hj->get_parameter(primary_mass);
	  sec_type    = hj->get_secondary_type();
	  sec_mass    = hj->get_parameter(secondary_mass);
	  if(prim_type==Main_Sequence && prim_mass>m_to)
	    p_bss = true;
	  else
	    p_bss = false;
	  if(sec_type==Main_Sequence && sec_mass>m_to)
	    s_bss = true;
	  else
	    s_bss = false;
	  if (prim_mass<sec_mass) {
	    ts = prim_type;
	    prim_type = sec_type;
	    sec_type = ts;
	    tb = p_bss;
	    p_bss=s_bss;
	    s_bss=tb;
	  }

	  if(prim_type>=0 && sec_type>=0) {
	    if(tpe_bin==Detached) {
	      N_detached++;
	      if (p_bss)
		blue_straggler_population[(int)sec_type]++; 
	      else if (s_bss)
		blue_straggler_population[(int)prim_type]++; 
	      else
		detached_population[(int)prim_type][(int)sec_type]++; 
	    }
	    else if(tpe_bin==Semi_Detached ||
		    tpe_bin==Contact) {
	      N_contact++;
	      if (p_bss)
		blue_straggler_population[(int)sec_type]++; 
	      else if (s_bss)
		blue_straggler_population[(int)prim_type]++; 
	      else
		contact_population[(int)prim_type][(int)sec_type]++; 
	    }
	    else {
	      N_single++;
	      single_population[(int)prim_type]++; 
	      if (p_bss)
		blue_straggler_population[Proto_Star]++; 
	      
	      if(tpe_bin==Disrupted) {
		N_single++;
		single_population[(int)sec_type]++; 
		if (p_bss)
		  blue_straggler_population[Disintegrated]++; 
	      }
	    }
	  }
	}


      delete hi;
      
      hi = next_hi;
      next_hi = get_history(hi, cin);
    }
    while (next_hi);

    PRL(binaries_read);
    if(normalize) {
	int p = cerr.precision(LOW_PRECISION);
	cerr << "Normalized to 100%"<<endl;
	binaries_read /= 100;

	for (int j = 0; j<no_of_stellar_type-1; j++)  {

	    single_population[j] /= binaries_read;
	    blue_straggler_population[j] /= binaries_read;

	    for (int i = 0; i < no_of_stellar_type-1; i++) {
		detached_population[i][j] /= binaries_read;
		contact_population[i][j] /= binaries_read;
	    }
	}
	cerr.precision(p);
    }

    if (N_detached>0) {
	cerr << "     Detached population"<<endl;
	print_binary_matric(detached_population);
        cerr << "\n     Summary:" << endl;
        print_binary_array(detached_population);
    }
    else
	cerr << "     ---No detached binaries---" <<endl;

    cerr <<  endl;

    if (N_contact>0) {
	cerr << "     Semi-detached/contact population"<<endl;
	print_binary_matric(contact_population);
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
    cerr << "\n     Blue straggler population (pl=merged, di=dissrupted)"
	 << endl;
    cerr << endl << "     ";
    for (int j = 0; j<no_of_stellar_type-1; j++) {
	cerr << " " << blue_straggler_population[j];
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

    cerr << "Time = " << snap_time << " [Myr]" << endl << endl;
    

    real binaries_read = extract_snapshot(hi, snap_time, normalize);

    cerr << "Total number of binaries read: " << binaries_read << endl;
}

#endif // endof: TOOLBOX



