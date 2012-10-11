
//// rdc_SeBa: read SeBa short dump
////           output requested binary parameters.
////          
//// Options:     none
//-----------------------------------------------------------------------------
//   version 1:  Sept 1998   Simon Portegies Zwart   spz@grape.c.u-tokyo.ac.jp
//                                                   University of Tokyo
//.............................................................................
//   non-local functions: 
//-----------------------------------------------------------------------------

#include "stdinc.h"
#include "node.h"
#include "double_star.h"
#include "main_sequence.h"
//
// SeBa_hist.h
//

#ifndef    _SeBa_HIST
#  define  _SeBa_HIST


enum binary_parameter {identity=0, bin_type, current_time,
		       primary_mass, primary_radius,
		       secondary_mass, secondary_radius,
		       semi_major_axis, eccentricity,
		       mass_ratio
                      };

/*-----------------------------------------------------------------------------
 * SeBa_hist  --  a linked list of SeBa histories.
 *-----------------------------------------------------------------------------
 */
class SeBa_hist {
    protected:

    int  number;
    real time;

    binary_type bin_tpe;
    mass_transfer_type mtrans_tpe;
    real semi;
    real ecc;

    //    char label_prim[3];
    int id_prim;
    stellar_type tpe_prim;
    real m_prim;
    real r_prim;
    real temp_prim;
    real mcore_prim;

    //char label_sec[3];
    int id_sec;
    stellar_type tpe_sec;
    real m_sec;
    real r_sec;
    real temp_sec;
    real mcore_sec;

    SeBa_hist * past;
    SeBa_hist * future;

    public:
       SeBa_hist(SeBa_hist* s=NULL) {
	   if (s) {
  	      past=s;
	      past->future=this;
	      future = NULL;
 	   }
	   else {
	      past=NULL;
	      future=NULL;
	  }

	   number = 0;
	   time   = 0;
	   bin_tpe = Detached;
	   mtrans_tpe = Unknown;
	   // label_prim = "1a";
	   //strcpy(label_prim, "1a");			// Steve, 20010907
	   // label_sec = "1b";
	   //strcpy(label_sec, "1b");
	   id_prim = 0;
	   id_sec = 1;
	   tpe_prim = tpe_sec = Main_Sequence;
	   m_prim=m_sec=r_prim=r_sec=temp_prim=temp_sec=mcore_prim=mcore_sec=0;
       }
       ~SeBa_hist(){
         if (future!=NULL) {
	     SeBa_hist *tmp = future;
	     future = NULL;
	     delete tmp;
	 }

	 if (past)
	    past->future = NULL;
       }

       SeBa_hist* get_past() {return past;}
       void set_past(SeBa_hist *sb) {past = sb;}
       SeBa_hist* get_future() {return future;}
       SeBa_hist* get_first() {
           if (past!=NULL)
              return past->get_first();
           else
              return this;
       }
       SeBa_hist* get_last() {
           if (future!=NULL)
              return future->get_last();
           else 
              return this;
       }

       void set_future(SeBa_hist* f) {future = f;}
       void set_last(SeBa_hist* f) {
            if(future!=NULL) 
              future->set_last(f);
            else
              future=f;
       }
       int get_number()             {return number;}
       void set_number(int n)        {number=n;} // (GN Nov 16 2001)

       real get_time()              {return time;}

       real set_stellar_radius(bool);
       void move_SeBa_hist_to(SeBa_hist*);
       bool read_SeBa_hist(istream&);

       void put_history(ostream&, bool);
       void put_single_reverse(ostream&);
       void put_first_formed_left(char*, real);

       bool binary_contains(char*, char *, binary_type);
       bool binary_limits(binary_parameter, real, real);
       real get_parameter(binary_parameter);
       binary_type get_binary_type() { return bin_tpe;}
       mass_transfer_type get_mass_transfer_type() { return mtrans_tpe;}
       stellar_type get_primary_type() { return tpe_prim;}
       stellar_type get_secondary_type() { return tpe_sec;}

       SeBa_hist* get_SeBa_hist_at_time(real);
				  
       friend ostream& operator<<(ostream& s, SeBa_hist&);
//       friend istream& operator>>(istream& s, cluster_table& table);
};

#define for_all_SeBa_hist(SeBa_hist, base, SeBa_hist_next)                    \
        for (SeBa_hist* SeBa_hist_next = base;                                \
	     SeBa_hist_next != NULL;                            \
	     SeBa_hist_next = SeBa_hist_next->get_future())

ostream& operator<<(ostream& s, SeBa_hist& hi) {
    // short dump output is without bin_type
    //s << "\t" << (int)hi.bin_tpe; 
    // now with bin_type and mt_type

    s << hi.number;
    s << "  " << (int)hi.bin_tpe; 
    s << "\t" << (int)hi.mtrans_tpe; 	      
    s << "\t" << hi.time; 
    s << "\t" << hi.semi << " " << hi.ecc; 
    //   s << "\t" << hi.id_prim << " " << hi.tpe_prim 
    s << "\t" << hi.tpe_prim 
      << " " << hi.m_prim << " " << hi.r_prim 
      << " " << hi.temp_prim << " " << hi.mcore_prim;
   //    s << "\t" << hi.id_sec << " " << hi.tpe_sec 
    s << "\t" << hi.tpe_sec 
      << " " << hi.m_sec << " " << hi.r_sec
      << " " << hi.temp_sec << " " << hi.mcore_sec;
    s << endl;

    return s;

}


void SeBa_hist::put_history(ostream& s, bool verbose=false) {

     if (verbose) {
	 s << "\nnumber= " << number <<"\tTime= "<< time;
	 s << "\t" <<  type_string(bin_tpe) 
	   << "\t a= " << semi << " e= " << ecc << endl;
	 s << "\t"  << type_string(tpe_prim) //<< "(" << label_prim << ")"
	   << "\tM= " << m_prim << " R= " << r_prim << endl;
	 s << "\t" << type_string(tpe_sec) //<< "(" << label_sec << ")"
	   << "\tm= " << m_sec << " r= " << r_sec << endl;
     }
     else {
	 s << *this;
     }

     if (future!=NULL) 
        future->put_history(s, verbose);
 }


void SeBa_hist::put_single_reverse(ostream& s) {

     s << number;
     //s << "\t" <<  type_string(bin_tpe);
     s << "  " << bin_tpe; 
     s << "\t" << mtrans_tpe;
     s << "\t" << time; 	       
 	       
     s << "\t" << semi << " " << ecc; 
     //     s << "\t" << id_sec << " " << tpe_sec << " " 
     s << "\t" << tpe_sec << " " 
       << m_sec << " " << r_sec << " " 
       << temp_sec << " " << mcore_sec;
     //     s << "\t" << id_prim << " " <<tpe_prim << " " 
     s << "\t" <<tpe_prim << " " 
       << m_prim << " " << r_prim << " "
       << temp_prim << " " << mcore_prim;
     s << endl;
 }

real SeBa_hist::get_parameter(binary_parameter param) {

     switch(param) {
     case identity:                          return number;
          break;
     case bin_type:                          return bin_tpe;
          break;
     case current_time:                      return time;
          break;
     case primary_mass:                      return m_prim;
          break;
     case primary_radius:                    return r_prim;
          break;
     case secondary_mass:                    return m_sec;
          break;
     case secondary_radius:                  return r_sec;
          break;
     case semi_major_axis:                   return semi;
          break;
     case eccentricity:                      return ecc;
          break;
     case mass_ratio:                        return m_sec/m_prim;
          break;
     }
}


real SeBa_hist::set_stellar_radius(bool primary = true) {

    real r = 1.e6;
    real q = m_prim/m_sec;
    real rl_r = r_prim;
    if (!primary) {
       rl_r = r_sec;
       q = 1/q;
    }
    
    if (rl_r >= 0) { // detached
      real q1_3 = pow(q, cnsts.mathematics(one_third));
      real q2_3 = pow(q1_3, 2);   
      
      r = semi*0.49*q2_3/(0.6*q2_3 + log(1 + q1_3)) - rl_r;
    }
    
    return r;
}

bool SeBa_hist::read_SeBa_hist(istream& s) {
    
    int tpe_1, tpe_2, btpe, mtpe;

    // (GN May 27 2004) adjust for new dump output

    s >> number >> btpe >> mtpe >> time 
      >> semi >> ecc 
//      >> tpe_1 >> m_prim >> r_prim
//      >> tpe_2 >> m_sec >> r_sec;
      >> id_prim >> tpe_1 >> m_prim >> r_prim >> temp_prim >> mcore_prim
      >> id_sec >> tpe_2 >> m_sec >> r_sec >> temp_sec >> mcore_sec;

    if(s.eof()) return false;


    tpe_prim = (stellar_type)tpe_1;
    tpe_sec = (stellar_type)tpe_2;

    bin_tpe = (binary_type)btpe;
    mtrans_tpe = (mass_transfer_type)mtpe;


#if 0 // (GN May 27 2004) obsolete
    // The safety factor semi<100000 for semi-detached and contact binaries
    // is a precision problem. If the orbital separation is very large and
    // the size fo the star is small, 
    // the size of the Roche-lobe - stellar radius =~ Roche-lobe.
    // This results in semi-detached binaries at very large separations. 
    if (ecc>=1) {
      bin_tpe = Disrupted;
    }
    else if (semi<=0) {
      bin_tpe = Merged;
    }
    else if (r_prim<=0 && r_sec<=0 && semi<100000) 
      bin_tpe = Contact;
    else if(r_prim<=0 || r_sec<=0 && semi<100000)  
      bin_tpe = Semi_Detached;
    
    if (bin_tpe != Merged && bin_tpe != Disrupted) {
      r_prim = set_stellar_radius(true);
      r_sec  = set_stellar_radius(false);
    } else {
      r_prim = -r_prim;
      r_sec  = -r_sec;
    }
#endif

    return true;

}

void SeBa_hist::move_SeBa_hist_to(SeBa_hist *next_hi) {
    
    next_hi = this;
    next_hi->future = NULL;
    next_hi->past   = NULL;
}

bool SeBa_hist::binary_contains(char *prim_string,
				char *sec_string,
				binary_type bt = Detached) {

  if (bin_tpe== bt && 
      (!strcmp(prim_string, type_string(tpe_prim))      ||
       !strcmp(prim_string, type_short_string(tpe_prim)) ||
       !strcmp(prim_string, type_string(summarize_stellar_type(tpe_prim))) ||
       !strcmp(prim_string, "any"))
                           && 
      (!strcmp(sec_string, type_string(tpe_sec))      ||
       !strcmp(sec_string, type_short_string(tpe_sec)) ||
       !strcmp(sec_string, type_string(summarize_stellar_type(tpe_sec))) ||
       !strcmp(sec_string, "any"))) {

    return true;
  }

#if 0
  PRC(number);
  PRC(prim_string);PRL(type_string(tpe_prim));
  PRC(!strcmp(prim_string, type_string(tpe_prim)));
  PRC(!strcmp(prim_string, type_short_string(tpe_prim)));
  PRL(!strcmp(prim_string, type_string(summarize_stellar_type(tpe_prim))));

  PRC(sec_string);PRL(type_string(tpe_sec));
  PRC(!strcmp(sec_string, type_string(tpe_sec)));
  PRC(!strcmp(sec_string, type_short_string(tpe_sec)));
  PRL(!strcmp(sec_string, type_string(summarize_stellar_type(tpe_sec))));
#endif
  
  return false;
       
}

void SeBa_hist::put_first_formed_left(char *string_type, real dt) {

// (GN Apr  6 1999)
//  for_all_SeBa_hist(SeBa_hist, this, hi) {
  for_all_SeBa_hist(SeBa_hist, get_first(), hi) {
//    cerr << get_time() << " " << hi->get_time() << endl;
//    cerr << string_type << " " 
//	 << type_string(summarize_stellar_type(hi->tpe_prim)) << endl;

    if (!strcmp(string_type, type_string(hi->tpe_prim))      ||
	!strcmp(string_type, type_short_string(hi->tpe_prim)) ||
	!strcmp(string_type,
	       type_string(summarize_stellar_type(hi->tpe_prim))))  {

      if ((get_time() - hi->get_time()) < dt) {
	cout << *get_first();
	// (GN Sep 27 1999) output also intermediate state
	cout << *hi;
	cout << *this;
      }
      break;
    }
    else if (!strcmp(string_type, type_string(hi->tpe_sec))      ||
	     !strcmp(string_type, type_short_string(hi->tpe_sec)) ||
	     !strcmp(string_type,
		    type_string(summarize_stellar_type(hi->tpe_sec))) &&
	((get_time() - hi->get_time()) < dt)) {

      if ((get_time() - hi->get_time()) < dt) {
	get_first()->put_single_reverse(cout);
	// (GN Sep 27 1999) output also intermediate state
	hi->put_single_reverse(cout);
	put_single_reverse(cout);
      }
      break;
    }
  }

}

bool SeBa_hist::binary_limits(binary_parameter param,
			      real lower, real upper) {


     if (get_parameter(param) >= lower &&
	 get_parameter(param) < upper)
       return true;

     return false;

}

SeBa_hist* SeBa_hist::get_SeBa_hist_at_time(real snap_time) {

  if (snap_time>get_last()->get_time())
    return NULL;

  for_all_SeBa_hist(SeBa_hist, get_future(), hi) {
    if (hi->get_past()->get_time()<= snap_time && 
	     snap_time<hi->get_time()) {
      return hi->get_past();
    }
  }

  err_exit("SeBa_hist::SeBa_hist_at_time(real time): time not found");
}

#endif //   _SeBa_HIST

#ifndef TOOLBOX
#else

local SeBa_hist* get_history(SeBa_hist *hi, istream& is) {

    while (hi->get_number()==hi->get_last()->get_number()) {

	if (is.eof())
	   return NULL;

	new SeBa_hist(hi->get_last());
	
	hi->get_last()->read_SeBa_hist(is);

    };

    SeBa_hist *next_hi;
    
    next_hi = hi->get_last();
    next_hi->get_past()->set_future(NULL);
    next_hi->set_past(NULL);

    return next_hi;

	   
}

local void extract_binary_on_stellar_types(SeBa_hist *ha,
					   char* prim_string,
					   char* sec_string,
					   binary_type bt,
					   real a_min, real a_max,
					   real mp_min, real mp_max,
					   real ms_min, real ms_max,
					   real e_min, real e_max,
					   real q_min, real q_max,
					   real dt,
					   bool first_occasion,
					   bool R_flag) {

  if (!strcmp(prim_string, sec_string)) {
    for_all_SeBa_hist(SeBa_hist, ha, hi) {
       if (hi->binary_contains(prim_string, sec_string, bt) &&
	   hi->binary_limits(semi_major_axis, a_min, a_max) &&
	   hi->binary_limits(primary_mass, mp_min, mp_max) &&
	   hi->binary_limits(secondary_mass, ms_min, ms_max) &&
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
	
	if ((hi->binary_contains(prim_string, sec_string, bt) &&
	     hi->binary_limits(mass_ratio, q_min, q_max) &&
	     hi->binary_limits(primary_mass, mp_min, mp_max) &&
	     hi->binary_limits(secondary_mass, ms_min, ms_max))) {
	  

	  if (R_flag) {
	    for_all_SeBa_hist(SeBa_hist, ha, ho) cout << *ho;	   
	  } 
	  else {
	    cout << *hi->get_first();
	    // (GN Oct  3 2000) temporary fot hewd
	    cout << *hi->get_past();
	    cout << *hi;
	  }

	  if (first_occasion) 
	    return;

	} else if ((hi->binary_contains(sec_string, prim_string, bt) &&
		    hi->binary_limits(mass_ratio, 1/q_max, 1/q_min) &&
		    hi->binary_limits(secondary_mass, mp_min, mp_max) &&
		    hi->binary_limits(primary_mass, ms_min, ms_max) )) {
	  // (GN Nov 14 2001) Beware: primary is now secondary and vice versa
	  
	  if (R_flag) {
	    for_all_SeBa_hist(SeBa_hist, ha, ho) cout << *ho;	   
	  } 
	  else {
	    hi->get_first()->put_single_reverse(cout);
	    // (GN Oct  3 2000) temporary fot hewd
	    hi->get_past()->put_single_reverse(cout);
	    hi->put_single_reverse(cout);
	  }

	  if (first_occasion) 
	    return;

	}	
      }
    }
  }
}


int **mkarray(int rows, int cols) {
    int ** ret = new int *[rows];
    for (int i=0; i<cols; i++)
	ret[i] = new int[cols];

    for (int i=0; i<rows; i++)
	for (int j=0; j<cols; j++)
	ret[i][j] = 0;

    return ret;
}

void print_binary_matric(int **population) {

    int N_row = 0; 
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

//#else

//-----------------------------------------------------------------------------
//  main  --  driver to reduce SeBa short dump data
//-----------------------------------------------------------------------------

int main(int argc, char ** argv) {

    bool  c_flag = false;
    bool  v_flag = false;
    bool  R_flag = false;
    bool  first_occasion = false;
    
    char * primary_type_string = "bla";
    char * secondary_type_string = "bla";

    binary_type bt = Detached;

    real a_max = VERY_LARGE_NUMBER;
    real a_min = 0;
    real mp_max = VERY_LARGE_NUMBER;
    real mp_min = 0;
    real ms_max = VERY_LARGE_NUMBER;
    real ms_min = 0;
    real e_max = 1;
    real e_min = 0;
    real q_max = VERY_LARGE_NUMBER;
    real q_min = 0;
    real dt = 1.37e4;

    real snap_time = -1;
    bool T_flag = false;

    int binaries_read = 0;

    char  *comment;
    check_help();

    extern char *poptarg;
    int c;
    char* param_string = "P:p:S:s:B:A:a:M:m:N:n:E:e:Q:q:fc:t:vRT:";

    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
	    case 'A': a_max = atof(poptarg);
		      break;
	    case 'a': a_min = atof(poptarg);
		      break;
	    case 'M': mp_max = atof(poptarg);
		      break;
	    case 'm': mp_min = atof(poptarg);
		      break;
	    case 'N': ms_max = atof(poptarg);
		      break;
	    case 'n': ms_min = atof(poptarg);
		      break;
	    case 'E': e_max = atof(poptarg);
		      break;
	    case 'e': e_min = atof(poptarg);
		      break;
	    case 'Q': q_max = atof(poptarg);
		      break;
	    case 'q': q_min = atof(poptarg);
		      break;
	    case 'f': first_occasion = true;
		      break;
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
            case 'p':
            case 'P': primary_type_string = poptarg;
		      break;
            case 's':
            case 'S': secondary_type_string = poptarg;
		      break;
            case 'B': bt = extract_binary_type_string(poptarg);
		      break;
            case 't': dt = atof(poptarg);
		      break;
            case 'T': T_flag = true;
		      snap_time = atof(poptarg);
		      break;
            case 'R': R_flag = true;
		      break;
	    case 'v': v_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      get_help();
		      exit(1);
	    }            

#if 0
    stellar_type ptype = extract_stellar_type_string(primary_type_string);
    if (P_flag) 
      goto P;
    else
      goto p;

P:  stellar_type prim_type = ptype;
    goto do_S;
p:  stellar_type_summary prim_type = summarize_stellar_type(ptype);
    goto do_S;

do_S: stellar_type stype = extract_stellar_type_string(secondary_type_string);
    if (S_flag) 
      goto S;
    else
      goto s;

S:  stellar_type sec_type = stype;
    goto start;
s:  stellar_type_summary sec_type= summarize_stellar_type(stype);
    goto start;

start:    

#endif    

    SeBa_hist* hi = new SeBa_hist;
    if (!hi->read_SeBa_hist(cin))
	   exit(-1);

    int** detached_population;
    int** contact_population;
    int* single_population;
    int* blue_straggler_population;
    bool p_bss, s_bss, tb;
    int N_detached=0, N_contact=0, N_single=0;
    binary_type tpe_bin;
    stellar_type prim_type, sec_type, ts;
    real prim_mass, sec_mass;
    real m_to;
    if(T_flag) {
      m_to = turn_off_mass(snap_time);

      detached_population = mkarray((int)no_of_stellar_type, 
				    (int)no_of_stellar_type);
      contact_population = mkarray((int)no_of_stellar_type, 
				   (int)no_of_stellar_type);
      single_population = new int[(int)no_of_stellar_type];
      blue_straggler_population = new int[(int)no_of_stellar_type];
      for(int i=0; i<(int)no_of_stellar_type; i++) {
	single_population[i]=0;
        blue_straggler_population[i] = 0;
      }

    }

    // (GN+SPZ Dec 13 1999) disrupted has e = 1
    if (bt == Disrupted) e_max = 2.;

    int mergers = 0;	
    SeBa_hist* next_hi = get_history(hi, cin);
    do {

      binaries_read++;

      if (hi->get_last()->get_binary_type() == Merged) 
	mergers++;
//	cerr << hi->get_parameter(identity) << endl;

      if(!T_flag) {
	extract_binary_on_stellar_types(hi, primary_type_string,
					secondary_type_string,
					bt,
					a_min, a_max, 
					mp_min, mp_max, ms_min, ms_max,
					e_min, e_max, q_min, q_max,
					dt, first_occasion,
					R_flag);
      }
      else {
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
      }

      if (v_flag)
	hi->put_history(cout, v_flag);
      

      //      if (R_flag) 
      //	for_all_SeBa_hist(SeBa_hist, hi, ho) cerr << *ho;

      delete hi;
      
      hi = next_hi;
      next_hi = get_history(hi, cin);
    }
    while (next_hi);

     extract_binary_on_stellar_types(hi->get_last(), primary_type_string,
				     secondary_type_string,
				     bt,
				     a_min, a_max, 
				     mp_min, mp_max, ms_min, ms_max,
				     e_min, e_max, q_min, q_max, 
				     dt, first_occasion,
				     R_flag);

     if (v_flag) {
       hi->put_history(cout, v_flag);
       cerr <<"Number of Mergers: " << mergers << endl;
     }


    if(T_flag) {

      if (N_detached>0) {
	cerr << "     Detached population"<<endl;
	print_binary_matric(detached_population);
      }
      else
	cerr << "     ---No detached binaries---" <<endl;
      if (N_contact>0) {
	cerr << "     Semi-detached/contact population"<<endl;
	print_binary_matric(contact_population);
      }
      else
	cerr << "     ---No contact Binaries---" <<endl;

      cerr << "     Single population"<<endl;
      cerr << "    ";
      for (int j=0; j<no_of_stellar_type-1; j++) {
	cerr << " " << type_short_string((stellar_type)j);
	if(j==(int)Main_Sequence || 
	   j==(int)Super_Giant   || 
	   j==(int)Thorn_Zytkow)
	  cerr << "  "; 
      }

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

    }
//     if (R_flag) 
//       for_all_SeBa_hist(SeBa_hist, hi, ho) cerr << *ho;

    cerr << "Total number of binaries read: " << binaries_read << endl;
}

#endif // endof: TOOLBOX



