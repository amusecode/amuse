#include "SeBa_hist.h"
 
bool scenarios_identical(SeBa_hist* hi, SeBa_hist* ha) {

  bool the_same = true;
  if(*hi->get_last() != *ha->get_last())
    the_same = false;
  else {
    SeBa_hist *SeBa_hist_hi, *SeBa_hist_hj;
    for (SeBa_hist_hi = hi->get_first(),
	 SeBa_hist_hj = ha->get_first();
	 (SeBa_hist_hi != NULL && SeBa_hist_hj != NULL);
	 SeBa_hist_hi = SeBa_hist_hi->get_future(),
	 SeBa_hist_hj = SeBa_hist_hj->get_future()) {

      if(*SeBa_hist_hi != *SeBa_hist_hj)
	the_same = false;
      
    }
  }

  return the_same;
}

bool SeBa_hist::operator == (SeBa_hist& ha) const {

  bool the_same = false;
  if(bin_tpe == ha.get_binary_type()) {
    if(bin_tpe == Merged && tpe_prim == ha.get_primary_type()) {
      
      the_same = true;
    }
    else if(tpe_prim == ha.get_primary_type() &&
	    tpe_sec == ha.get_secondary_type()) {
      
      the_same = true;
    }
  }

  return the_same;
}

bool SeBa_hist::operator != (SeBa_hist& ha) const {

  bool different = true;
  if (*this == ha)
    different = false;

  return different;
}


ostream& operator<<(ostream& s, SeBa_hist& hi) {

  s << hi.number << " "<< hi.bin_tpe << " " << hi.mttype
    << " " << hi.time;
    // short dump output is without bin_type
    //<< " " << (int)hi.bin_tpe; 
    s << "\t" << hi.semi << " " << hi.ecc; 
    s << "\t" << hi.label_prim << " " << hi.tpe_prim 
//    s << "\t" << hi.tpe_prim 
              << " " << hi.m_prim << " " << hi.r_prim;
    s << "\t" << hi.label_sec << " " << hi.tpe_sec 
//    s << "\t" << hi.tpe_sec 
              << " " << hi.m_sec << " " << hi.r_sec;
    s << endl;

    return s;

}


void SeBa_hist::put_history(ostream& s, bool verbose=false) {

     if (verbose) {
	 s << "\nnumber= " << number <<"\tTime= "<< time;
	 s << "\t" <<  type_string(bin_tpe) 
	   << "\t" <<  type_string(mttype) 
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

  s << number << " "<< bin_tpe << " " << mttype << " " << time;
     s << "\t" << semi << " " << ecc; 
     s << "\t" << label_sec << " " << tpe_sec << " " 
//     s << "\t" << tpe_sec << " " 
               << m_sec << " " << r_sec;
     s << "\t" << label_prim << " " <<tpe_prim << " " 
//     s << "\t" <<tpe_prim << " " 
               << m_prim << " " << r_prim;
     s << endl;
 }

real SeBa_hist::get_parameter(binary_parameter param) {

     switch(param) {
     case identity:                          return number;
          break;
     case bin_type:                          return bin_tpe;
          break;
     case mtrtype:                           return mttype;
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
    
    int tpe_1, tpe_2;
    real lT_prim, lT_sec;
    real mc_prim, mc_sec;

    int ibtp, imttp;
    s >> number >> ibtp >> imttp >> time >> semi >> ecc 
//      >> tpe_1 >> m_prim >> r_prim
//      >> tpe_2 >> m_sec >> r_sec;
      >> label_prim >> tpe_1 >> m_prim >> r_prim >> lT_prim >> mc_prim
      >> label_sec >> tpe_2 >> m_sec >> r_sec >> lT_sec >> mc_sec;

    bin_tpe = (binary_type)ibtp;
    mttype = (mass_transfer_type)imttp;
    if(s.eof()) return false;


    tpe_prim = (stellar_type)tpe_1;
    tpe_sec = (stellar_type)tpe_2;

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
      if(m_prim<=0) {
	m_prim = m_sec;
	tpe_prim = tpe_sec;
	r_prim = r_sec;
      }
    }
    else if (r_prim<=0 && r_sec<=0 && semi<100000) 
      bin_tpe = Contact;
    else if(r_prim<=0 || r_sec<=0 && semi<100000)  
      bin_tpe = Semi_Detached;
    
//    if (bin_tpe != Merged && bin_tpe != Disrupted) {
//      r_prim = set_stellar_radius(true);
//      r_sec  = set_stellar_radius(false);
//    } else {
//      r_prim = -r_prim;
//      r_sec  = -r_sec;
//    }
    
    return true;

}

void SeBa_hist::move_SeBa_hist_to(SeBa_hist *next_hi) {
    
    next_hi = this;
    next_hi->future = NULL;
    next_hi->past   = NULL;
}

bool SeBa_hist::binary_contains(char *prim_string,
				char *sec_string,
				binary_type bt = Detached,
				mass_transfer_type mt = Unknown) {


  if (bin_tpe== bt && (mttype == mt || mt == Unknown) &&
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
  return 0;
}


SeBa_hist* get_history(SeBa_hist *hi, istream& is) {

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

void SeBa_hist::add_to_SeBa_hist(SeBa_hist *next_hi) {

  cerr << "add: "; 
  next_hi->put_state(cerr);
  cerr << "\nto ";
  put_state(cerr);
  cerr << endl;

  SeBa_hist *last = get_last();
  last->set_future(next_hi);
  next_hi->set_past(last);
}

/*-----------------------------------------------------------------------------
 *  put_state --
 *-----------------------------------------------------------------------------
 */

void SeBa_hist::put_state(ostream & s) {

  switch (bin_tpe) {
  case Merged:
    s << "{" << type_short_string(tpe_prim) << "}";
    break;
  case Disrupted:
    s << ")";
    s << type_short_string(tpe_prim) << ",";
    s << type_short_string(tpe_sec);
    s << "(";
    break;
  case Detached:
    s << "(";
    s << type_short_string(tpe_prim) << ",";
    s << type_short_string(tpe_sec);
    s << ")";
    break;
  case Contact:
    s << "[";
    s << type_short_string(tpe_prim) << ",";
    s << type_short_string(tpe_sec);
    s << "]";
    break;
  case Semi_Detached:
    if (r_prim<=0 && semi<100000)  {
      s << "[";
      s << type_short_string(tpe_prim) << ",";
      s << type_short_string(tpe_sec);
      s << ")";
    }
    else {
      s << "(";
      s << type_short_string(tpe_prim) << ",";
      s << type_short_string(tpe_sec);
      s << "]";
    }
    break;
  default:
    cerr << "Unknown binary type in SeBa_hist->put_state() "<< endl;
  }
  s << ";";
}

void put_state(SeBa_hist * hi, ostream & s) {

  for_all_SeBa_hist(SeBa_hist, hi, ha) {
    ha->put_state(s);
  }

  s << endl;
}

