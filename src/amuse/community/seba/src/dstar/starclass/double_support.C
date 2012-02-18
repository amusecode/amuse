//
// double_support: Helper functions to aid in setting up and manipulating
//

#include "double_support.h"
#include "double_star.h"

binary_type extract_binary_type_string(const char *type) {

  if (!strcmp("synchronized", type))          return Synchronized;
  else if (!strcmp("strong_encounter", type)) return Strong_Encounter;
  else if (!strcmp("detached", type))         return Detached;
  else if(!strcmp("semi_detached", type))     return Semi_Detached;
  else if(!strcmp("contact", type))           return Contact;
  else if(!strcmp("common_envelope", type))   return Common_Envelope;
  else if(!strcmp("double_spiral_in", type))  return Double_Spiral_In;
  else if(!strcmp("merged", type))            return Merged;
  else if(!strcmp("disrupted", type))         return Disrupted;
  else if(!strcmp("spiral_in", type))         return Spiral_In;
  
  return Unknown_Binary_Type;

}

const char* type_string(binary_type tpe) {
   
  switch(tpe) {
  case Strong_Encounter:	return "strong_encounter";
  case Synchronized:		return "synchronized";
  case Detached:		return "detached";
  case Semi_Detached:		return "semi_detached";
  case Contact:			return "contact";
  case Common_Envelope:		return "common_envelope";
  case Double_Spiral_In:	return "double_spiral_in";
  case Merged:			return "merged";
  case Disrupted:		return "disrupted";
  case Spiral_In:		return "spiral_in";
  default:			return "unknown_binary_type";
  }
}

//
// double_hist: helper history function.
//

void double_hist::put_double_hist() {

  cout << binary_age<< " " << semi <<" "<< eccentricity<<endl;
}


//
// double_init: helper history function.
//

void double_init::read_element() {

  cin >> mass_prim
      >> q
      >> semi
      >> eccentricity;
}

void double_init::put_element() {

  ofstream outfile("binev.data", ios::app|ios::out);
  if(!outfile) cerr << "error: couldn't create file binev.data"<<endl;

  outfile << mass_prim << " "
	  << q << " "
	  << semi << " "
	  << eccentricity << endl;
}

void double_init::dump(ostream & s) {

  s << " " << mass_prim
    << " " << mass_prim*q
    << " " << semi
    << " " << eccentricity
    << endl;

}

void double_init::dump(const char* filename) {

  ofstream s(filename, ios::app|ios::out);
  if(!s) cerr << "error: couldn't create file "<<filename<<endl;

  s << " " << mass_prim
    << " " << mass_prim*q
    << " " << semi
    << " " << eccentricity
    << endl;
}


double_star * get_new_binary(double_init& init, const int id) {

  //        double_star * binary = new double_star;

  //        binary->initialize(init, id); 

  //        return binary;
  return NULL;
}

void pptime(real time,
	    ostream & s,		// default = cerr
	    const char *t) {		// default = "time"

  if (time<1.e-3) 
    ppperiod(365.25*time*1.e+6, s, t);
  else {
    s << "   " << t << " = ";
    if (time < 1)
      s << time*1000 << " kyr";
    else if (time < 1000)
      s << time << " Myr";
    else if (time < 1.e+6)
      s << time/1000 << " Gyr";
    else if (time < 1.e+9)
      s << time/1.e+6 << " Tyr";
    else if (time < 1.e+12)
      s << time/1.e+9 << " Pyr";
    else 
      s << time/1.e+12 << " Eyr";
  }
}

void ppperiod(real period,
	      ostream & s,		// default = cerr
	      const char *p) {		// default ="Porb"

  if (period>1000*365.25)
    pptime(period*1.e-6/365.25, s, p);
  else {
    s << "   " << p << " = ";
    if (period < 6.94e-4) 
      s << period*1440 << " seconds";
    else if (period < 0.0417) 
      s << period*1440 << " minutes";
    else if (period < 1) 
      s << period*24 << " hours";
    else if (period < 365.25) 
      s << period << " days";
    else
      s << period/365.25 << " yr";
  }

}

/*-----------------------------------------------------------------------------
 *  put_state --
 *-----------------------------------------------------------------------------
 */

void put_state(double_state d,
	       ostream & s) {

  switch (d.type) {
  case Merged:
    s << "{";
    d.primary.put_star_state(s);
    s << "}";
    break;
  case Disrupted:
    s << ")";
    d.primary.put_star_state(s);
    s << ", ";
    d.secondary.put_star_state(s);
    s << "(";
    break;
  default:
    if (d.primary.class_spec[Rl_filling])
      s << "[";
    else
      s << "(";
    d.primary.put_star_state(s);
    s << ", ";
    d.secondary.put_star_state(s);
    if (d.secondary.class_spec[Rl_filling])
      s << "]";
    else
      s << ")";
    real period = 2*PI*sqrt(pow(d.semi
				*cnsts.parameters(solar_radius), 3.) 
			    / (cnsts.physics(G)*cnsts.parameters(solar_mass)
			       *(d.primary.mass+d.secondary.mass)))
      / cnsts.physics(days);

    ppperiod(period, s, "Porb");
  }
  s << endl;
}

/*-----------------------------------------------------------------------------
 *  make_profile --
 *-----------------------------------------------------------------------------
 */
void  make_profile(int id, real start_time,
                   double_profile& binary, double_init& init) {

  cerr<<"double_support make_profile"<<endl;
  real dt = (init.end_time - start_time)/init.n_steps; 

  star_state primary, secondary;
  double_star * b = get_new_binary(init, id);

  primary.init_star_state((star*)b->get_primary());
  secondary.init_star_state((star*)b->get_secondary());
  binary.init_double_profile((double_star*) b, primary, secondary);
  //	      put_initial_conditions(init);

  //              Evolve binary.
  for (int j = 0; j<init.n_steps; j++)
    b->evolve_element(start_time+dt*(j+1));

//              Create binary profile
  if (primary.identity == b->get_primary()->get_identity()) {
    primary.make_star_state((star*)b->get_primary());
    secondary.make_star_state((star*)b->get_secondary());

    //         put_state(b, primary, secondary);
    binary.enhance_double_profile((double_star*) b, primary, secondary);
  }
  else {
    secondary.make_star_state((star*)b->get_primary());
    primary.make_star_state((star*)b->get_secondary());

    //         put_state(b, secondary, primary);
    binary.enhance_double_profile((double_star*) b, secondary, primary);
  }

}

void double_profile::init_double_profile(double_star* b,
                                         star_state& p, star_state& s) {

  init.identity = b->get_identity();
  init.time = b->get_current_time();
  init.type = b->get_bin_type();
  init.total_mass = b->get_total_mass();
  init.primary = p;
  init.secondary = s;
}

void double_profile::enhance_double_profile(double_star* b,
                                            star_state& p, star_state& s) {

  final.identity = b->get_identity();
  final.time = b->get_current_time();
  final.type = b->get_bin_type();
  final.total_mass = b->get_total_mass();
  final.primary = p;
  final.secondary = s;

  mdot = init.total_mass - final.total_mass;
}

void double_profile::init_double_profile(double_star* b) {

  init.identity = b->get_identity();
  init.time = b->get_current_time();
  init.type = b->get_bin_type();
  init.semi = b->get_semi();
  init.ecc  = b->get_eccentricity();
  init.velocity = b->get_velocity();
  init.total_mass = b->get_total_mass();

  star_state pa, sa;
  if(b->get_use_hdyn()) {
    pa.make_star_state(dynamic_cast(star*, b->get_primary()));
    sa.make_star_state(dynamic_cast(star*, b->get_secondary()));
  }
  else {
    pa.make_star_state(dynamic_cast(star*, b->get_initial_primary()));
    sa.make_star_state(dynamic_cast(star*, b->get_initial_secondary()));
  }

  init.primary = pa;
  init.secondary = sa;
}

void double_profile::enhance_double_profile(double_star* b) {

  final.identity = b->get_identity();
  final.time = b->get_current_time();
  final.type = b->get_bin_type();
  final.semi = b->get_semi();
  final.ecc  = b->get_eccentricity();
  final.velocity = b->get_velocity();
  final.total_mass = b->get_total_mass();

  star_state pa, sa;
  if(b->get_use_hdyn()) {
    pa.make_star_state(dynamic_cast(star*, b->get_primary()));
    sa.make_star_state(dynamic_cast(star*, b->get_secondary()));
  }
  else {
    pa.make_star_state(dynamic_cast(star*, b->get_initial_primary()));
    sa.make_star_state(dynamic_cast(star*, b->get_initial_secondary()));
  }

  final.primary = pa;
  final.secondary = sa;
}

void put_profile(double_profile& d) {

  cerr <<"Initial: " << endl;
  put_state(d.init, cerr);
  cerr << "final: " << endl;
  put_state(d.final, cerr);

}

#if 0
double_star * triple_star(double_init& inner,
                          double_init& outer, int id) {

      single_star * p = new single_star();
      p->initialize(inner.mass_prim, id, inner.start_time);
      single_star * s = new single_star();
      s->initialize(inner.mass_prim*inner.q, id+1, inner.start_time);

//              Inner binary
      double_star * b = new double_star(p, s, inner);

      single_star * t = new single_star();
      t->initialize(outer.mass_prim*outer.q, id+2, outer.start_time);

//              Outer binary
      double_star * triple = new double_star(b, t, outer);

      return triple;
   }
#endif

double_state make_state(double_star* b) {

  double_state dbl;
  dbl.identity = b->get_identity();
  dbl.time = b->get_current_time();
  //dbl.type = b->get_bin_type(); // substituted by (SPZ:2/1998)
  dbl.type = b->obtain_binary_type();
  dbl.semi = b->get_semi();
  dbl.ecc  = b->get_eccentricity();
  dbl.velocity = b->get_velocity();
  dbl.total_mass = b->get_total_mass();

  star_state pa, sa;
  if(b->get_use_hdyn()) {
    pa.make_star_state(dynamic_cast(star*, b->get_primary()));
    sa.make_star_state(dynamic_cast(star*, b->get_secondary()));
  }
  else {
    pa.make_star_state(dynamic_cast(star*, b->get_initial_primary()));
    sa.make_star_state(dynamic_cast(star*, b->get_initial_secondary()));
  }
  
  dbl.primary = pa;
  dbl.secondary = sa;

  return dbl;
}

