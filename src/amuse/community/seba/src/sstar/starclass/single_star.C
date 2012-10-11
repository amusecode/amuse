// 
// single_star.C
//
#include <fstream>
#include <string>
#include "single_star.h"
#include "proto_star.h"
//#include "main_sequence.h"

// Only to make SPZDCH star known to this file 
#include "SPZDCH_star.h"

#define REPORT_MASS_TRANSFER_TIMESCALE false

single_star * new_single_star(stellar_type type,	// All defaults are
			      int  id,			// now specified in
			      real z,
			      real t_cur,
			      real t_rel,
			      real m_rel,
			      real m_tot,
			      real m_core,
			      real co_core,
			      real p_rot,
			      real b_fld,
			      node* n)
{ single_star* element = 0;
  switch(type) {
  case SPZDCH_Star: element = new SPZDCH_star(n);
    break;
  case Proto_Star: element = new proto_star(n);
    break;
  case Planet:
  case Brown_Dwarf: element = new brown_dwarf(n);
    break;
  case Main_Sequence: element = new main_sequence(n);
    break;
  case Hertzsprung_Gap: element = new hertzsprung_gap(n);
    break;
  case Sub_Giant: element = new sub_giant(n);
    break;
  case Horizontal_Branch: element = new horizontal_branch(n);
    break;
  case Super_Giant: element = new super_giant(n);
    break;
  case Carbon_Star:  
  case Helium_Star:  
  case Helium_Giant: element = new helium_star(n);
    break;
  case Hyper_Giant: element = new hyper_giant(n);
    break;
  case Helium_Dwarf:
  case Carbon_Dwarf:
  case Oxygen_Dwarf: element = new white_dwarf(n);
    break;
  case Thorn_Zytkow: element = new thorne_zytkow(n);
    break;
  case Xray_Pulsar:
  case Radio_Pulsar:
  case Neutron_Star: element = new neutron_star(n);
    element->set_rotation_period(p_rot);
    element->set_magnetic_field(b_fld);
    break;
  case Black_Hole: element = new black_hole(n);
    break;
  case Disintegrated: element = new disintegrated(n);
    break;
  default: cerr << "No default stellar type in new_single_star()"<<endl;
           exit(1);
		    
  }
  
  element->initialize(id, z, t_cur, t_rel, m_rel, m_tot, m_core, co_core);

  element->post_constructor();

  //element->get_seba_counters()->add_sstar++;

  return element;
}

single_star::single_star(node* n) : star(n) {

    star_type = NAS;
    for (int i=NAC; i<no_of_spec_type; i++)
      spec_type[i] = NAC;

    metalicity = 0;
    current_time=relative_age=0;
    last_update_age = next_update_age=0;
    relative_mass = accreted_mass=0;
    envelope_mass=core_mass=0;
    core_radius=effective_radius=radius=0;
    COcore_mass = 0;
    luminosity=0;
    velocity=0;
    wind_constant=0;
    magnetic_field = rotation_period = 0;
    birth_mass=0;
}


single_star::single_star(single_star & rv) : star(rv) {

  identity      = rv.identity;

  for (int i=NAC; i<no_of_spec_type; i++)
    spec_type[i] = rv.spec_type[i];

  //		copy star
  metalicity  = rv.metalicity;

  current_time  = rv.current_time;
  relative_age  = rv.relative_age;
  relative_mass = rv.relative_mass;
  envelope_mass = rv.envelope_mass;
  core_mass     = rv.core_mass;
  COcore_mass   = rv.COcore_mass;

  core_radius   = rv.core_radius;
  radius        = rv.radius;
  luminosity    = rv.luminosity;
  effective_radius = rv.effective_radius;
  last_update_age = rv.last_update_age;
  next_update_age = rv.next_update_age;
  velocity        = rv.velocity;
  wind_constant   = rv.wind_constant;
  accreted_mass   = rv.accreted_mass;
  magnetic_field  = rv.magnetic_field;
  rotation_period = rv.rotation_period;
  birth_mass      = rv.birth_mass;
			
  //              copy stellar history.
  previous.current_time    = rv.previous.current_time;
  previous.last_update_age = rv.previous.last_update_age;
  previous.next_update_age = rv.previous.next_update_age;
  previous.relative_age    = rv.previous.relative_age;
  previous.relative_mass   = rv.previous.relative_mass;
  previous.envelope_mass   = rv.previous.envelope_mass;
  previous.core_mass       = rv.previous.core_mass;
  previous.COcore_mass     = rv.previous.COcore_mass;

  previous.radius          = rv.previous.radius;
  previous.effective_radius= rv.previous.effective_radius;
  previous.magnetic_field  = rv.previous.magnetic_field;
  previous.rotation_period = rv.previous.rotation_period;
  previous.birth_mass      = rv.previous.birth_mass;

}

#if 0
void single_star::post_constructor() {
    //(GN+SPZ Apr 28 1999) stars with M < 8 Msun have wind per phase:
  update_wind_constant();

    // MEmory refrsh to prevent helium star to become white dwarf with
    // zero-mass core.
    // Not clear yet if companion and binary must also be updated?
    // (SPZ+GN, 10 Nov 1998)
      if (is_binary_component())
          get_binary()->refresh_memory();
      else
          refresh_memory();


    if (is_binary_component() && get_binary()->get_bin_type()==Detached)
      get_binary()->set_first_contact(false);

    if (is_binary_component()){
        get_binary()->dump("SeBa.data", true);
    }
    else {
    }
// (GN May 12 1999) skrews up output
//    {
//      dump("SeBa.data", true);
//      cerr << endl;
//    }


    if (remnant()) {
        if (is_binary_component()){
	  get_binary()->dump("binev.data", false);
        }
	else
	  dump("binev.data", false);
    }
}

#endif

void single_star::post_constructor() {

  //(GN+SPZ Apr 28 1999) stars with M < 8 Msun have wind per phase:
  update_wind_constant();

  // MEmory refrsh to prevent helium star to become white dwarf with
  // zero-mass core.
  // Not clear yet if companion and binary must also be updated?
  // (SPZ+GN, 10 Nov 1998)
  if (is_binary_component())
    get_binary()->refresh_memory();
  else
    refresh_memory();
  //    refresh_memory();
    
  //    if (is_binary_component() && get_binary()->get_bin_type()==Detached)
  if (is_binary_component() && 
      get_binary()->get_bin_type() != Merged &&
      get_binary()->get_bin_type() != Disrupted) {

    get_binary()->set_first_contact(false);

    if (remnant() || !hydrogen_envelope_star()) {
      get_binary()->set_bin_type(Detached);
      get_binary()->set_current_mass_transfer_type(Unknown);
      
    }
  }

  if (is_binary_component() && 
      get_binary()->roche_radius(this) > radius) {
    get_binary()->dump("SeBa.data", true);
  }

  if (remnant()) {
    if (is_binary_component()) {
      get_binary()->dump("binev.data", false);
    }
   else
      dump("binev.data", false);
  }

}


void single_star::initialize(int id, real z, real t_cur,
			     real t_rel, real m_rel, real m_tot,
			     real m_core, real co_core) {

  real m_env = m_tot - m_core;
//  if (m_rel<=0 || m_rel>cnsts.parameters(maximum_main_sequence)) {
//    cerr << "Mass of initial star is prefered to be "
//	 << "\nwithin reasonable limits <0, "
//	 << cnsts.parameters(maximum_main_sequence)
//	 << "] \nbut found (mass = " << m_rel << ")" << endl;
//    if (m_rel<=0)
//      exit (-1);
//  }

  identity = id;
  luminosity = 0.01;
  metalicity = z;
  current_time = t_cur;
  relative_age = t_rel;
  relative_mass = max(m_rel, m_env+m_core);
  previous.envelope_mass = envelope_mass = m_env;
  core_mass = m_core;
  COcore_mass = co_core;

  instantaneous_element();

  // (SPZ+GN: 26 Jul 2000) 
  // update_wind_constant() must be called before
  update_wind_constant();
  adjust_next_update_age();

  // Was previously done after adjust_next_update_age, which may be wrong?
  instantaneous_element();

  update();
}

#if 0
// Determine size of stellar core from fits to
// temperature and luminosity from figure of
// Iben, I.Jr., and Tutukov, A.V., 1985, ApJSS 58, 661.
real single_star::helium_core_radius() {

  // Safety. Should ne be required.
  if (core_mass==0) {
    cerr << "WARNING single_star::helium_core_radius():" << endl;
    cerr << "        core_radius == 0" << endl;
    cerr << "        action: return 0;" << endl;
    return 0;
  }

  real r_he_core;
  real tmp = log10(core_mass)
            * (0.4509 - 0.1085*log10(core_mass)) + 4.7143;
  real lum  = 8.33*tmp - 36.8;
  tmp = 1.e-3*pow(10., tmp);
  lum  = pow(10., lum);

  // helium core is smaller than helium star of same mass.
  // reasoning: after spiral-in class helium_star is supposed to
  // handle further mass transfer.
  //  real fraction = 0.5; // (SPZ+GN: 6 Oct 1998)
  real fraction =1.; // But gives unrealistic mass transfer after spiral-in
  r_he_core = fraction*33.45*sqrt(lum)/pow(tmp, 2);
 
  return min(radius, r_he_core);
}

#endif


real single_star::nucleair_evolution_time(const real mass, const real mass_tot,
					  const real z) {

    real t_tagb = TAGB_time(mass, mass_tot, z);
    return t_tagb;
}


real single_star::nucleair_evolution_time() {
    cerr<<"single_star::nucleair_evolution_time is currently not used (used to be part of wind_constant)"<<endl;
    real t_nuc = nucleair_evolution_time(relative_mass, get_total_mass(),
                       metalicity);
    return t_nuc;
}

real single_star::temperature() {
    // Effective radius is the radius of the star as it really is.
    // radius is the equilibrium radius of the star.
   //return 1000*pow(1130.*luminosity/(radius*radius), 0.25); // in [1000K]


    real T_eff = cnsts.parameters(Tsun)
             * pow(luminosity/pow(effective_radius, 2), 0.25); // in [K]
 
    // (GN+PimvO Jun  8 2012) safety for mergers
    if (effective_radius <= 0) T_eff = 1e-99;

    return T_eff;
}

real single_star::magnitude() {

  return -2.5*log10(luminosity) + 4.75 - bolometric_correction();
}

// real function bolometric_correction()
// calculates bolometric correction for giant stars
// input: stellar effective surface temperature in kKelvin.
// output: bolometric correction for giant star in cgs.
real single_star::bolometric_correction() {

  // temperature() is defined in Kelvin.
  // here we should use old 10^3K implementation 
  // (SPZ+GN: 1 Oct 1998)
  real temp_in_kK = 0.001 * temperature();
  real bc;

  if (temp_in_kK < 4.195)
    bc = 2.5*log10((1.724e-7*pow(temp_in_kK,11.) + 1.925e-2)
		   / (1. + 1.884e-9*pow(temp_in_kK,14.)));
  else if (temp_in_kK >=  4.195 &&
	   temp_in_kK <= 10.89)
    bc = 2.5*log10((7.56e-2*pow(temp_in_kK,1.5))
		   / (1. + 6.358e-5*pow(temp_in_kK, 4.5)));
  else
    bc = 2.5*log10((2728/pow(temp_in_kK,3.5) + 1.878e-2*temp_in_kK)
		   /(1. + 5.362e-5*pow(temp_in_kK,3.5)));

  return bc;
}

// wind_velocity(): Calculates stellar wind velocoty.
// Steller wind velocity is 2.5 times stellar escape velocity
real single_star::wind_velocity() {

  real v_esc2 = cnsts.physics(G)
              * cnsts.parameters(solar_mass)*get_total_mass()
              / (effective_radius*cnsts.parameters(solar_radius));
  real v_wind = 2.5*sqrt(v_esc2)/cnsts.physics(km_per_s);

  return v_wind;
}

real single_star::kelvin_helmholds_timescale() {

  // Equilibrium radius is 'single_star::radius'
  // Effective radius may be puffed up by accretion or subtraction of mass.

  return 31.56*pow(relative_mass,2.)/(radius*luminosity); // [Myr]
}

real single_star::nucleair_evolution_timescale() {
  // overloaded for main_sequence:: as
  // t_nuc = 10^10 [years] Msun/Lsun.
  // Assumed that 0.1 Msun is thermalized.

  // but for giants etc. we use the following, based on
  // dimensional grounds and empirical comparison with
  // Pols, 1994, A&A 290, 119
  // (SPZ+GN:30 Sep 1998)
  
//  real fused_mass = 0.1*relative_mass;
//
//  real t_nuc = cnsts.parameters(energy_to_mass_in_internal_units)
//             * fused_mass/luminosity;
//  
//  real t_kh = kelvin_helmholds_timescale();
//
//  return sqrt(t_nuc * t_kh);

// (GN+SPZ Apr 28 1999) giant lifetime is ~ 10% of ms life time
//  return 0.1*main_sequence_time();


// (GN+SilT August 5 2010) t_nuc ~ delta t * R/ delta R
// old prescription gave long timescales which destables the mass transfer
     real t_nuc = radius * (current_time - previous.current_time)/(radius -
previous.radius);

    if (t_nuc < 0) t_nuc = 0.1*main_sequence_time();
    return t_nuc;
}

real single_star::dynamic_timescale() {

  return 5.08e-11*sqrt(pow(radius, 3.)/relative_mass);    // [Myr]
}



void single_star::read_element() {

  real total_mass, log_temp, log_lum;
  real bol_corr, temp, mv;
  int type = (int)star_type;
  cin >> type >> relative_age >> relative_mass
      >> total_mass >> core_mass >> radius >> log_temp
      >> log_lum >> bol_corr >> mv >> velocity
      >> magnetic_field >>  rotation_period;

  envelope_mass = total_mass - core_mass;
  temp = pow(log_temp, 10);
  luminosity = pow(log_lum, 10);
}

void single_star::put_element() {

  ofstream outfile("binev.data", ios::app|ios::out);
  if (!outfile) cerr << "error: couldn't create file binev.data"<<endl;

  real bol_corr = bolometric_correction();
  real mv = magnitude();
  outfile << star_type << " " << relative_age << " "
	  << relative_mass << " " << get_total_mass() << " "
	  << core_mass << " " << radius << " "
	  << log10(temperature()) << " "
	  << log10(luminosity) << " " << bol_corr << " "
	  << mv << " " << " " << velocity << " "
	  << magnetic_field <<  rotation_period << endl;

  outfile.close();
}

void single_star::dump(ostream & s, bool brief) {
    s.precision(HIGH_PRECISION);
  if (brief) {
    
    s << get_node()->format_label() << " "
      << get_current_time() << " "
      << "0 1   ";
    s << get_node()->format_label() << " "	  
      << get_element_type() << " "
      << get_total_mass() << " "
      << radius << "    ";
    s << "0 0 0 0" << endl;
  }
  else
  {   
      real bol_corr = bolometric_correction();
      int s1, s2, s3, s4, s5, s6;
      s1 = get_spec_type(Emission);
      s2 = get_spec_type(Blue_Straggler);
      s3 = get_spec_type(Barium);
      s4 = get_spec_type(Rl_filling) + get_spec_type(Accreting);
      s5 = get_spec_type(Runaway);
      s6 = get_spec_type(Merger);

      s << "1\n ";
      s << get_element_type() << " " << s1 << " " << s2 << " "
        << s3 << " " << s4 << " " << s5 << " " << s6 << endl;
      
      s << " " << identity
	<< " " << metalicity
	<< " " << current_time
	<< " " << relative_age 
	<< " " << relative_mass
	<< " " << get_total_mass()
	<< " " << core_mass
	<< " " << radius
	<< " " << effective_radius
	<< " " << log10(temperature())
	<< " " << log10(luminosity)
	<< " " << bol_corr
	<< " " << magnitude()
	<< " " << velocity
	<< " " << magnetic_field
	<< " " << rotation_period
	<< " " << core_radius
	<< endl;

      s.precision(STD_PRECISION);
#if 0
      PRC(identity);
	PRC(metalicity);
	PRC(current_time);
	PRC(relative_age);
	PRC(relative_mass);
	PRC(get_total_mass());
	PRC(core_mass);
	PRC(radius);
	PRC(effective_radius);
	PRC(log10(temperature()));
	PRC(log10(luminosity));
	PRC(bol_corr);
	PRC(magnitude());
	PRC(velocity);
	PRC(magnetic_field);
	PRL(rotation_period);
#endif
    }
}

void single_star::dump(char * filename, bool brief) {

  ofstream s(filename, ios::app|ios::out);
  if (!s) cerr << "error: couldn't create file "<<filename<<endl;

  if (brief) {
    
    s << get_node()->format_label() << " "
      << get_current_time() << " "
      << "0 1   ";
    s << get_node()->format_label() << " "	  
      << get_element_type() << " "
      << get_total_mass() << " "
      << radius << " ";
    s << "0 0 0 0" << endl;

  }
  else
    {

      real bol_corr = bolometric_correction();
      int s1, s2, s3, s4, s5, s6;
      s1 = get_spec_type(Emission);
      s2 = get_spec_type(Blue_Straggler);
      s3 = get_spec_type(Barium);
      s4 = get_spec_type(Rl_filling) + get_spec_type(Accreting);
      s5 = get_spec_type(Runaway);
      s6 = get_spec_type(Merger);

      s << "1\n ";
      s << get_element_type() << " " << s1 << " " << s2 << " "
        << s3 << " " << s4 << " " << s5 << " " << s6 << endl;
      s <<" " << identity
	<< " " << metalicity
	<< " " << current_time
	<< " " << relative_age
	<< " " << relative_mass
	<< " " << get_total_mass()
	<< " " << core_mass
	<< " " << radius
	<< " " << effective_radius
	<< " " << log10(temperature())
	<< " " << log10(luminosity)
	<< " " << bol_corr
	<< " " << magnitude()
	<< " " << velocity
	<< " " << magnetic_field
	<< " " << rotation_period
	<< endl;

    }
  s.close();
}

void single_star::put_state() {

  star_state str;
  str.init_star_state(dynamic_cast(single_star*, this));
  str.make_star_state(dynamic_cast(single_star*, this));

  str.put_star_state();
}

void single_star::put_hrd(ostream & s) {

  int s1, s2, s3, s4, s5, s6;
  s1 = get_spec_type(Emission);
  s2 = get_spec_type(Blue_Straggler);
  s3 = get_spec_type(Barium);
  s4 = get_spec_type(Rl_filling) + get_spec_type(Accreting);
  s5 = get_spec_type(Runaway);
  s6 = get_spec_type(Merger);

  s << "1\n ";
  s << s1 << " " << s2 << " " << s3 << " " << s4 << " "
    << s5 << " " << s6 << " ";
  s << get_total_mass()
    << " " << log10(temperature()) 
    << " " << log10(luminosity) 
    << endl;
}

void single_star::refresh_memory() {

  previous.star_type = get_element_type();
  previous.metalicity = metalicity;
  previous.current_time = current_time;
  previous.last_update_age = last_update_age;
  previous.next_update_age = next_update_age;
  previous.relative_age = relative_age;
  previous.relative_mass = relative_mass;
  previous.envelope_mass = envelope_mass;
  previous.core_mass = core_mass;
  previous.COcore_mass = COcore_mass;

  previous.radius = radius;
  previous.effective_radius = effective_radius;

  //             Pulsars
  previous.magnetic_field  = magnetic_field;
  previous.rotation_period = rotation_period;
  previous.birth_mass      = birth_mass;
}

void single_star::recall_memory() {

  star_type        = previous.star_type;
  metalicity       = previous.metalicity;
  current_time     = previous.current_time;
  last_update_age  = previous.last_update_age;
  next_update_age  = previous.next_update_age;
  relative_age     = previous.relative_age;
  relative_mass    = previous.relative_mass;
  envelope_mass    = previous.envelope_mass;
  core_mass        = previous.core_mass;
  COcore_mass      = previous.COcore_mass;

  radius           = previous.radius;
  effective_radius = previous.effective_radius;

  // Pulsars
  magnetic_field   = previous.magnetic_field;
  rotation_period  = previous.rotation_period;
  birth_mass       = previous.birth_mass;
}

// Mass transfer timescale is checked and updated at the moment the
// mass-ratio is reversed.

// in version 3.3, the mass transfer timescale is updated each
// timestep in ::recursive(), therefore this function is not required.
// It hangs the code
// because mass loss by stellar wind occurs after
// ::mass_ratio_mdot_limit().
// therefore the mass ratio > 1 after ::recursive step.
// (SPZ+GN:11 Oct 1998)
real single_star::mass_ratio_mdot_limit(real mdot) {

    
    // No limit for the moment.
    return mdot;
    
    real accretor_mass = 0;

    if (is_binary_component()) 
      get_companion()->get_total_mass();

    if (accretor_mass<get_total_mass()) {
	real mdot_max = get_total_mass() - accretor_mass;
	if (mdot>mdot_max) 
	mdot = mdot_max;
    }

    int p = cerr.precision(HIGH_PRECISION);
    PRC(accretor_mass);PRL(get_total_mass());
    cerr.precision(p);
  
    return mdot;
}

// Computes expansion of acceptor due to mass accretion.
// At moment accretor fills own Roche-lobe mass transfer becomes
// inconservative.
real single_star::accretion_limit(const real mdot, const real dt) {

  if (dt < 0) return mdot;
  // Conservative mass transfer.
  //   return mdot;

  // Non-conservative mass transfer.
  // Based on Pols & Marinus,1994, A&A,288, 475
  real r_rl = effective_radius;
  if (is_binary_component())
    r_rl = get_binary()->roche_radius(this);
  real mdot_kh = dt*relative_mass/kelvin_helmholds_timescale();
  real accretion = log10(r_rl/effective_radius)
                 / pow(10, expansionB(relative_mass));

  accretion = max(accretion, 0.);
  real mdot_max = mdot_kh*pow(accretion, 1./expansionA(relative_mass));
  mdot_max = max(mdot_max, 0.);	
  return min(mdot, mdot_max);


  // (SPZ+GN: 26 Jul 2000) Test Kelvin Helmholtz accretion
  // (GN Jul 28 1999) test non conservative Algol evolution
  // return min(mdot, mdot_kh);

}


real single_star::accretion_limit_eddington(const real mdot, const real dt) {


  if (dt < 0) return mdot;
  real eddington = 1.5e-08*cnsts.parameters(solar_radius)*radius*dt;

  if(mdot>=eddington) return eddington;

  return mdot;
}


void single_star::adjust_donor_radius(const real delta_m) {


  real zeta_th = zeta_thermal();

  // -1 because delta_m is defined positive and star loses mass.
  real delta_r = -1 * effective_radius * zeta_th * delta_m/relative_mass;

  effective_radius += delta_r;
  
  // Effective radius return to equilibrium radius when mass transfer
  // is on nuclear timescale
  // (SPZ+GN:30 Sep 1998)
  // or the aml timescale
  // (SilT: 22 Apr 2012)

  if (is_binary_component() &&
      (get_binary()->get_current_mass_transfer_type() == Nuclear | 
      get_binary()->get_current_mass_transfer_type() == AML_driven)) {
    effective_radius = radius;
  }


}

// Adding mass to a star causes it to expand.
void single_star::adjust_accretor_radius(const real mdot, const real dt) {
// Based on Pols & Marinus,1994, A&A,288, 475

  // function returns directly: effective radius is used to determine RLOF.
  // Allowing accretor to grow in size results in contact systems.
  // do not allow bloating
  // (SPZ+GN:28 Sep 1998)
//(GN+SPZ Apr 28 1999) star do bloat however... bloating on again
//  return;
  

  //cerr<<"void star::adjust_accretor_radius()"<<endl;
  if (mdot>0) {
    real mdot_kh = relative_mass*dt/kelvin_helmholds_timescale();

    real r_fr = expansionA(relative_mass)*log10(mdot/mdot_kh)
      + expansionB(relative_mass);
    if (r_fr>0.5) // (GN+SPZ Apr 28 1999) was 1 
      r_fr = 0.5;
// (GN+SPZ Apr 28 1999) radius is equilibrium radius
//    effective_radaius = radius = radius*pow(10, pow(10, r_fr));

    
    if (is_binary_component())
      effective_radius= min(effective_radius, get_binary()->roche_radius(this));
    effective_radius = max(effective_radius, radius*pow(10, pow(10, r_fr)));
  
  } 
}


real single_star::expansionA(const real m) {

  // Lineair interpolation.
  real rc, inc;
  if (m<=3) {
    rc = 0.120;  	//base on: 	m1 = 2; A1 = 0.599;
    inc = 0.359; 	//	  	m2 = 3; A2 = 0.719;
  }
  else if (m<=5) {
    rc = 0.144;	        //based on:	m1 = 3; A1 = 0.719;
    inc = 0.289;	//		m2 = 5; A2 = 1.006;
  }
  else if (m<=10) {
    rc = 0.123;	        //based on:	m1 = 5; A1 = 1.006;
    inc = 0.393;	//		m2 = 10; A2 = 1.619;
  }
  else {
    rc = 0.085;	        //based on:	m1 = 10; A1 = 1.619;
    inc = 0.772;	//		m2 = 17; A2 = 2.212;
  }
     
  return inc + rc*m;
}

real single_star::expansionB(const real m) {

//              Lineair interpolation.
  real rc, inc;
  if (m<=3) {
    rc = -0.273;     //base on:      m1 = 2; B1 = -1.374;
    inc = -0.828;    //              m2 = 3; B2 = -1.647;
  }
  else if (m<=5) {
    rc = -0.192;     //based on:     m1 = 3; B1 = -1.647;
    inc = -1.071;    //              m2 = 5; B2 = -2.030;
  }
  else if (m<=10) {
    rc = -0.106;    //based on:     m1 = 5; B1 = -2.030;
    inc = -1.50;    //              m2 = 10; B2 = -2.560;
  }
  else {
    rc = -7.08e-2;   //based on:     m1 = 10; B1 = -2.560;
    inc = -1.852;    //              m2 = 17; B2 = -3.056;
  }

  real value = inc + rc*m;

  return min(1., value);
}

// Merges cores and envelopes of two stars.
// Star is updated.
star* single_star::merge_elements(star* str) {

  cerr << "Merge element single star" << endl;

  star* merged_star = this;

  real m_conserved = get_total_mass() + str->get_total_mass();


  if (str->get_element_type()!=Main_Sequence) {
    // adding two cores of giants together should not result in
    // rejuvenation.
    // previous method appeared to make mergers too old.
    // (SPZ+GN:11 Oct 1998)
    //add_mass_to_core(str->get_core_mass());
    // (GN Oct 25 2010) proper adding of core mass as non-hydrogen via add_mass_to_accretor
    add_mass_to_accretor(str->get_core_mass(), false);

      if ((str->get_element_type()==Neutron_Star ||
	     str->get_element_type()==Black_Hole)   &&
	     core_mass < cnsts.parameters(helium2neutron_star)) {
          real dm = cnsts.parameters(helium2neutron_star) - core_mass;
	  add_mass_to_accretor(dm, false);

          if (envelope_mass<dm){
	    // (GN Oct 26 2010) mass not conserved in this case!
	        m_conserved += dm-envelope_mass;
                envelope_mass = 0;
          }
          else{
                envelope_mass -= dm;
          }
      }

            //		What to do here is put in SPH!
      if (str->get_envelope_mass()>1.e-11) {

            add_mass_to_accretor(0.5*str->get_envelope_mass(), str->hydrogen_envelope_star());
      }

    }
    else {
        add_mass_to_accretor(str->get_total_mass(), str->hydrogen_envelope_star());
    }

    if (get_total_mass()-m_conserved > 1.e-11) {
        cerr.precision(HIGH_PRECISION);
        cerr << "ERROR: Mass is not conserved in single_star::merge elements()"
        << endl;

        PRC(get_total_mass());PRC(m_conserved);
        PRL(get_total_mass()-m_conserved);

        exit(-1);
    }

    spec_type[Merger]=Merger;
    adjust_next_update_age();
    


  // Redundant: in core mass addition star is no rejuvenated but aged.
  // addition of envelope mass causes rejuvenation.
  // (SPZ+GN:27 Sep 1998)
  //if (str->get_element_type()!=Main_Sequence) {
  //real aged = 0.01*(core_mass/agb_core_mass(relative_mass));
  //relative_age = (1+aged)*next_update_age;
  //}


   if (get_relative_mass() >= cnsts.parameters(massive_star_mass_limit) &&
       hydrogen_envelope_star()){ 
       merged_star =  reduce_mass(envelope_mass);
   }
   else {
       instantaneous_element();
       // calling evolve element may cause segmentation fault if star changes.
      //evolve_element(current_time);
   }

    cerr << "Merged star: "<<endl;
    merged_star->dump(cerr, false);

    return merged_star;
}



// Determine mass transfer stability according to
// `zeta' discription.
// New version // (SilT + GN) 10 Feb 2011 z_l depends on timescale  
// calculate two zeta_l's and compare to appropriate zeta_star
// use timescale of type of mass transfer to check if this specific kind of mass transfer is stable 
//    used to be one z_l (with md_timescale = previous md_timescale);
real single_star::mass_transfer_timescale(mass_transfer_type &type) {
    
  type = Unknown;

  real z_l_th = 0, z_l_ad = 0, z_star = 0;

  if (is_binary_component()) {

    z_l_ad = get_binary()->zeta(this, get_companion(), kelvin_helmholds_timescale());
    z_l_th = get_binary()->zeta(this, get_companion(), nucleair_evolution_timescale());
  }

  real z_ad = zeta_adiabatic();
  real z_th = zeta_thermal();
//  PRC(z_ad);PRC(z_th);PRC(z_l_th);PRL(z_l_ad);

  real mtt;
  if (z_ad > z_l_th && z_th > z_l_th) {

    mtt = nucleair_evolution_timescale();
    type = Nuclear;
    z_star = z_th;
  }
  else if (z_ad >= z_l_ad) {

    mtt = kelvin_helmholds_timescale();
    type = Thermal;
    z_star = z_ad;
  }
  else if (z_l_ad > z_ad) {

    mtt = sqrt(kelvin_helmholds_timescale()*dynamic_timescale());
    type = Dynamic;
    z_star = z_ad;    
  }
  else {
    cerr << "No clear indication for mass transfer timescale: "
	 << "Kelvin-Helmholds time-scale assumed."<<endl;

    mtt = kelvin_helmholds_timescale();
    type = Unknown;
    z_star = z_ad;    
  }

  if (low_mass_star()) {

    real mdot = 0;
    if (is_binary_component())
      mdot = get_binary()
              ->mdot_according_to_roche_radius_change(this, get_companion(), z_star);
    if (mdot > 0) {

      real mtt_rl = get_relative_mass()/mdot;
      if(mtt>mtt_rl) {

	mtt = mtt_rl;
	type = AML_driven;
      }      
    }
  }

  if (REPORT_MASS_TRANSFER_TIMESCALE) {
    cerr << "single_star::mass_transfer_timescale()" << endl;
    cerr << "    star id = " << identity
	 << "  Zeta (lobe_ad, lobe_th, ad, th) = ("
	 << z_l_ad <<", "<< z_l_th<<", "<<z_ad<<", "<<z_th<<") : " << endl;
    cerr << type_string(type);
    cerr<<":    dm/dt=" <<get_relative_mass()/(mtt*1.e+6)
	<< " [Msun/yr]" << endl;
	
  }

  return mtt;
}


//		Stellar stability functions.
real single_star::zeta_adiabatic() {
// (GN+SPZ Apr 28 1999) this is used by sub_giant: all stars with HW87
// all stars should have their own actually...(when we have time)

// (GN+SPZ Apr 28 1999) fit from Lev Yungelson private communication
// for giants with not too convective envelope = radiative envelope

  real r_dconv = 2.4*pow(relative_mass,1.56);
  if (relative_mass > 10 )
    r_dconv = 5.24*pow(relative_mass,1.32);
  else if (relative_mass > 5)
    r_dconv = 1.33*pow(relative_mass,1.93);
    
  if (radius < r_dconv) {
    return 12.25;
//    cerr << "radius < r_dconv" << endl;
  }
  else {
//   		Hjellming and Webbink 1987 ApJ, 318, 804
    real x = core_mass/get_total_mass();
    real A = -0.220823;
    real B = -2.84699;
    real C = 32.0344;
    real D = -75.6863;
    real E = 57.8109;

// (GN+SPZ Apr 28 1999) not for (sub) giants    
//  if (low_mass_star())
//    z = -cnsts.mathematics(one_third);
//  else 
    return A + x*(B + x*(C + x*(D + x*E)));

  }

}

// Values of zeta are changed (SPZ+GN:28 Sep 1998)
real single_star::zeta_thermal() {
  //cerr<<"real single_star::zeta_thermal()"<<endl;

  real z;
  if (low_mass_star())
    z = 0;
  else 
    z = 0; // (GN+SPZ Apr 28 1999) radius determined by core only (was -1) 

  return z;
}

// not used anymore, use add_mass_to_accretor instead
// add two cores. Is performed in ::merge_elements();
//void single_star::add_mass_to_core(const real mdot) {    
//  if (mdot<=0) {
//    cerr << "single_star::add_mass_to_core(mdot=" << mdot << ")"<<endl;
//    cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
//    
//  }
//  else {
//      cerr<<"adjust_accretor_age bool = FALSE, why exactly?"<<endl;
//    adjust_accretor_age(mdot, false);
//    core_mass += mdot;
//    accreted_mass += mdot;
//    
//    if (relative_mass<get_total_mass()) 
//      update_relative_mass(get_total_mass());
//    
//  }
//
//}



real single_star::add_mass_to_accretor(real mdot, bool hydrogen, const real dt) {
    cerr<<"currently single_star::add_mass_to_accretor should not be used"<<endl;
    PRL(dt);

    exit(-1);
    
  if (mdot<0) {
    cerr << "single_star::add_mass_to_accretor(mdot=" << mdot 
	 << ", dt=" << dt << ")"<<endl;
    cerr << "mdot (" << mdot << ") smaller than zero!" << endl;
    return 0;
  }

    
    if(hydrogen){
      //hydrogen accretion
      mdot = accretion_limit(mdot, dt);

      // For now, no rejuvenation of SG, CHeB, AGB or He giant accretor   
      // adjust_accretor_age(mdot);
      envelope_mass += mdot;
      accreted_mass += mdot;

    //  if (relative_mass<get_total_mass()) 
    //    update_relative_mass(get_total_mass());

      adjust_accretor_radius(mdot, dt);
      
      set_spec_type(Accreting);
    }
    else{
        //for the moment assume helium accretion
        cerr<<"helium accretion"<<endl;        
    }
  return mdot;
}

// Post-evolution stellar update.
// Makes sure age and radius are updated.
void single_star::update() {

  // New core mass determination occurs in ::evolve_element.
  // (SPZ+GN:09/1998)
  // real m_tot = get_total_mass();
  // core_mass = helium_core_mass();
  // envelope_mass = m_tot - core_mass;

  core_radius = helium_core_radius();

  // (GN+SPZ Apr 28 1999)
  // effective_radius can be larger than  radius
  // (SPZ:  8 Jul 2001) 
  // except for horizontal branch stars
  effective_radius = max(radius,effective_radius);


  // last update age is set after stellar expansion timescale is set.
// (GN+SPZ May  4 1999) last_update_age now used as time of last type change
//  last_update_age = relative_age;

  detect_spectral_features();
  //dump(cerr, false);

//    real (single_star::*fptr)() = &single_star::nucleair_evolution_timescale;        
//    real result = (this->*fptr)();
//    real test = linear_function_inversion(fptr, relative_mass, core_mass);     
}


void single_star::detect_spectral_features() {

  //		Clean old spectral specifications.
  set_spec_type(Emission, false);
  set_spec_type(Blue_Straggler, false);
  set_spec_type(Barium, false);

  // set new
  if (is_binary_component())
    if (!remnant() && 
	!get_companion()->hydrogen_envelope_star() &&
	accreted_mass >= cnsts.parameters(Barium_star_mass_limit)
	               * relative_mass)
      set_spec_type(Barium);
}

// In the case of a binary, the companion star might accrete from a
// stellar wind or post-AGB mass-shell.
// Bondi, H., and Hoyle, F., 1944, MNRAS 104, 273 (wind accretion.
// Livio, M., Warner, B., 1984, The Observatory 104, 152.
real single_star::accrete_from_stellar_wind(const real mdot, const real dt) {

  real alpha_wind = 0.5;
  real v_wind = get_companion()->wind_velocity();

  real acc_radius = pow(cnsts.physics(G)*cnsts.parameters(solar_mass)
			* get_total_mass()
			/ pow(v_wind*cnsts.physics(km_per_s), 2),2)
                  / cnsts.parameters(solar_radius);
  real wind_acc = alpha_wind/(sqrt(1-pow(get_binary()->get_eccentricity(), 2))
			      * pow(cnsts.parameters(solar_radius)
				    * get_binary()->get_semi(),2));
  real v_factor = 1/pow(1+pow(velocity/v_wind, 2), 3./2.);

  real mass_fraction = acc_radius*wind_acc*v_factor;
// (GN+SPZ May  4 1999) Do not multiply with dt*cnsts.physics(Myear)!

//  PRC(v_wind);PRC(acc_radius);PRC(wind_acc);PRL(v_factor);
//  PRL(mass_fraction);PRL(mdot);

  mass_fraction = min(0.9, mass_fraction);
  return add_mass_to_accretor(mass_fraction*mdot, get_companion()->hydrogen_envelope_star(), dt);
}

// Tidal energy dissipation during close encounter.
// Portegies Zwart, SF & Meinen AT, 1992, AA 280, 174.
// for polytrope n=1.5
real single_star::tf2_energy_diss(const real eta) {

  const real coeff2[] = {-0.397, 1.678, 1.277, -12.42, 9.446, -5.550};

  real y = log10(eta);
  real logT = ((((coeff2[5]*y + coeff2[4])*y + coeff2[3])*y 
		+ coeff2[2])*y + coeff2[1])*y + coeff2[0];

  return pow(10., logT);
}

// Tidal energy dissipation during close encounter.
// Portegies Zwart, SF & Meinen AT, 1992, AA 280, 174.
// for polytrope n=1.5
real single_star::tf3_energy_diss(const real eta) {

  const real coeff3[] = {-0.909, 1.574, 12.37, -57.40, 80.10, -46.43};

  real y = log10(eta);
  real logT = ((((coeff3[5]*y + coeff3[4])*y + coeff3[3])*y  
		+ coeff3[2])*y + coeff3[1])*y + coeff3[0];

  return pow(10., logT);
}

//	pretty-print
void single_star::print_status() {
  
  star_state str;
  str.init_star_state((single_star*)this);
  str.make_star_state((single_star*)this);
  
  cout << " [M = " << get_total_mass() 
       << ", R = " << effective_radius 
       << ", v = " << velocity << "] "
       << type_string(get_element_type()) << " ";
  str.put_star_state();

}

//	print data for EPJ Roche program.
void single_star::print_roche() {

  real r = effective_radius;
  if (effective_radius<radius) r = 1.e6;

  cerr << get_total_mass() << " " << r << " ";

}

void  single_star::set_spec_type(star_type_spec s,
				 bool on) { 	// default =true

  if (on) spec_type[s] = s;
  else    spec_type[s] = NAC;
}


// Angular momentum of homogeneous sphere.
real single_star::angular_momentum() {
       
  real m = get_total_mass()*cnsts.parameters(solar_mass);
  real r = effective_radius*cnsts.parameters(solar_radius);

// (GN+SPZ May  5 1999) effective_radius may increase when rl shrinks
  if (is_binary_component()) {
    r = min(r, 
	    get_binary()->roche_radius(this)*cnsts.parameters(solar_radius));
  }

  real o = 0;                            // default rotation.
  if(rotation_period>0)
    o = 2*PI/rotation_period;
  else if (is_binary_component()) 
    o = 2*PI/(get_binary()->get_period()*cnsts.physics(days));

  return gyration_radius_sq()*m*o*pow(r, 2);
	
}

real single_star::rejuvenation_fraction(const real mdot_fr) {

      real rejuvenation = (1-pow(mdot_fr,
				 cnsts.parameters(rejuvenation_exponent)));

//      // no rejuvenation if companion has no hydrogen envelope.
//      if(is_binary_component() &&
//         !get_companion()->hydrogen_envelope_star()) {
//          cerr<<"rejuvenation_fraction without hydrogen accretion should not happen anymore"<<endl;
//          cerr<<"maybe for helium stars & giants"<<endl;
//        //rejuvenation = 1;
//      }
      return rejuvenation;
}

void single_star::stellar_wind(const real dt) {
#if 0
        // (GN+SPZ Apr 28 1999) wind for low mass stars per phase
    real end_time = next_update_age - last_update_age;
//    real prev_rel_time = max(0.,previous.relative_age - last_update_age);
//    real relative_time = min(relative_age - last_update_age, end_time);
    real relative_time = relative_age - last_update_age;

// for high mass stars over whole evolution
// (GN May 12 1999)
// except for stars more massive than 85 that become WR after ms immediately
    if (relative_mass >= cnsts.parameters(super_giant2neutron_star) &&
	relative_mass < 85.) {
      end_time = nucleair_evolution_time();
      relative_time = relative_age;
    }

    real wind_mass = wind_constant 
                   * (pow(relative_time/end_time,
			cnsts.parameters(massive_star_mass_loss_law))
	           -  pow((relative_time-dt)/end_time,
			cnsts.parameters(massive_star_mass_loss_law)));

    // Previous second term according to GN.
    //	           -  pow((prev_rel_time)/end_time,
    //			cnsts.parameters(massive_star_mass_loss_law)));
#endif
    
    // (ST 17 Sep 2009) 
    // wind_constant in solar masses per year
    update_wind_constant();
    real wind_mass = wind_constant * dt * 1E6;  
    
    //PRC(wind_constant);PRC(dt);PRL(wind_mass);

    if (wind_mass > 0){    
        // (GN+SPZ Apr 28 1999) wind induced helium star formation can happen
        // because stellar_wind is last function in evolve_element
        if (wind_mass>=envelope_mass) {
            wind_mass = envelope_mass;
            radius = core_radius;
        }
      
        if (is_binary_component())
            get_binary()->adjust_binary_after_wind_loss(this, wind_mass, dt);
        else {
            // (GN Oct 11 1999) for single stars: previous used for stellar wind! (?)
            // refresh_memory();
            reduce_mass(wind_mass);
        }
    }
//    else 
//        cerr<<"No wind mass loss"<<endl;

    return;
}

void single_star::update_relative_mass(const real new_relative_mass) {
  relative_mass = new_relative_mass;
  adjust_next_update_age();
  update_wind_constant();

}

void single_star::lose_envelope_decent() {

  //  cerr << "single_star::lose_envelope_decent" << endl;
  //  dump(cerr, false);
  if (envelope_mass>0) {
    if (is_binary_component()) {
      get_binary()->adjust_binary_after_wind_loss(
		    this, envelope_mass, POST_AGB_TIME);
    } 
    else 
      reduce_mass(envelope_mass);
  }

  if (is_binary_component()) {
    get_companion()->set_effective_radius(get_companion()->get_radius());
    get_companion()->set_spec_type(Accreting, false);
  }
}

// wind_constant is fraction of envelope lost in lifetime
// of stars. Should be updated after mass accretion
// (SPZ+GN: 3 Oct 1998)
void single_star::update_wind_constant() {
  
#if 0
  if (relative_mass >= cnsts.parameters(massive_star_mass_limit))
    wind_constant = (get_relative_mass()-core_mass)
                  * cnsts.parameters(massive_star_envelope_fraction_lost);
  else 
    wind_constant = get_relative_mass()
                  * cnsts.parameters(non_massive_star_envelope_fraction_lost);
    
#endif
#if 0 
    // (GN+SPZ Apr 28 1999) new fits to Maeder, de Koter and common sense
//  cerr << "update_wind_constant"<<endl;

  if (relative_mass >= cnsts.parameters(super_giant2neutron_star)) {

    real meader_fit_dm = 0.01*pow(relative_mass,2.);
    
    if (relative_mass < 85)
      wind_cons tant = meader_fit_dm;
    else {// constant
      real final_mass = 43; // final mass after ms
      wind_constant = relative_mass - final_mass;
    }

  } 
  else { // no wind for low mass ms stars
    wind_constant = 0;
  }

  wind_constant = max(wind_constant, 0.0);
#endif
    
    //single_star::update_wind_constant is not used
    // only by WD, NS, BH who don't have stellar_wind
    wind_constant = 0.0;
}


real single_star::potential_energy() {
  
     real GM2_R = cnsts.physics(G)*pow(cnsts.parameters(solar_mass), 2)
                / cnsts.parameters(solar_radius);
     real p = GM2_R*get_total_mass()*get_total_mass()
            / get_effective_radius();
     
     return -p;
}

real single_star::kinetic_energy() {
  
     real Mkm_s2 = cnsts.parameters(solar_mass)
                 * pow(cnsts.physics(km_per_s), 2);
     real k = 0.5*Mkm_s2*get_total_mass()*pow(velocity, 2);
     
     return k;
}

real single_star::total_energy() {
     return kinetic_energy() + potential_energy();
}

real single_star::get_evolve_timestep() {
// (GN+SPZ Apr 28 1999) was a bit too small
//  return max(next_update_age - relative_age
//	     -0.5*cnsts.safety(minimum_timestep),
//	     cnsts.safety(minimum_timestep));

// (GN+SPZ May  5 1999) type change time must be small because of rapid
// growth of giants at end phase 0.0001 seems to be OK (?)
//  return max(next_update_age - relative_age - (0.5*0.001), 0.001);


// (GN + SilT Nov 23 2009) go to end of phase in a maximum of 2 steps
// in stead of going to 90% of phase until the timestep is smaller than minimum_timestep
//  return max(next_update_age - relative_age, 0.0001);

  real timestep = min((next_update_age - last_update_age )/ cnsts.safety(number_of_steps), 
                      next_update_age - relative_age - 0.5 * cnsts.safety(minimum_timestep));   

    //temper LBV massloss rate
//    real timestep_lbv = timestep;
//    real x_lbv = 1.0E-5*radius*sqrt(luminosity);
//    if(hydrogen_envelope_star() && luminosity > 6.0E5 && x_lbv > 1.0){// dm 1% of m_env
//        timestep_lbv = 0.1* envelope_mass *pow(x_lbv -1.0, -3.0) / (luminosity/6.0E5 -1.0) /1.0E6;
//    }
    
//   timestep = min(timestep, timestep_lbv);             

    //temper WR massloss rate
//    real timestep_wr = timestep;
//    real mu = (get_total_mass()-core_mass)/get_total_mass() * min(5.0,max(1.2, pow(luminosity/7.E4,-0.5)));
//    if ( mu < 1.){ //dm 1% of m_env
//        timestep_wr = envelope_mass / (1.-mu) / 1.38 * pow(get_total_mass(), -2.87);
//    }
//    timestep = min(timestep, timestep_wr);             
    
   return max(timestep, cnsts.safety(minimum_timestep));
    
}



// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// New metalicity dependencies from thesis of Hurley, J., 2000
// starts here.
//

// This zeta is the metalicity parameter and has no relation to the 
// zeta used to determine the mass transfer timescale
real single_star::get_zeta(const real z) {

  real zeta = log10(z/cnsts.parameters(solar_metalicity));
  return zeta;
}

// Eq.2
real single_star::helium_flash_mass(const real z) {

  real zeta = get_zeta(z);
  real m_HeF = 1.995 + zeta*(0.087*zeta + 0.25);

  return m_HeF;
}

// Eq.3
real single_star::helium_ignition_mass(const real z) {

  real zsun = z/cnsts.parameters(solar_metalicity);
  real m_FGB = 13.048 * pow(zsun, 0.06)
            / (1 + 0.0012 * pow(1./zsun, 1.27));

  return m_FGB;
}

// Eq.49
real single_star::helium_ignition_luminosity(const real mass, 
					     const real z) {
  real b10 = smc.b(10,z);

  real l_HeI; 
  real m_HeF = helium_flash_mass(z);
  if(mass < m_HeF) {
    real l_HeF = (smc.b(11, z) + smc.b(12, z)*pow(m_HeF, 3.8))
             / (smc.b(13, z) + pow(m_HeF, 2));
    real a = (smc.b(9, z)*pow(m_HeF, b10) - l_HeF)/l_HeF;
    l_HeI = smc.b(9, z)*pow(mass, smc.b(10, z))
          / (1 + a*exp(15*(mass - m_HeF)));
  }
  else {

      l_HeI = (smc.b(11, z) + smc.b(12, z)*pow(mass, 3.8))
             / (smc.b(13, z) + pow(mass, 2));
  }

  return l_HeI;
}

// Eq.50
real single_star::helium_ignition_radius(const real mass, const real mass_tot, const real z) {
    // ST: watch out when changing this function, because
    // for mass>max(12.0, m_FGB) r_mhe = min(r_mhe, r_agb(L_HeI))
    // for m_FGB<mass<12.0 unclear what we should do with r_agb(L_HeI)
    // maybe r_agb(L_HeI)> r_mHe for normal metallicities
    
  real r_hi;
  real m_FGB = helium_ignition_mass(z);
  real r_mHe = minimum_blue_loop_radius(mass, mass_tot, z);
  if (mass<=m_FGB) { 
    real l_HeI = helium_ignition_luminosity(mass, z);
    r_hi = giant_branch_radius(l_HeI, mass_tot, z);
  }
  else if(mass>=max(m_FGB,12.)) {
    r_hi = r_mHe;
  }
  else {
    real mu = log10(mass/12)/log10(m_FGB/12);
    real l_HeI = helium_ignition_luminosity(mass, z);
    real r_gb = giant_branch_radius(l_HeI, mass_tot, z);
    r_hi = r_mHe*pow(r_gb/r_mHe, mu);
  }

  return r_hi;
}

bool single_star::low_mass_star(const real mass,
				const real z) {

    bool is_low_mass_star = false;

    if (!remnant()) {
        if(relative_mass <= cnsts.parameters(low_mass_star_mass_limit)) 
            is_low_mass_star = true;
    }
    else if (get_total_mass() <=
             cnsts.parameters(low_mass_star_mass_limit)) {
        is_low_mass_star = true;
    }
    
    return is_low_mass_star;
}
 
bool single_star::low_mass_star() {

  return low_mass_star(relative_mass, metalicity);
}

bool single_star::intermediate_mass_star(const real mass, 
					 const real z) {

  return (!low_mass_star(mass, z) && 
	  !high_mass_star(mass, z))?true:false;
}

bool single_star::medium_mass_star() {

  return intermediate_mass_star(relative_mass, metalicity);
}

bool single_star::high_mass_star(const real mass,
				 const real z) {

    if(remnant())
        return (get_total_mass()>cnsts.parameters(medium_mass_star_mass_limit))
        ?true :false;
    else
        return (get_relative_mass()>cnsts.parameters(medium_mass_star_mass_limit))
        ?true :false;
    
}

bool single_star::high_mass_star() {

  return high_mass_star(relative_mass, metalicity);
}


// Eq.4
// supersedes real single_star::base_giant_branch_time(const real mass,
//                                                     const real t_ms);
real single_star::base_giant_branch_time(const real mass,
					 const real z) {

    real t_bgb;
    real pow_mass_7 = pow(mass, 7);
    real teller = smc.a(1, z) +smc.a(2, z)*pow(mass, 4) 
                      +smc.a(3, z)*pow(mass, 5.5)  
                      +       pow_mass_7;
    real noemer =  smc.a(4, z)*pow(mass, 2) +smc.a(5, z)*pow_mass_7; 
    t_bgb = teller/noemer;
  
    return t_bgb;
}

//real single_star::base_giant_branch_time(const real mass) {
//
//  return base_giant_branch_time(mass, metalicity);
//}

real single_star::base_giant_branch_time() {

  return base_giant_branch_time(relative_mass, metalicity);
}

// This function uses internal parameters explicitely!
real single_star::convective_envelope_mass(const real z) {

  real m_ce;
  if(low_mass_star(relative_mass, z)) {
    m_ce = 0.4 * envelope_mass;
  }
  else {
    m_ce = cnsts.mathematics(one_third) * envelope_mass;
  }

  return m_ce;
}


// Eq.10
real single_star::base_giant_branch_luminosity(const real mass, 
                                               const real z) {

  real c2 = 9.301992;
  real c3 = 4.637345;
  real a29 = smc.a(29, z);
  real l_bgb = (smc.a(27, z)*pow(mass, smc.a(31, z)) 
		+ smc.a(28, z)*pow(mass, c2))
             / (a29 + smc.a(30, z)*pow(mass, c3) 
	                     + pow(mass, smc.a(32, z)));

  return l_bgb;
}

real single_star::base_giant_branch_luminosity(const real mass) { 
  
 return base_giant_branch_luminosity(mass, metalicity);
}

real single_star::base_giant_branch_luminosity() {

  return base_giant_branch_luminosity(relative_mass, metalicity);
}

// Eq.3
real single_star::FGB_mass(const real z) {

  real zsun = z/cnsts.parameters(solar_metalicity);
  real m_FGB = 13.048 * pow(zsun, 0.06)
            / (1 + 0.0012 * pow(1./zsun, 1.27));

  return m_FGB;
}

real single_star::get_hydrogen_fraction(const real z) {

  real X = 0.76 - 3*z;
  return X;
}




//++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//Main_Sequence:


// Eq.1
real single_star::main_sequence_hook_mass(const real z) {

  real zeta = get_zeta(z);
  real m_hook = 1.0185 + zeta * (zeta * 0.0892 + 0.16015);

  return m_hook;
}

//Eq. 5
real single_star::main_sequence_time(const real mass, const real z) {

  real t_ms = base_giant_branch_time(mass, z)
            * max(stars_with_main_sequence_hook(mass, z),
                  stars_without_main_sequence_hook(mass, z)); 
 
  return t_ms;
}

real single_star::main_sequence_time() {

  return main_sequence_time(relative_mass, metalicity);
}

real single_star::main_sequence_hook_time(const real mass,
					    const real z) { 

  real mu = stars_with_main_sequence_hook(mass, z);
  real t_bgb = base_giant_branch_time(mass, z);
  real t_hook = mu * t_bgb;

  return t_hook;
}

// Eq.6 identified as 'x' by Hurley
real single_star::stars_without_main_sequence_hook(const real mass, 
						    const real z) {
  
  real zeta = get_zeta(z);
  real x = max(0.95, min(0.99,
                         0.95 - 0.03*(zeta + 0.30103)));
  
  return x;
} 

// Eq.7 identified as 'mu' by Hurley
real single_star::stars_with_main_sequence_hook(const real mass, 
						const real z) {
  
  real mu = max(0.5, 
		1.0 - 0.01*max(smc.a(6, z)/pow(mass, smc.a(7, z)),
			       smc.a(8, z) + smc.a(9, z)
			       /pow(mass, smc.a(10, z))));

  return mu;
}


// Eq.8
real single_star::terminal_main_sequence_luminosity(const real mass, 
                                                    const real z) {

  real teller = pow(mass, 3)*(smc.a(11, z) + mass*smc.a(12, z)) 
                                     +smc.a(13, z)*pow(mass,smc.a(16, z)+1.8);
  real noemer = smc.a(14, z) + smc.a(15, z)*pow(mass, 5) + pow(mass, smc.a(16, z));

  real l_tms = teller/noemer;

  return l_tms;
}

// Eq. 9a and Eq.9b
real single_star::terminal_main_sequence_radius(const real mass, 
                                                const real z) {

  real a17 = smc.a(17, z);
  real c1 = -0.08672073; 
  real r_tams;
  if (mass<=a17) { 
     r_tams = (smc.a(18, z) + smc.a(19, z)*pow(mass, smc.a(21, z)))
            / (smc.a(20, z) + pow(mass, smc.a(22, z)));
  }
  else if (mass>=a17+0.1) { 
     r_tams = (c1*pow(mass, 3) + smc.a(23, z)*pow(mass, smc.a(26, z))
                                  + smc.a(24, z)*pow(mass, smc.a(26, z)+1.5))
            / (smc.a(25, z) + pow(mass, 5));
  }
  else {
    real r1_tms = terminal_main_sequence_radius(a17, z); 
    real r2_tms = terminal_main_sequence_radius(a17+0.1, z); 
     r_tams = lineair_interpolation(mass, a17, a17+0.1, r1_tms, r2_tms);
  }

  if (mass<0.5) // not depending on metalicity
    r_tams = max(r_tams, 1.5*base_main_sequence_radius(mass, z));

  return r_tams;
}


// Tout 1996 fits
real single_star::base_main_sequence_radius(const real mass, const real z) { 
    real mx = pow(mass, 0.5);
    real teller = (smc.c(8,z)*pow(mass,2) + smc.c(9,z)*pow(mass,6))*mx + smc.c(10,z)*pow(mass,11) +(smc.c(11,z) + 
                                  smc.c(12,z)*mx)*pow(mass,19);
    real noemer = smc.c(13,z) + smc.c(14,z)*pow(mass,2) + (smc.c(15,z)*pow(mass,8) + pow(mass,18) + 
                                                smc.c(16,z)*pow(mass,19))*mx;
    return teller/noemer;
}




//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Hertzsprung gap stars
//

// based on Eq. 28
real single_star::initial_hertzsprung_gap_core_mass(const real mass, const real z){
    real m5_25 = pow(mass, 5.25);
    real mc_ehg = terminal_hertzsprung_gap_core_mass(mass, z);
    
    real rho = (1.586 + m5_25) / (2.434 + 1.02*m5_25);
        
    return mc_ehg*rho; 
}


//Eq.28
real single_star::terminal_hertzsprung_gap_core_mass(const real mass, 
                                                         const real z) {
    
    real m_core;
    real m_HeF = helium_flash_mass(z);
    real m_FGB = helium_ignition_mass(z);
    if (mass < m_HeF) {
        real l_bgb = base_giant_branch_luminosity(mass, z);
        m_core = FGB_core_mass_luminosity_relation(l_bgb, mass, z);
    }
    else if (mass >= m_FGB) {
        real mc_HeI = helium_ignition_core_mass(mass, z);//sect.5.3: Eq.67
        m_core = mc_HeI;
    }
    else {
        real mc_bgb = base_giant_branch_core_mass(mass, z); //Eq.44    
        m_core = mc_bgb;
        
    }
    return m_core;
}

//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Sub giant branch
//

// Eq.46
real single_star::giant_branch_radius(const real l_hb,
				      const real mass_tot, 
				      const real z) {

  real b1 = smc.b(1,z);
  real b2 = smc.b(2,z);
  real b3 = smc.b(3,z);
  real b4 = smc.b(4,z);
  real b5 = smc.b(5,z);
  real b6 = smc.b(6,z);
  real b7 = smc.b(7,z);

  real A = min(b4/pow(mass_tot, b5), b6/pow(mass_tot, b7)); 

  real r_hb = A * (pow(l_hb, b1) + b2 * pow(l_hb, b3));   

  return r_hb;
}


// Eq.51
real single_star::minimum_horizontal_branch_luminosity(const real mass, 
						       const real z) {
  real m_FGB = helium_ignition_mass(z);
  real c = smc.b(17, z)/pow(m_FGB, 0.1) 
         + (smc.b(16, z)*smc.b(17, z) - smc.b(14, z))
         / pow(m_FGB, smc.b(15, z)+0.1);
  real l_HeI = helium_ignition_luminosity(mass, z);
  real l_min_hb = l_HeI * (smc.b(14,z) + c*pow(mass, smc.b(15,z)+0.1))
    / (smc.b(16,z) + pow(mass, smc.b(15,z)));

  return l_min_hb;
}

// Eq.53
real single_star::base_horizontal_branch_luminosity(const real mass, 
						    const real z) {

  real b18 = smc.b(18,z);
  real b19 = smc.b(19,z);
  real b20 = smc.b(20,z);
  real m_HeF = helium_flash_mass(z);
  // mc represents here the core mass of a star of mass m at the beginning of the 
  // core_helium_burning_fase
  real l_HeI = helium_ignition_luminosity(mass,z);
  real mc = FGB_core_mass_luminosity_relation(l_HeI, mass, z);
  real mu = (mass-mc)/(m_HeF-mc);
  if (mu < 0) mu = 0;
  
  real l_zHe = helium_star_luminosity_for_solar_metalicity(mc);
  real l_min_He = minimum_horizontal_branch_luminosity(m_HeF, z);
  real a = (b18 + l_zHe - l_min_He) / (l_min_He - l_zHe);
  
  real l_bhb = l_zHe + (1+b20)/(1+b20*pow(mu, 1.6479))
             * b18*pow(mu, b19)
             / (1 + a*exp(15*(mass-m_HeF)));
                
  return l_bhb;
}

// Eq.55
real single_star::minimum_blue_loop_radius(const real mass, 
                                           const real mass_tot, const real z) {
    // watch out when changing this function, because
    // for mass>max(12.0, m_FGB) r_mhe = min(r_mhe, r_agb(L_HeI))
    // for m_FGB<mass<12.0 unclear what we should do with r_agb(L_HeI)
    // maybe r_agb(L_HeI)> r_mHe for normal metallicities
    
    real b24 = smc.b(24,z);
    real b25 = smc.b(25,z);
    real b26 = smc.b(26,z);
    real b27 = smc.b(27,z);
    real b28 = smc.b(28,z);
    
    
    real r_mHe = (b24*mass_tot + pow(b25*mass_tot, b26)*pow(mass_tot, b28))
    / (b27 + pow(mass_tot, b28));
    
    real m_HeF = helium_flash_mass(z); 
    if(mass<m_HeF) { 
        real mu = mass_tot/m_HeF;
        real l_bhb = base_horizontal_branch_luminosity(mass, z);
        real r_hb = giant_branch_radius(l_bhb, mass_tot, z);
        real l_HeF = base_horizontal_branch_luminosity(m_HeF, z);
        real r_HeF = giant_branch_radius(l_HeF, m_HeF, z);
        real r_mHe_HeF = (b24*m_HeF + pow(b25*m_HeF, b26)*pow(m_HeF, b28))
        / (b27 + pow(m_HeF, b28)); 
        r_mHe = r_hb * pow(r_mHe_HeF/r_HeF, mu);
    }
    //these lines are not in the HPT2000 article, but they are in the HPT2000 code
    //in case a massive star skips the blue loop phase, 
    // the stellar radius should continue smoothly
    if (mass >= max(12.0, helium_ignition_mass(z))){
        // In this mass range r_x = r_HeI = r_mHe.
        // In case the star skips the blue loop phase, 
        // r_mHe should be the minimum of r_mHe and r_agb.
        real l_HeI = helium_ignition_luminosity(mass, z);
        real r_agb = AGB_radius(l_HeI, mass, mass_tot, z);
        r_mHe = min(r_agb, r_mHe);
    }
    if (mass > helium_ignition_mass(z) && mass <12.0){
        real l_HeI = helium_ignition_luminosity(mass, z);
        real r_agb = AGB_radius(l_HeI, mass, mass_tot, z);
        if (r_agb < r_mHe) {
            PRC(r_agb);PRC(r_mHe);
            cerr<<"WARNING in single_star::minimum_blue_loop_radius: R_AGB(L_HeI) < R_mHe, skipping blue loop?"<<endl;
        }
    }
    
    return r_mHe;
}


// Eq.54 R_ZAHB
real single_star::base_horizontal_branch_radius(const real mass, 
						const real mass_tot, const real z) {
  real b21 = smc.b(21,z);
  real b22 = smc.b(22,z);
  real b23 = smc.b(23,z);
  real m_HeF = helium_flash_mass(z);
  
  // mc represents here the core mass of a star of mass m at the beginning of the 
  // core_helium_burning_fase
  real l_HeI = helium_ignition_luminosity(mass, z);
  real mc = FGB_core_mass_luminosity_relation(l_HeI, mass, z);
    
  real mu = (mass_tot-mc)/(m_HeF-mc);
  real f = (1 + b21)*pow(mu, b22)
         / (1 + b21*pow(mu, b23));
  real l_bhb = base_horizontal_branch_luminosity(mass, z);
  real r_gb = giant_branch_radius(l_bhb, mass_tot, z);
  real r_czahe = helium_star_radius_for_solar_metalicity(mc);
  real r_bhb = (1-f)*r_czahe + f*r_gb;
    
  return r_bhb;
}

// Shell helium burining
// Eq.56
real single_star::base_AGB_luminosity(const real mass, const real z) {

    real b29 = smc.b(29,z);
    real b30 = smc.b(30,z);
    real b31 = smc.b(31,z);
    real b32 = smc.b(32,z);
    real b33 = smc.b(33,z);
    real b34 = smc.b(34,z);
    
    real m_HeF = helium_flash_mass(z);
    real l_bagb;
    if (mass >= m_HeF){ 
        l_bagb = (b31 + b32*pow(mass, b33+1.8))
            / (b34 + pow(mass, b33));
    }
    else {
        real l_bagb_mHeF = (b31 + b32*pow(m_HeF, b33+1.8))
              / (b34 + pow(m_HeF, b33));
        real a = (b29*pow(m_HeF, b30) - l_bagb_mHeF)/l_bagb_mHeF;
        l_bagb = b29*pow(mass, b30) / (1 + a * exp(15*(mass-m_HeF)));
    }
  return l_bagb;
}


//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Core helium burning with envelope (sub giant)
// Hurley et al Sect 5.2 First Giant branch
//

//Eq.33
real single_star::hydrogen_rate_constant(const real mass) { 

  real log_Ah = -4.84;     // Classical value.
  
  log_Ah = max(-4.8, min(-5.7 + 0.8*mass, -4.1 + 0.14*mass));

  real Ah = pow(10., log_Ah);
  return Ah;
}

real single_star::sub_giant_Ah_estimator(const real mass) { 
  
  return hydrogen_rate_constant(mass);
}

real single_star::sub_giant_B_factor(const real mass) { 

  real B = max(30000., 500 + 17500*pow(mass, 0.6));

  return B;
}

real single_star::sub_giant_D_factor(const real mass) { 
  cerr<<"Using sub_giant_D_factor without metalicity as function parameter"<<endl;
  cerr<<"Check if this should be the case"<<endl;
  return single_star::sub_giant_D_factor(mass, metalicity);
}

//Eq.38
real single_star::sub_giant_D_factor(const real mass, 
				     const real z) {

  real zeta = log10(z/cnsts.parameters(solar_metalicity));
  real D0 = 5.37  + 0.135*zeta;

  real m_hef = helium_flash_mass(z);
  real D;
  if (mass<=m_hef)
    D = pow(10., D0);
  else {
    real log_D = max(-1., max(0.975*D0 - 0.18*mass, 0.5*D0 - 0.06*mass));
    if (mass<2.5) {
      real D1 = log_D;
      log_D = lineair_interpolation(mass, m_hef, 2.5, D0, D1);
    }
    D = pow(10., log_D);
  }

  return D;
}

real single_star::sub_giant_p_parameter(const real mass, 
					const real z) {
  
  real m_hef = helium_flash_mass(z);
  real p; 
  if (mass<=m_hef)
    p = 6;
  else if(mass>=2.5)
    p = 5;
  else 
    p = lineair_interpolation(mass, m_hef, 2.5, 6., 5.);

  return p;
}

real single_star::sub_giant_q_parameter(const real mass, 
					const real z) {

  real m_hef = helium_flash_mass(z);
  real q; 
  if (mass<=m_hef)
    q = 3;
  else if(mass>=2.5)
    q = 2;
  else 
    q = lineair_interpolation(mass, m_hef, 2.5, 3., 2.);

  return q;
}
/*
//Eq.34
real single_star::FGB_core_mass_luminosity_relation(const real time, 
						    const real mass, 
						    const real z) {
  cerr << "FGB_core_mass_luminosity_relation(...)"<<endl;

  //  real lum = get_luminosity();
  //  real mc = FGB_core_mass_luminosity_relation(lum, mass, z);
  //  PRL(mc);

  real A_H = sub_giant_Ah_estimator(mass);
  real p = sub_giant_p_parameter(mass, z);
  real D = sub_giant_D_factor(mass, z);
  real t_bgb = base_giant_branch_time(mass, z);
  real l_bgb = base_giant_branch_luminosity(mass, z);

  real t_inf = specific_time_limit(A_H, t_bgb,
				   D, l_bgb, p);

  PRC(time);PRL(t_inf);
  real arg = (p-1)*A_H*D*(t_inf-time);
  real m_core = pow(arg, 1/(1-p));

  return m_core;
}*/

//Eq. 37
real single_star::FGB_core_mass_luminosity_relation(const real lum, 
						    const real mass,
                                                    const real z) {

    real l_x = FGB_x_luminosity(mass,z);
    real mc;
    if (lum <= l_x){
        real D = sub_giant_D_factor(mass, z); 
        real p = sub_giant_p_parameter(mass, z);
        mc=pow(lum/D,1./p);
    }
    else {
        real B = sub_giant_B_factor(mass);
        real q = sub_giant_q_parameter(mass, z);
        mc=pow(lum/B, 1./q);
    }
  return mc;
}

//Eq. 37 for Mx
real single_star::FGB_x_luminosity(const real mass, const real z){
    real D = sub_giant_D_factor(mass, z);
    real p = sub_giant_p_parameter(mass, z);
    real mx = FGB_x_mass(mass, z);
    return D*pow(mx, p);
}

//Eq. 38
real single_star::FGB_x_mass(const real mass, const real z){
    real B = sub_giant_B_factor(mass);
    real D = sub_giant_D_factor(mass, z);
    real p = sub_giant_p_parameter(mass, z); 
    real q = sub_giant_q_parameter(mass, z);
    return pow(B/D, 1./(p-q)); 
}

//Eq.34 
real single_star::FGB_luminosity_core_mass_relation(const real time, 
						    const real mass, 
						    const real z) {
  real A_H = sub_giant_Ah_estimator(mass);
  real t_bgb = base_giant_branch_time(mass, z);
  real l_bgb = base_giant_branch_luminosity(mass, z);
  real D = sub_giant_D_factor(mass, z);
  real p = sub_giant_p_parameter(mass, z); 
  real l_x = FGB_x_luminosity(mass,z);
  real t_x = specific_time_boundary(mass, A_H, t_bgb, l_bgb, D, p, l_x);
  real l_sg; 
  
  if (time  <= t_x){
        
        real t_inf1 = specific_time_limit(A_H, t_bgb,
				   D, l_bgb, p);
        real arg = (p-1)*A_H*D*(t_inf1-time);   
        l_sg = D * pow(arg, p/(1-p));      
  }
  else {
      real q = sub_giant_q_parameter(mass, z);
      real B = sub_giant_B_factor(mass);
      real t_inf2 = specific_time_limit(A_H, t_x,
                                        B, l_x, q);
      real arg = (q-1)*A_H*B*(t_inf2-time);   
      l_sg = B * pow(arg, q/(1-q));      
  }  
  return l_sg; 
}

real single_star::FGB_luminosity_core_mass_relation(const real m_core, 
						    const real mass) {

  real A_H = sub_giant_Ah_estimator(mass);
  cerr<<"Using sub_giant_D_factor without metalicity as function parameter"<<endl;
  real D = sub_giant_D_factor(mass);
  real lum = D*m_core/A_H;
  
  return lum;
}

// Eq.39 see also Eq.34
//real single_star::sub_giant_core_mass(const real time,
//				      const real mass, 
//				      const real z) {
//  real t_bgb = base_giant_branch_time(mass, z);
//  real p = sub_giant_p_parameter(mass, z);
//  real D = sub_giant_D_factor(mass, z);
//  real l_bgb = base_giant_branch_luminosity(mass, z);
//  real A_H = sub_giant_Ah_estimator(mass);
//  real l_x = FGB_x_luminosity(mass,z);
//  real t_x = specific_time_boundary(mass, z, A_H, t_bgb, l_bgb, D, p, l_x);
//
//  real m_core;
//  if(time<=t_x) {
//    real t_inf1 = specific_time_limit(A_H, t_bgb,
//				      D, l_bgb, p);
//    m_core = pow((p-1)*A_H*D*(t_inf1-time), 1./(1-p));
//  }
//  else {
//    real B = sub_giant_B_factor(mass);
//    real q = sub_giant_q_parameter(mass, z);
//    real t_inf2 = specific_time_limit(A_H, t_x,
//				      B, l_x, q);
//    m_core= pow((q-1)*A_H*B*(t_inf2-time), 1./(1-q));
//  }
//  return m_core;
//}

// Eq.39 see also Eq.34
real single_star::determine_core_mass(const real time,
				      const real mass, 
				      const real z, 
				      const real A,  
				      const real t_b,  
				      const real l_b) {
  
  real D = sub_giant_D_factor(mass, z);
  real p = sub_giant_p_parameter(mass, z);
  real l_x = FGB_x_luminosity(mass,z);
  real t_x = specific_time_boundary(mass, A, t_b, l_b, D, p, l_x);
  real m_core;
  if(time<=t_x) {
    real t_inf1 = specific_time_limit(A, t_b,
				      D, l_b, p);
    m_core = pow((p-1)*A*D*(t_inf1-time), 1./(1-p));
   }
  else {
    real B = sub_giant_B_factor(mass);
    real q = sub_giant_q_parameter(mass, z);
    real t_inf2 = specific_time_limit(A, t_x,
				      B, l_x, q);
    m_core= pow((q-1)*A*B*(t_inf2-time), 1./(1-q));
   }
  return m_core;
}


real single_star::determine_age(const real m_core, const real mass,
                                const real z, const real A, 
                                const real t_b, const real l_b){
    real age;  
    real p = sub_giant_p_parameter(mass, z);
    real D = sub_giant_D_factor(mass, z);
    
    real mx = FGB_x_mass(mass, z);
    if (m_core <= mx){        
        real t_inf1 = specific_time_limit(A, t_b,
                                          D, l_b, p);
        age = t_inf1 - pow(m_core, 1.-p)/A/D/(p-1.);
    }
    else{
        real q = sub_giant_q_parameter(mass, z);
        real B = sub_giant_B_factor(mass);
        real l_x = FGB_x_luminosity(mass, z);
        real t_x = specific_time_boundary(mass, A, t_b, l_b, D, p, l_x);
        real t_inf2 = specific_time_limit(A, t_x,
                                          B, l_x, q);
        age = t_inf2 - pow(m_core, 1.-q)/A/B/(q-1.);
    }
    return age;
}




//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Core helium burning (horizontal branch)
//

//Eq.37
//real single_star::sub_giant_luminosity(const real time,
//				       const real mass, 
//				       const real z) {
//
//  real A_H = sub_giant_Ah_estimator(mass);
//  real D = sub_giant_D_factor(mass, z);
//  real p = sub_giant_p_parameter(mass, z);
//  real t_bgb = base_giant_branch_time(mass, z);
//  real l_bgb = base_giant_branch_luminosity(mass, z);
//  real t_inf = specific_time_limit(A_H, t_bgb,
//				   D, l_bgb, p);
//
//  real l_sg = D * pow((p-1)*A_H*D*(t_inf-time), p/(1-p));
//  return l_sg;
//}

// Eq.41
real single_star::specific_time_boundary(const real mass,
                                         const real A,
                                         const real t_b,
                                         const real l_b,
                                         const real D, 
                                         const real p,
                                         const real l_x) { 
    
    real t_inf1 = specific_time_limit(A, t_b,
                                      D, l_b, p);
//    real l_x = FGB_x_luminosity(mass, z);
    real t_x = t_inf1 - (t_inf1-t_b)*pow(l_b/l_x, (p-1.)/p);
    return t_x;
}


// Eq.43
real single_star::helium_ignition_time(const real mass, 
				       const real z) {

  real t_HeI;
  real t_bgb = base_giant_branch_time(mass, z);
  if( mass > helium_ignition_mass(z)){
        t_HeI = t_bgb;
  }
  else{
      real l_HeI = helium_ignition_luminosity(mass, z);
      real l_x = FGB_x_luminosity(mass, z);

      real A_H = sub_giant_Ah_estimator(mass);
      real l_bgb = base_giant_branch_luminosity(mass, z);
      real D = sub_giant_D_factor(mass, z);
      real p = sub_giant_p_parameter(mass, z);
        
      if (l_x>=l_HeI) {
        real t_inf1 = specific_time_limit(A_H, t_bgb,
                          D, l_bgb, p);
        
        t_HeI = t_inf1 - pow(D/l_HeI, (p-1.)/p) / ((p-1.)*A_H*D);
      }
      else {
        real B = sub_giant_B_factor(mass);
        real q = sub_giant_q_parameter(mass, z);
        real t_x = specific_time_boundary(mass, A_H, t_bgb, l_bgb, D, p, l_x);
        real t_inf2 = specific_time_limit(A_H, t_x,
                          B, l_x, q);


        t_HeI = t_inf2 - pow(B/l_HeI, (q-1.)/q) / ((q-1.)*A_H*B);
      }
  }  
  return t_HeI;
}

real single_star::helium_ignition_time() {

  return helium_ignition_time(relative_mass, metalicity);
} 
  

// Eq.44
real single_star::base_giant_branch_core_mass(const real mass, real z) {

    real m_HeF = helium_flash_mass(z);
    real l_bgb = base_giant_branch_luminosity(m_HeF, z);
    real mc = FGB_core_mass_luminosity_relation(l_bgb, m_HeF,z);
  
    real mc_bagb = base_AGB_core_mass(mass, z);
    real c1 = 9.20925e-5;
    real c2 = 5.402216;
    real C = pow(mc, 4) - c1*pow(m_HeF, c2);
    real mc_bgb = min(0.95*mc_bagb, pow(C + c1*pow(mass, c2), 0.25));
    return mc_bgb;
}

// Eq.57
real single_star::core_helium_burning_timescale(const real mass, 
						const real z) {

  real m_HeF = helium_flash_mass(z);
  real t_He;
  if (mass<m_HeF) {
      real b39 = smc.b(39,z);
      real b40 = smc.b(40,z);
      
      // mc represents here the core mass of a star of mass m at the beginning of the 
      // core_helium_burning_fase
      real l_HeI = helium_ignition_luminosity(mass, z); 
      real mc = FGB_core_mass_luminosity_relation(l_HeI, mass, z);
      real t_Hems = helium_main_sequence_time_for_solar_metalicity(mc);
      real a = (core_helium_burning_timescale(m_HeF,z) - b39)/b39;
      real mu = (mass-mc)/(m_HeF-mc);
      t_He = (b39 + (t_Hems - b39)*pow(1-mu, b40))
         * (1 + a*exp(15*(mass-m_HeF)));
  }
  else {
      real b41 = smc.b(41,z);
      real b42 = smc.b(42,z);
      real b43 = smc.b(43,z);
      real b44 = smc.b(44,z);
      
      real t_bgb = base_giant_branch_time(mass, z);
      t_He = t_bgb*(b41*pow(mass, b42) + b43*pow(mass, 5))
         / (b44 + pow(mass, 5));
  }
  return t_He;
}

    
real single_star::core_helium_burning_timescale() {

  real t_cHe = core_helium_burning_timescale(relative_mass,
					     metalicity);

  return t_cHe;
}

// Eq.66
real single_star::base_AGB_core_mass(const real mass, const real z) {

  real b36 = smc.b(36,z);
  real b37 = smc.b(37,z);
  real b38 = smc.b(38,z);

  real m_core = pow(b36*pow(mass, b37) + b38, 0.25);

  return m_core;
}

// Eq.66 inverted
real single_star::base_AGB_relative_mass(const real m_core, const real z) {
    
    real b36 = smc.b(36,z);
    real b37 = smc.b(37,z);
    real b38 = smc.b(38,z);
    
    real m_rel = pow((pow(m_core, 4) - b38) / b36, 1/b37); 
    
    return m_rel;
}


//sect.5.3: Eq.66-
real single_star::helium_ignition_core_mass(const real mass, 
					    const real z) {
    
  real mc_HeI;
  real m_HeF = helium_flash_mass(z);
  if (mass < m_HeF){ 
    real l_HeI = helium_ignition_luminosity(mass, z);
    mc_HeI = FGB_core_mass_luminosity_relation(l_HeI, mass, z);}
  else {
    //mc_HeI = base_giant_branch_core_mass(mass, z); //Eq.44    
    // should be Eq. 44 but with l_bgb replaced by LHeI 

    real l_HeI = helium_ignition_luminosity(m_HeF, z);
    real mc = FGB_core_mass_luminosity_relation(l_HeI, m_HeF,z);
  
    real mc_bagb = base_AGB_core_mass(mass, z);
    real c1 = 9.20925e-5;
    real c2 = 5.402216;
    real C = pow(mc, 4) - c1*pow(m_HeF, c2);
    mc_HeI = min(0.95*mc_bagb, pow(C + c1*pow(mass, c2), 0.25));
  }

  return mc_HeI;
}

// function needed for possible track hydrogen accreting helium star can turn into horizontal branch star, not implemented currently

// Eq.67
//real single_star::core_helium_burning_core_mass_from_helium_star(const real time,
//                                                const real mass, 
//                                                const real z,
//                                                const real mc) {

//real single_star::core_helium_burning_core_mass_from_helium_star(const real mass, 
//                                                const real z) {
//    
//    real mc_bagb = base_AGB_core_mass(mass, z);
//    real mc_HeI = helium_ignition_core_mass(mass, z);
    
// //    real t_hems = helium_main_sequence_time_for_solar_metalicity(mc);//helium_core_mass
// //    real tau = time / t_hems;    
// //cerr<<"evt relative age setten, zodat tau niet nodig, en dan later relative age setten als fractie van cheb time"<<endl;

//    real tau = relative_age;
//    PRL(tau);
//
//    real m_core = (1-tau)*mc_HeI + tau*mc_bagb;
//    return m_core;
//    
//}



//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Shell helium burning and exhausted Hydrogen core (AGB)
//

// Eq.36, Eq.42, Eq.70 
real single_star::specific_time_limit(const real A,
				      const real t_min,
				      const real DB,
				      const real l_min,
				      const real pq) {

  real t_inf = t_min + pow(DB/l_min, (pq-1)/pq)/((pq-1)*DB*A);
  return t_inf;
}

//Eq.68
real single_star::AGB_A_He_estimator() {

    real A_He = 7.66E-5;
    //A_He=8.E-5;
    return A_He;
}


// Eq.74
real single_star::AGB_radius(const real lum, const real mass,
			     const real mass_tot, const real z) {
  real b1 = smc.b(1,z);
  real b2 = smc.b(2,z);

  real A, b50;
  real m_HeF = helium_flash_mass(z);
  if(mass <= m_HeF-0.2) {

    b50 = smc.b(3,z);
    A = smc.b(56,z) + smc.b(57,z)*mass_tot;
  }
  else if(mass >= m_HeF) {
    real b51 = smc.b(51,z);
    real b52 = smc.b(52,z);
    real b53 = smc.b(53,z);
    real b54 = smc.b(54,z);
    b50 = smc.b(55,z)*smc.b(3,z);
    A = min(b51/pow(mass_tot, b52), b53/pow(mass_tot, b54));
  }
  else{
    b50 = lineair_interpolation(mass_tot, m_HeF-0.2, m_HeF,smc.b(3,z), smc.b(55,z)*smc.b(3,z) );
    real b51 = smc.b(51,z);
    real b52 = smc.b(52,z);
    real b53 = smc.b(53,z);
    real b54 = smc.b(54,z);
    A = lineair_interpolation(mass_tot, m_HeF-0.2, m_HeF,smc.b(56,z) + smc.b(57,z)*(m_HeF-0.2), min(b51/pow(m_HeF, b52), b53/pow(m_HeF, b54)) );
  }
  real r_agb = A*(pow(lum, b1) + b2*pow(lum, b50)); 
  return r_agb;
}


real single_star::base_AGB_time(const real mass, const real z) {

    real t_HeI = helium_ignition_time(mass, z); 
    real t_He = core_helium_burning_timescale(mass, z);
    real t_bagb = t_HeI + t_He;
    return t_bagb;
    
}


real single_star::TAGB_time(const real mass, const real mass_tot, const real z) {
    real t_tagb;
    real mc_max = maximum_AGB_core_mass(mass,z);
    
    //core should minimally grow 5% on the AGB
    real A_He = AGB_A_He_estimator();
    real t_bagb = base_AGB_time(mass, z);
    real l_bagb = base_AGB_luminosity(mass, z);
    real mc_bagb = determine_core_mass(t_bagb, mass, z, 
                                       A_He, t_bagb, l_bagb); //co core mass
    mc_max = max(mc_max, 1.05*mc_bagb);
    
    real mc_du = dredge_up_core_mass(mass, z);
    real m_x = FGB_x_mass(mass,z);
    real D = sub_giant_D_factor(mass, z); 
    real p = sub_giant_p_parameter(mass, z);
    
    
    if (mc_max <= mc_du) {
        real t_bagb = base_AGB_time(mass, z);
        real A_He = AGB_A_He_estimator();
        real l_bagb = base_AGB_luminosity(mass, z);
        if (mc_max <= m_x){
            real t_inf1 = specific_time_limit(A_He, t_bagb,
                                              D, l_bagb, p);
            t_tagb = t_inf1 - pow(mc_max,1-p)/(p-1)/A_He/D; 
        }
        else{
            real B = sub_giant_B_factor(mass);
            real q = sub_giant_q_parameter(mass, z);
            real l_x = FGB_x_luminosity(mass,z);
            real t_x = specific_time_boundary(mass, A_He, t_bagb, l_bagb, D, p, l_x);
            real t_inf2 = specific_time_limit(A_He, t_x,
                                              B, l_x, q);
            
            t_tagb = t_inf2 - pow(mc_max, 1-q)/(q-1)/A_He/B; 
        }
    }
    else{
        //stop evolving star if the core mass reaches
        // m_Ch = cnsts.parameters(Chandrasekar_mass) or
        // mc_ignite_CO = 0.773*mc_bagb - 0.35 or
        // get_total_mass()
        mc_max = min(mc_max, mass_tot);

        real lambda =  min(0.9, 0.3+0.001*pow(mass, 5)); // Eq.73
        mc_max = (mc_max - lambda* mc_du) / (1.0 - lambda);

        real AH_He = TPAGB_AH_He_estimator();
        real t_du = dredge_up_time(mass, z);
        real l_du = AGB_luminosity(mc_du, mass, z); 
        real l_x = FGB_x_luminosity(mass,z);
        
        if (l_du <= l_x) {
            if (mc_max <= m_x){
                real t_inf1 = specific_time_limit(AH_He, t_du,
                                                  D, l_du, p);
                t_tagb = t_inf1 - pow(mc_max,1-p)/(p-1)/AH_He/D;  
            }
            else{
                real B = sub_giant_B_factor(mass);
                real q = sub_giant_q_parameter(mass, z);
                
                real l_x = FGB_x_luminosity(mass,z);            
                real t_x = specific_time_boundary(mass, AH_He, t_du, l_du, D, p, l_x);
                real t_inf2 = specific_time_limit(AH_He, t_x,
                                                  B, l_x, q);
                t_tagb = t_inf2 - pow(mc_max, 1-q)/(q-1)/AH_He/B; 
            }
        }
        else {
            real B = sub_giant_B_factor(mass);
            real q = sub_giant_q_parameter(mass, z);
            real t_inf2 = specific_time_limit(AH_He, t_du,
                                              B, l_du, q);
            t_tagb = t_inf2 - pow(mc_max, 1-q)/(q-1)/AH_He/B;
            
            if(mc_max<=m_x) {
                cerr<<"ERROR in single_star::TAGB_time"<<endl;
                cerr<<"mc_max <= m_x need tinf1"<<endl;
                PRC(mc_bagb);PRC(mc_du);PRL(m_x);
                PRC(l_du); PRL(l_x);
                PRC(mc_max);PRC( maximum_AGB_core_mass(mass,z));
                PRC(relative_mass);PRC(get_total_mass());PRC(core_mass);PRL(COcore_mass);
                exit(-1);
            }
            
        }
    }
    return max(t_bagb,t_tagb);
}


// Eq.69
real single_star::dredge_up_core_mass(const real mass, 
                                      const real z) {
    
    real mc_bagb = base_AGB_core_mass(mass, z);
    real mc_du;
    if(mc_bagb > 0.8 && mc_bagb < 2.25){
        mc_du = 0.44*mc_bagb + 0.448;
    }
    else{
        mc_du = mc_bagb;
    }
    return mc_du;
}


real single_star::TPAGB_AH_He_estimator() {
    
    real AH_He = 1.27E-5;
    
    return AH_He;
}


// Eq.75
real single_star::maximum_AGB_core_mass(const real mass,
                                        const real z) {
    
    real mc_bagb = base_AGB_core_mass(mass, z);
    real m_Ch = cnsts.parameters(Chandrasekar_mass);
    real mc_max = max(m_Ch, 0.773*mc_bagb - 0.35);
    return mc_max;
}

//Eq.37AGB
real single_star::AGB_luminosity(const real m_core,
                                 const real mass, 
                                 const real z) {
    real B = sub_giant_B_factor(mass);
    real D = sub_giant_D_factor(mass, z);
    real p = sub_giant_p_parameter(mass, z);
    real q = sub_giant_q_parameter(mass, z);
    real l_agb = min(B*pow(m_core, q), D*pow(m_core, p));
    
    return l_agb;
}


// Eq.70
real single_star::dredge_up_time(const real mass, const real z) {
    
    real t_du;   
    real l_du = dredge_up_luminosity(mass, z);
    real l_x = FGB_x_luminosity(mass, z);
    real t_bagb = base_AGB_time(mass, z);
    real l_bagb = base_AGB_luminosity(mass, z);
    real A_He = AGB_A_He_estimator();    // [Msun/(Lsun Myear)]
    real p = sub_giant_p_parameter(mass, z);
    real D = sub_giant_D_factor(mass, z);
    
    if(l_du<=l_x) {
        real t_inf1 = specific_time_limit(A_He, t_bagb,
                                          D, l_bagb, p);
        t_du = t_inf1 - pow(D/l_du, (p-1)/p)/((p-1)*D*A_He);
    }
    else {
        real B = sub_giant_B_factor(mass);
        real q = sub_giant_q_parameter(mass, z);
        real t_x=specific_time_boundary(mass, A_He, t_bagb, l_bagb, D, p, l_x);
        real t_inf2 = specific_time_limit(A_He, t_x, B, l_x, q);
        t_du = t_inf2 - pow(B/l_du, (q-1)/q)/((q-1)*B*A_He);
    }
    return t_du;
}

//Eq.70
real single_star::dredge_up_luminosity(const real mass, const real z) {
    
    real mc_du = dredge_up_core_mass(mass, z);
    real l_du = AGB_luminosity(mc_du, mass, z);
    
    return l_du;
}



//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Core helium burning without envelope (helium star & giant)
//

// Eq.77
real single_star::helium_star_luminosity_for_solar_metalicity(const real mass) {
//He ZAMS
  real l_zahe = 15262 * pow(mass, 10.25)
              / (pow(mass, 9) + 29.54*pow(mass, 7.5) 
                              + 31.18*pow(mass, 6) + 0.0469);   
  
  return l_zahe;
}

// Eq.78
real single_star::helium_star_radius_for_solar_metalicity(const real mass_tot) {
//He ZAMS
  real r_zahe = 0.2391 * pow(mass_tot, 4.6)
              / (pow(mass_tot, 4) + 0.162*pow(mass_tot, 3)+ 0.0065);
  
  return r_zahe;
}

// Eq.79
real single_star::helium_main_sequence_time_for_solar_metalicity(const real mass) {

  real t_Hems = (0.4129 + 18.81*pow(mass, 4) + 1.853*pow(mass, 6))
              / pow(mass, 6.5);
  
    return t_Hems;
}
       
//Eq.80
real single_star::helium_main_sequence_luminosity(const real time, const real mass){
    real L_hezams = helium_star_luminosity_for_solar_metalicity(mass);
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(mass);
    real tau = time / t_Hems;
    real alpha = max(0.0, 0.85-0.08*mass);
    return L_hezams * (1.0 + 0.45*tau + alpha * pow(tau,2));
}
    
    
    
//Eq.81
real single_star::helium_main_sequence_radius(const real time, const real mass, const real mass_tot){
    real R_hezams = helium_star_radius_for_solar_metalicity(mass_tot);
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(mass);
    real tau = time / t_Hems;
    real beta = max(0.0, 0.4-0.22*log10(mass_tot));
    return R_hezams*(1.0+beta*tau-beta*pow(tau,6));
        
}

    
real single_star::terminal_helium_main_sequence_luminosity(const real mass){
    real L_hezams = helium_star_luminosity_for_solar_metalicity(mass);
    real alpha = max(0.0, 0.85-0.08*mass);
    return L_hezams * (1.45 + alpha);
}


//Eq. 84
real single_star::helium_giant_luminosity_from_core_mass(const real m_core, const real mass, const real z){
    real l_He; 
    real m_x = helium_giant_x_mass(mass);
    if (m_core  <= m_x){
        real p = helium_giant_p_parameter();
        real D = helium_giant_D_factor(mass);
        l_He = D * pow(m_core, p);      
    }
    else {        
        real q = helium_giant_q_parameter();
        real B = helium_giant_B_factor();
        l_He = B * pow(m_core, q);
    }  
    return l_He; 
}


//Eq.85
real single_star::helium_giant_radius(const real lum, const real mass, const real mass_tot, const real z){
    
    real lambda = 500*(2.0+pow(mass_tot,5))/pow(mass_tot, 2.5);
    real Rzhe = helium_star_radius_for_solar_metalicity(mass_tot);
    real L_tHems = terminal_helium_main_sequence_luminosity(mass); 
    real R1 = Rzhe * pow(lum/L_tHems, 0.2) + 0.02*(exp(lum/lambda)-exp(L_tHems/lambda));//He Hg
    real R2 = 0.08*pow(lum, 0.75); //He GB
    return min(R1, R2);    
}
    
    

//Eq. 38
real single_star::helium_giant_x_mass(const real mass){
    real p = helium_giant_p_parameter();
    real D = helium_giant_D_factor(mass);
    real q = helium_giant_q_parameter();
    real B = helium_giant_B_factor();
    return pow(B/D, 1./(p-q)); 
} 

//Eq. 37 for Mx
real single_star::helium_giant_x_luminosity(const real mass){
    real p = helium_giant_p_parameter();
    real D = helium_giant_D_factor(mass);
    real mx = helium_giant_x_mass(mass);
    return D*pow(mx, p);
}


real single_star::helium_giant_B_factor() { 
    real B = 4.1E4;
    return B;
}

real single_star::helium_giant_D_factor(const real mass) { 
    real D = 5.5E4 / (1.0 + 0.4*pow(mass, 4));
    return D;
}


real single_star::helium_giant_p_parameter() {
    real p = 5.0;    
    return p;
}

real single_star::helium_giant_q_parameter() {
    real q = 3.0;    
    return q;
}    

//related to Eq. 84
real single_star::helium_giant_initial_core_mass(const real mass, const real z ){
    // z is not used, needed for conformity with linear_function_inversion
    real A_He = AGB_A_He_estimator();
    real t_Hems = helium_main_sequence_time_for_solar_metalicity(mass);
    real time = t_Hems; //only difference with helium_giant::helium_giant_core_mass
    real l_tHems = terminal_helium_main_sequence_luminosity(mass);
    real p = helium_giant_p_parameter();
    real D = helium_giant_D_factor(mass);
    real l_x = helium_giant_x_luminosity(mass);
    real t_x = specific_time_boundary(mass, A_He, t_Hems, l_tHems, D, p, l_x);
    real m_core;
    
    if (time  <= t_x){
        real t_inf1 = specific_time_limit(A_He, t_Hems,
                                          D, l_tHems, p);
        real arg = (p-1)*A_He*D*(t_inf1-time);   
        m_core = pow(arg, 1.0/(1-p));     
    }
    else {
        real q = helium_giant_q_parameter();
        real B = helium_giant_B_factor();
        real t_inf2 = specific_time_limit(A_He, t_x,
                                          B, l_x, q);
        real arg = (q-1)*A_He*B*(t_inf2-time);   
        m_core = pow(arg, 1.0/(1-q));      
    }  
    return m_core;     
}

    
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Small envelope behaviour 
//
    
void single_star::small_envelope_perturbation(){
    real mu = small_envelope_mu(luminosity, get_total_mass(), core_mass);
    if(mu < 1.){
        real lum_c = small_envelope_core_luminosity();
        luminosity = perturb_luminosity(luminosity, lum_c, get_total_mass(), core_mass, mu);
        real rad_c = small_envelope_core_radius();
        if(rad_c < radius){
            radius = perturb_radius(radius, rad_c, get_total_mass(), core_mass, mu);
        }
    }
}
    
    
//Eq.97
real single_star::small_envelope_mu(const real lum, const real mass_tot, const real m_core){
    real kappa = -0.5;
    real lum_0 = 7.0E4;
    real mu = (mass_tot - m_core) / mass_tot * min(5.0, max(1.2, pow(lum/lum_0, kappa)));
    return mu;
    }
    
    
//Eq.99
real single_star::perturb_luminosity(const real lum, const real lum_c, const real mass_tot, const real m_core, const real mu){
    real s = s_helper(mass_tot, mu);
    real lum_perturb = lum_c * pow(lum/lum_c, s);
    return lum_perturb;    
}

//Eq.100
real single_star::perturb_radius(const real rad, const real rad_c, const real mass_tot, const real m_core, const real mu){
    real r = r_helper(rad, rad_c, mass_tot, mu);
    real rad_perturb = rad_c * pow(rad/rad_c, r);
    return rad_perturb;
}
    
    
//Eq.101     
real single_star::s_helper(const real mass_tot, const real mu){
    real b = 0.002 * max(1.0, 2.5/mass_tot);
    real s = (1.0 + pow(b, 3)) * pow(mu/b, 3) / (1. + pow(mu/b,3));
    return s;
    
}

//Eq.102    
real single_star::r_helper(const real rad, const real rad_c, const real mass_tot, const real mu){
    real c = 0.006 * max(1.0, 2.5/mass_tot); 
    real q = log(rad/rad_c);
    real r = (1. + pow(c, 3)) * pow(mu/c, 3) * pow(mu, 0.1/q) / (1. + pow(mu/c, 3));
    return r;
}
    
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// White dwarf
//


// (SPZ+GN: 27 Jul 2000) Gijs sais is better: Nelemans et al 2000
real single_star::white_dwarf_radius(real mass, real time) {
    //safety
    if (time < 1) time = 1;
 
    real r;
    if (mass < 0.8) { 
        real a,b;
        if (mass < 0.4) {
            if (mass < 0.2) mass = 0.2; // radius for M < 0.2 equal to M = 0.2
            
            a = lineair_interpolation(mass , 0.2, 0.4, 0.1, 0.03);
            b = lineair_interpolation(mass , 0.2, 0.4, 0.0175, 0.0044);
        }
        else if (mass < 0.6) {
            
            a = lineair_interpolation(mass, 0.4, 0.6, 0.03, 0.017);
            b = lineair_interpolation(mass, 0.4, 0.6, 0.0044, 0.001);
        }
        else {
            
            a = lineair_interpolation(mass, 0.6, 0.8, 0.017, 0.011);
            b = lineair_interpolation(mass, 0.6, 0.8, 0.001, 0.0005);
        }
        r = a - b*log10(time);
    }
    else { 
        //		Nauenberg, M, 1972, Apj 175, 417
        //		mu=4.039 is the mean molecular weight.
        real mu = 2.;
        real m_rel = min(0.99, mass / cnsts.parameters(Chandrasekar_mass));
        r = (0.0225/mu) 
        * sqrt(1./pow(m_rel, cnsts.mathematics(two_third))
               - pow(m_rel, cnsts.mathematics(two_third)));
    }
    return min(0.1, r);
}

    
//
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Miscellaneous
//
    
real single_star::update_core_and_envelope_mass(const real m_core) {
  
  bool successful_update = false;
  real dm_core = m_core-core_mass;
    
  if (m_core > get_total_mass()){
    cerr << "single_star::update_core_and_envelope_mass limits new core mass to total mass." << endl;
    real m_tot = get_total_mass();
    envelope_mass = 0;
    core_mass = m_tot;
    successful_update = true;
  }    
  else if (dm_core/core_mass >= 0.0 - cnsts.safety(tiny)){
    envelope_mass -= dm_core; 
    core_mass += dm_core;
    successful_update = true;
  }
  else {
      cerr << "WARNING: void update_core_and_envelope_mass(m_core)"<< endl;
      cerr << "    old core mass exceeds new core mass." << endl;
      dump(cerr, false);
      exit(-1);
  }
  return successful_update;
}

real single_star::update_core_and_envelope_mass_TPAGB(const real m_core) {
    //difference with update_core_and_envelope_mass is in the fact that on the TPAGB
    // the core mass can decrease during the second dredge_up 
    bool successful_update = false;
    
    if (m_core > get_total_mass()){
        cerr << "single_star::update_core_and_envelope_mass limits new core mass to total mass." << endl;
        real m_tot = get_total_mass();
        envelope_mass = 0;
        core_mass = m_tot;
    }
    else{
        real dm_core = m_core - core_mass;
        envelope_mass -= dm_core; 
        core_mass += dm_core;
    }
    successful_update = true;
    return successful_update;
}

real single_star::update_COcore_mass(const real mco_core) {
    
    bool successful_update = false;    
    if(mco_core>=COcore_mass- cnsts.safety(tiny) &&
       mco_core<=get_total_mass() + cnsts.safety(tiny)) {
        COcore_mass = mco_core;
        successful_update = true;
    }
    else {
        cerr << "WARNING: void update_COcore(mco_core)"<< endl;
        if(mco_core < COcore_mass- cnsts.safety(tiny)) {
            cerr << "    old COcore mass exceeds new COcore mass." << endl;
            PRC(mco_core);PRL(COcore_mass);
            dump(cerr, false);
            exit(-1);
        }
        else {
            cerr.precision(HIGH_PRECISION);
            cerr << "    new COcore mass exceeds total mass." << endl;
            cerr.precision(STD_PRECISION);
            PRC(mco_core);PRL(get_total_mass());
            dump(cerr, false);
            exit(-1);
        }
    }
    
    return successful_update;
}


// SPZ & SilT Feb 4 2010
// Added routine for inverting relatively smooth functions
// Based on linear iterative procedure until satisfied 
// based on Secant method
real single_star::linear_function_inversion(real (single_star::*fptr)(real, const real), const real x_guess, const real y_value, const real z, const real xmin, const real xmax) {
//    real gx = x_guess; // x_guess is currently not used
    real gx = xmin;
    real gy = (this->*fptr)(gx, z);//gx;

    real gx_old = xmax; 
    real gy_old = (this->*fptr)(gx_old, z);//gx;

    real gx_new;
    real gy_new;

    real ymin = (this->*fptr)(xmin, z);
    real ymax = (this->*fptr)(xmax, z);
//    real mu = (xmax-xmin)/(ymax-ymin);
    

    if (ymin < ymax){
        if (y_value > ymax || y_value < ymin){
            PRC(y_value);PRC(ymin);PRL(ymax);
            cerr << "WARNING!!! single_star::linear_function_inversion is unable to find a solution" << endl;
    	    exit(-1);
        }
    }
    else{
        if (y_value < ymax || y_value > ymin){
            PRC(y_value);PRC(ymin);PRL(ymax);
            cerr << "WARNING!!! single_star::linear_function_inversion is unable to find a solution" << endl;
    	    exit(-1);
        }
    }
    
    real dy = 0, dy_old = 0;
    bool within_range = 1; 
    do {
        gx_new = gx - (gy-y_value) * (gx-gx_old)/ (gy-gy_old);  
        
        if (gx_new > xmax || gx_new < xmin){
            within_range = 0;
            gx_new = min(max(gx_new, xmin), xmax);
        }
        
        gy_new = (this->*fptr)(gx_new, z);
//        PRC(gx_new);PRL(gy_new);
        
        dy_old = dy;
        dy = y_value - gy_new;

        if (!within_range && dy_old * dy > 0){
            // with secant method this should not be possible
            // gx lies out of range and the solution does not lie between gx and the previous one     
            PRC(within_range);PRC(dy_old);PRL(dy);
            cerr << "WARNING!!! single_star::linear_function_inversion is unable to find a solution" << endl;
	        exit(-1);
        }
        
        gy_old = gy;
        gy = gy_new;
        gx_old = gx;
        gx = gx_new; 
        within_range = 1; 
    }
    while(abs(dy/y_value)>cnsts.safety(minimum_inversion_precision));
    
//    PRC(x_guess);PRL(gx);
//    PRC(y_value);PRL((this->*fptr)(gx, z));

    return gx;
}
