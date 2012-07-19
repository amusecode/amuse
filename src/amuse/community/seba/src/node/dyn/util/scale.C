
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~
 
//// (Re)scale an N-body system to specified M, Q (=T/U), and E.  Note
//// that if the required Q implies E > 0, no additional energy scaling
//// is applied.
////
//// Note: Only top-level nodes are considered in scale_virial and 
//// scale_energy.  Also, only internal positions and velocities
//// are modified -- the root pos and vel are left untouched.
//// Option "-s" is equivalent to "-M 1 -R 1 -Q 0.5".
//// For systems in which stellar data have been specified, the
//// physical mass, length, and time scales are already defined.
//// Normally those will correspond to the cluster mass, length,
//// and time units, so the "-s" option is probably most appropriate
//// in most cases.  This tool completes the connection between
//// properties and physical units.  In all cases, command-line
//// parameters are specified in N-body units.  As of 7/01, systems
//// with embedded tidal or other external fields are also properly
//// scaled (assuming that the Jacobi radius scales with the virial
//// radius).  Note that the virial radius is defined in terms of the
//// *internal* potential energy only.  For now, only allow the energy
//// to be specified if there are no external (non-tidal) fields.
////
//// Usage: scale [OPTIONS] < input > output
////
//// Options:    
////              -c    zero the center of mass position and velocity [no]
////              -d    debug mode [off]
////              -e    specify softening parameter [0]
////              -E    specify total energy [don't scale]
////              -m/M  specify total mass [don't scale]
////              -q/Q  specify virial ratio [don't scale]
////              -r/R  specify virial radius [don't scale]
////              -s/S  scale to "standard" units (-m 1 -r 1 -q 0.5) [not set]
////
//// Scale writes the following quantities to the root log story if
//// system_time = 0:
////
////      initial_mass
////      initial_total_energy
////      initial_rvirial
////      total_energy
////
//// If no arguments are provided and system_time = 0, scale will compute
//// whatever quantities are needed to set these values.
////
//// If a leaf dyn story contains R_eff, it is scaled to follow r_virial.
////
//// Written by Steve McMillan.
////
//// Report bugs by starlab@sns.ias.edu.

// Steve McMillan, April 1993
//
// Significant changes to both command-line options and internal operation
// by SMcM, July 2001 -- added external fields and redefined virial radius.

// ***  MUST make sure that the definitions of virial radius and virial  ***
// ***  equilibrium are consistent between scale and sys_stats.          ***

#include "dyn.h"

#ifndef TOOLBOX

real get_mass(dyn *b)
{
    real mass = 0;
    for_all_daughters(dyn, b, bb)
	mass += bb->get_mass();
    return mass;
}

void scale_mass(dyn* b, real mscale)
{
    for_all_nodes(dyn, b, bb)			// N.B. all nodes
	bb->scale_mass(mscale);
}

// Note: scale_pos and scale_vel *never* change the root pos and vel.

void scale_pos(dyn* b, real rscale,
	       vec com_pos)			// default = 0
{
    cerr << "scale: "; PRL(rscale);
    cerr << "scale: "; PRL(com_pos);
    for_all_daughters(dyn, b, bb) {		// N.B. all daughters
	bb->set_pos(com_pos
		     + rscale*(bb->get_pos()-com_pos));
	if (bb->is_leaf()) {
	    story *s = bb->get_dyn_story();
	    if (s) {
		real rad = getrq(s, "R_eff");
		if (rad > 0)
		    putrq(s, "R_eff", rad*rscale);

	    }
	}
    }
}

void scale_vel(dyn* b, real vscale,
	       vec com_vel)			// default = 0
{
     for_all_daughters (dyn, b, bb)		// N.B. all daughters
	bb->set_vel(com_vel
		     + vscale*(bb->get_vel()-com_vel));
}

// Note that kinetic energies are all defined in the "root" frame
// -- i.e. they do not include any bulk motion associated with the
// root node itself.  (Usually what we want in this function.)

real get_top_level_kinetic_energy(dyn *b)	// top-level nodes only
{
    real kinetic_energy = 0;

    for_all_daughters(dyn, b, bb)
	kinetic_energy += bb->get_mass() * square(bb->get_vel());

    kinetic_energy *= 0.5;
    return kinetic_energy;
}

real get_kinetic_energy(dyn *b)			// all nodes
{
    real kinetic_energy = 0;

    for_all_leaves(dyn, b, bb)
	kinetic_energy += bb->get_mass()
	    * square(something_relative_to_root(bb, &dyn::get_vel));

    kinetic_energy *= 0.5;
    return kinetic_energy;
}

// NOTE: get_top_level_energies does *not* resolve binaries, and
// does not include any external fields.

void get_top_level_energies(dyn *b, real eps2,
			    real& potential_energy,
			    real& kinetic_energy)
{
    real total_energy;
    calculate_energies(b, eps2,
		       potential_energy, kinetic_energy, total_energy,
		       true);		// ==> CM approximation
}

void scale_virial(dyn *b, real q,
		  real potential_energy,		// internal
		  real& kinetic_energy,			// internal
		  vec com_vel)			// default = 0
{
    // Set the virial ratio by scaling the velocities.
    // Also rescale the kinetic energy.

    if (q > 0) q = -q;  	// only need specify |Q|

    real vscale = sqrt(q*potential_energy/kinetic_energy);
    scale_vel(b, vscale, com_vel);
    kinetic_energy = q*potential_energy;
}

real scale_energy(dyn * b,
		  real e, real& energy,			// internal
		  vec com_pos,			// default = 0
		  vec com_vel)			// default = 0
{
    // Set the energy by scaling positions and velocities, keeping
    // the virial ratio fixed.  Note that eps = 0 is implicit.

    if (energy >= 0) return 1;
    if (e > 0) e = -e;  // only need specify |E| and only E < 0 makes sense!

    real xscale = energy/e;
    real vscale = sqrt(1./xscale);

    scale_pos(b, xscale, com_pos);
    scale_vel(b, vscale, com_vel);

    energy = e;
    return xscale;
}

void scale(dyn *b, real eps,
	   bool c_flag,
	   bool e_flag, real e,
	   bool m_flag, real m,
	   bool q_flag, real q,
	   bool r_flag, real r,
	   bool debug,				     // default = false
	   void (*top_level_energies)(dyn*, real,    // default =
				      real&, real&)) // get_top_level_energies()
{
    // Another consistency check:

    if (q_flag && find_qmatch(b->get_log_story(), "alpha3_over_alpha1")) {
	warning("scale: can't reset virial ratio for aniso_king model");
	q_flag = false;
    }

    bool no_flags = !(c_flag || e_flag || m_flag || q_flag || r_flag);

    // Check for physical parameters.

    real phys_mass, phys_length, phys_time;
    bool phys = get_physical_scales(b, phys_mass, phys_length, phys_time);

    // If phys, we are using physical units.  Phys_mass, length, and time
    // are the equivalents of 1 N-body unit, in Msun, pc, and Myr.

    if (phys) {

	// Flag a non-standard choice of units if physical units have
	// been specified (see note in --help text).

	if ((m_flag && !twiddles(m, 1))
	     || (r_flag && !twiddles(r, 1))
	     || (q_flag && !twiddles(q, 0.5))
	     || (e_flag && !twiddles(abs(e), 0.25)))
	    cerr << "scale: non-standard choice of units for system"
		 << " with physical data" << endl;
    }

    // Optionally transform to the center of mass frame.  In all cases,
    // make sure that the root node is at the center of mass of its
    // daughters.

    if (c_flag) {
	cerr << "scale: transforming to com frame" << endl;
	b->to_com();
    } else
	b->reset_com();		// root node absorbs any numerical errors

    // See if we are dealing with test particles only (still can set
    // mass and radius, but velocity scaling will ignore the internal
    // potential).

    check_set_ignore_internal(b);

    // Define various relevant quantities.  Trust the data in the input
    // snapshot, if current.

    // Always need to know the mass.

    real mass = 0;

    // Believe initial_mass if system_time = 0.

    if ((real)b->get_system_time() == 0
	&& find_qmatch(b->get_log_story(), "initial_mass"))
	mass = getrq(b->get_log_story(), "initial_mass");
    else
	mass = get_mass(b);

    if (debug) {
	cerr << "debug: "; PRC(mass); PRL(get_mass(b));
    }

    // Need to know some energies if the E, Q, or R flags are set.
    // Note that the definition of r_virial now includes *only* the
    // internal potential energy.

    real r_virial = 0, pot_int = 0, pot_ext = 0, kin = 0;

    if (no_flags || e_flag || r_flag || q_flag) {
	if ((real)b->get_system_time() == 0
	    && find_qmatch(b->get_log_story(), "initial_rvirial")) {

	    // Avoid N^2 calculation if possible.

	    r_virial = getrq(b->get_log_story(), "initial_rvirial");
	    pot_int = -0.5*mass*mass/r_virial;
	    kin = get_top_level_kinetic_energy(b);

	    if (debug) {

		// Basic checks:

		real pe_int, ke_tot;
		top_level_energies(b, eps*eps, pe_int, ke_tot);
		cerr << "debug: "; PRC(pot_int); PRC(pe_int);
		cerr << "debug: "; PRC(kin); PRL(ke_tot);
	    }

	} else {

	    // Compute the potential energy.

	    top_level_energies(b, eps*eps, pot_int, kin);
	    r_virial = -0.5*mass*mass/pot_int;
	}

	// Get external potential parameters (including any tidal field)
	// from the input data and compute the external potential energy
	// (excluding any tidal field).

	check_set_external(b);
	pot_ext = get_external_pot(b);

	if (debug) {
	    cerr << "debug: "; PRL(pot_ext);
	}

	// Variable kira_initial_jacobi_radius is set by check_set_external,
	// but it will be rescaled, so just remove it.

	rmq(b->get_log_story(), "kira_initial_jacobi_radius");

	// If an external (non-tidal) field exists, we probably won't
	// naturally want to specify the energy.  In addition, the
	// procedure for doing this is complicated (iterative, and may
	// not converge).  For now, at least, only allow e to be set
	// in the case of no external fields.

	if (e_flag && pot_ext != 0)
	    err_exit("Can't set energy in case of external field.");
    }

    // The relevant kinetic energy for scaling should be in the center
    // of mass frame, which should be the frame of the root node.
    // Compute the CM velocity here and correct if necessary.

    vec com_pos, com_vel;
    compute_com(b, com_pos, com_vel);

    com_pos -= b->get_pos();		// com includes root node; energies
    com_vel -= b->get_vel();		// computed here are relative to the
					// root node

    real com_kin = 0.5*mass*square(com_vel);	// ideally, should be zero

//    if (debug) {
	cerr << "scale: internal "; PRC(com_vel); PRL(com_kin);
//    }

    // First do the mass scaling, if any.  NOTE: Scaling the total mass
    // changes both the total energy and the virial ratio.  No attempt is
    // made to preserve either in cases where they are not specified on
    // the command line.

    if (m_flag) {


        real mfac = m/mass;
	// PRL(mfac);

	scale_mass(b, mfac);

	mass = m;
	pot_int *= mfac*mfac;
	pot_ext *= mfac;
	kin *= mfac;
	com_kin *= mfac;

	cerr << "scale: "; PRL(mass);
	if (debug) {
	    cerr << "debug: "; PRL(com_kin);
	}
    }

    // Simplest choice now is r_flag; e_flag is more complicated,
    // particularly in the presence of an external field.

    if (r_flag) {

	// Rescale all positions to force the virial radius to the
	// specified value.

	real oldvir = pot_int + get_external_virial(b);	  // denominator of
							  // virial ratio
	if (debug) {
	    cerr << "debug: "; PRL((kin-com_kin)/oldvir);
	}

	real rfac =  r/r_virial;
	if (debug) {
	    cerr << "debug: "; PRL(rfac);
	}
	scale_pos(b, rfac, com_pos);

	r_virial = r;
	pot_int /= rfac;

	pot_ext = get_external_pot(b);
	if (debug) {
	    cerr << "debug: "; PRL(pot_ext);
	}

	cerr << "scale: "; PRL(r_virial);

	// Rescale all internal velocities to preserve the virial
	// ratio.  Scaling is already OK in case of a tidal field,
	// so get_external_virial() includes only non-tidal terms.

	if (debug) {
	    cerr << "debug: "; PRC(pot_int); PRL(get_external_virial(b));
	}
	real vfac = (pot_int + get_external_virial(b))/oldvir;
	if (vfac <= 0)
	    cerr << "scale: unable to preserve virial ratio" << endl;
	else {
	    vfac = sqrt(vfac);
	    if (debug) {
		cerr << "debug: "; PRL(vfac);
	    }
	    scale_vel(b, vfac, com_vel);
	    kin = com_kin + vfac*vfac*(kin - com_kin);
	}

	if (debug) {

	    // Check:

	    real pe_int, ke;
	    top_level_energies(b, eps*eps, pe_int, ke);
	    real p_ext = get_external_pot(b);
	    real vir = get_external_virial(b);

	    cerr << "debug: "; PRC(pot_int); PRL(pe_int);
	    cerr << "debug: "; PRC(kin); PRL(ke);
	    cerr << "debug: "; PRC(p_ext); PRL((ke-com_kin)/(pe_int+vir));
	    cerr << "debug: "; PRL(-0.5*mass*mass/pe_int);
	}
    }

    if (e_flag) {

	// Attempt to set the energy while holding the virial ratio
	// constant, or else set both e and q.  Assume no external
	// fields (see above).  Procedure is OK as is in the presence
	// of a tidal field.

	if (q_flag) {
	    kin -= com_kin;
	    scale_virial(b, q, pot_int, kin, com_vel);	// scales kin
	    kin += com_kin;

	    // NOTE: Changing the virial ratio changes the total kinetic
	    // energy of the system, but preserves the potential energy
	    // and total mass.

	    q_flag = false;

	    cerr << "scale: "; PRL(q);
	}

	// NOTE: The energy scaling is done on the assumption that eps = 0,
	//	 and preserves the value of the virial ratio in that case.
	//	 The total mass is always left unchanged.

	real ee = kin - com_kin + pot_int;		// internal energy
	real fac = scale_energy(b, e, ee, com_pos, com_vel);
	ee += com_kin;					// total_energy

	// Update all relevant quantities.

	kin = com_kin + (kin - com_kin)/fac;
	pot_int /= fac;
	r_virial *= fac;

	// For Roche-lobe filling systems, this also rescales
	// alpha1 and alpha3 by fac^{-3}.

	cerr << "scale: "; PRL(e);
    }

    if (q_flag) {

	// Scale the velocities to set q and preserve r_virial.

	real vir = get_external_virial(b);
	cerr << "scale: "; PRL(vir);

	real denominator = vir;
	if (!b->get_ignore_internal()) denominator += pot_int;

	real qvir = -(kin - com_kin) / denominator;
	real vfac = sqrt(q/qvir);
	cerr << "scale: "; PRL(vfac);

	scale_vel(b, vfac, com_vel);
	kin = com_kin + vfac*vfac*(kin - com_kin);

	cerr << "scale: "; PRL(q);
	if (debug) {
	    cerr << "debug: "; PRL((kin-com_kin)/denominator);
	}
    }

    // Update the root log story -- probably best to do this only if
    // system_time = 0.

    if ((real)b->get_system_time() == 0) {
	putrq(b->get_log_story(), "initial_mass", mass, HIGH_PRECISION);
	putrq(b->get_log_story(), "initial_total_energy", kin+pot_int+pot_ext);
	putrq(b->get_log_story(), "initial_rvirial", r_virial);
	putrq(b->get_dyn_story(), "total_energy", kin+pot_int+pot_ext);
    }
}

bool parse_scale_main(int argc, char *argv[],
		      real& eps, bool& c_flag,
		      bool& e_flag, real& e,
		      bool& m_flag, real& m,
		      bool& q_flag, real& q,
		      bool& r_flag, real& r,
		      bool& debug,
		      const char *cvs_id, const char *source)
{
    // Defaults:

    debug = false;
    c_flag = false;
    e_flag = false;
    m_flag = false;
    q_flag = false;
    r_flag = false;
    bool eps_flag = false;
    bool s_flag = false;

    m = 0, q = -1, e = 0, r = 0;
    eps = 0;

    extern char *poptarg;
    int c;
    const char *param_string = "0cdm:M:q:Q:e:E:r:R:sS";

    while ((c = pgetopt(argc, argv, param_string, cvs_id, source)) != -1)
	switch(c) {

	    case '0': break;			// for hdyn compatibility

	    case 'c': c_flag = true;
		      break;
	    case 'd': debug = true;
		      break;
	    case 'E': e_flag = true;
	    	      r_flag = false;
		      e = atof(poptarg);
		      break;
	    case 'e': eps_flag = true;
		      eps = atof(poptarg);
		      break;
	    case 'M':
	    case 'm': m_flag = true;
		      m = atof(poptarg);
		      break;
	    case 'Q':
	    case 'q': q_flag = true;
		      q = atof(poptarg);
		      break;
	    case 'R':
	    case 'r': r_flag = true;
	    	      e_flag = false;
	    	      r = atof(poptarg);
		      break;
	    case 'S':
	    case 's': s_flag = true;
		      m_flag = true;
		      r_flag = true;
		      q_flag = true;
		      m = 1;
	    	      q = 0.5;
		      r = 1;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
                      return false;
	}

    if (e_flag && e >= 0) warning("scale: specified energy >= 0");
    if (m_flag && m <= 0) warning("scale: specified mass <= 0");
    if (q_flag && q <  0) warning("scale: specified virial ratio < 0");
    if (r_flag && r <= 0) warning("scale: specified virial radius <= 0");

    return true;
}

#else

main(int argc, char ** argv)
{
    check_help();

    real m = 0, q = -1, e = 0, r = 0;
    real eps = 0;

    bool debug = false;
    bool c_flag = false;
    bool e_flag = false;
    bool m_flag = false;
    bool q_flag = false;
    bool r_flag = false;

    if (!parse_scale_main(argc, argv, eps,
			  c_flag,
			  e_flag, e, m_flag, m,
			  q_flag, q, r_flag, r,
			  debug,
			  "$Revision: 1.27 $", _SRC_)) {
	get_help();
	exit(1);
    }

    dyn *b;

    while (b  = get_dyn()) {
	b->log_history(argc, argv);
	scale(b, eps, c_flag, e_flag, e, m_flag, m, q_flag, q, r_flag, r, debug);
	put_dyn(b);
	rmtree(b);
    }
}

#endif
