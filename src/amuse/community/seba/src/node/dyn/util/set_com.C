
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Modify positions and velocities to set the center-of-mass position
//// and velocity.  Write com_pos and com_vel to the root dyn story.  If an
//// external field has been specified, the velocity is taken to be in units
//// of the circular orbit speed at the specified location.  Otherwise, the
//// velocity is taken as is.  Positions and velocities are in physical units
//// (parsecs and km/s, if relevant) if physical stellar parameters are known,
//// and N-body units otherwise (or if -n is specified).
////
//// Usage: set_com [OPTIONS] < input > output
////
//// Options:     
////		  -c    add a comment to the output snapshot [false]
////              -n    force interpretation of r and v in N-body units [no]
////              -r    specify center of mass position [(0,0,0)]
////              -v    specify center of mass velocity [(0,0,0)]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.


//   version 1:  Aug/Sep 2001   Steve McMillan

#include "dyn.h"

#ifndef TOOLBOX

void dyn::set_com(vec set_pos, vec set_vel)	// defaults = 0, 0
{
    // Place the root node at the specified pos and vel, and at the
    // center of mass of all top-level nodes.

    vec com_pos;
    vec com_vel;

    compute_com(this, com_pos, com_vel); 	// including 'this' pos and vel

    // Force daughters to have zero CM quantities relative to root.

    com_pos -= pos;
    com_vel -= vel;

    for_all_daughters(dyn, this, bb) {
	bb->inc_pos(-com_pos);
	bb->inc_vel(-com_vel);
    }

    // All new CM information resides in the root node.

    pos = set_pos;
    vel = set_vel;

    // Correct entries in the dyn story.
    // Use INT_PRECISION here because these quentities may be used
    // in detailed calculations elsewhere.

    putvq(get_dyn_story(), "com_pos", set_pos,
	  INT_PRECISION);
    putvq(get_dyn_story(), "com_vel", set_vel,
	  INT_PRECISION);

    // Note: Do not modify the "center" of any external field.
    // This function adjusts only the N-body center of mass.
}

void dyn::reset_com()
{
    // Place the root node at the center of mass of the system, keeping the
    // absolute coordinates and velocities of all top-level nodes unchanged.
    // Note: Do not modify the "center" of any external field.

    vec com_pos;
    vec com_vel;

    compute_com(this, com_pos, com_vel); 	// including 'this' pos and vel

    // Force daughters to have zero CM quantities relative to root.

    com_pos -= pos;
    com_vel -= vel;

    for_all_daughters(dyn, this, bb) {
	bb->inc_pos(-com_pos);
	bb->inc_vel(-com_vel);
    }

    // All new CM information resides in the root node.

    pos += com_pos;
    vel += com_vel;

    // Correct entries in the dyn story.

    putvq(get_dyn_story(), "com_pos", pos,
	  INT_PRECISION);
    putvq(get_dyn_story(), "com_vel", vel,
	  INT_PRECISION);
}

void dyn::offset_com()
{
    // Place the root node at the origin, keeping the absolute coordinates
    // and velocities of all top-level nodes unchanged.  This transformation
    // moves from the standard representation to that used internally by kira.
    // Note: Do not modify the "center" of any external field.

    for_all_daughters(dyn, this, bb) {
	bb->inc_pos(pos);
	bb->inc_vel(vel);
    }

    pos = vel = 0;
}

#else

main(int argc, char ** argv)
{
    bool  c_flag = FALSE;
    char  *comment;
    vec r = 0;
    vec v = 0;

    bool n_flag = false;

    check_help();

    extern char *poptarg;
    extern char *poparr[];
    int c;
    const char *param_string = "c:nr:::v:::";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.10 $", _SRC_)) != -1)
	switch(c) {
	    case 'c':	c_flag = TRUE;
	    		comment = poptarg;
	    		break;
	    case 'n':	n_flag = true;
			break;

	    case 'r':	r = vec(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	    		break;
	    case 'v':	v = vec(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	    		break;

            case '?':	params_to_usage(cerr, argv[0], param_string);
	    		get_help();
	    		exit(1);
        }            

    dyn *b;

    while (b = get_dyn()) {

        if (c_flag == TRUE)
            b->log_comment(comment);
        b->log_history(argc, argv);

	// Check if we have to reinterpret r and v in light of possible
	// external fields and physical parameters.

	real mass, length, time;
	bool phys = get_physical_scales(b, mass, length, time);

	// If phys, we are using physical units.  Mass, length, and time are
	// the physical equivalents of 1 N-body unit, in Msun, pc, and Myr.

	check_set_external(b, false);
	bool ext = (b->get_external_field() != 0);

	// If ext, we have an external field, assumed attractive.  Note
	// that this function doesn't make much sense if ext = false...
	// If ext, the external field has already been converted to
	// N-body units, regardless of phys.

	// Possibilities:
	//
	// no ext field, no physical units:	set r, v as is in N-body units
	// no ext field, physical units:	convert r, v to N-body and set
	// ext field, no physical units:	use N-body units, but use v/vc
	// ext field, physical units:		convert to N-body and use v/vc

	if (phys && !n_flag) {

	    cerr << "set_com:  converting input physical parameters"
		 << " to N-body units" << endl;
	    cerr << "          r = (" << r << ") pc,  v = ("
		 << v << ") v_circ" << endl;
	    cerr << "          N-body length scale = " << length << " pc"
		 << endl;

	    // Convert r and v to N-body units.

	    r /= length;
	    if (!ext) {

		// Input v is in km/s; length/time is the velocity unit
		// in pc/Myr = 3.086e13/3.156e13 km/s = 0.978 km/s.

		real vunit = 0.978*length/time;
		v /= vunit;
	    }

	} else

	    cerr << "set_com:  interpreting input parameters"
		 << " as N-body units" << endl;

	real vc = 0;

	if (ext) {

	    // Convert v from units of vcirc into N-body units.

	    vc = vcirc(b, r);

	    if (vc > 0) {
		v *= vc;
		cerr << "          circular orbit speed = " << vc;
		if (phys && !n_flag) cerr << " (N-body)";
		cerr << ",  period = " << 2*M_PI*abs(r)/vc
		     << endl;
	    }
	}

	// Now r and v are in N-body units, appropriately scaled.

        b->set_com(r, v);
	cerr << "          r = (" << r << "),  v = (" << v << ")"
	     << endl;

	// Add data to the root log story.
	putvq(b->get_log_story(), "initial_CM_pos", r);
	putvq(b->get_log_story(), "initial_CM_vel", v);
	if (ext && vc > 0)
	    putrq(b->get_log_story(), "initial_CM_period", 2*M_PI*abs(r)/vc);

	put_dyn(b);
	rmtree(b);
    }
}

#endif
