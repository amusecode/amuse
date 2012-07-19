
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Add tidal parameters to an input snapshot and write it out 
//// again.  Do not change the input data.  Use before scaling.
////
//// Usage: add_tidal [OPTIONS] < input > output
////
//// Options:      
////               -c    add comment [none]
////               -F    specify tidal_field_type [1]
////                         1) point mass
////                         2) isothermal halo
////                         3) disk (solar neighborhood)
////               -J    specify Jacobi radius scaling factor [no default]
////                     (scales Jacobi radius if set, otherwise virial radius)
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

//   Code moved here from kira_init.C

#ifdef TOOLBOX

#include "dyn.h"

main(int argc, char *argv[])
{
    bool c_flag = false;
    char *comment;		// comment string

    bool J_flag = false;
    real input_r_jacobi = -1;

    // Use of this function means that kira's "-Q" flag is implicit,
    // and the default tidal field type is 1 (point-mass).

    bool F_flag = false;
    int  tidal_field_type = 0;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:F:J:";

    // Parse the argument list:

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.10 $", _SRC_)) != -1) {
	switch (c) {
	    case 'c':	c_flag = TRUE;
			comment = poptarg;
			break;
	    case 'F':	F_flag = true;
	    		tidal_field_type = atoi(poptarg);
			break;
	    case 'J':	J_flag = true;
	    		input_r_jacobi = atof(poptarg);
			break;
	    default:
	    case '?':	params_to_usage(cerr, argv[0], param_string);
			get_help();
			return false;
	}
    }

    dyn *b = get_dyn();
    if (b == NULL) err_exit("Can't read input snapshot");

    b->log_history(argc, argv);
    if (c_flag)	b->log_comment(comment);

    if (tidal_field_type < 1) tidal_field_type = 1;
    if (tidal_field_type > 4) tidal_field_type = 4;

    // Action of add_tidal depends on how the input data were created.
    //
    // 1. Some initialization functions (e.g. makeplummer, makesphere, etc.)
    //    don't know about tidal radii.  They write nothing to the log story
    //    unless they are asked to scale, in which case they write initial_mass
    //    and initial_rvirial.
    //
    // 2. Some (e.g. makeking) have limited knowledge, and write initial_mass,
    //    initial_rvirial, and initial_rtidal_over_rvirial to the log story.
    //
    // 3. Others (e.g. make_aniso_king) know a lot, and write initial_mass,
    //    initial_rvirial, initial_rtidal_over_rvirial, alpha3_over_alpha1,
    //    and kira_tidal_field_type to the log story.
    //
    // In case 1, simply add the ratio initial_rtidal_over_rvirial to the
    // log story.  If initial_rvirial is unknown, it will have to be
    // computed (by scale) in any case before the data are usable.  When
    // the data are used, initial_rtidal_over_rvirial will translate into
    // r_jacobi_over_r_virial.  The tidal field type is taken from the
    // command line.
    //
    // In case 2, initial_rtidal_over_rvirial is already contained in the
    // log story.  In keeping with the older conventions in kira, use
    // the 'J' input variable to scale this quantity (allows "underfilling"
    // King models).  The tidal field type is taken from the command line.
    //
    // In case 3, there is nothing to be done, as the tidal properties
    // are completely determined once the virial radius is set (or reset
    // by scale).  Flag a warning on any attempt to change parameters.

    bool print = false;
    real RtRv;

    if (!find_qmatch(b->get_log_story(), "initial_rtidal_over_rvirial")) {

	// Case 1:  use r_jacobi to set the tidal radius.

	if (input_r_jacobi <= 0)
	    err_exit("add_tidal: insufficient data to set tidal radius.");

	putrq(b->get_log_story(), "initial_rtidal_over_rvirial",
	      input_r_jacobi);
	putiq(b->get_log_story(), "kira_tidal_field_type", tidal_field_type);

	print = true;
	RtRv = input_r_jacobi;

    } else {

	if (find_qmatch(b->get_log_story(), "alpha3_over_alpha1")
	    || find_qmatch(b->get_log_story(), "kira_tidal_field_type")) {

	    // Case 3:  no action, but flag any attempt to reset.

	    if (F_flag) warning("add_tidal: ignored -F flag");
	    if (J_flag) warning("add_tidal: ignored -J flag");

	} else {

	    // Case 2:  scale the existing tidal radius by input_r_jacobi,
	    //          if set.

	    if (J_flag) {
		if (input_r_jacobi <= 0)
		    err_exit("add_tidal: input Jacobi radius < 0");
	    }

	    real r_tidal_over_rvirial
		= getrq(b->get_log_story(), "initial_rtidal_over_rvirial");

	    if (r_tidal_over_rvirial <= 0)
	      err_exit("add_tidal: error reading initial_rtidal_over_rvirial");

	    if (input_r_jacobi > 0)
		r_tidal_over_rvirial *= input_r_jacobi;
	    else
		warning("add_tidal: ignored input_r_jacobi < 0");

	    putrq(b->get_log_story(), "initial_rtidal_over_rvirial",
		  r_tidal_over_rvirial);
	    putiq(b->get_log_story(), "kira_tidal_field_type",
		  tidal_field_type);

	    print = true;
	    RtRv = r_tidal_over_rvirial;
	}
    }

    if (print)
	cerr << "add_tidal:  Rt/Rv = " << RtRv << ", field type = "
	     << tidal_field_type << endl;

    put_dyn(b);
}
#endif
