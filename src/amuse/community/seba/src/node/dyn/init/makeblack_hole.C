
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Replace the specified star, or the star closest to the system center
//// of mass, by a black hole.  Only top-level nodes are considered.
////                
//// Usage:  makeblack_hole [OPTIONS] < input > output
////
//// Options:
////      -c    add comment to snapshot
////      -i    specify black hole index (integer)
////      -M    specify black hole mass; if < 0 specify -fraction of total
////      -r    specify black hole position
////      -s    skip specified mass as black hole candidate
////      -v    specify black hole velocity [relative to circular speed]
////
//// Written by Simon Portegies Zwart.
////
//// Report bugs to starlab@sns.ias.edu.

//                 Simon Portegies Zwart, MIT November 2000
//                Steve McMillan, March 2006

#include "dyn.h"

#define  SEED_STRING_LENGTH  256
char  tmp_string[SEED_STRING_LENGTH];

#ifdef TOOLBOX

#define TOL 1.e-6

local void makeblack_hole(dyn* b, int id, real m_bh, real m_skip,
			  bool r_flag, vec r_bh, bool v_flag, vec v_bh) {

    b->to_com();

    dyn *bh = NULL;

    if (id > 0) {
	for_all_daughters(dyn, b, bi) {
	    if(bi->get_index() == id) {
		bh = bi;
		break;
	    }
	}
	if (!bh) {
	    cerr << "Couldn't find index " << id
		 << "; using central particle." << endl;
	    id = -1;
	} else if (twiddles(bh->get_mass(), m_skip, TOL)) {
	    cerr << "Skipping central particle " << bh->format_label()
		 << endl;
	    id = -1;
	}
    }

    if (id < 0) {
	real r_min = VERY_LARGE_NUMBER;
	for_all_daughters(dyn, b, bi) {
	    if (abs(bi->get_pos()) < r_min
	       && !twiddles(bi->get_mass(), m_skip, TOL)) {
		r_min = abs(bi->get_pos());
		bh = bi;
	    }
	}
    }

    if (!bh) err_exit("makeblack_hole: no replacement candidate found");
    cerr << "makeblack_hole: replacing particle " << bh->format_label()
	 << " (mass = " << bh->get_mass() << ")" << endl;

    if(m_bh <= 1) {
	cerr << "makeblack_hole: fractional black hole mass" << endl;
	m_bh *= -(b->get_mass() - bh->get_mass());
	PRL(m_bh);
    }

    cerr << "makeblack_hole: black hole mass = " << m_bh << endl;
    cerr << "makeblack_hole: black hole pos = " << r_bh << endl;
    cerr << "makeblack_hole: black hole vel = " << v_bh << endl;

    real prev_mass = bh->get_mass();
    bh->set_mass(m_bh);

    if(r_flag) 
	bh->set_pos(r_bh);
    if(v_flag) 
	bh->set_vel(v_bh);
    else
	bh->set_vel(bh->get_vel() * sqrt(prev_mass/m_bh));

    putiq(bh->get_log_story(), "black_hole", 1);

    real m_sum = 0;
    for_all_leaves(dyn, b, bi) {
	m_sum += bi->get_mass();
    }

    b->set_mass(m_sum);

    sprintf(tmp_string,
	    "         black hole added, total mass = %8.2f", m_sum); 
    b->log_comment(tmp_string);
}

local inline real mass(dyn *b, real r) 
{

  real m_sum = 0;
  for_all_leaves(dyn, b, bi) {
    if(abs(bi->get_pos())<=r)
      m_sum += bi->get_mass();
  }
  return m_sum;
}

local inline real vcirc2(dyn *b, real r)	// circular orbit speed:
						// recall vc^2 = r d(phi)/dr
{
    if (r > 0)
	return mass(b, r)/r;
    else
	return 0;
}

int main(int argc, char ** argv) {

    bool  c_flag = FALSE;
    char  *comment;

    real m_bh = 0;
    int id = -1;		// default: most central particle is selected
    bool r_flag = false;
    vec r_bh = 0;
    bool v_flag = false;
    vec v_bh = 0;
    real m_skip = 0;

    check_help();

    extern char *poptarg;
    extern char *poparr[];

    int c;
    const char *param_string = "c:i:M:r:::s:v:::";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.16 $", _SRC_)) != -1)
	switch(c) {
	    case 'c':	c_flag = TRUE;
	    		comment = poptarg;
	    		break;
	    case 'M': 	m_bh = atof(poptarg);
		      	break;
	    case 'r':	r_bh = vec(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	                r_flag = true;
	    		break;
	    case 's': 	m_skip = atof(poptarg);
		      	break;
	    case 'v':	v_bh = vec(atof(poparr[0]),
				   atof(poparr[1]),
				   atof(poparr[2]));
	                v_flag = true;
	    		break;
	    case 'i': 	id = atoi(poptarg);
		      	break;
            case '?': 	params_to_usage(cerr, argv[0], param_string);
		      	exit(1);
	}

    dyn* b;
    while (b = get_dyn()) {

      if (c_flag == TRUE)
	  b->log_comment(comment);
      b->log_history(argc, argv);

      // Check if we have to reinterpret r and v in light of possible
      // external fields and physical parameters.

      //    real mass, length, time;
      //    bool phys = get_physical_scales(b, mass, length, time);

//      if (id>b->n_leaves())
//	err_exit("selected id exceeds particle number");

      if(v_flag) {
	  real v_disp = sqrt(2.);
	  v_bh = v_bh * v_disp * sqrt(vcirc2(b, abs(r_bh)));
      }

      makeblack_hole(b, id, m_bh, m_skip, r_flag, r_bh, v_flag, v_bh);

      real initial_mass = getrq(b->get_log_story(), "initial_mass");

      if (initial_mass > -VERY_LARGE_NUMBER)
	  putrq(b->get_log_story(), "initial_mass", b->get_mass(),
		HIGH_PRECISION);

      real m_sum = b->get_mass();
      real old_mtot = b->get_starbase()->conv_m_dyn_to_star(1);
      if(old_mtot!=m_sum) {
	  real old_r_vir= b->get_starbase()->conv_r_star_to_dyn(1);
	  real old_t_vir= b->get_starbase()->conv_t_star_to_dyn(1);
	  b->get_starbase()->set_stellar_evolution_scaling(m_sum,
							   old_r_vir,
							   old_t_vir);
      }

      put_dyn(b);
      rmtree(b);
    }
}
#endif
