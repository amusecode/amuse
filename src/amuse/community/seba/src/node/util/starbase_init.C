
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

#include "starbase.h"
#include "node.h"

real starbase::m_conv_star_to_dyn = -1; // mass conversion factor
real starbase::r_conv_star_to_dyn = -1; // length conversion factor
real starbase::t_conv_star_to_dyn = -1; // time conversion factor

bool starbase::use_hdyn = true;
 
#define Rsun_pc 2.255e-8		// R_sun/1 parsec = 6.960e+10/3.086e+18;

#ifndef TOOLBOX

void starbase::set_stellar_evolution_scaling(real m_max, real r_hm, real t_hc)
{
    // Set scaling factors from input physical parameters.

    m_conv_star_to_dyn = 1./m_max;  	// m_max = total mass in Msun
    r_conv_star_to_dyn = Rsun_pc/r_hm;  // r_hm = virial radius in pc
    t_conv_star_to_dyn = 1./t_hc;  	// t_hc = crossing time in Myr
}

void starbase::print_stellar_evolution_scaling(ostream & s)
{
    // Print current scaling factors.  Note that the mass unit is
    // artificial, as it is the initial value of N<m>...

    cerr << "    [m]: " << 1. / m_conv_star_to_dyn  << " M_sun\n";
    cerr << "    [R]: " << Rsun_pc / r_conv_star_to_dyn << " pc\n";
    cerr << "    [T]: " << 1. / t_conv_star_to_dyn << " Myr\n";
}

bool starbase::get_stellar_evolution_scaling()
{
    // Attempt to read and set scaling factors.
    // If they appeared in the input stream, the factors should
    // have been set by scan_star_story.  Nevertheless, if they
    // not set, check the star story, just in case.

    // Return true iff successful.

    if (m_conv_star_to_dyn > 0
	&& r_conv_star_to_dyn > 0
	&& t_conv_star_to_dyn > 0) return true;

    if (star_story == NULL) return false;

    if (!find_qmatch(star_story, "mass_scale")) return false;
    if (!find_qmatch(star_story, "size_scale")) return false;
    if (!find_qmatch(star_story, "time_scale")) return false;

    real m = getrq(star_story, "mass_scale");
    real l = getrq(star_story, "size_scale");
    real t = getrq(star_story, "time_scale");

    // Accept scale factors only if all are acceptable.

    if (m == -VERY_LARGE_NUMBER) return false;
    if (l == -VERY_LARGE_NUMBER) return false;
    if (t == -VERY_LARGE_NUMBER) return false;

    m_conv_star_to_dyn = m;		// 1/total mass in Msun
    r_conv_star_to_dyn = l;		// 1/virial radius in Rsun
    t_conv_star_to_dyn = t;		// 1/crossing time in Myr

    return true;
}

starbase::starbase(node* n)
{
    if (n) {
	the_node = n;
	if (the_node->get_star_story()) {
	    star_story = the_node->get_star_story();
	    // cerr << "using existing star story" << endl;
	} else {
	    star_story = mk_story_chapter();
	    // cerr << "creating new star story" << endl;
	}

	// cerr << "starbase:  node = " << n << " story = "
	//	<< star_story << endl;

//	cerr << "About to delete the_node->get_starbase() --1"<<flush <<endl;
//	cerr << the_node->get_starbase() << flush << endl;
	if (the_node->get_starbase())
           delete the_node->get_starbase();
	else 
	    cerr << "No the_node->get_starbase() in starbase::starbase(node*)" << endl;

	the_node->set_starbase(this);
    }
    else {
	the_node = NULL;
	star_story = mk_story_chapter();
    }

    // use_hdyn = false;
}

starbase::starbase(starbase& sb)
{
    the_node   = sb.the_node;
    star_story = sb.star_story;

    // use_hdyn   = sb.use_hdyn;

    sb.star_story = NULL;
    the_node->set_starbase(this);

    // Note that sb is likely to be left inaccessible if not
    // immediately deleted on return...

}


#else

void main()
{
    cerr << "Seems to work..." << endl;
}
  
#endif

