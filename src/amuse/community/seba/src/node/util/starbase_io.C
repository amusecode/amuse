
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// starbase_io.C
//

#include "starbase.h"
#include "node.h"
#include "util_io.h"

istream & starbase::scan_star_story(istream& s, int level)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    while (get_line(s,input_line), !matchbracket(END_STAR, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

	if (val) {

	    // NOTE: The code below is significantly different from the old
	    // version, which used is_root() and appears to have been wrong.
	    // The problem is that, because of the logic of get_node_recursive,
	    // is_root can't be used because the parent pointers haven't yet
	    // been set.  Fix is to include level as an additional parameter.

	    // cerr << "node: " << the_node->format_label() << "  "; PRL(level);

	    if (level == 0 && !strcmp("mass_scale", keyword)) {
		m_conv_star_to_dyn = strtod(val, NULL);

	    } else if (level == 0 && !strcmp("size_scale", keyword)) {
		r_conv_star_to_dyn = strtod(val, NULL);

	    } else if (level == 0 && !strcmp("time_scale", keyword)) {
		t_conv_star_to_dyn = strtod(val, NULL);

	    } else {
//		cerr << "Adding " << input_line << " to star story for "
//		     << get_node()->format_label() << endl;
		add_story_line(star_story, input_line);
	    }
	}
    }
    return s;
}

ostream& starbase::print_star_story(ostream& s,
				    int short_output)	// default = 0
{
    put_story_header(s, STAR_ID);

    if (the_node->is_root()) {
       put_real_number(s, "  mass_scale     =  ", m_conv_star_to_dyn);
       put_real_number(s, "  size_scale     =  ", r_conv_star_to_dyn);
       put_real_number(s, "  time_scale     =  ", t_conv_star_to_dyn);
    }

    // Note from Steve (5/01): It seems that, if this virtual function
    // is actually being used, then we have a starbase without a "star"
    // type, and hence no evolving stellar properties.  In that case,
    // we probably have no need for short_format output here, since
    // that is relevant only to the evolution code, so suppress it.

    if (star_story && !short_output)
        put_story_contents(s, *star_story);

    put_story_footer(s, STAR_ID);
    
    return s;
}


//		Function calls for stellar evolution link.
bool starbase::get_use_hdyn()       {return use_hdyn;}
void starbase::set_use_hdyn(bool u) {use_hdyn = u;}

//seba_counters* starbase::get_seba_counters() {return sbc;} 
//void starbase::set_seba_counters(seba_counters *sb) {sbc = sb;}

 

void starbase::dump(ostream&, bool) {} 
real starbase::get_total_mass() {return 0;}
real starbase::get_effective_radius() {return 0;}
real starbase::get_current_time() {return 0;}
real starbase::get_relative_age() {return 0;}
real starbase::get_evolve_timestep() {return 0;}

real starbase::temperature() {return 0;}
real starbase::get_luminosity() {return 0;}

vec starbase::get_anomal_velocity() {
                 vec v; return v;}
void starbase::set_anomal_velocity(const vec v) {}
void starbase::evolve_element(const real) {}
star* starbase::merge_elements(star*) { return NULL; }    // HELP HELP

real starbase::get_semi() {return 0;}
void starbase::set_semi(real a) {}
real starbase::get_eccentricity()    {return 0;}
void starbase::set_eccentricity(real e)    {}
binary_type starbase::get_bin_type() {return Unknown_Binary_Type;}

//     -----  Scaling:  -----

real starbase::conv_m_star_to_dyn(real ms)	// input:  mass (solar)
{return ms*m_conv_star_to_dyn;}			// return: mass (code)

real starbase::conv_r_star_to_dyn(real rs)	// input:  length (solar)
{return rs*r_conv_star_to_dyn;}			// return: length (code)

real starbase::conv_t_star_to_dyn(real ts)	// input:  time (Myr)
{return ts*t_conv_star_to_dyn;}			// return: time (code)

real starbase::conv_m_dyn_to_star(real md)	// input:  mass (code)
{return md/m_conv_star_to_dyn;}			// return: mass (solar)

real starbase::conv_r_dyn_to_star(real rd)	// input:  length (code)
{return rd/r_conv_star_to_dyn;}			// return: length (solar)

real starbase::conv_t_dyn_to_star(real td)	// input:  time (code)
{return td/t_conv_star_to_dyn;}			// return: time (Myr)

// Needed for elegant stellar creation, but makes things less clean...
stellar_type starbase::get_element_type() {return NAS;}

real starbase::sudden_mass_loss() {return 0;}
