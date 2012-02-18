//
// single_star_io.C
//

//#include "hdyn.h"
#include "single_star.h"
#include "util_io.h"

istream & single_star::scan_star_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];

    cerr << "\n\nUsing single_star_io version of scan_star_story()...\n\n";

    while(get_line(s,input_line), strcmp(END_STAR, input_line)){
	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

	if (val) {
	    real number;
	    if(0){   // trick to keep the if() statement, for the else below
	    }else if(!strcmp("Type",keyword)){
		char str_tpe[MAX_INPUT_LINE_LENGTH];
		sscanf(val,"%s",str_tpe);

	    }else if(!strcmp("Class",keyword)){
		char str_cls[MAX_INPUT_LINE_LENGTH];
		sscanf(val,"%s",str_cls);

	    }else if(!strcmp("T_cur",keyword)){
		set_current_time( strtod(val,NULL) );

	    }else if(!strcmp("T_rel",keyword)){
		set_relative_age( strtod(val,NULL) );

	    }else if(!strcmp("M_rel",keyword)){
		set_relative_mass( strtod(val,NULL) );

	    }else if(!strcmp("M_env",keyword)){
		set_envelope_mass( strtod(val,NULL) );

	    }else if(!strcmp("M_core",keyword)){
		set_core_mass( strtod(val,NULL) );

	    }else if(!strcmp("M_COcore",keyword)){
		set_COcore_mass( strtod(val,NULL) );

	    }else if(!strcmp("T_eff",keyword)){
		number = strtod(val,NULL);
		/* XXX ignore value?? */

	    }else if(!strcmp("L_eff",keyword)){
		set_luminosity( strtod(val,NULL) );

	    }else if(!strcmp("P_rot",keyword)){           // experimental extra 
	  					          // information for 
		set_rotation_period( strtod(val,NULL) );  // radio pulsars.

	    } else if(!strcmp("B_fld",keyword)){
		set_magnetic_field( strtod(val,NULL) );

	    }else{
		add_story_line(star_story, input_line);
	    }
	}
    }

    //if(current_time<=0) current_time = relative_age;
    return s;
}

ostream& single_star::print_star_story(ostream& s,
				       int short_output)  // default = 0
{
    put_story_header(s, STAR_ID);

    if (short_output) {

	// If we were to handle short output in the Star part of
	// the output stream, this is where it should be done.
	// However, for now it is simpler to handle all stellar
	// quantities in the Dyn output (in hdyn_io).

#if 0
	put_string(s,      "  S  =  ", type_short_string(get_element_type()));
	put_real_number(s, "  T  =  ", temperature());
	put_real_number(s, "  L  =  ", get_luminosity());
#endif

    } else {

	put_string(s,      "  Type   =  ", type_string(get_element_type()));
	if (get_spec_type(Accreting)==Accreting)
	    put_string(s,  "  Class   =  ",
		       type_string(get_spec_type(Accreting)));
	put_real_number(s, "  T_cur  =  ", get_current_time());
	if (get_current_time()!=get_relative_age())
	    put_real_number(s, "  T_rel  =  ", get_relative_age());
	put_real_number(s, "  M_rel  =  ", get_relative_mass());
	put_real_number(s, "  M_env  =  ", get_envelope_mass());
	put_real_number(s, "  M_core =  ", get_core_mass());
	put_real_number(s, "  T_eff  =  ", temperature());
	put_real_number(s, "  L_eff  =  ", get_luminosity());

	// Extra output for stars with CO cores
	if (star_with_COcore()) {
	    put_real_number(s, "  M_COcore  =  ", get_COcore_mass());
	}

	// Extra output for radio- and X-ray pulsars.
	if (get_element_type()==Xray_Pulsar ||
	    get_element_type()==Radio_Pulsar ||
	    get_element_type()==Neutron_Star) {
	    put_real_number(s, "  P_rot  =  ", get_rotation_period());
	    put_real_number(s, "  B_fld  =  ", get_magnetic_field());
	}

	if (star_story)
	    put_story_contents(s, *star_story);
    }

    put_story_footer(s, STAR_ID);

    return s;
}

void extract_line_text(stellar_type& type, real& t_cur, real& t_rel, 
		       real& m_rel, real& m_env, real& m_core, real& co_core,
		       real& t_eff, real& l_eff,
		       real& p_rot, real& b_fld,
		       story& s)
    {

        char line [MAX_INPUT_LINE_LENGTH];
	char type_string[MAX_INPUT_LINE_LENGTH];

        strcpy(line, s.get_text());
        char keyword[MAX_INPUT_LINE_LENGTH];
        const char *val = getequals(line, keyword);

	if (val) {
	    real number;
	    if(0){   // trick to keep the if() statement, for the else below
	    }else if(!strcmp("Type",keyword)){
		char str_tpe[MAX_INPUT_LINE_LENGTH];
		sscanf(val,"%s",type_string); 
		type = extract_stellar_type_string(type_string);
	    }else if(!strcmp("T_cur",keyword)){
		t_cur = strtod(val,NULL);
	    }else if(!strcmp("T_rel",keyword)){
		t_rel = strtod(val,NULL);
	    }else if(!strcmp("M_rel",keyword)){
		m_rel = strtod(val,NULL);
	    }else if(!strcmp("M_env",keyword)){
		m_env = strtod(val,NULL);
	    }else if(!strcmp("M_core",keyword)){
		m_core = strtod(val,NULL);
	    }else if(!strcmp("M_COcore",keyword)){
		co_core = strtod(val,NULL);
	    }else if(!strcmp("T_eff",keyword)){
		t_eff = strtod(val,NULL);
	    }else if(!strcmp("L_eff",keyword)){
		l_eff = strtod(val,NULL);
	    }else if(!strcmp("P_rot",keyword)){       // experimental extra 
		p_rot = strtod(val,NULL);     	      // information for 
						      // radio pulsars
	    }else if(!strcmp("B_fld",keyword)){       // Only for neutron stars
		b_fld = strtod(val,NULL);
	    }
	}
    }

//! currently this does not print metallicity Feb 19th 2009 ST
void extract_story_chapter(stellar_type& type, real& z, real& t_cur, real& t_rel, 
                           real& m_rel, real& m_env, real& m_core, 
			   real& co_core,
                           real& T_eff, real& L_eff,
			   real& p_rot, real& b_fld, story& s)
{
    
    if (!s.get_chapter_flag())
        {
        cerr << "extract_story_chapter: not a story\n";
        exit(1);
        }

    t_rel=-1;
      
    for (story * d = s.get_first_daughter_node();
	 d != NULL;
	 d = d->get_next_story_node()) {

      if (d->get_chapter_flag()) 
	extract_story_chapter(type, z,t_cur, t_rel,
			      m_rel, m_env, m_core, co_core,
			      T_eff, L_eff, p_rot, b_fld, *d);
      else 
	extract_line_text(type, t_cur, t_rel,
			  m_rel, m_env, m_core, co_core,
			  T_eff, L_eff, p_rot, b_fld, *d);
    }
    
    if (t_rel<0)
      t_rel = t_cur;
}

void single_star::star_transformation_story(stellar_type new_type)
{
    char info_line[MAX_STORY_LINE_LENGTH];
    stellar_type old_type = get_element_type();
    real time = get_current_time();

    sprintf(info_line,
	    "%s_to_%s_at_time = %6.2f",
	    type_string(old_type), type_string(new_type), time);

    add_story_line(get_node()->get_log_story(), info_line);

    // cerr output for debugging!

    cerr << endl << get_node()->format_label() << " " << info_line
	 << " Myr (old mass = " << get_total_mass() << ")" << endl;

    if(is_binary_component()) {
#if 0 // This should be fixed some time SPZ:26/02/03
      if(get_binary()->obtain_binary_type() == Merged) {
	cerr << "binary is merged." << endl;
      }
      else if(get_binary()->obtain_binary_type() == Disrupted) {
	cerr << "binary is disrupted, other component: " << endl
	     << type_string(get_companion()->get_element_type()) << endl;
      }
      else {
	cerr << "binary companion: " 
	     << type_string(get_companion()->get_element_type()) << endl;
      }
#endif
      cerr << "binary companion: " 
	   << type_string(get_companion()->get_element_type()) << endl;
      cerr << "parent: " << get_node()->get_parent()->format_label() 
	   << endl;
    }
}

void single_star::post_supernova_story()
{
    char info_line[MAX_STORY_LINE_LENGTH];

    sprintf(info_line, "v_kick = %6.2f, %6.2f, %6.2f",
	    anomal_velocity[0], anomal_velocity[1], anomal_velocity[2]);

    add_story_line(get_node()->get_log_story(), info_line);
}

void single_star::first_roche_lobe_contact_story(stellar_type accretor)
{
    char info_line[MAX_STORY_LINE_LENGTH];
    real time = get_current_time();
     
    sprintf(info_line,
	    "%s_%s_%s_RLOF_to_%s_at_time = %6.2f",
	    type_string(get_element_type()),
	    type_string(get_spectral_class(temperature())),
	    type_string(get_luminosity_class(temperature(), luminosity)),
	    type_string(accretor), time);

    add_story_line(get_node()->get_log_story(), info_line);

}


/* commented out by ST Feb 19th 2009 
void single_star::merge_two_stars_story(stellar_type that_type)
{
    char info_line[MAX_STORY_LINE_LENGTH];
    stellar_type this_type = get_element_type();
    real time = get_current_time();

    sprintf(info_line,
	    "merge_%s_and_%s_at_time = %6.2f",
	    type_string(this_type), type_string(that_type), time);

    add_story_line(get_node()->get_log_story(), info_line);

    // cerr output for debugging!

    //    cerr << endl << get_node()->format_label() << " " << info_line
    //	 << " Myr (old mass = " << get_total_mass() << ")" << endl;

    //    if(is_binary_component()) {
    //      cerr << "binary companion: " 
    //	   << type_string(get_companion()->get_element_type()) << endl;
    //      cerr << "parent: " << get_node()->get_parent()->format_label() 
    //	   << endl;
    //    }
}*/
