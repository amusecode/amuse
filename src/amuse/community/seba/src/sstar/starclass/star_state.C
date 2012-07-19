//
// star_support: Helper functions to aid in setting up and manipulating
//

#include "star.h"
#include "star_state.h"

#include "main_sequence.h"

inline  single_star * new_single_star(stellar_type type,  // All defaults are
				      int  id,		  // now specified in
				      real z,             // star_state.h
				      real t_cur,	  
				      real t_rel,	  
				      real m_rel,
				      real m_tot,
				      real m_core,
				      real co_core,
				      real p_rot,
				      real b_fld,
				      node* n)
{

     single_star* element = 0;
      switch(type) {
         case Planet:
         case Brown_Dwarf: element = new brown_dwarf(n);
                             break;
         case Main_Sequence: element = new main_sequence(n);
                            break;
         case Hyper_Giant: element = new hyper_giant(n);
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
         default: element = new main_sequence(n);
      }

      element->initialize(id, z, t_cur, t_rel, m_rel, m_tot, m_core, co_core);
      
      return element;
   }

//
// star_hist: helper history function.
//

void star_hist::put_star_hist() {

     printf("\n %d %f %f %f %f",
             star_type, relative_age, relative_mass,
             core_mass, radius);

}

//
// star_state: state op star.
//

star_state::star_state() {
 
  identity = 0;
  type = NAS;
  class_tpe = O5; //O0_O9;
  for (int i=Emission; i<no_of_spec_type; i++) 
    class_spec[i] = NAC;
  mass = radius = velocity = mdot = 0;
  lclass=no_luminosity_class;
}

star_state::star_state(star* str) {

  identity = 0;
  type = NAS;
  class_tpe = O5; //O0_O9;
  for (int i=Emission; i<no_of_spec_type; i++) 
    class_spec[i] = NAC;
  mass = radius = velocity = mdot = 0;
  lclass=no_luminosity_class;

  make_star_state(str);

}

bool remnant(stellar_type type) {

  if ( //type==Brown_Dwarf || type==Planet ||
      type==Carbon_Star || type==Helium_Star || type==Helium_Giant ||
      type==Helium_Dwarf || type==Carbon_Dwarf || type==Oxygen_Dwarf ||
      type==Xray_Pulsar || type==Radio_Pulsar || type==Neutron_Star ||
      type==Black_Hole ||  type==Disintegrated)

    return true;

  else
    return false;
}


void star_state::put_star_state(ostream & s) {

     if(type!=Disintegrated)
       s << type_string(class_tpe);
     
     if (!remnant(type) && type!=Double) {
       if (lclass!=no_luminosity_class)
	 s << type_string(lclass);
       s << type_short_string(type);
     }
     for (int i=Emission; i<no_of_spec_type; i++)
         if (class_spec[i])
            s << "-" << type_short_string((star_type_spec)i);
}

void star_state::init_star_state(star* str) {
    
     identity = str->get_identity();

     mass = str->get_total_mass();
     velocity = str->get_velocity();
     mdot = 0;

     type = str->get_element_type();

     if (!str->remnant()) {
        class_tpe = get_spectral_class(str->temperature()); 
                  //get_spectral_class(mass);
     }
     else if(type==Carbon_Star || Helium_Star || Helium_Giant)
        class_tpe = he;
     else if(type==Helium_Dwarf || type==Carbon_Dwarf ||type==Oxygen_Dwarf)
        class_tpe = wd;
     else if(type==Xray_Pulsar || type==Radio_Pulsar || type==Neutron_Star)
        class_tpe = ns;
     else if(type==Black_Hole)
        class_tpe = bh;
     else if(type==Brown_Dwarf || type==Planet)
        class_tpe = bd;
     else if(type==Disintegrated)
        class_tpe = di;
     else if(type==Double)
        class_tpe = bin;
     else
        class_tpe = no_spectral_class;

     for (int i=Emission; i<no_of_spec_type; i++)
         class_spec[i] = str->get_spec_type(
		         dynamic_cast(star_type_spec, i));

     if(str->is_binary_component() &&
        str->get_effective_radius() >= str->get_binary()->roche_radius(str)) {
       class_spec[Rl_filling] = 1;
     }
}

void star_state::make_star_state(star* str) {
  //  cerr<<"star_state::make_star_state"<<endl;
  //  str->dump(cerr);
  
     type = str->get_element_type();
     lclass = get_luminosity_class(str->temperature(),
				   str->get_luminosity());
     if (identity) // is star initialized?
        mdot = str->get_total_mass() - mass;
     velocity = str->get_velocity();
     mass = str->get_total_mass();
     // Changed the effective radius to real radius so that
     // the double star can see that it is indeed filling its Roche-lobe.
     // (SPZ:2/1998)
     //     radius = str->get_effective_radius();

     if (!str->remnant()){
        class_tpe = get_spectral_class(str->temperature());
		    //get_spectral_class(str->get_total_mass());
     }
     else if(type==Carbon_Star || type==Helium_Star || type==Helium_Giant)
        class_tpe = he;
     else if(type==Helium_Dwarf || type==Carbon_Dwarf || type==Oxygen_Dwarf)
        class_tpe = wd;
     else if(type==Xray_Pulsar || type==Radio_Pulsar || type==Neutron_Star)
        class_tpe = ns;
     else if(type==Black_Hole)
        class_tpe = bh;
     else if(type==Brown_Dwarf || type==Planet)
        class_tpe = bd;
     else if(type==Disintegrated)
        class_tpe = di;
     else if(type==Double)
        class_tpe = bin;
     else
        class_tpe = no_spectral_class;

     for (int i=NAC; i<no_of_spec_type; i++)
         class_spec[i] = str->get_spec_type(
		         dynamic_cast(star_type_spec, i));

     if(str->is_binary_component() &&
        str->get_effective_radius() >= str->get_binary()->roche_radius(str)) {
         class_spec[Rl_filling] = 1;
       }
   }

bool star_state::special() {

  bool special_star = false;

  if (!class_spec[Dsntgr] &&
      (class_spec[Emission] ||
       class_spec[Merger] ||
       class_spec[Blue_Straggler] ||
       class_spec[Barium] ||
       (remnant(type) && class_spec[Accreting])))
    special_star = true;

  return special_star;
}

void put_state(star_state st, ostream & s) {
  //  cerr<<"put_state: "<<endl;
  //cerr << st.type<<" "<<st.lclass<<" "
  //   <<st.mass<<" "<<st.radius<<" "<<st.velocity<<endl;

     if (st.type!=Disintegrated)
       s << type_string(st.class_tpe);
     
     if (!remnant(st.type) && st.type!=Double) {
       if (st.lclass!=no_luminosity_class)
	 s << type_string(st.lclass);
       s << type_short_string(st.type);
     }
     for (int i=Emission; i<no_of_spec_type; i++)
         if (st.class_spec[i])
            s << "-" << type_short_string((star_type_spec)i);
}

void put_short_state(star_state st, ostream & s) {
  
     s << type_string(st.class_tpe);
     if (!remnant(st.type) && st.type!=Double) {
       if (st.lclass!=no_luminosity_class)
	 s << type_string(st.lclass);
       s << type_short_string(st.type);
     }
}

char* type_dominant_state(star_state sts) {

  if (sts.class_spec[Blue_Straggler])
    return type_string(Blue_Straggler);
  
  else if (sts.class_spec[Barium])
    return type_string(Barium);
  
  else if (sts.class_spec[Emission])
    return type_string(Emission);
  
  else if (sts.class_spec[Accreting])
    return type_string(Accreting);
  
  else if (sts.class_spec[Merger])
    return type_string(Merger);
  
  return "?";
  
}

void print_star(starbase* b, ostream & s) {
  
  //  if (has_dstar(dynamic_cast(hdyn*, b->get_node));
  if (b->get_element_type()==Double) {

    // SPZ+SMcM: Temporarily Removed at July 21, 1998
    cerr << "SPZ+SMcM: Temporarily Removed at July 21, 1998"<<endl;
    cerr << "in print_star" << endl;
    //    put_state(make_state(dynamic_cast(double_star*, b)), s);
    
  //  else if (has_sstar(dynamic_cast(hdyn*, b->get_node));
  }  else if (b->get_element_type()) {
    star_state ss(dynamic_cast(star*, b));
    put_state(ss, s);
  }
  else {
    s << "No star attached to starbase of: ";
    b->get_node()->pretty_print_tree(s);
  }
}

void pretty_print_star(starbase* b, int depth_level, ostream & s) {
    int  k = depth_level;
    while (k--)
	s << "  ";
    print_star(b, s);
    s << endl;
    if (b->get_node()->is_parent())
	for_all_daughters(node, b->get_node(), n)
	    pretty_print_star(n->get_starbase(), depth_level + 1, s);
    }

void pretty_print_star(starbase* b, ostream & s) {

      pretty_print_star(b, 0, s);
    }



