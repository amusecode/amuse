//
// star.C
//
// derived class from starbase

#include "star.h"


//#include "seba.h"
//seba_counters star::get_seba_counters() {return sbc;}
//void star::set_seba_counters(seba_counters sb) {sbc=sb;}

// MEMO (JM, 17-Aug-1996)
//
// In stellar evolution package, each node should tell
// if it is a binary component or not, by checking
// the star part of its parent, and not by checking
// if it has kepler or not... I changed this part, which
// might cause some other troubles...
//

star::star(node* n) : starbase(n) { }

star::star(star& st) : starbase(st) { }

bool star::is_binary() {
    return (the_node->is_parent() && !the_node->is_root());
}

bool star::is_binary_component() {

  if (the_node->is_low_level_leaf()){
	if(the_node->get_parent()->
	   get_starbase()->get_element_type()==Double) {
	    return true;
	}
    }
    return false;
}

bool star::is_star_in_binary() {
  
    return (!the_node->is_parent() &&
	    is_binary_component());
}

star* star::get_binary()
{

  if (is_binary_component())
    return (star*)the_node->get_parent()->get_starbase();
  else
    err_exit("star* get_binary_node: no binary node.");
}

star* star::get_companion()
{

  if (is_binary_component())
    return (star*)the_node->get_binary_sister()
                          ->get_starbase();
  else
    err_exit("star* get_companion: no companion.");
}

star* star::get_companion(star* str)
{
  
  if (str->is_binary_component())
    return str->get_companion();
  else
    err_exit("star* get_companion(star*): no companion.");
}

star* star::get_primary()
{

  if (is_binary_component())
    if (get_total_mass()>get_companion()->get_total_mass())
      return this;
    else
      return get_companion();
  else
    if (the_node->get_oldest_daughter()->get_binary_sister()
	        ->get_starbase()
	        ->get_total_mass()
	>=
	the_node->get_oldest_daughter()
	        ->get_starbase()
	        ->get_total_mass())
      
      return (star*)the_node->get_oldest_daughter()
	                    ->get_binary_sister()
	                    ->get_starbase();
  
    else
      return (star*)the_node->get_oldest_daughter()
	                    ->get_starbase();
}

star* star::get_secondary() {

    if (is_binary_component())
	if (get_total_mass() < get_companion()->get_total_mass())
	    return this;
	else
	    return get_companion();
    else
	if (the_node->get_oldest_daughter()->get_binary_sister()
	    	    ->get_starbase()
	    	    ->get_total_mass() <
	    the_node->get_oldest_daughter()
	    	    ->get_starbase()
	    	    ->get_total_mass())
	    return (star*)the_node->get_oldest_daughter()
				  ->get_binary_sister()->get_starbase();
	else
	    return (star*)the_node->get_oldest_daughter()->get_starbase();
}

star* star::get_initial_primary() {

  if (get_primary()->get_identity() < get_secondary()->get_identity())
    return get_primary();
  else
    return get_secondary();
}

star* star::get_initial_secondary() {
  
  return get_companion(get_initial_primary());

}
  

