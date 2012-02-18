#include "node.h"
#include "double_star.h"

bool has_dstar(node * bi) {

    if(bi->is_low_level_node() && (bi->get_elder_sister() == NULL) && 
       bi->get_parent()->get_starbase()->get_element_type()==Double) {
	return true;
    }else{
	return false;
    }
}

