
// dyn_ev.C
// member functions for the dyn class that are related to orbit integration.
//

#include "dyn.h"

local void  accumulate_acceleration(dyn * bj,    // n-body system pointer
				    dyn * bi,    // the particle to calculate
			                         // force on
				    real eps_squared)  // softening length
                                                       //squared
{
    if(bj->get_oldest_daughter() != NULL){
	for(dyn * bb = bj->get_oldest_daughter(); bb != NULL;
	    bb = bb->get_younger_sister()){
	    accumulate_acceleration(bb, bi, eps_squared);
	}
    }else{
	if(bi != bj){
	    vec d_pos = bi->get_pos() - bj->get_pos();
	    real soft_d_pos_squared = d_pos * d_pos + eps_squared;
	    real inverse_d_pos_cubed =
	    1 / ( soft_d_pos_squared * sqrt( soft_d_pos_squared ));      
	    
	    bi->inc_acc(-inverse_d_pos_cubed * bj->get_mass() * d_pos);
	}
    }
}

void dyn::calculate_acceleration(dyn * b,
				 real eps_squared)  // softening length squared
{

    if (get_oldest_daughter() != NULL) {
	for_all_daughters(dyn, b, bb) 
	    bb->calculate_acceleration(b, eps_squared);
    } else {
	clear_acc();
	accumulate_acceleration(b, this, eps_squared);
    }
}
