#include "octgrav.h"

double octgrav::CUDA_evaluate_gravity() 
{
  
  return host_evaluate_gravity(1.0/opening_angle,
			       softening*softening,
			       
			       n_bodies,
			       dev.bodies_pos,
			       dev.bodies_grav,

			       children_list.size(),
			       dev.children,

			       node_list.size(),
			       tree.get_pos(),
			       tree.get_com(),
			       dev.node_pos,
			       dev.node_com,
			       dev.node_Qu,
			       dev.node_Qd,
			       dev.Oct1,
			       dev.Oct2,
			       dev.Oct3,
			       dev.n_in_node,
			       dev.node_bodies_offset,

			       cell_list.size(),
			       dev.cell_pos,
			       dev.cell_com,
			       dev.n_in_cell,
			       dev.cell_bodies_offset);
 
}

