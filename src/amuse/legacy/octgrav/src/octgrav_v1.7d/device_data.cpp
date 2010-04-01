#include "octgrav.h"

void octgrav::prepare_data_for_device() {

  for (size_t i = 0; i < children_list.size(); i++) {
    hst.children[i] = children_list[i];
  }
  
  for (size_t i = 0; i < node_list.size(); i++) {
    if (node_list[i] != NULL) {
      hst.node_com[i]  = node_list[i]->get_com();
      hst.node_pos[i]  = node_list[i]->get_boundary();
//       hst.node_pos[i]  = node_list[i]->get_pos();

      hst.node_Qu[i]   = node_list[i]->get_Qu();
      hst.node_Qd[i]   = node_list[i]->get_Qd();

      hst.Oct1[i]   = node_list[i]->get_Oct1();
      hst.Oct2[i]   = node_list[i]->get_Oct2();
      hst.Oct3[i]   = node_list[i]->get_Oct3();

      hst.n_in_node[i] = node_list[i]->get_n_bodies();
      hst.node_bodies_offset[i] = node_list[i]->offset; 
    }
  }

//   for (size_t i = 0; i < leaf_list.size(); i++) {
//     hst.leaf_com[i]  = leaf_list[i]->get_com();
// //     hst.leaf_pos[i]  = leaf_list[i]->get_pos();
//     hst.leaf_pos[i]  = leaf_list[i]->get_boundary();
//     hst.n_in_leaf[i] = leaf_list[i]->get_n_bodies();
//     hst.leaf_bodies_offset[i] = leaf_list[i]->offset; 
//  }

  for (size_t i = 0; i < cell_list.size(); i++) {
    hst.cell_com[i]  = cell_list[i]->get_com();
//     hst.cell_pos[i]  = cell_list[i]->get_pos();
    hst.cell_pos[i]  = cell_list[i]->get_boundary();
    hst.n_in_cell[i] = cell_list[i]->get_total_bodies();
    hst.cell_bodies_offset[i] = cell_list[i]->offset; 
 }
  
}

void octgrav::copy_data_to_device() {
  copyArrayToDevice(dev.bodies_pos,  hst.bodies_pos,  n_bodies * sizeof(float4));
//   copyArrayToDevice(dev.bodies_grav, hst.bodies_grav, n_bodies * sizeof(float4));

  copyArrayToDevice(dev.children, hst.children, children_list.size() * sizeof(float4));

  copyArrayToDevice(dev.node_Qu,       hst.node_Qu,       node_list.size() * sizeof(float4));
  copyArrayToDevice(dev.node_Qd,       hst.node_Qd,       node_list.size() * sizeof(float4));

  copyArrayToDevice(dev.Oct1,       hst.Oct1,       node_list.size() * sizeof(float4));
  copyArrayToDevice(dev.Oct2,       hst.Oct2,       node_list.size() * sizeof(float4));
  copyArrayToDevice(dev.Oct3,       hst.Oct3,       node_list.size() * sizeof(float2));


  copyArrayToDevice(dev.node_com,      hst.node_com,      node_list.size() * sizeof(float4));
  copyArrayToDevice(dev.node_pos,      hst.node_pos,      node_list.size() * sizeof(float4));
  copyArrayToDevice(dev.n_in_node,     hst.n_in_node,     node_list.size() * sizeof(int));
  copyArrayToDevice(dev.node_bodies_offset, hst.node_bodies_offset, node_list.size() * sizeof(int));

//   copyArrayToDevice(dev.leaf_com,  hst.leaf_com,  leaf_list.size() * sizeof(float4));
//   copyArrayToDevice(dev.leaf_pos,  hst.leaf_pos,  leaf_list.size() * sizeof(float4));
//   copyArrayToDevice(dev.n_in_leaf, hst.n_in_leaf, leaf_list.size() * sizeof(int));
//   copyArrayToDevice(dev.leaf_bodies_offset, hst.leaf_bodies_offset, leaf_list.size() * sizeof(int));

  copyArrayToDevice(dev.cell_com,  hst.cell_com,  cell_list.size() * sizeof(float4));
  copyArrayToDevice(dev.cell_pos,  hst.cell_pos,  cell_list.size() * sizeof(float4));
  copyArrayToDevice(dev.n_in_cell, hst.n_in_cell, cell_list.size() * sizeof(int));
  copyArrayToDevice(dev.cell_bodies_offset, hst.cell_bodies_offset, cell_list.size() * sizeof(int));
}

void octgrav::copy_data_from_device() {
  copyArrayFromDevice(hst.bodies_grav, dev.bodies_grav, n_bodies * sizeof(float4));
}
