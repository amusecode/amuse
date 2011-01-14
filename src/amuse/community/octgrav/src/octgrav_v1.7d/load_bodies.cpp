#include "octgrav.h"

float4 octgrav::load_data(vector<float4>& bodies_pos_in) {

  n_bodies = bodies_pos_in.size();
  bodies_pos.resize(n_bodies);
  bodies_id.resize(n_bodies);
  bodies_id_orig.resize(n_bodies);
  
  double4 com = {0,0,0,0};
  for (int i = 0; i < n_bodies; i++) {
    float4 pos = bodies_pos_in[i];
    bodies_pos[i] = pos;
    bodies_id[i]  = i;

    com.x += pos.w * pos.x;
    com.y += pos.w * pos.y;
    com.z += pos.w * pos.z;
    com.w += pos.w;
  }
  com.x *= 1.0/com.w;
  com.y *= 1.0/com.w;
  com.z *= 1.0/com.w;

  return (float4){com.x, com.y, com.z, com.w};
}
