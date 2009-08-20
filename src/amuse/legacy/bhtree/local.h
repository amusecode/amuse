
// Accessors used only for communication between the C++ main driver
// and the local dynamics module:

void set_dt_dia(double dt);
void set_timestep(double dt);
void set_eps2_for_gravity(double eps2);
void set_theta_for_tree(double theta);
void set_use_self_gravity(int use);
void set_ncrit_for_tree(int ncrit);

