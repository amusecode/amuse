
int get_mass(int index_of_the_particle, double * mass);

int commit_particles();

int get_time(double * time);

int get_theta_for_tree(double * theta_for_tree);

int set_mass(int index_of_the_particle, double mass);

int evolve(double time);

int get_index_of_first_particle(int * index_of_the_particle);

int get_dt_dia(double * dt_dia);

int get_total_radius(double * radius);

int get_potential_at_point(double eps, double x, double y, double z, double * phi);

int new_particle(int * index_of_the_particle, double mass, double radius, double x, double y, double z, double vx, double vy, double vz);

int get_total_mass(double * mass);

int reinitialize_particles();

int set_eps2(double epsilon_squared);

int set_theta_for_tree(double theta_for_tree);

int get_eps2(double * epsilon_squared);

int set_dt_dia(double dt_dia);

int get_index_of_next_particle(int index_of_the_particle, int * index_of_the_next_particle);

int delete_particle(int index_of_the_particle);

int get_potential(int index_of_the_particle, double * potential);

int synchronize_model();

int set_state(int index_of_the_particle, double mass, double radius, double x, double y, double z, double vx, double vy, double vz);

int get_state(int index_of_the_particle, double * mass, double * radius, double * x, double * y, double * z, double * vx, double * vy, double * vz);

int get_time_step(double * time_step);

int set_use_self_gravity(int use_self_gravity);

int recommit_particles();

int get_kinetic_energy(double * kinetic_energy);

int get_number_of_particles(int * number_of_particles);

int get_epsilon_squared(double * epsilon_squared);

int set_acceleration(int index_of_the_particle, double ax, double ay, double az);

int get_indices_of_colliding_particles(int * index_of_particle1, int * index_of_particle2);

int get_center_of_mass_position(double * x, double * y, double * z);

int set_time_step(double timestep);

int set_epsilon_squared(double epsilon_squared);

int get_center_of_mass_velocity(double * vx, double * vy, double * vz);

int get_radius(int index_of_the_particle, double * radius);

int set_ncrit_for_tree(int ncrit_for_tree);

int set_radius(int index_of_the_particle, double radius);

int cleanup_code();

int recommit_parameters();

int initialize_code();

int get_use_self_gravity(int * use_self_gravity);

int get_potential_energy(double * potential_energy);

int get_gravity_at_point(double eps, double x, double y, double z, double * forcex, double * forcey, double * forcez);

int get_velocity(int index_of_the_particle, double * vx, double * vy, double * vz);

int get_ncrit_for_tree(int * ncrit_for_tree);

int get_position(int index_of_the_particle, double * x, double * y, double * z);

int set_position(int index_of_the_particle, double x, double y, double z);

int get_acceleration(int index_of_the_particle, double * ax, double * ay, double * az);

int commit_parameters();

int set_velocity(int index_of_the_particle, double vx, double vy, double vz);

