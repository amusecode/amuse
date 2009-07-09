int setup_module();

int cleanup_module();

void new_id(int id);

void delete_id(int id);

void set_mass(int id, double value);

void set_radius(int id, double radius);

void set_position(int id, double p1, double p2, double p3);

void set_velocity(int id, double v1, double v1, double v3);

void evolve(double dt);

