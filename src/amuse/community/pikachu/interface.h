typedef struct {
    double mass;                                        /// mass
    double radius;                                      /// radius
    double x, y, z;                                     /// position
    double vx, vy, vz;                                  /// velocity
} dynamics_state;

bool found_particle(int particle_identifier, int *index);
void get_identifier_of_particle_with_index(int index, int *particle_identifier);
