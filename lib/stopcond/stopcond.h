
// public methods
int enable_stopping_condition(int type);
int get_number_of_stopping_conditions_set(int * result);
int is_stopping_condition_set(int type, int * result);
int is_stopping_condition_enabled(int type, int * result);
int disable_stopping_condition(int type);
int has_stopping_condition(int type, int * result);
int get_stopping_condition_info(int index, int * index_of_the_condition);
int get_stopping_condition_particle_index(int index, int index_in_the_condition, int * index_of_particle);


// private methods
#define MAX_NUMBER_OF_SIMULTANIOS_CONDITIONS_SET 100
#define MAX_NUMBER_OF_PARTICLES_PER_INDEX        10

#define COLLISION_DETECTION  0
#define PAIR_DETECTION       1
#define ESCAPER_DETECTION    2


#define COLLISION_DETECTION_BITMAP  1
#define PAIR_DETECTION_BITMAP       2
#define ESCAPER_DETECTION_BITMAP    4

extern long enabled_conditions;
extern long supported_conditions;
extern long set_conditions;

int reset_stopping_conditions();
int next_index();
int set_stopping_condition_info(int index, int type);
int set_stopping_condition_particle_index(int index, int index_in_the_condition, int index_of_particle);
