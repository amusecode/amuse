#include "stopcond.h"

#include <iostream>

int type_of_stopping_condition_set[MAX_NUMBER_OF_SIMULTANIOS_CONDITIONS_SET];
int index_of_particle_in_stopping_condition[MAX_NUMBER_OF_SIMULTANIOS_CONDITIONS_SET][MAX_NUMBER_OF_PARTICLES_PER_INDEX];

long enabled_conditions = 0;
long set_conditions = 0;
long number_of_stopping_conditions_set = 0;

int enable_stopping_condition(int type)
{
    if(type > 32) {
        return -1;
    }
    enabled_conditions |= 1l << type;
    
    //std::cerr<< enabled_conditions << std::endl;
    //std::cerr<< type<< std::endl;
    return 0;
}

int get_number_of_stopping_conditions_set(int * result)
{
    *result = number_of_stopping_conditions_set;
    
    return 0;
}

int is_stopping_condition_set(int type, int * result)
{
    long mask =  1l << type;
    *result =  (mask & set_conditions) > 0;
    return 0;
}

int is_stopping_condition_enabled(int type, int * result)
{   
    long mask =  1l << type;
    *result = (mask & enabled_conditions) > 0;
    return 0;
}

int disable_stopping_condition(int type)
{
    long mask =  ~(1l << type);
    enabled_conditions &= mask;
    return 0;
}

int has_stopping_condition(int type, int * result)
{
    long mask =  1l << type;
    *result = (mask & supported_conditions) > 0;
    return 0;
}

int get_stopping_condition_info(int index, int * type)
{
    if(index >= number_of_stopping_conditions_set) {
        return -1;
    }
    
    *type = type_of_stopping_condition_set[index];

    return 0;
}

int get_stopping_condition_particle_index(
    int index, 
    int index_in_the_condition, 
    int * index_of_particle)
{
    if(index >= number_of_stopping_conditions_set) {
        return -1;
    }
    if(index_in_the_condition >= MAX_NUMBER_OF_PARTICLES_PER_INDEX) {
        return -1;
    }
    *index_of_particle = index_of_particle_in_stopping_condition[index][index_in_the_condition];

    return 0;
}


int reset_stopping_conditions() {
    number_of_stopping_conditions_set = 0;
    set_conditions = 0;
}

int next_index() {
    number_of_stopping_conditions_set++;
    return number_of_stopping_conditions_set - 1;
}

int set_stopping_condition_info(int index, int type) {
    if(index >= number_of_stopping_conditions_set) {
        return -1;
    }
    long mask =  1l << type;
    set_conditions |= mask;
    
    type_of_stopping_condition_set[index] = type;
    return 0;

}

int set_stopping_condition_particle_index(int index, int index_in_the_condition, int index_of_the_particle) {
    if(index >= number_of_stopping_conditions_set) {
        return -1;
    }
    if(index_in_the_condition >= MAX_NUMBER_OF_PARTICLES_PER_INDEX -1 ) {
        return -1;
    }
    index_of_particle_in_stopping_condition[index][index_in_the_condition] = index_of_the_particle;
    index_of_particle_in_stopping_condition[index][index_in_the_condition + 1] = -1;
    return 0;
}

