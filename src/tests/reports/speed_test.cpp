#include <stdio.h>
#include <mpi.h>
#include <new>
#include <iostream>

extern "C" {

struct data {
    double x;
    double y;
    double z;
};

int number_of_points_in_one_dimension = 0;
data * model = 0;

int set_data(int index, double vx, double vy, double vz)
{
    if(!model)
    {
        return -1;
    }
    
    if(index > (number_of_points_in_one_dimension * number_of_points_in_one_dimension * number_of_points_in_one_dimension))
    {
        return -2;
    }
    data & m = model[index];
    m.x = vx;
    m.y = vy;
    m.z = vz;
    return 0;
}


int get_data(int index, double * vx, double  * vy, double  * vz)
{
    double data_in[6], data_out[6];
    int status_in,status_out;
    if(!model)
    {
        return -1;
    }
    
    if(index > (number_of_points_in_one_dimension * number_of_points_in_one_dimension * number_of_points_in_one_dimension))
    {
        return -2;
    }
    data & m = model[index];
    *vx = m.x;
    *vy = m.y;
    *vz = m.z;
    data_in[0] = data_in[1] = data_in[2] = 0.0;
    data_in[3] = data_in[4] = data_in[5] = 0.0;
    data_out[0] = data_out[1] = data_out[2] = 0.0;
    data_out[3] = data_out[4] = data_out[5] = 0.0;
    status_in = status_out = 0;
    /*
    MPI::COMM_WORLD.Allreduce(data_in, data_out, 6, MPI::DOUBLE,MPI::SUM);
    MPI::COMM_WORLD.Barrier();
    MPI::COMM_WORLD.Allreduce(&status_in, &status_out, 1, MPI::DOUBLE,MPI::SUM);
    */
    return 0;
}
  
int step() 
{
    if(!model) {
        return -1;
    }
    for(int xindex ; xindex < number_of_points_in_one_dimension; xindex++)
    {
        for(int yindex ; yindex < number_of_points_in_one_dimension; yindex++)
        {
            for(int zindex ; zindex < number_of_points_in_one_dimension; zindex++)
            {
                int index = xindex * number_of_points_in_one_dimension * number_of_points_in_one_dimension;
                index += yindex * number_of_points_in_one_dimension;
                index += zindex;
                
                model[index].x = index;
                model[index].y = model[index].x / (1.0 + model[index].y);
                model[index].z = model[index].x  * model[index].y / (model[index].z + 1e-7);
            }
        }
    }
    return 0;
}    
  
int set_number_of_points_in_one_dimension(int value)
{
    if(model) {
        delete model;
    }
    
    try {
        model = new data[value*value*value];
    } catch (std::bad_alloc &e) {
        number_of_points_in_one_dimension = 0;
        return -1;
    }
    number_of_points_in_one_dimension = value;
    
    return 0;
    
}

int set_data_to_same(int n, double vx, double vy, double vz) {
    for(int i = 0; i < n; i++) {
        set_data(i, vx, vy, vz);
    }
    return 0;
}

int reset()
{
    if(model) {
        delete model;
    }
    model = 0;
    return 0;
}

}

