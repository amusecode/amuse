from amuse.support.units import units
from amuse.support.data.core import Particles

from amuse.legacy import *

class EVtwin(LegacyInterface): 
    include_headers=['muse_worker.h']
    
    def __init__(self):
        LegacyInterface.__init__(self)    
    

    @legacy_function   
    def wrapper_swap_in_star():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        return function


    @legacy_function   
    def wrapper_swap_out_star():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        return function


    @legacy_function   
    def wrapper_select_star():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        return function


    @legacy_function   
    def wrapper_flush_star():
        return RemoteFunction()        
        
    
    @legacy_function   
    def wrapper_get_luminosity():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    @legacy_function   
    def wrapper_get_mass():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    
    @legacy_function   
    def wrapper_get_radius():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    @legacy_function   
    def wrapper_get_temperature():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    @legacy_function   
    def wrapper_get_age():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    @legacy_function   
    def wrapper_get_stellar_type():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def wrapper_get_stellar_type():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def wrapper_load_zams_star():
        function = RemoteFunction()  
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.addParameter('age_tag', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def wrapper_twin_evolve():
        function = RemoteFunction()
        function.result_type = 'i'
        return function
    
    @legacy_function
    def wrapper_set_init_dat_name():
        function = RemoteFunction()
        function.addParameter('new_init_dat_name', dtype='uint8', direction=function.IN)
        return function
        
    @legacy_function
    def wrapper_set_init_run_name():
        function = RemoteFunction()
        function.addParameter('new_init_run_name', dtype='uint8', direction=function.IN)
        return function

        
    @legacy_function
    def wrapper_initialise_twin():
        function = RemoteFunction()
        function.addParameter('path', dtype='uint8', direction=function.IN)
        function.addParameter('nstars', dtype='i', direction=function.IN)
        function.addParameter('z', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

#extern double *wrapper_allocate_stellar_model_data(void);
#extern void wrapper_get_meshpoint_data(int id, int meshpoint_index, double * masscoordinate, double * radiuscoordinate, double * pressure, double * density, double *
#extern void wrapper_free_stellar_model_data(double *data);
#extern int wrapper_import_stellar_merger(int nmesh, int nvar, double *data, double age);
#extern void wrapper_export_stellar_model(int jstar, double *data);
#extern int wrapper_get_number_of_meshpoints(void);
