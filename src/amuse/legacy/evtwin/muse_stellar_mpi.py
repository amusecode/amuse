from amuse.support.units import units
from amuse.support.data.core import Particles

from amuse.legacy import *

class EVtwinC(LegacyInterface): 
    include_headers=['muse_worker.h']
    
    def __init__(self):
        LegacyInterface.__init__(self,name_of_the_worker = 'muse_worker_c')    
    

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
        function.addParameter('new_init_dat_name', dtype='string', direction=function.IN)
        return function
        
    @legacy_function
    def wrapper_set_init_run_name():
        function = RemoteFunction()
        function.addParameter('new_init_run_name', dtype='string', direction=function.IN)
        return function

        
    @legacy_function
    def wrapper_initialise_twin():
        function = RemoteFunction()
        function.addParameter('path', dtype='string', direction=function.IN)
        function.addParameter('nstars', dtype='i', direction=function.IN)
        function.addParameter('z', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
class EVtwin(LegacyInterface): 
    use_modules=['twin_library']
    
    def __init__(self):
        LegacyInterface.__init__(self)    
    

    @legacy_function   
    def swap_in_star():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        return function


    @legacy_function   
    def swap_out_star():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        return function


    @legacy_function   
    def select_star():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        return function


    @legacy_function   
    def flush_star():
        return RemoteFunction()        
        
    
    @legacy_function   
    def get_luminosity():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    @legacy_function   
    def get_mass():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    
    @legacy_function   
    def get_radius():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    @legacy_function   
    def get_temperature():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    @legacy_function   
    def get_age():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'd'
        return function
        
    @legacy_function   
    def get_stellar_type():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def get_stellar_type():
        function = RemoteFunction()  
        function.addParameter('id', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def load_zams_star():
        function = RemoteFunction()  
        function.addParameter('mass', dtype='d', direction=function.IN)
        function.addParameter('age_tag', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function
        
    @legacy_function   
    def twin_evolve():
        function = RemoteFunction()
        function.result_type = 'i'
        return function
    
    @legacy_function
    def set_init_dat_name():
        function = RemoteFunction()
        function.addParameter('new_init_dat_name', dtype='string', direction=function.IN)
        return function
        
    @legacy_function
    def set_init_run_name():
        function = RemoteFunction()
        function.addParameter('new_init_run_name', dtype='string', direction=function.IN)
        return function

        
    @legacy_function
    def initialise_twin():
        function = RemoteFunction()
        function.addParameter('path', dtype='string', direction=function.IN)
        function.addParameter('nstars', dtype='i', direction=function.IN)
        function.addParameter('z', dtype='string', direction=function.IN)
        function.result_type = 'i'
        return function

#extern double *allocate_stellar_model_data(void);
#extern void get_meshpoint_data(int id, int meshpoint_index, double * masscoordinate, double * radiuscoordinate, double * pressure, double * density, double *
#extern void free_stellar_model_data(double *data);
#extern int import_stellar_merger(int nmesh, int nvar, double *data, double age);
#extern void export_stellar_model(int jstar, double *data);
#extern int get_number_of_meshpoints(void);
