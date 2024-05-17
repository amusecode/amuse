from amuse.rfi.core import *


class TestCode(CodeInterface):
    include_headers = ['c_interface.h']
    
    def __init__(self, exefile):
        CodeInterface.__init__(self, exefile)
         
    @legacy_function
    def set_number_of_points_in_one_dimension():
        """
        Set the set number of points in one dimension (N), the total model
        size will be qubed (N*N*N)
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('value',
            dtype='int32',
            direction=function.IN,
            description =  
                "The number of points in one direction")
        function.result_type = 'int32'
        return function  
        
    @legacy_function
    def step():
        """
        Do one step over the N * N * N grid
        """
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def reset():
        """
        Restore the model to its original state
        """
        function = LegacyFunctionSpecification()  
        function.result_type = 'int32'
        return function  
        
    @legacy_function
    def set_data():
        """
        set example vector data
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index',
            dtype='int32',
            direction=function.IN,
            description =  
                "index in the array in range 0 <= index < (N*3)")
        function.addParameter('vx',
            dtype='float64',
            direction=function.IN,
            description =  
                "x component of the vector")
        function.addParameter('vy',
            dtype='float64',
            direction=function.IN,
            description =  
                "y component of the vector")
        function.addParameter('vz',
            dtype='float64',
            direction=function.IN,
            description =  
                "z component of the vector")
        function.can_handle_array = True
        function.result_type = 'int32'
        return function 
        
    @legacy_function
    def set_data_to_same():
        """
        set all vector data to same value
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('max',
            dtype='int32',
            direction=function.IN,
            description =  
                "index in the array in range 0 <= index < (N*3)")
        function.addParameter('vx',
            dtype='float64',
            direction=function.IN,
            description =  
                "x component of the vector")
        function.addParameter('vy',
            dtype='float64',
            direction=function.IN,
            description =  
                "y component of the vector")
        function.addParameter('vz',
            dtype='float64',
            direction=function.IN,
            description =  
                "z component of the vector")
        function.can_handle_array = True
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_data():
        """
        retrieve example vector data
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('index',
            dtype='int32',
            direction=function.IN,
            description =  
                "index in the array in range 0 <= index < (N*3)")
        function.addParameter('vx',
            dtype='float64',
            direction=function.OUT,
            description =  
                "x component of the vector")
        function.addParameter('vy',
            dtype='float64',
            direction=function.OUT,
            description =  
                "y component of the vector")
        function.addParameter('vz',
            dtype='float64',
            direction=function.OUT,
            description =  
                "z component of the vector")
        function.can_handle_array = True
        function.result_type = 'int32'
        return function  
