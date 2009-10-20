"""
Example function for documentation purposess
"""

from amuse.legacy.support.core import legacy_function, RemoteFunction

class ExampleInterface(object):
    @legacy_function   
    def example_function():
        """
        Example template for the other functions defined in this 
        specification. All functions should follow this example..
        """
        function = RemoteFunction()  
        function.addParameter('input', dtype='int32', direction=function.IN
            , description="Typical input parameter, the argument is passed by value to the function.")
        function.addParameter('output', dtype='float64', direction=function.OUT
            , description="""
            Typical output parameter, the argument is passed by reference.
            The argument should point to a valid memory location.
             """)
        function.addParameter('inout', dtype='float64', direction=function.INOUT
            , description="""
            Some arguments can be both input as well as output. The function will
            update the value of the passed argument.
             """)
        function.result_type = 'int32'
        function.result_doc = "Function will return an error code."
        return function
        
    @legacy_function
    def get_example_parameter():
        """
        Retrieve the current value of the parameter. Note, values can be any
        of the supported types.
        """
        function = RemoteFunction()  
        function.addParameter('value', dtype='float64', direction=function.OUT,
            description = "The current value of the parameter.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value was retrieved
        -1 - ERROR
            The code does not have support for this parameter, use this when
            a code does not support a parameter pre-defined in a physical domain
        """
        return function
        
    
    @legacy_function
    def set_example_parameter():
        """
        Update the value of the parameter. The type of the new value argument
        must be the same as the ''get_example_parameter'' function.
        """
        function = RemoteFunction()  
        function.addParameter('value', dtype='float64', direction=function.IN,
            description = "The new value of the parameter.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New value of the parameter was set
        -1 - ERROR
            The code does not have support for this parameter
        """
        return function   
        
    
    @legacy_function
    def initialize_code():
        """
        Let the code perform initialization actions after all parameters have been set.
        Should be called once per running code instance.
        """
        function = RemoteFunction()  
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Code is initialized
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention 
        """
        return function  

        