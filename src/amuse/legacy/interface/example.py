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

        