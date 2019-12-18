'''
This is the interface file for the HOP group-finding algorithm. HOP starts by determining 
the densest neighbor of each particle (which can be the particle proper). The algorithm 
then groups particles together by jumping from densest neighbor to densest neighbor until 
a density peak is reached.

core HOP algorithm functions:
calculate_densities()
    - calculates densities
do_hop()
    - determines densest neighbors then groups the particles and calculates boundaries
      between groups
'''


from amuse.units import generic_unit_system
from amuse.community import *
from amuse.community.interface.common import CommonCodeInterface, CommonCode

class HopInterface(CodeInterface, CommonCodeInterface, LiteratureReferencesMixIn):
    '''
        .. [#] Eisenstein, DJ, Hut, P, HOP: A new group-finding algorithm for N-body simulations, ApJ 498 (1998)
    '''
    include_headers = ['worker_code.h']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="hop_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
    
    @legacy_function
    def new_particle():
        '''
        Add a new particle with the specified position.
        Returns the particle index of the new particle as an integer (note zero offset).
        '''
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT,
          description='index of the particle')
        function.addParameter('mass', dtype='float64', direction=function.IN,
          description='mass of the particle')
        function.addParameter('x', dtype='float64', direction=function.IN,
          description='particle position on x-axis')
        function.addParameter('y', dtype='float64', direction=function.IN,
          description='particle position on y-axis')
        function.addParameter('z', dtype='float64', direction=function.IN,
          description='particle position on z-axis')
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def delete_particle():
        '''
        Remove a particle.
        '''
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
          description='index of the particle to be removed')
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found and removed
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_position():
        '''
        Retrieve the position of a particle.
        '''
        function = LegacyFunctionSpecification()
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
          description='index of a particle')
        function.addParameter('x', dtype='float64', direction=function.OUT,
          description = 'particle position on x-axis')
        function.addParameter('y', dtype='float64', direction=function.OUT,
          description = 'particle position on y-axis')
        function.addParameter('z', dtype='float64', direction=function.OUT,
          description = 'particle position on z-axis')
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found and its position returned
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_mass():
        '''
        Retrieve the position of a particle.
        '''
        function = LegacyFunctionSpecification()
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
          description='index of a particle')
        function.addParameter('mass', dtype='float64', direction=function.OUT,
          description = 'particle position on x-axis')
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found and its position returned
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_position():
        '''
        Set the position of a particle.
        '''
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
          description='index of a particle')
        function.addParameter('x', dtype='float64', direction=function.IN,
          description = 'particle position on x-axis')
        function.addParameter('y', dtype='float64', direction=function.IN,
          description = 'particle position on y-axis')
        function.addParameter('z', dtype='float64', direction=function.IN,
          description = 'particle position on z-axis')
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found and its position altered
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_density():
        '''
        Retrieve the density of a particle.
        '''
        function = LegacyFunctionSpecification()
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
          description='index of a particle')
        function.addParameter('density', dtype='float64', direction=function.OUT,
          description='the density of the specified particle')
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found and its density returned
        -1 - ERROR
            particle could not be found, unphysical density encountered
            or density was not yet calculated
        """
        return function

    @legacy_function
    def set_density():
        '''
        Set the density of a particle.
        '''
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
          description='index of a particle')
        function.addParameter('density', dtype='float64', direction=function.IN,
          description='the density of the specified particle')
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found and its density altered
        -1 - ERROR
            particle could not be found
        """
        return function
    
    @legacy_function
    def get_densest_neighbor():
        '''
        Retrieve the particle index of the densest neighbor of a particle.
        '''
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
          description='index of a particle')
        function.addParameter('index_of_densest_neighbor', dtype='int32', direction=function.OUT,
          description='particle index of the densest neighbor of the specified particle')
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found and its densest neighbor index returned
        -1 - ERROR
            particle could not be found or its densest neighbor was not yet determined
        """
        return function
    
    @legacy_function
    def get_group_id():
        '''
        Retrieve the group ID of a particle (-1 for ungrouped particles).
        '''
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
          description='index of a particle')
        function.addParameter('group_id', dtype='int32', direction=function.OUT,
          description='group id of the group which the specified particle is a part of (-1 for ungrouped particles)')
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found and its group id returned
        -1 - ERROR
            particle could not be found or it was not part of a group
        """
        return function

    @legacy_function
    def get_number_of_particles():
        '''
        Retrieve the number of particles in the HopInstance.
        '''
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True 
        function.addParameter('value', dtype='int32', direction=function.OUT,
          description='the total number of particles')
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def calculate_densities():
        '''
        Calculate the densities of the particles in the HopInstance.
        Use set_density_method() to specify the density calculation 
        method (DEFAULT: gather-scatter cubic spline).
        '''
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            calculating densities successful
        -1 - ERROR
            no particles or parameter values too high (higher than number of partcles)
        '''
        return function
    
    @legacy_function
    def do_hop():
        '''
        Find the densest neighbors, then group the particles. do_hop() assigns 
        group ids to particles, starting with 0 for the largest group.
        '''
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            hop successful
        -1 - ERROR
            no particles or unphysical density encountered when reading densities or
            parameter values too high (higher than number of particles)
        '''
        return function
        
    @legacy_function
    def get_densest_particle_in_group():
        '''
        Retrieve the particle index of the densest particle in a group.
        (Requires do_hop() to be done first.)
        '''
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('group_id', dtype='int32', direction=function.IN,
          description='group identification number')
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT,
          description='particle index of the particle with highest density in specified group')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            group was found and index of densest particle returned
        -1 - ERROR
            group could not be found
        '''
        return function
        
    @legacy_function
    def get_number_of_particles_in_group():
        '''
        Retrieve the number of particles in a group.
        '''
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('group_id', dtype='int32', direction=function.IN,
          description='group identification number')
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT,
          description='number of particles in the specified group')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            group was found and the number of particles in the group returned
        -1 - ERROR
            group could not be found
        '''
        return function
    
    @legacy_function
    def get_average_boundary_density_of_groups():
        '''
        Retrieve the average density of the boundary between two groups.
        (This is O(n_groups) so it might take a while to retrieve.)
        (Requires do_hop() to be done first.)
        '''
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('first_group_id', dtype='int32', direction=function.IN,
          description='identification number of the first group')
        function.addParameter('second_group_id', dtype='int32', direction=function.IN,
          description='identification number of the second group')
        function.addParameter('boundary_density', dtype='float64', direction=function.OUT,
          description='average density at the boundary between specified groups')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            both groups were found and the average density of their boundary returned
        -1 - ERROR
            one of the groups could not be found
        -2 - ERROR
             the groups did not share a boundary
        '''
        return function
    
    @legacy_function
    def get_number_of_groups():
        '''
        Retrieve the number of groups in the HopInstance.
        (Requires do_hop() to be done first.)
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_groups', dtype='int32', direction=function.OUT,
          description='the number of groups found by do_hop()')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            number of groups found by hop returned
        -1 - ERROR
            do_hop was not run first
        '''
        return function
    
    @legacy_function
    def get_number_of_particles_outside_groups():
        '''
        Retrieve the number of particles that do not belong to a group.
        (Requires do_hop() to be done first.)
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT,
          description='the number of particles that were not assigned a group id by do_hop()')
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_nBucket():
        '''
        Set nBucket, to tune the performance of the kd-tree search.
        DEFAULT: 16
        '''        
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.IN,
          description='the value of nBucket affects the performance of the kd-tree search (DEFAULT: 16)')
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_nBucket():
        '''
        Return the bucket parameter to tune the performance of the kd-tree search.
        DEFAULT: 16
        '''        
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.OUT,
          description='the value of nBucket affects the performance of the kd-tree search (DEFAULT: 16)')
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_nDens():
        '''
        Set nDens, the number of particles to smooth over when calculating densities.
        DEFAULT: 64
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.IN,
          description='the value of nDens, the number of particles to smooth over when calculating densities (DEFAULT: 64)')
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def get_nDens():
        '''
        Return the number of particles to smooth over when calculating densities.
        DEFAULT: 64
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.OUT,
          description='the value of nDens, the number of particles to smooth over when calculating densities (DEFAULT: 64)')
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_nHop():
        '''
        Set nHop, the number of particles over which to look for density maximum.
        DEFAULT: 64
        (minimum allowed value of nHop is nMerge+1)
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.IN,
          description='the value of nHop, the number of particles over which to look for density maximum (DEFAULT: 64)')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            nHop was set
        -1 - ERROR
            requested value for nHop was too low
        '''
        return function
    @legacy_function
    def get_nHop():
        '''
        Return the number of particles over which to look for density maximum.
        DEFAULT: 64
        (minimum allowed value of nHop is nMerge+1)
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.OUT,
          description='the value of nHop, the number of particles over which to look for density maximum (DEFAULT: 64)')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            nHop was set
        -1 - ERROR
            requested value for nHop was too low
        '''
        return function
    
    
    @legacy_function
    def set_fDensThresh():
        '''
        Set fDensThresh, the density below which particles are not assigned to any group.
        DEFAULT: -1.0 (no threshold)
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'float64', direction=function.IN,
            description='value of fDensThresh, the minimum density of grouped particles (DEFAULT: no minimum)',
            unit = generic_unit_system.density)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_fDensThresh():
        '''
        Return the density below which particles are not assigned to any group.
        DEFAULT: -1.0 (no threshold)
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'float64', direction=function.OUT,
            description='value of fDensThresh, the minimum density of grouped particles (DEFAULT: no minimum)',
            unit = generic_unit_system.density)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_saddle_densthresh():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'float64', direction=function.IN,
            description='For two groups to merge, the density at their boundary must '
            'exceed this threshold (DEFAULT: 1.75 * outer_density_threshold)',
            unit = generic_unit_system.density)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_saddle_densthresh():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'float64', direction=function.OUT,
            description='For two groups to merge, the density at their boundary must '
            'exceed this threshold (DEFAULT: 1.75 * outer_density_threshold)',
            unit = generic_unit_system.density)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_peak_densthresh():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'float64', direction=function.IN,
            description='Groups with density below this threshold (fringe groups) are '
            'attached to other (proper) groups, or dropped (DEFAULT: max(d_saddle, 2.0 * d_outer))',
            unit = generic_unit_system.density)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_peak_densthresh():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'float64', direction=function.OUT,
            description='Groups with density below this threshold (fringe groups) are '
            'attached to other (proper) groups, or dropped (DEFAULT: max(d_saddle, 2.0 * d_outer))',
            unit = generic_unit_system.density)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def set_saddle_density_threshold_factor():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'float64', direction=function.IN,
            description='For two groups to merge, the density at their boundary must '
            'exceed this factor times the lowest of the two peak densities (set relative_saddle_density_threshold to True)',
            unit = NO_UNIT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_saddle_density_threshold_factor():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'float64', direction=function.OUT,
            description='For two groups to merge, the density at their boundary must '
            'exceed this factor times the lowest of the two peak densities (set relative_saddle_density_threshold to True)',
            unit = NO_UNIT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def set_relative_saddle_density_threshold():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.IN,
            description="Flag to use a saddle-density-threshold relative to the lowest peak density, "
                "instead of the absolute saddle_density_threshold",
            unit = NO_UNIT)
        function.result_type = 'int32'
        return function
    @legacy_function
    def get_relative_saddle_density_threshold():
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.OUT,
            description="Flag to use a saddle-density-threshold relative to the lowest peak density, "
                "instead of the absolute saddle_density_threshold",
            unit = NO_UNIT)
        function.result_type = 'int32'
        return function
    
    
    @legacy_function
    def set_fPeriod():
        '''
        Set the x, y and z periodicity of the simulation box.
        DEFAULT: infinite
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('x', dtype = 'float64', direction=function.IN,
          description='periodicity in x direction')
        function.addParameter('y', dtype = 'float64', direction=function.IN,
          description='periodicity in y direction')
        function.addParameter('z', dtype = 'float64', direction=function.IN,
          description='periodicity in z direction')
        function.result_type = 'int32'
        return function

    @legacy_function
    def get_fPeriod():
        '''
        Get the x, y and z periodicity of the simulation box.
        DEFAULT: infinite
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('x', dtype = 'float64', direction=function.OUT,
          description='periodicity in x direction')
        function.addParameter('y', dtype = 'float64', direction=function.OUT,
          description='periodicity in y direction')
        function.addParameter('z', dtype = 'float64', direction=function.OUT,
          description='periodicity in z direction')
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def set_nMerge():
        '''
        Set nMerge, the number of particles to catalogue group boundaries.
        DEFAULT: 4
        (maximum allowed value of nMerge is nHop-1)
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.IN,
          description='the value of nMerge, the numbder of particles catalogue group boundaries')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            nMerge was set
        -1 - ERROR
            requested value for nMerge was too high
        '''
        return function

    @legacy_function
    def get_nMerge():
        '''
        Return the number of particles to catalogue group boundaries.
        DEFAULT: 4
        (maximum allowed value of nMerge is nHop-1)
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.OUT,
          description='the value of nMerge, the numbder of particles catalogue group boundaries')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            nMerge was set
        -1 - ERROR
            requested value for nMerge was too high
        '''
        return function
        
    @legacy_function
    def set_density_method():
        '''
        Set the density calculation method used by calculate_densities().
        0 - gather-scatter cubic spline kernel (DEFAULT)
        1 - gather-only cubic spline kernal
        2 - tophat kernal
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.IN,
          description='value representing the density calculation method (DEFAULT: 0; gather scatter cubic spline)')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            density method set
        -1 - ERROR
            invalid density method requested
        '''
        return function

    @legacy_function
    def get_density_method():
        '''
        Get the density calculation method used by calculate_densities().
        0 - gather-scatter cubic spline kernel (DEFAULT)
        1 - gather-only cubic spline kernal
        2 - tophat kernal
        '''
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype = 'int32', direction=function.OUT,
          description='value representing the density calculation method (DEFAULT: 0; gather scatter cubic spline)')
        function.result_type = 'int32'
        function.result_doc = '''
        0 - OK
            density method set
        -1 - ERROR
            invalid density method requested
        '''
        return function
        
    @legacy_function
    def show_parameters():
        '''
        Print all parameters and their values.
        '''
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    def commit_particles(self):
        pass
    
    def recommit_particles(self):
        pass
    


class Hop(CommonCode):

    def __init__(self, unit_converter = None, **options):
    
        self.unit_converter = unit_converter
        
        InCodeComponentImplementation.__init__(self,  HopInterface(**options), **options)
    
    def define_converter(self, handler):
        if self.unit_converter is None:
            return
        handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())
    
    def define_errorcodes(self, handler):
        handler.add_errorcode(-1, 'Something went wrong...')
        handler.add_errorcode(-2, 'Not implemented.')
        handler.add_errorcode(-3, 'A particle with the given index was not found.')
        handler.add_errorcode(-4, 'Negative density encountered.')
        handler.add_errorcode(-5, 'Too few particles.')
    
    def define_methods(self, builder):
        
        builder.add_method(
            'new_particle',
            (generic_unit_system.mass, generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
            (builder.INDEX, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'delete_particle',
            (builder.INDEX,),
            (builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_position',
            (builder.INDEX,),
            (generic_unit_system.length, generic_unit_system.length, generic_unit_system.length, builder.ERROR_CODE),
            public_name = 'get_position'
        )
        
        builder.add_method(
            'set_position',
            (builder.INDEX, generic_unit_system.length, generic_unit_system.length, generic_unit_system.length,),
            (builder.ERROR_CODE),
            public_name = 'set_position'
        )
        
        builder.add_method(
            'get_density',
            (builder.INDEX,),
            (generic_unit_system.density, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_mass',
            (builder.INDEX,),
            (generic_unit_system.mass, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'set_density',
            (builder.INDEX, generic_unit_system.density,),
            (builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_densest_neighbor',
            (builder.INDEX,),
            (builder.INDEX, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_group_id',
            (builder.INDEX,),
            (builder.NO_UNIT, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_number_of_particles',
            (),
            (builder.INDEX, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'calculate_densities',
            (),
            (builder.ERROR_CODE)
        )
        
        builder.add_method(
            'do_hop',
            (),
            (builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_densest_particle_in_group',
            (builder.INDEX,),
            (builder.INDEX, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_number_of_particles_in_group',
            (builder.INDEX,),
            (builder.INDEX, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_average_boundary_density_of_groups',
            (builder.INDEX, builder.INDEX,),
            (generic_unit_system.density, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_number_of_groups',
            (),
            (builder.NO_UNIT, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'get_number_of_particles_outside_groups',
            (),
            (builder.NO_UNIT, builder.ERROR_CODE)
        )
        
        builder.add_method(
            'show_parameters',
            (),
            (builder.ERROR_CODE)
        )
        
        builder.add_method(
            "get_nHop",
            (),
            (builder.NO_UNIT, builder.ERROR_CODE,)
        )
        
        builder.add_method(
            "set_nHop",
            (builder.NO_UNIT, ),
            (builder.ERROR_CODE,)
        )
        
        builder.add_method(
            "get_nDens",
            (),
            (builder.NO_UNIT, builder.ERROR_CODE,)
        )
        
        builder.add_method(
            "set_nDens",
            (builder.NO_UNIT, ),
            (builder.ERROR_CODE,)
        )
        
        builder.add_method(
            "get_nBucket",
            (),
            (builder.NO_UNIT, builder.ERROR_CODE,)
        )
        
        builder.add_method(
            "set_nBucket",
            (builder.NO_UNIT, ),
            (builder.ERROR_CODE,)
        )
        
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_fDensThresh", 
            "set_fDensThresh",
            "outer_density_threshold", 
            "the density below which particles are not assigned to any group", 
            default_value = -1.0 | generic_unit_system.density
        )
        handler.add_method_parameter(
            "get_saddle_densthresh", 
            "set_saddle_densthresh",
            "saddle_density_threshold", 
            "For two groups to merge, the density at their boundary must exceed this "
            "threshold (DEFAULT: 1.75 * outer_density_threshold)",
            default_value = -1.0 | generic_unit_system.density
        )
        handler.add_method_parameter(
            "get_peak_densthresh", 
            "set_peak_densthresh",
            "peak_density_threshold", 
            "Groups with density below this threshold (fringe groups) are attached to "
            "other (proper) groups, or dropped (DEFAULT: max(saddle_density_threshold, 2.0 * outer_density_threshold))",
            default_value = -1.0 | generic_unit_system.density
        )
        handler.add_method_parameter(
            "get_saddle_density_threshold_factor", 
            "set_saddle_density_threshold_factor",
            "saddle_density_threshold_factor", 
            "For two groups to merge, the density at their boundary must exceed this "
            "factor times the lowest of the two peak densities (set relative_saddle_density_threshold to True)",
            default_value = 0.80
        )
        handler.add_boolean_parameter(
            "get_relative_saddle_density_threshold",
            "set_relative_saddle_density_threshold",
            "relative_saddle_density_threshold",
            "Flag to use a saddle-density-threshold relative to the lowest peak density, "
            "instead of the absolute saddle_density_threshold",
            False
        )
        
        handler.add_method_parameter(
            "get_nHop", 
            "set_nHop",
            "number_of_neighbors_for_hop", 
            "The number of neighbors to search for local density maximum (search is performed iteratively until a true maximum is found)", 
            default_value = 64
        )

        handler.add_method_parameter(
            "get_density_method", 
            "set_density_method",
            "density_method", 
            "method used for density computation (0,1,2 = gather-scatter spline, gather spline, tophat)", 
            default_value = 0
        )


        handler.add_method_parameter(
            "get_nDens", 
            "set_nDens",
            "number_of_neighbors_for_local_density", 
            "Return the number of particles to smooth over when calculating densities.", 
            default_value = 64
        )
        
        handler.add_method_parameter(
            "get_nBucket", 
            "set_nBucket",
            "number_of_buckets", 
            "Return the bucket parameter to tune the performance of the kd-tree search.", 
            default_value = 16
        )
        
        handler.add_method_parameter(
            "get_nMerge", 
            "set_nMerge",
            "number_of_particles_per_group_pair_boundary", 
            "The number of (densest) particles per boundary between each pair of groups, for merging.", 
            default_value = 4
        )
  
    def define_particle_sets(self, builder):
        builder.define_set('particles', 'index_of_the_particle')
        builder.set_new('particles', 'new_particle')
        builder.set_delete('particles', 'delete_particle')
        builder.add_setter('particles', 'set_position')
        builder.add_setter('particles', 'set_density', names=('density',))
        builder.add_getter('particles', 'get_position')
        builder.add_getter('particles', 'get_density', names=('density',))
        builder.add_getter('particles', 'get_mass', names=('mass',))
        #builder.add_getter('particles', 'get_densest_neighbor')
        builder.add_getter('particles', 'get_group_id', names=('group_id',))
    
    def define_state(self, handler):
        CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        handler.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        handler.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('EDIT', 'before_get_parameter')
        handler.add_method('UPDATE','before_get_parameter')
        
        handler.add_method('EDIT', 'new_particle')
        handler.add_method('EDIT', 'delete_particle')
        handler.add_method('UPDATE', 'new_particle')
        handler.add_method('UPDATE', 'delete_particle')
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('RUN', 'UPDATE', 'new_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        handler.add_method('RUN', 'calculate_densities')
        handler.add_method('RUN', 'do_hop')
        handler.add_method('RUN', 'get_mass')
        handler.add_method('RUN', 'get_position')
        handler.add_method('RUN', 'get_density')
        handler.add_method('RUN', 'get_group_id')
    
    def no_group(self):
        return self.particles.select(lambda group_id: group_id == -1, ["group_id"])
    
    def groups(self):
        number_of_groups = self.get_number_of_groups()
        group_id = self.particles.group_id
        
        for index in range(number_of_groups):
            result = self.particles[group_id == index]
            result.add_function_attribute('id_of_group', lambda particles: particles[0].group_id)
            #result.add_function_attribute('get_density_of_group', lambda particles: self.get_group_density(particles.id_of_group()))
            yield result










  
