from amuse.community import *

class CraterInterface(CodeInterface, LiteratureReferencesMixIn):
    """
    Calculating the size of a crater as a function of impact and target objects.
        .. [#] ... Collins, G. S.; Melosh, H. J.; Marcus, R. A. 2005, M&PS 40, 718 [2005M&PS...40..817C]
    """
    
    use_modules = ['CraterRadius']
    
    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker="crater_worker", **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)
        

    @legacy_function
    def set_target_type():
        """
        Set the target_type value. 
        """
        function = LegacyFunctionSpecification()
        function.addParameter('target_type', dtype='int32',
                              direction=function.IN,
                              description="target type [1, 2 or 3]")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_target_type():
        """
        Get the target type value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('target_type', dtype='int32',
                              direction=function.OUT,
                              description="value of the target type")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def set_target_gravity():
        """
        Set the target_gravity value
        """
        function = LegacyFunctionSpecification()
        function.addParameter('target_gravity', dtype='float64',
                              direction=function.IN,
                              description="target gravity")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_target_gravity():
        """
        Get the target_gravity value
        """
        function = LegacyFunctionSpecification()
        function.addParameter('target_gravity', dtype='float64',
                              direction=function.OUT,
                              description="target gravity")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def set_projectile_type():
        """
        Set the projectile_type value. 
        """
        function = LegacyFunctionSpecification()
        function.addParameter('projectile_type', dtype='int32',
                              direction=function.IN,
                              description="projectile type [1, 2 or 3]")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_projectile_type():
        """
        Get the projectile type value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('projectile_type', dtype='int32',
                              direction=function.OUT,
                              description="value of the projectile type")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    
    @legacy_function
    def set_projectile_density():
        """
        Set the target diatemter value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('projectile_density', dtype='float64',
                              direction=function.IN,
                              description="Projectile density")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_projectile_density():
        """
        Get the projectile_density value
        """
        function = LegacyFunctionSpecification()
        function.addParameter('projectile_density', dtype='float64',
                              direction=function.OUT,
                              description="target gravity")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_target_density():
        """
        Get the target density value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('target_density', dtype='float64',
                              direction=function.OUT,
                              description="value of the target density")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def set_projectile_diameter():
        """
        Set the target diatemter value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('projectile_diameter', dtype='float64',
                              direction=function.IN,
                              description="Projectile diameter")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_projectile_diameter():
        """
        Get the projectile diameter value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('projectile_diameter', dtype='float64',
                              direction=function.OUT,
                              description="value of the projectile diameter")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def set_target_density():
        """
        Set the target diatemter value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('target_density', dtype='float64',
                              direction=function.IN,
                              description="Target density")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def set_impact_angle():
        """
        Set the impact_angle value in degrees.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('impact_angle', dtype='float64',
                              direction=function.IN,
                              description="impact angle in degrees")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_impact_angle():
        """
        Get the projectile diameter value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('impact_angle', dtype='float64',
                              direction=function.OUT,
                              description="value of the projectile diameter")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    
    @legacy_function
    def set_impact_velocity():
        """
        Set the impact_velocity value in degrees.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('impact_velocity', dtype='float64',
                              direction=function.IN,
                              description="impact velocity in degrees")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    @legacy_function
    def get_impact_velocity():
        """
        Get the impact velocity
        """
        function = LegacyFunctionSpecification()
        function.addParameter('impact_velocity', dtype='float64',
                              direction=function.OUT,
                              description="value of the impact velocity")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_crater_diameter():
        """
        Get the crater diatmeter value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('crater_diameter', dtype='float64',
                              direction=function.OUT,
                              description="value of the crater diameter")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_crater_type():
        """
        Get the crater diatmeter value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('crater_type', dtype='int32',
                              direction=function.OUT,
                              description="value of the crater type")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_formation_time():
        """
        Get the projectile diameter value.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('formation_time', dtype='float64',
                              direction=function.OUT,
                              description="time-scale for crater formation")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function
    
    
    @legacy_function
    def crater_radius():
        function = LegacyFunctionSpecification()  
        function.addParameter('rhoproj', dtype='float64', direction=function.IN)
        function.addParameter('L', dtype='float64', direction=function.IN)
        function.addParameter('v', dtype='float64', direction=function.IN)
        function.addParameter('theta', dtype='float64', direction=function.IN)
        function.addParameter('rhotarget', dtype='float64', direction=function.IN)
        function.addParameter('g', dtype='float64', direction=function.IN)

        function.addParameter('targtype', dtype='int32', direction=function.IN)

        function.addParameter('cratertype', dtype='int32', direction=function.OUT)
        function.addParameter('Dyield', dtype='float64', direction=function.OUT)
        function.addParameter('Dgault', dtype='float64', direction=function.OUT)
        function.addParameter('Tform', dtype='float64', direction=function.OUT)
        function.addParameter('crater_diameter', dtype='float64', direction=function.OUT)
        function.addParameter('Lyield', dtype='float64', direction=function.OUT)
        function.addParameter('Lgault', dtype='float64', direction=function.OUT)

        function.result_type = 'int32'
        function.can_handle_array = True
        return function
        

class Crater(InCodeComponentImplementation):

    def __init__(self, **options):
        InCodeComponentImplementation.__init__(self,  CraterInterface(**options), **options)

    def initialize_code(self):
        print "ahhh."

    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_target_density",
            "set_target_density",
            "target_density",
            "initialize the density of the target",
            default_value = 1.0 | units.kg/units.m**3
            )
        handler.add_method_parameter(
            "get_target_type",
            "set_target_type",
            "target_type",
            "initialize the type of target",
            default_value = 1
            )
        handler.add_method_parameter(
            "get_target_gravity",
            "set_target_gravity",
            "target_gravity",
            "initialize the gravity of the target",
            default_value = 9.8 | units.m/units.s**2
            )
        handler.add_method_parameter(
            "get_projectile_density",
            "set_projectile_density",
            "projectile_density",
            "initialize the density of the target",
            default_value = 1
            )
        handler.add_method_parameter(
            "get_projectile_diameter",
            "set_projectile_diameter",
            "projectile_diameter",
            "initialize the diameter of the target",
            default_value = 1
            )
        handler.add_method_parameter(
            "get_impact_angle",
            "set_impact_angle",
            "impact_angle",
            "initialize the impact angle",
            default_value = 1
            )
        handler.add_method_parameter(
            "get_impact_velocity",
            "set_impact_velocity",
            "impact_velocity",
            "initialize the impact velocity",
            default_value = 1
            )
        handler.add_method_parameter(
            "get_crater_diameter",
            "set_crater_diameter",
            "crater_diameter",
            "initialize the crater diameter",
            default_value = 1
            )
        handler.add_method_parameter(
            "get_crater_type",
            "set_crater_type",
            "crater_type",
            "initialize the crater type",
            default_value = 1
            )
        handler.add_method_parameter(
            "get_formation_time",
            "set_formation_time",
            "formation_time",
            "initialize the formation time",
            default_value = 1
            )
        
    def define_methods(self, handler):
        handler.add_method(
            "crater_radius",
            (units.kg/units.m**3, units.km, units.kms, units.deg, units.kg/units.m**3, units.m/units.s**2, handler.NO_UNIT),
            (handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, handler.NO_UNIT, units.km, handler.NO_UNIT, handler.NO_UNIT, handler.ERROR_CODE,)
            )
        handler.add_method(
            "set_target_density",
            (units.kg/units.m**3, ),
            (handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_target_density",
            ( ),
            (units.kg/units.m**3, handler.ERROR_CODE,)
            )
        handler.add_method(
            "set_target_type",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_target_type",
            ( ),
            (handler.NO_UNIT, handler.ERROR_CODE,)
            )
        handler.add_method(
            "set_target_gravity",
            (units.m/units.s**2),
            (handler.ERROR_CODE,)
            )
        handler.add_method(
            "set_projectile_type",
            (handler.NO_UNIT, ),
            (handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_projectile_type",
            ( ),
            (handler.NO_UNIT, handler.ERROR_CODE,)
            )
        handler.add_method(
            "set_projectile_density",
            (units.kg/units.m**3),
            (handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_projectile_density",
            ( ),
            (units.kg/units.m**3, handler.ERROR_CODE,)
            )
        handler.add_method(
            "set_projectile_diameter",
            (units.m),
            (handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_projectile_diameter",
            ( ),
            (units.m, handler.ERROR_CODE,)
            )
        handler.add_method(
            "set_impact_angle",
            (units.deg),
            (handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_impact_angle",
            ( ),
            (units.deg, handler.ERROR_CODE,)
            )
        handler.add_method(
            "set_impact_velocity",
            (units.kms),
            (handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_impact_velocity",
            ( ),
            (units.kms, handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_crater_diameter",
            ( ),
            (units.m, handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_crater_type",
            ( ),
            (handler.NO_UNIT, handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_formation_time",
            ( ),
            (units.s, handler.ERROR_CODE,)
            )
        handler.add_method(
            "get_target_gravity",
            ( ),
            (units.m/units.s**2, handler.ERROR_CODE,)
            )

    def evolve_model(self, end_time=None):
        L = self.get_projectile_diameter()
        g = self.get_target_gravity()
        print "g=", g
        rhoproj = self.get_projectile_density()
        rhotarget = self.get_target_density()
        #target_diameter = self.get_target_diameter()
        targtype = self.get_target_type()
        theta = self.get_impact_angle()
        v = self.get_impact_velocity()
        print "xx=", L, g, rhoproj, rhotarget
        print "yy=", targtype, theta, v

        self.crater_radius(rhoproj, L, v, theta, rhotarget, g, targtype)
        evolve_success = 1
        return evolve_success
