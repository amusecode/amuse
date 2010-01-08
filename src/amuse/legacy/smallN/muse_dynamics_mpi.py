from amuse.legacy import *
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data.core import Particle

class SmallN(LegacyInterface):
    include_headers = ['amuse_interface.h']

    parameter_definitions = [
        parameters.ModuleMethodParameterDefinition(
            "get_eps", "set_eps",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.0 | nbody_system.length * nbody_system.length
        ),
    ]


    def __init__(self, convert_nbody = None):
        LegacyInterface.__init__(self, name_of_the_worker='muse_worker_starlab')
        if convert_nbody is None:
            convert_nbody = nbody_system.nbody_to_si.get_default()
        self.convert_nbody = convert_nbody
        self.parameters = parameters.Parameters(self.parameter_definitions, self)
        self.has_run = False
        self.eps2 = 0.0
        self.time = 0.0

    @legacy_function
    def report_multiples():
        function = LegacyFunctionSpecification()
        function.addParameter('level', 'i', function.IN)
        return function;

    @legacy_function   
    def get_total_energy():
        function = LegacyFunctionSpecification()  
        function.result_type = 'd'
        return function;

    @legacy_function
    def add_to_interaction():
        function = LegacyFunctionSpecification()
        function.addParameter('i', 'i', function.IN)    # "i" is an ID, not an index
        for x in ['m', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, 'd', function.IN)
        return function;

    @legacy_function
    def get_particle_result():
        function = LegacyFunctionSpecification()
        function.addParameter('k', 'i', function.IN)
        function.addParameter('id', 'i', function.OUT)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, 'd', function.OUT)
        return function;

    @legacy_function
    def get_particle_original():
        function = LegacyFunctionSpecification()
        function.addParameter('k', 'i', function.IN)
        function.addParameter('id', 'i', function.OUT)
        for x in ['mass', 'x', 'y', 'z', 'vx', 'vy', 'vz']:
            function.addParameter(x, 'd', function.OUT)
        return function;


    @legacy_function
    def integrate_multiple():
        function = LegacyFunctionSpecification()
        function.addParameter('end_time', 'd', function.OUT)
        function.addParameter('verbose', 'i', function.IN)
        function.addParameter('eps2', 'd', function.IN)
        return function;

    def get_eps2(self):
        return self.eps2

    def set_eps2(self, new_eps2_value):
        self.eps2 = new_eps2_value

    def get_time(self):
        return self.convert_nbody.to_si(self.time | nbody_system.time)

    @legacy_function
    def add_binary():
        function = LegacyFunctionSpecification()
        function.addParameter('id1', 'i', function.IN)
        function.addParameter('id2', 'i', function.IN)
        function.addParameter('mass1', 'd', function.IN)
        function.addParameter('mass2', 'd', function.IN)
        function.addParameter('period', 'd', function.IN)
        function.addParameter('eccentricity', 'd', function.IN)
        return function;

    @legacy_function
    def clear_multiple():
        function = LegacyFunctionSpecification()
        return function;


    def evolve(self, verbose=False, super_verbose=False):
        """ Evolves the system until a stable configuration is reached. """
        self.has_run = True
        verbose_int = 0
        if verbose:
            verbose_int = 1
        if super_verbose:
            verbose_int = 100
        end_time = self.integrate_multiple(verbose_int, self.eps2)
        self.time = end_time

    def add_particle(self, particle):
        """ Adds the specified Particle to the simulation. """
        assert isinstance(particle, Particle), "You must pass in a Particle object."
        mass = self.convert_nbody.to_nbody(particle.mass).number
        x = self.convert_nbody.to_nbody(particle.x).number
        y = self.convert_nbody.to_nbody(particle.y).number
        z = self.convert_nbody.to_nbody(particle.z).number
        vx = self.convert_nbody.to_nbody(particle.vx).number
        vy = self.convert_nbody.to_nbody(particle.vy).number
        vz = self.convert_nbody.to_nbody(particle.vz).number
        self.add_to_interaction(particle.key, mass, x, y, z, vx, vy, vz)

    def new_particle(self, mass, radius, x, y, z, vx, vy, vz):
        particle_index = -1
        #TODO

    def get_particle_by_index(self, index):
        """ Returns a Particle by index.  """
        if self.has_run:
            (key, mass, x, y, z, vx, vy, vz) = self.get_particle_result(index)
        else:
            (key, mass, x, y, z, vx, vy, vz) = self.get_particle_original(index)

        p = Particle(key)
        p.mass = self.convert_nbody.to_si(mass | nbody_system.mass)
        p.x = self.convert_nbody.to_si(x | nbody_system.length)
        p.y = self.convert_nbody.to_si(y | nbody_system.length)
        p.z = self.convert_nbody.to_si(z | nbody_system.length)
        nbody_speed = nbody_system.length / nbody_system.time
        p.vx = self.convert_nbody.to_si(vx | nbody_speed)
        p.vy = self.convert_nbody.to_si(vy | nbody_speed)
        p.vz = self.convert_nbody.to_si(vz | nbody_speed)
        return p

    def reset_close_encounter(self):
        """ Resets the internal variables so that a new close encounter can be
        run.  This method should be called once before every close encounter to
        ensure that no data remains from the previous close encounter. """
        self.clear_multiple()
        self.time = 0.0

class BinaryStar(Particle):
    def __init__(self, key, component_particle1, component_particle2, period, \
                 eccentricity, particles_set=None, **keyword_arguments):
        Particle.__init__(self, key, particles_set, *keyword_arguments)
        self.component_particle1 = component_particle1
        self.component_particle2 = component_particle2
        self.period = period
        self.eccentricity = eccentricity

