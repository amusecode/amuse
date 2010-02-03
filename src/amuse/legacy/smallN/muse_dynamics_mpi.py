from amuse.legacy import *
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.support.data.core import Particle, Particles, ParticlesWithUnitsConverted
from amuse.support.data import binding
from amuse.legacy.interface.gd import NBodyGravitationalDynamicsBinding

class SmallNInterface(LegacyInterface):
    """
        Interface to the Kira Small-N Integrator and Kepler modules from
        Starlab.  http://www.ids.ias.edu/~starlab/
        
        You will need to download Starlab from the above site, make it, install
        it, and then set the STARLAB_INSTALL_PATH variable to be equal to the
        installation directory (typically something like ~/starlab/usr).

        Starlab is available under the GNU General Public Licence (version 2),
        and is developed by:
            * Piet Hut
            * Steve McMillan
            * Jun Makino
            * Simon Portegies Zwart
        Other Starlab Contributors:
            * Douglas Heggie
            * Kimberly Engle
            * Peter Teuben 
            
    """
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
        function.addParameter('id', 'i', function.OUT)
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
        return function

    @legacy_function
    def integrate_multiple():
        function = LegacyFunctionSpecification()
        function.addParameter('start_time', 'd', function.IN)
        function.addParameter('end_time', 'd', function.OUT)
        function.addParameter('verbose', 'i', function.IN)
        function.addParameter('eps2', 'd', function.IN)
        return function;

    def delete_particle(self, id):
        return -1   # -1 == Not yet implemented.

    def get_time(self):
        return self.convert_nbody.to_si(self.time | nbody_system.time)

    def set_time(self, new_time_value):
        self.time = self.convert_nbody.to_nbody(new_time_value).number

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

    @legacy_function
    def get_number_of_particles():
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', 'i', function.OUT)
        return function;

    def evolve(self, verbose=False, super_verbose=False):
        """ Evolves the system until a stable configuration is reached. """
        self.has_run = True
        verbose_int = 0
        if verbose:
            verbose_int = 1
        if super_verbose:
            verbose_int = 100
        end_time = self.integrate_multiple(self.time, verbose_int, self.eps2)
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
        return self.add_to_interaction(mass, x, y, z, vx, vy, vz)

    def new_particle(self, mass, radius, x, y, z, vx, vy, vz):
        newp = Particle()
        newp.mass = mass | nbody_system.kg
        newp.radius = radius | nbody_system.length
        newp.x = x | nbody_system.length
        newp.y = y | nbody_system.length
        newp.z = z | nbody_system.length
        newp.vx = vx | nbody_system.speed
        newp.vy = vy | nbody_system.speed
        newp.vz = vz | nbody_system.speed
        return self.add_particle(newp)

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

    def get_state(self, index):
        if self.has_run:
            (key, mass, x, y, z, vx, vy, vz) = self.get_particle_result(index)
        else:
            (key, mass, x, y, z, vx, vy, vz) = self.get_particle_original(index)
        radius = 0.0
        return (mass, radius, x, y, z, vx, vy, vz), 0   # The last zero is the OK code.

class SmallNInCodeAttributeStorage(InCodeAttributeStorage):
    new_particle_method = binding.NewParticleMethod(
        "new_particle", 
        (
            ("mass", "mass", nbody_system.mass),
            ("radius", "radius", nbody_system.length),
            ("x", "x", nbody_system.length),
            ("y", "y", nbody_system.length),
            ("z", "z", nbody_system.length),
            ("vx", "vx", nbody_system.speed),
            ("vy", "vy", nbody_system.speed),
            ("vz", "vz", nbody_system.speed),
        )
    )

    getters = (
        binding.ParticleGetAttributesMethod(
            "get_state",
            (
                ("mass", "mass", nbody_system.mass),
                ("radius", "radius", nbody_system.length),
                ("x", "x", nbody_system.length),
                ("y", "y", nbody_system.length),
                ("z", "z", nbody_system.length),
                ("vx", "vx", nbody_system.speed),
                ("vy", "vy", nbody_system.speed),
                ("vz", "vz", nbody_system.speed),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_mass",
            (
                ("mass", "mass", nbody_system.mass),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_radius",
            (
                ("radius", "radius", nbody_system.length),
            )
        ),
        binding.ParticleGetAttributesMethod(
            "get_position",
            (
                ("x", "x", nbody_system.length),
                ("y", "y", nbody_system.length),
                ("z", "z", nbody_system.length),
            )
        ),
    )

class SmallNBinding(NBodyGravitationalDynamicsBinding):
    parameter_definitions = [
         parameters.ModuleAttributeParameterDefinition(
            "eps2",
            "epsilon_squared", 
            "smoothing parameter for gravity calculations", 
            nbody_system.length * nbody_system.length, 
            0.0 | nbody_system.length * nbody_system.length
        ),
        parameters.ModuleAttributeParameterDefinition(
            "number_of_particles",
            "number_of_particles", 
            "The number of particles being managed by the SmallN module", 
            units.none, 
            0 | units.none
        )
    ]
   
    attribute_definitions = [ ]

    def __init__(self, convert_nbody = None):
        NBodyGravitationalDynamicsBinding.__init__(self, convert_nbody)
       
        self.nbody_particles = Particles(storage =  SmallNInCodeAttributeStorage(self))
        self.particles = ParticlesWithUnitsConverted(self.nbody_particles, self.convert_nbody.as_converter_from_si_to_nbody())
   
    def current_model_time(self):
        return self.convert_nbody.to_si( self.t | nbody_system.time)
           
   
    def evolve_model(self, time_end):
        result = self.evolve(self.convert_nbody.to_nbody(time_end).value_in(nbody_system.time))
        return result

class SmallN(SmallNInterface, SmallNBinding):
    def __init__(self, convert_nbody = None):
        SmallNInterface.__init__(self)
        SmallNBinding.__init__(self, convert_nbody)

