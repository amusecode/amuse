import numpy
import os.path
from amuse.community import *
from amuse.datamodel import Particles
from amuse import datamodel
from amuse.units import nbody_system

class FractalClusterInterface(CodeInterface,  LiteratureReferencesMixIn):
    """
    makes fractal of nstar particles of dimension fdim, using ndiv 
    subunits forcing the number of cells if force=.true.
    
    reference:
        .. [#] Simon Goodwin & Ant Whitworth (2004, A&A, 413, 929)

    """

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'fractal_worker'

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('nstar', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def generator():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def get_fractal_dimension():
        function = LegacyFunctionSpecification()
        function.addParameter('fdim', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function
    def set_fractal_dimension():
        function = LegacyFunctionSpecification()
        function.addParameter('fdim', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_random_seed():
        function = LegacyFunctionSpecification()
        function.addParameter('seed', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function
    def set_random_seed():
        function = LegacyFunctionSpecification()
        function.addParameter('seed', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_nstar():
        function = LegacyFunctionSpecification()
        function.addParameter('nstar', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function
    def set_nstar():
        function = LegacyFunctionSpecification()
        function.addParameter('nstar', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

class MakeFractalCluster(object):
    def __init__(self,N=None,masses=None,convert_nbody = None,do_scale=False,
                     random_seed=None,fractal_dimension=1.6,virial_ratio=0.5,verbose=False):
        if N is None and masses is not None:
          self.N=len(masses)
        if masses is not None and len(masses)!=N:
          print "warning: provided mass array not equal to masses"          
        self.N=N
        self.convert_nbody=convert_nbody
        self.do_scale=do_scale
        self.random_seed=random_seed
        self.masses=masses
        self.fractal_dimension=fractal_dimension
        self.virial_ratio=virial_ratio
        self.verbose=verbose

    def new_model(self):
        if self.verbose:
          frac=FractalClusterInterface(redirection="none")
        else:  
          frac=FractalClusterInterface()
        frac.set_nstar(self.N)
        frac.set_fractal_dimension(self.fractal_dimension)
        if self.random_seed is not None:
            frac.set_random_seed(self.random_seed)
        frac.generator()
        
        x,y,z,vx,vy,vz,err=frac.get_state(range(1,self.N+1))

        if self.masses is None:
            masses=numpy.array([0.]*len(x))+1./self.N
        else:
            masses=self.masses/self.masses.sum()
        positions =  numpy.hstack((x,y,z))
        velocities =  numpy.hstack((vx,vy,vz))
        
        return masses,positions,velocities

    @property
    def result(self):
        masses, positions, velocities = self.new_model()
        result = datamodel.Particles(self.N)
        result.mass = nbody_system.mass.new_quantity(numpy.hstack(masses))
        result.position = nbody_system.length.new_quantity(positions)
        result.velocity = nbody_system.speed.new_quantity( velocities)
        result.radius = 0 | nbody_system.length

        result.move_to_center()

        ep=result.potential_energy(G=nbody_system.G)
        ek=result.kinetic_energy()
        fac=(self.virial_ratio*abs(ep)/ek)**0.5
        result.velocity*=fac

        if self.do_scale:
            result.scale_to_standard()

        if not self.convert_nbody is None:
            result = datamodel.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_generic())
            result = result.copy()

        return result

def new_fractal_cluster_model(**keyword_arguments):
    """
    create a fractal stellar distribution

    :argument N: Number of particles, if None then masses must be provided [None]
    :argument masses: optional masses, if not provided masses will be 1/N [None]
    :argument fractal_dimension: fractal dimension of distribution, between 1.6 and 3 [1.6]
    :argument random_seed: random seed [None]
    :argument convert_nbody:  When given will convert the resulting set to SI units
    :argument do_scale: scale the result to exact nbody units (M=1, K=0.25, U=-0.5)
    :argument virial_ratio: ratio of kinetic to potential energy [0.5]
    """
    uc = MakeFractalCluster(**keyword_arguments)
    return uc.result
