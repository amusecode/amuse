import numpy
from amuse.units import constants, units
from amuse.datamodel import Grid


class SPH2StellarModel(object):
    """
    Converts a set of SPH particles to a 1D stellar evolution model. If the SPH 
    model of the star included a (non-SPH) 'core' particle, supply it via the 
    optional core_particle keyword argument.
    
    SPH particles are sorted using the pressure
    Useful for continuing the stellar evolution of merged stars.
    
    :argument sph_particles: The SPH particles to be converted to a stellar model
    :argument core_particle: Gravitational particle representing the stellar core (optional)
    """
    
    def __init__(self, sph_particles, core_particle = None, particles_per_zone=1):
        self.sph_particles = sph_particles
        self.core_particle = core_particle
        self.particles_per_zone = particles_per_zone
    
    def derive_stellar_structure(self):
        sorted = self.sph_particles.pressure.argsort()[::-1]
        binned = sorted.reshape((-1, self.particles_per_zone))
        
        stellar_model = Grid(binned.shape[0])
        stellar_model.dmass = self.sph_particles.mass[binned].sum(axis=1)
        stellar_model.mass = stellar_model.dmass.accumulate()
        stellar_model.pressure= self.sph_particles.pressure[binned].sum(axis=1)
        stellar_model.rho = stellar_model.dmass / (self.sph_particles.mass / self.sph_particles.density)[binned].sum(axis=1)
        stellar_model.radius = ((3 / (4 * numpy.pi)) * stellar_model.dmass / stellar_model.rho).accumulate()**(1.0/3.0) * 1
        stellar_model.temperature = ((self.sph_particles.mass * self.sph_particles.u * self.sph_particles.mu)[binned].sum(axis=1) / 
            (1.5 * constants.kB * stellar_model.dmass)).as_quantity_in(units.K)
        zeros = numpy.zeros(len(stellar_model.dmass))
        stellar_model.luminosity = zeros - 1 | units.LSun
        
        attribute_names = self.sph_particles.get_attribute_names_defined_in_store()
        for attribute, name in [("h1", "X_H"), ("he4", "X_He"), ("c12", "X_C"), ("n14", "X_N"), 
                ("o16", "X_O"), ("ne20", "X_Ne"), ("mg24", "X_Mg"), ("si28", "X_Si"), ("fe56", "X_Fe")]:
            if attribute in attribute_names:
                setattr(stellar_model, name, (self.sph_particles.mass * getattr(self.sph_particles, attribute))[binned].sum(axis=1) / stellar_model.dmass)
            else:
                setattr(stellar_model, name, zeros)
        
        return stellar_model
    

def convert_SPH_to_stellar_model(sph_particles, **keyword_arguments):
    """
    Converts a set of SPH particles to a 1D stellar evolution model. If the SPH 
    model of the star included a (non-SPH) 'core' particle, supply it via the 
    optional core_particle keyword argument (Not yet supported).
    
    SPH particles are sorted using the pressure.
    Useful for continuing the stellar evolution of merged stars.
    
    :argument sph_particles: The SPH particles to be converted to a stellar model
    :argument core_particle: Gravitational particle representing the stellar core (optional)
    :argument particles_per_zone: The number of sph particles within each mesh cell (default=1)
    """
    converter = SPH2StellarModel(sph_particles, **keyword_arguments)
    return converter.derive_stellar_structure()

