import math
import numpy

from operator import itemgetter

from amuse.community.mesa.interface import MESA
from amuse.support.exceptions import AmuseException
from amuse.support.exceptions import AmuseWarning
from amuse.units import constants
from amuse.units import units

from amuse.datamodel import Particles
from amuse.datamodel.particles import ParticlesSubset
from amuse.datamodel.particles import ParticlesSuperset
class RocheLobeOverflow(object):
    """
    Applies Roche-lobe overflow to a set of stars
    Currently only MESA particles can be added.
    """
    def __init__(self, accretion_efficiency = 1.0, dynamics_code = None):
        self.accretion_efficiency = accretion_efficiency
        self.particles  = []
        self.companions = []
        self.overflow_radii = [] | units.RSun
        self.dynamics_code = dynamics_code
    
    def add_particle(self, particle, overflow_radius, companion = None):
        self.particles.append(particle)
        self.companions.append(companion)
        self.overflow_radii.append(overflow_radius)
    
    def add_particles(self, particles, overflow_radii, companions = None):
        if companions is None:
            companions = [None]*len(particles)
        for particle, companion, overflow_radius in zip(particles, companions, overflow_radii):
            self.add_particle(particle, overflow_radius, companion)
    
    def update_roche_radii(self):
        if self.dynamics_code is None:
            raise AmuseException("Specifying a dynamics_code is required for automatic updating of the overflow radii. Otherwise use set_roche_radii.")
        if None in self.companions:
            raise AmuseException("Specifying companion stars is required for automatic updating of the overflow radii. Otherwise use set_roche_radii.")
        for i, (part, comp) in enumerate(zip(self.particles, self.companions)):
            dyn_particles = ParticlesSuperset([part.as_set(), comp.as_set()]).get_intersecting_subset_in(self.dynamics_code.particles)
            separation = (dyn_particles.position[0] - dyn_particles.position[1]).length()
            cube_root_q = (part.mass / comp.mass).value_in(units.none)**(1/3.0)
            # Eggleton 1983:
            self.overflow_radii[i] = separation * 0.49 * cube_root_q**2 / \
                (0.6 * cube_root_q**2 + math.log(1 + cube_root_q))
    
    def set_roche_radii(self, particles, radii):
        for r, particle in zip(radii, particles):
            for i, local_part in enumerate(self.particles):
                if particle.key == local_part.key:
                    self.overflow_radii[i] = r
                    break
            else:
                raise AmuseException("A particle (key: {0}) wasn't found.".format(particle.key))
    
    def do_roche_lobe_overflow(self):
        result = [] | units.MSun
        for part, comp, r in zip(self.particles, self.companions, self.overflow_radii):
            mass_lost = do_roche_lobe_overflow(part, r)
            if comp:
                orig_comp_mass = comp.mass
                new_comp_mass  = orig_comp_mass + mass_lost * self.accretion_efficiency
                comp.mass = new_comp_mass
                if hasattr(comp, "set_density_profile"):
                    # rescale densities
                    mass_scaling = new_comp_mass / orig_comp_mass
                    rho_profile = comp.get_density_profile()
                    rho_profile *= mass_scaling
                    comp.set_density_profile(rho_profile)
            result.append(mass_lost)
        if self.dynamics_code:
            if None in self.companions:
                raise AmuseException("Specifying companion stars is required for automatic updating of the dynamics_code.particles masses. Otherwise use dynamics_code = None.")
            for particle, companion in zip(self.particles, self.companions):
                dyn_particles = ParticlesSuperset([part.as_set(), comp.as_set()]).get_intersecting_subset_in(self.dynamics_code.particles)
                dyn_particles.mass = [part.mass, comp.mass]
        return result
    

class LowMassXrayBinary(object):
    """
    Simulates a Low Mass X-ray Binary, using a single Stellar Evolution star model 
    and a recipe for Roche-lobe overflow. 
    """
    def __init__(self, mass_primary, mass_secondary, accretion_efficiency = 1.0):
        self.binary = Particles(2)
        self.primary = self.binary[0]
        self.secondary = self.binary[1]
        self.binary.mass = [mass_primary, mass_secondary]
        self.stellar_evolution_code = MESA
    
    def initialize(self, initial_age):
        self.stellar_evolution = self.stellar_evolution_code()
        self.stellar_evolution.initialize_module_with_default_parameters() 
        self.stellar_evolution.particles.add_particle(self.secondary)
        self.stellar_evolution.commit_particles()
        self.stellar_evolution.evolve_model(initial_age)
        self.secondary = self.stellar_evolution.particles[0]
    
    def roche_lobe_overflow_evolve(self, overflow_radius):
        d_mass = do_roche_lobe_overflow(self.secondary, overflow_radius)
        self.secondary.evolve_one_step()
        return self.primary.mass, self.secondary.mass, self.secondary.age
    

def new_low_mass_xray_binary(mass_primary, mass_secondary, initial_age, **keyword_arguments):
    """
    Creates a low mass X-ray binary, of which the secondary is evolved using 
    a Stellar Evolution legacy code. Mass is peeled off the star during 
    evolution to mimic Roche lobe overflow to the primary (neutron star or
    black hole).
    
    :argument mass_primary: Mass of the primary star (assumed to be a neutron star or black hole)
    :argument mass_secondary: Mass of the secondary star (assumed to be a low mass star filling its Roche lobe)
    :argument initial_age: Age of the secondary star when Roche lobe oveflow is assumed to start
    """
    low_mass_xray_binary = LowMassXrayBinary(mass_primary, mass_secondary, **keyword_arguments)
    low_mass_xray_binary.initialize(initial_age)
    return low_mass_xray_binary


def determine_RLOF_mass_excess(star, rmax):
    """
    Determine how much mass sticks out of radius rmax.
    This procedure is primarily intended for usage by do_roche_lobe_overflow,
    and returns additional info of the star so they don't need to be gathered
    twice. Returns: mass_excess, star-info-tuple
    (star-info-tuple holds current mass, radius profile, and mass-transfer flag)
    """
    #obtain stellar structure (radius and mass)
    r_profile = star.get_radius_profile()
    mf_profile = star.get_mass_profile() # really the mass fraction, i.e. sums to 1 instead of star.mass
    current_mass = star.mass
    #identify the shell where the radius exceeds rmax
    index = numpy.searchsorted(r_profile, rmax)
    if index == len(r_profile):
        return 0 | units.MSun, (None, None, False)
    # For the innermost mass-losing shell calculate the fraction of mass lost:
    # First determine the volume of the shell that sticks out rmax, and relate 
    # it to the total volume of the shell.
    # (4/3 pi) factor is ignored since it cancels out for fraction_V.
    V_shell = r_profile[index]**3 - r_profile[index-1]**3
    V_shell_out = r_profile[index]**3 - rmax**3
    fraction_V = V_shell_out / V_shell
    dmf_shell = mf_profile[index] * fraction_V 
    # calculate the mass fraction of the star with rstar>rmax
    delta_f_mass = dmf_shell + mf_profile[index+1:].sum()
    mass_transfer_flag = True if delta_f_mass > 0.0 else False
    return delta_f_mass * current_mass, (current_mass, r_profile, mass_transfer_flag)

def do_roche_lobe_overflow(star, roche_radius):
    """
    Apply Roche-lobe overflow to an existing star particle (currently only
    possible with MESA star particles). All mass of the star outside 
    roche_radius is lost. Returns the amount of mass lost.
    
    :argument star:         star particle (in MESA)
    :argument roche_radius: Roche-lobe radius
    """
    d_mass, (current_mass, r_profile, mass_transfer_flag) = \
        determine_RLOF_mass_excess(star, roche_radius)
    if not mass_transfer_flag:
        return 0 | units.MSun
    mass_scaling = (current_mass - d_mass) / current_mass
    star.mass = current_mass - d_mass

    # rescale stellar radius
    radial_scaling = roche_radius/r_profile[-1]
    r_profile *= radial_scaling
    star.set_radius_profile(r_profile)
    
    # rescale densities
    rho_profile = star.get_density_profile()
    rho_profile *= mass_scaling / radial_scaling**3
    star.set_density_profile(rho_profile)

    return d_mass




