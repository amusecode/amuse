import sys
import numpy

from amuse.units import units
from amuse.datamodel import Particles
from amuse.support.exceptions import AmuseException


class StickySpheres(object):
    """
    Resolves collisions between particles by treating them as "sticky spheres", 
    i.e. perfectly inelastic collisions. Mass and momentum are conserved, while 
    all energy is lost in the center of mass frame.
    
    Optionally a mass_loss fraction (between 0 and 1) can be given. In this case, 
    a fraction of the total mass escapes the system during the collision, 
    carrying away momentum and energy. The velocity of the center of mass is 
    still conserved.
    """
    
    stellar_evolution_code_required = False
    gravity_code_required = False
    
    def __init__(self, mass_loss=0):
        if 0 <= mass_loss < 1:
            self.mass_loss = mass_loss
        else:
            raise AmuseException("Mass-loss fraction must be in the range [0, 1)")
    
    def handle_collision(self, primary, secondary):
        colliders = primary + secondary
        result = Particles(1)
        result.mass = colliders.total_mass() * (1 - self.mass_loss)
        result.position = colliders.center_of_mass()
        result.velocity = colliders.center_of_mass_velocity()
        if hasattr(colliders, "radius"):
            result.radius = colliders.radius.amax()
        return result

