from amuse.units import constants
from amuse.datamodel import Particles
from amuse.units.quantities import zero
from amuse.units.quantities import AdaptingVectorQuantity

class directsum(object):
    def __init__(self, systems, G=constants.G):
        self.systems=systems
        self.G=G
    
    def get_gravity_at_point(self,radius,x,y,z):
        part=Particles(0)
        for s in self.systems:
          part.add_particles(s.particles)
        ax=AdaptingVectorQuantity()
        ay=AdaptingVectorQuantity()
        az=AdaptingVectorQuantity()
        for rr,xx,yy,zz in zip(radius,x,y,z):
          dr2=((part.x-xx)**2+(part.y-yy)**2+(part.z-zz)**2+rr**2+part.radius**2)
          ax.append( (self.G*part.mass*(part.x-xx)/dr2**1.5).sum() )
          ay.append( (self.G*part.mass*(part.y-yy)/dr2**1.5).sum() )
          az.append( (self.G*part.mass*(part.z-zz)/dr2**1.5).sum() )
        return ax,ay,az
                          
    def get_potential_at_point(self,radius,x,y,z):
        part=Particles(0)
        for s in self.systems:
          part.add_particles(s.particles)
        phi=AdaptingVectorQuantity()
        for rr,xx,yy,zz in zip(radius,x,y,z):
          dr2=((part.x-xx)**2+(part.y-yy)**2+(part.z-zz)**2+rr**2+part.radius**2)
          phi.append( (-self.G*part.mass/dr2**0.5).sum() )
        return phi  
