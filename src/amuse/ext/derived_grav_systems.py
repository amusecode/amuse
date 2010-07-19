class center_of_mass(object):
  """
  com=center_of_mass(grav_instance)
  derived system, returns center of mass as skeleton grav system
  provides: get_gravity_at_point, get_potential_at_point
  """

  def __init__(self,baseclass):
    self.baseclass=baseclass

  def get_gravity_at_point(self,radius,x,y,z):
    mass=self.baseclass.total_mass
    xx,yy,zz=self.baseclass.get_center_of_mass_position()
    
    eps2=self.baseclass.parameters.epsilon_squared
    
    dr2=((xx-x)**2+(yy-y)**2+(zz-z)**2+eps2)
    
    ax=constants.G*mass*(xx-x)/dr2**1.5
    ay=constants.G*mass*(yy-y)/dr2**1.5
    az=constants.G*mass*(zz-z)/dr2**1.5
    
    return ax,ay,az

  def get_potential_at_point(self,radius,x,y,z):
    mass=self.baseclass.total_mass
    xx,yy,zz=self.baseclass.get_center_of_mass_position()
    
    eps2=self.baseclass.parameters.epsilon_squared
    dr2=((xx-x)**2+(yy-y)**2+(zz-z)**2+eps2)
    
    phi=-constants.G*mass/dr2**0.5
    
    return phi

class copycat(object):
  """
  copy=copycat(base_class,grav_instance, converter)
  derived system, returns copy of grav instance with
  get_gravity_at_point, get_potential_at_point reimplemented in 
  base_class
  """
  def __init__(self,baseclass, system,converter):
    self.baseclass=baseclass
    self.system=system
    self.converter=converter
        
  def get_gravity_at_point(self,radius,x,y,z):
    instance=self.baseclass(self.converter)

    instance.initialize_code()
    instance.parameters.epsilon_squared = self.system.parameters.epsilon_squared
    parts=self.system.particles.copy()
    instance.particles.add_particles(parts)

    ax,ay,az=instance.get_gravity_at_point(radius,x,y,z)
    
    instance.stop()
    return ax,ay,az

  def get_potential_at_point(self,radius,x,y,z):
    instance=self.baseclass(self.converter)

    instance.initialize_code()
    instance.parameters.epsilon_squared = self.system.parameters.epsilon_squared
    parts=self.system.particles.copy()
    instance.particles.add_particles(parts)

    phi=instance.get_potential_at_point(radius,x,y,z)
    
    instance.stop()
    return phi
