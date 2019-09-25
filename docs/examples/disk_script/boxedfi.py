
from amuse.units import nbody_system
from amuse.units import units

import amuse.datamodel as core
from amuse.community.fi.interface import Fi
from amuse.ext.gasplummer import MakePlummerGasModel



class BoxedFi(Fi):
  def __init__(self, *args, **kargs):
    Fi.__init__(self, *args, **kargs)
    self.escapers=core.Particles(0)
  
  def evolve_model(self, *args, **kargs):
    self.stopping_conditions.out_of_box_detection.enable()
    outofbox=0.9*self.parameters.periodic_box_size/2
    self.parameters.stopping_conditions_out_of_box_size = outofbox
#    Fi.evolve_model(self,*args,**kargs)
    self.overridden().evolve_model(*args,**kargs)
    while self.stopping_conditions.out_of_box_detection.is_set():
      escapers=self.particles.select_array(
        lambda x,y,z: (x**2+y**2+z**2 > outofbox**2), ["x","y","z"])
      print "***", len(escapers)
      if len(escapers)>0:
        self.escapers.add_particles(escapers)
        self.particles.remove_particles(escapers)
#      Fi.evolve_model(self,*args, **kargs)
      self.overridden().evolve_model(*args,**kargs)
      

if __name__=="__main__":
  Ngas=1000
  conv = nbody_system.nbody_to_si(100 | units.MSun, 1 | units.parsec)
  dt=conv.to_si(1|nbody_system.time)/100
  print dt.in_(units.Myr)
  parts=MakePlummerGasModel(Ngas,convert_nbody=conv).result
  parts.h_smooth=0 | units.parsec
  outofbox=0.9*10. | units.parsec
  escapers=parts.select_array(
        lambda x,y,z: (x**2+y**2+z**2 > outofbox**2), ["x","y","z"])
  print "**",len(escapers),outofbox.in_(units.parsec)

  parts.remove_particles(escapers)
  print len(parts)

  sph=BoxedFi(convert_nbody=conv,use_gl=True)
  sph.parameters.periodic_box_size=20. | units.parsec
  sph.parameters.timestep=dt
  sph.parameters.self_gravity_flag=False

  sph.gas_particles.add_particles(parts)
  
  sph.start_viewer()
  
  sph.evolve_model(dt*1000)
    
  print len(sph.gas_particles)  
  print len(sph.particles)

