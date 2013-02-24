from amuse.test import amusetest
import numpy


class harmonic_oscillator(object):
  def __init__(self,x,v,method=None):
    self.x=x
    self.v=v
    self.model_time=0
    self.method=method
  def kick(self,dt):
    self.v+=-self.x*dt
  def drift(self,dt):
    self.x+=self.v*dt
  def total_energy(self):
    return (self.v**2+self.x**2)/2  
  def evolve_step(self,dt):
    self.method(self.kick,self.drift,dt)
    self.model_time+=dt

def run_harmonic_oscillator(dt,method):
  h=harmonic_oscillator(0.,1.,method=method)
  tend=100*2*numpy.pi

  data=dict()
  data['x']=[]
  data['time']=[]
  data['de']=[]
  E0=h.total_energy()

  while h.model_time<tend-dt/2:
    h.evolve_step(dt)
    data['x'].append(h.x)
    data['time'].append(h.model_time)
    E=h.total_energy()
    data['de'].append(abs(E0-E)/E0)
  return data


class TestSymplecticCompositions(amusetest.TestCase):
  def test1(self):
    from amuse.ext.composition_methods import LEAPFROG
    dt1=.1
    data=run_harmonic_oscillator(dt1,method=LEAPFROG)  
    de1=max(data['de'])    
    dt2=0.01
    data=run_harmonic_oscillator(dt2,method=LEAPFROG)  
    de2=max(data['de'])
    order=int(numpy.log(de2/de1)/numpy.log(dt2/dt1)+0.5)
    self.assertEqual(order,2)
  def test2(self):
    from amuse.ext.composition_methods import SPLIT_4TH_S_M6
    dt1=.5
    data=run_harmonic_oscillator(dt1,method=SPLIT_4TH_S_M6)  
    de1=max(data['de'])    
    dt2=0.05
    data=run_harmonic_oscillator(dt2,method=SPLIT_4TH_S_M6)  
    de2=max(data['de'])
    order=int(numpy.log(de2/de1)/numpy.log(dt2/dt1)+0.5)
    self.assertEqual(order, 4)
  def test3(self):
    from amuse.ext.composition_methods import SPLIT_4TH_S_M5
    dt1=.5
    data=run_harmonic_oscillator(dt1,method=SPLIT_4TH_S_M5)  
    de1=max(data['de'])    
    dt2=0.05
    data=run_harmonic_oscillator(dt2,method=SPLIT_4TH_S_M5)  
    de2=max(data['de'])
    order=int(numpy.log(de2/de1)/numpy.log(dt2/dt1)+0.5)
    self.assertEqual(order, 4)        
  def test4(self):
    from amuse.ext.composition_methods import SPLIT_4TH_S_M4
    dt1=.5
    data=run_harmonic_oscillator(dt1,method=SPLIT_4TH_S_M4)  
    de1=max(data['de'])    
    dt2=0.05
    data=run_harmonic_oscillator(dt2,method=SPLIT_4TH_S_M4)  
    de2=max(data['de'])
    order=int(numpy.log(de2/de1)/numpy.log(dt2/dt1)+0.5)
    self.assertEqual(order, 4)
  def test5(self):
    from amuse.ext.composition_methods import SPLIT_6TH_SS_M11
    dt1=.5
    data=run_harmonic_oscillator(dt1,method=SPLIT_6TH_SS_M11)  
    de1=max(data['de'])    
    dt2=0.05
    data=run_harmonic_oscillator(dt2,method=SPLIT_6TH_SS_M11)  
    de2=max(data['de'])
    order=int(numpy.log(de2/de1)/numpy.log(dt2/dt1)+0.5)
    self.assertEqual(order, 6)
  def test6(self):
    from amuse.ext.composition_methods import SPLIT_6TH_SS_M13
    dt1=.5
    data=run_harmonic_oscillator(dt1,method=SPLIT_6TH_SS_M13)  
    de1=max(data['de'])    
    dt2=0.05
    data=run_harmonic_oscillator(dt2,method=SPLIT_6TH_SS_M13)  
    de2=max(data['de'])
    order=int(numpy.log(de2/de1)/numpy.log(dt2/dt1)+0.5)
    self.assertEqual(order, 6)                
  def test7(self):
    from amuse.ext.composition_methods import SPLIT_8TH_SS_M21
    dt1=1.
    data=run_harmonic_oscillator(dt1,method=SPLIT_8TH_SS_M21)  
    de1=max(data['de'])    
    dt2=0.25
    data=run_harmonic_oscillator(dt2,method=SPLIT_8TH_SS_M21)  
    de2=max(data['de'])
    order=int(numpy.log(de2/de1)/numpy.log(dt2/dt1)+0.5)
    self.assertEqual(order, 8)
  def test8(self):
    from amuse.ext.composition_methods import SPLIT_10TH_SS_M35
    dt1=1.
    data=run_harmonic_oscillator(dt1,method=SPLIT_10TH_SS_M35)  
    de1=max(data['de'])    
    dt2=0.5
    data=run_harmonic_oscillator(dt2,method=SPLIT_10TH_SS_M35)  
    de2=max(data['de'])
    order=int(numpy.log(de2/de1)/numpy.log(dt2/dt1)+0.5)
    self.assertEqual(order, 10)        
