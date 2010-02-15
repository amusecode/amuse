import os
import sys

from amuse.legacy.fi.interface import fi
from amuse.ext.evrard_test import MakeEvrardTest

from legacy_support import TestWithMPI

import numpy

from amuse.legacy.support.channel import MessageChannel

MessageChannel.no_redirection()


class testMPIInterface(TestWithMPI):

  def test0(self):
    instance=fi()  
    instance.setup_module()
    del instance

  def test1(self):
    instance=fi()  
    instance.setup_module()

    for x,l in [('usesph',0),('radiate',1),('starform',1),('cosmo',1),
                ('sqrttstp',1),('acc_tstp',0),('freetstp',1),('usequad',1),
                ('directsum',1),('selfgrav',0),('fixthalo',1),
                ('adaptive_eps',1),('gdgop',0),('smoothinput',1),
                ('consph',0),('sphinit',0),('uentropy',1),('isotherm',1),
                ('eps_is_h',0)]:
      result,err=eval("instance.get_"+x)()
      self.assertEquals( (x,result),(x,l))
      err=eval("instance.set_"+x)(1)
      result,err=eval("instance.get_"+x)()
      self.assertEquals( (x,result),(x,1))
      err=eval("instance.set_"+x)(0)
      result,err=eval("instance.get_"+x)()
      self.assertEquals((x,result),(x,0))

    for x,i in [ ('firstsnap',0),('stepout',5),('steplog',5),('max_tbin',4096),
                 ('minppbin',1),('targetnn',32),('verbosity',0),('nsmooth',64)]:
      result,err=eval("instance.get_"+x)()
      self.assertEquals( (x,result),(x,i))
      err=eval("instance.set_"+x)(12345)
      result,err=eval("instance.get_"+x)()
      self.assertEquals((x,result),(x,12345))

    for x,r in [ ('pboxsize',300.),('unitl_in_kpc',1.),('unitm_in_msun',1.e9),('dtime',1.),
      ('tstepcrit',1.),('tstpcr2',0.25),('freev',0.5),('freea',0.35),('freevexp',0.),
      ('freeaexp',-1.),('bh_tol',0.5),('gdgtol',0.01),('nn_tol',0.1),
      ('epsgas',0.005),('gamma',1.6666667),('alpha',0.5),('beta',1.0),('epssph',0.01),
      ('courant',0.3),('removgas',0.25),('consthsm',0.2),('nsmtol',0.1),
      ('graineff',0.05),('crionrate',3.6),('heat_par1',0.),('heat_par2',0.),
      ('cool_par',1.),('optdepth',0.),('tcollfac',1.),('masscrit',1.e5),
      ('sfeff',0.25),('tbubble',3.e7),('sne_eff',0.),('tsnbeg',3.e6),
      ('rhomax',100.),('eps',0.)]:
      result,err=eval("instance.get_"+x)()
      self.assertAlmostEquals(result,r,7)
      err=eval("instance.set_"+x)(0.)
      result,err=eval("instance.get_"+x)()
      self.assertEquals(result,0.)
      err=eval("instance.set_"+x)(0.12345)
      result,err=eval("instance.get_"+x)()
      self.assertAlmostEquals(result,0.12345,7)

    for x,s in [('halofile','none'),('feedback','fuv'),('sfmode','gerritsen'),
      ('hupdatemethod','mass'),('sph_visc','sph')]:
      result,err=eval("instance.get_"+x)()
      self.assertEquals((x,result),(x,s))
      err=eval("instance.set_"+x)("123")
      result,err=eval("instance.get_"+x)()
      self.assertEquals((x,result),(x,"123"))

    del instance

  def test2(self):
    instance=fi()  
    instance.setup_module()
    instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    retrieved_state = instance.get_state(1)
    self.assertEquals(11.0,  retrieved_state['mass'])
    self.assertEquals(2.0, retrieved_state['radius'])
    self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
    instance.cleanup_module()
    del instance

  def test3(self):
    instance=fi()
    instance.setup_module()
    instance.set_eps(0.001)
    instance.set_directsum(0)
    instance.new_particle( 
       [1.0,1.0,1.0],
       [0.0,0.0,0.0],
       [1.0,0.0,-1.0],
       [0.0,0.0,0.0],
       [0.0,0.0,0.0],
       [0.0,1.0,0.0],
       [0.0,0.0,0.0],
       [0.0,0.0,0.0] )
    instance.initialize_particles(0.0) 
    instance.synchronize_model()        
    Ep=instance.get_potential_energy()['potential_energy']
    Ek=instance.get_kinetic_energy()['kinetic_energy']

    self.assertAlmostEqual( Ek, 0.5,10)
    self.assertAlmostEqual( Ep, -2.5,10)
    instance.delete_particle(2)
    instance.reinitialize_particles() 
    instance.synchronize_model()        
    n=instance.get_number_of_particles()['number_of_particles']
    Ep=instance.get_potential_energy()['potential_energy']
    Ek=instance.get_kinetic_energy()['kinetic_energy']
    self.assertEqual( n, 2)
    self.assertAlmostEqual( Ek, 0.,10)
    self.assertAlmostEqual( Ep, -0.5,10)        

    instance.cleanup_module()
    del instance

class testEvrard(TestWithMPI):

  def test0(self):
    evrard=MakeEvrardTest(1000,grid=False)
    mass,x,y,z,vx,vy,vz,u=evrard.new_model()
    smooth=numpy.zeros_like(mass)
    nb = fi()
    nb.setup_module()

    nb.set_stepout(99999)
    nb.set_steplog(99999)
    nb.set_usesph(0)
    nb.set_radiate(1)
    nb.set_dtime(0.05)
    nb.set_gdgop(0)
    nb.set_uentropy(0)
    nb.set_verbosity(0)
    
    ids,error = nb.new_sph_particle(mass,smooth,x,y,z,vx,vy,vz,u)
    if filter(lambda x: x != 0, error) != []: raise Exception
 
    nb.initialize_particles(0.0)

    nb.synchronize_model()
    Ek,ret=nb.get_kinetic_energy()
    Ep,ret=nb.get_potential_energy()
    Eth,ret=nb.get_thermal_energy()

    self.assertAlmostEqual( Ek, 0.,3)
    self.assertAlmostEqual( Ep, -0.6611,3)        
    self.assertAlmostEqual( Eth, 0.05,3)        

    nb.evolve(0.5)
    nb.synchronize_model()
    Ek,ret=nb.get_kinetic_energy()
    Ep,ret=nb.get_potential_energy()
    Eth,ret=nb.get_thermal_energy()

    self.assertAlmostEqual( Ek, 0.129577,4)
    self.assertAlmostEqual( Ep, -0.831976,4)        
    self.assertAlmostEqual( Eth,0.08567999 ,4)        

    del evrard
    del nb
