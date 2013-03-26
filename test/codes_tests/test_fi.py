import os
import sys
import numpy
import time
from amuse.test.amusetest import TestWithMPI

from amuse.community.fi.interface import FiInterface, Fi
from amuse.ext.evrard_test import new_evrard_gas_sphere
from amuse.ext.spherical_model import new_uniform_spherical_particle_distribution
from amuse.units import nbody_system
from amuse.units import units, constants
from amuse import datamodel
from amuse.rfi import channel
from amuse.ic.plummer import new_plummer_model

class TestFiInterface(TestWithMPI):

    def test1(self):
        instance=FiInterface()
        instance.initialize_code()
        instance.stop()
    
    def test2(self):
        instance=FiInterface()
        instance.initialize_code()
        
        for x,l in [('use_hydro',1),('radiate',0),('starform',0),('cosmo',1),
                    ('sqrttstp',0),('acc_tstp',1),('freetstp',0),('usequad',0),
                    ('directsum',0),('selfgrav',1),('fixthalo',0),
                    ('adaptive_eps',0),('gdgop',1),('smoothinput',0),
                    ('consph',1),('sphinit',1),('uentropy',1),('isotherm',0),
                    ('eps_is_h',1),('balsara',0)]:
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals( (x,result),(x,l))
            err=getattr(instance, 'set_'+x)(1)
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals( (x,result),(x,1))
            err=getattr(instance, 'set_'+x)(0)
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals((x,result),(x,0))
        
        for x,i in [ ('firstsnap',0),('stepout',5),('steplog',5),('max_tbin',4096),
                     ('minppbin',1),('targetnn',32),('verbosity',0),('nsmooth',64)]:
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals( (x,result),(x,i))
            err=getattr(instance, 'set_'+x)(12345)
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals((x,result),(x,12345))
        
        for x,r in [ ('pboxsize',10000.),('unitl_in_kpc',1.),('unitm_in_msun',1.e9),('dtime',1.),
                     ('tstepcrit',1.),('tstpcr2',0.25),('freev',0.5),('freea',0.35),('freevexp',0.),
                     ('freeaexp',-1.),('bh_tol',0.5),('gdgtol',0.01),('nn_tol',0.1),
                     ('epsgas',0.005),('gamma',1.6666667),('alpha',0.5),('beta',1.0),('epssph',0.01),
                     ('courant',0.3),('removgas',0.25),('consthsm',0.2),('nsmtol',0.1),
                     ('graineff',0.05),('crionrate',3.6),('heat_par1',0.),('heat_par2',0.),
                     ('cool_par',1.),('optdepth',0.),('tcollfac',1.),('masscrit',1.e5),
                     ('sfeff',0.25),('tbubble',3.e7),('sne_eff',0.),('tsnbeg',3.e6),
                     ('rhomax',100.),('eps',0.)]:
            result,err=getattr(instance, 'get_'+x)()
            self.assertAlmostEquals(result,r,7)
            err=getattr(instance, 'set_'+x)(0.)
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals(result,0.)
            err=getattr(instance, 'set_'+x)(0.12345)
            result,err=getattr(instance, 'get_'+x)()
            self.assertAlmostEquals(result,0.12345,7)
        
        for x,s in [('halofile','none'),('feedback','fuv'),('sfmode','gerritsen'),
                    ('hupdatemethod','mass'),('sph_visc','sph')]:
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals((x,result),(x,s))
            err=getattr(instance, 'set_'+x)("123")
            result,err=getattr(instance, 'get_'+x)()
            self.assertEquals((x,result),(x,"123"))
        
        instance.stop()
    
    def test3(self):
        instance=FiInterface()
        instance.initialize_code()
        instance.commit_parameters()        
        instance.new_particle(11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        instance.cleanup_code()
        instance.stop()

  
    def test4(self):
        instance=FiInterface()
        instance.initialize_code()
        instance.commit_parameters()
        instance.new_particle(11.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0)
        instance.set_time_step(2.0)
        retrieved_state = instance.get_time_step()
        self.assertEquals(2.0, retrieved_state['time_step'])
        instance.cleanup_code()
        instance.stop()
    
    def test5(self):
        instance=FiInterface()
        instance.initialize_code()
        instance.set_eps(0.001)
        instance.set_directsum(1)
        instance.commit_parameters()
        instance.new_particle( 
           [1.0,1.0,1.0],
           [0.0,0.0,0.0],
           [1.0,0.0,-1.0],
           [0.0,0.0,0.0],
           [0.0,0.0,0.0],
           [0.0,1.0,0.0],
           [0.0,0.0,0.0],
           [0.0,0.0,0.0] )
        instance.commit_particles()
        self.assertEqual(instance.get_number_of_particles()['number_of_particles'], 3)
        instance.synchronize_model()
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
        
        self.assertAlmostEqual( Ek, 0.5,10)
        self.assertAlmostEqual( Ep, -2.5,10)
        instance.delete_particle(2)
        instance.recommit_particles()
        instance.synchronize_model()
        n=instance.get_number_of_particles()['number_of_particles']
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
        self.assertEqual( n, 2)
        self.assertAlmostEqual( Ek, 0.,10)
        self.assertAlmostEqual( Ep, -0.5,10)
        
        instance.cleanup_code()
        instance.stop()

    def test5a(self):
        instance=FiInterface()
        self.assertEquals(0,instance.initialize_code())
        self.assertEquals(0,instance.set_eps(0.001))
        self.assertEquals(0,instance.set_directsum(1))
        instance.new_particle( 
           [1.0,1.0],
           [0.0,0.0],
           [1.0,-1.0],
           [0.0,0.0],
           [0.0,0.0],
           [0.0,0.0],
           [0.0,0.0],
           [0.0,0.0] )
        self.assertEquals(0,instance.commit_particles())
        self.assertEquals(0, instance.commit_parameters())
        self.assertAlmostEqual(-0.500 , instance.get_potential(1)['potential'], places=1)
        print instance.get_potential([1])
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test6(self):
        print "Testing FiInterface get_hydro_state_at_point"
        number_sph_particles = 10000
        length = 1.0 | units.kpc
        mass = 1.0e9 | units.MSun
        
        gas = new_uniform_spherical_particle_distribution(number_sph_particles, length, mass)
        mass = [1.0 / number_sph_particles] * number_sph_particles
        h    = [0.01] * number_sph_particles
        x, y, z = gas.x.value_in(units.kpc), gas.y.value_in(units.kpc), gas.z.value_in(units.kpc)
        vx, vy, vz = [[0.0] * number_sph_particles] * 3
        u = [0.05] * number_sph_particles
        indices = range(1,number_sph_particles+1)
        
        instance=FiInterface()
        instance.initialize_code()
        instance.set_unitl_in_kpc(1.0)
        instance.set_unitm_in_msun(1.e9)
        instance.set_nsmooth(64)
        instance.set_nsmtol(0.2)
        instance.set_uentropy(1)
        instance.commit_parameters()
        instance.new_sph_particle(mass,x, y, z, vx, vy, vz, u, h)
        instance.commit_particles()
        instance.synchronize_model()

        h = instance.get_smoothing_length(indices)['h_smooth']
        self.assertIsOfOrder((instance.get_nsmooth()['nsmooth']*1.0 / number_sph_particles)**(1.0/3), h)
        
        hydrostate = instance.get_hydro_state_at_point(0, 0, 0, 0, 0, 0)
        density = 1.0 / (4.0/3.0 * numpy.pi * 1.0**3)
        self.assertAlmostEqual(hydrostate['rho'],   density, places=3)
        self.assertAlmostEqual(hydrostate['rhovx'],       0, places=3)
        self.assertAlmostEqual(hydrostate['rhovy'],       0, places=3)
        self.assertAlmostEqual(hydrostate['rhovz'],       0, places=3)
        self.assertAlmostEqual(hydrostate['rhoe'], density*u[0], places=3)
        
        instance.stop()

    



    def test7(self):
        instance=FiInterface(mode=FiInterface.MODE_PERIODIC_BOUNDARIES)
        instance.initialize_code()
        #instance.set_periodic(1)
        instance.set_use_hydro(0)
        instance.set_selfgrav(0)
        print instance.get_pboxsize()
        instance.set_pboxsize(2.)
        print instance.get_pboxsize()
        instance.set_dtime(0.1)
        instance.commit_parameters()
        ids,err=instance.new_particle( 
           [1.0,1.0,1.0],
           [0.5,0.0,0.0],
           [0.0,-0.5,0.0],
           [0.0,0.0,0.5],
           [1.0,0.0,0.0],
           [0.0,-1.0,0.0],
           [0.0,0.0,1.0])
        instance.commit_particles()
        
        instance.evolve_model(0.1)
        m,x,y,z,vx,vy,vz,r,err=instance.get_state(ids)
        self.assertAlmostEqual(x, [0.6,0.,0.], places=7)
        self.assertAlmostEqual(y, [0.,-0.6,0.], places=7)
        self.assertAlmostEqual(z, [0.,0.,0.6], places=7)
        instance.evolve_model(1.0)
        m,x,y,z,vx,vy,vz,r,err=instance.get_state(ids)
        self.assertAlmostEqual(x, [-0.5,0.,0.], places=7)
        self.assertAlmostEqual(y, [0.,0.5,0.], places=7)
        self.assertAlmostEqual(z, [0.,0.,-0.5], places=7)
        instance.cleanup_code()
        instance.stop()
        
    
    def test8(self):
        instance=FiInterface(redirection="none")
        instance.initialize_code()
        instance.set_use_hydro(0)
        instance.set_selfgrav(1)
        instance.set_pboxsize(1.0)
        instance.set_dtime(0.2)
        instance.commit_parameters()
        index,err=instance.new_sph_particle( 
           0.1, #mass
           0.5, #x,y,z
           0.0,
           0.0,
           1.0, #vx,vy,vz
           0.0,
           0.0,
           0.0
        )
        index2,err=instance.new_sph_particle( 
           0.001, #mass
           0.1, #x,y,z
           0.0,
           0.0,
           1.0, #vx,vy,vz
           0.0,
           0.0,
           0.0
        )
        instance.commit_particles()
        instance.evolve_model(0.1)
        m,x,y,z,vx,vy,vz,r,err=instance.get_state(index)
        self.assertAlmostRelativeEquals(x, 0.5)
        nremoved, error = instance.get_number_of_sph_particles_removed()
        self.assertEquals(error,0)
        self.assertEquals(nremoved,0)
        instance.evolve_model(0.15)
        nremoved, error = instance.get_number_of_sph_particles_removed()
        self.assertEquals(error,0)
        self.assertEquals(nremoved,1)
        idremoved, error = instance.get_id_of_removed_sph_particle(0)
        self.assertEquals(idremoved,index)
        
    
class TestEvrard(TestWithMPI):

    def xtest0(self):
        evrard=MakeEvrardTest(1000)
        mass,x,y,z,vx,vy,vz,u=evrard.new_model()
        smooth=numpy.zeros_like(mass)
        nb = FiInterface()
        nb.initialize_code()
        
        nb.set_stepout(99999)
        nb.set_steplog(99999)
        nb.set_use_hydro(1)
        nb.set_radiate(0)
        nb.set_dtime(0.05)
        nb.set_gdgop(1)
        nb.set_uentropy(1)
        nb.set_verbosity(0)
        nb.commit_parameters()
        
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
        
        nb.evolve_model(0.5)
        nb.synchronize_model()
        Ek,ret=nb.get_kinetic_energy()
        Ep,ret=nb.get_potential_energy()
        Eth,ret=nb.get_thermal_energy()
        self.assertAlmostEqual( Ek, 0.129577,3)
        self.assertAlmostEqual( Ep, -0.831976,3)
        self.assertAlmostEqual( Eth,0.08567999,3)
        
        del evrard
        nb.stop()
        

class TestFi(TestWithMPI):

    def test0(self):
        instance=Fi(nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.parameters.timestep = 0.5 | units.day
        self.assertEquals(instance.parameters.timestep, 0.5 | units.day)
        instance.commit_parameters()
        instance.commit_particles()
        self.assertEquals(instance.parameters.timestep, 0.5 | units.day)
        instance.stop()
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = Fi(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00000001 | units.AU**2
        instance.parameters.timestep = 0.5 | units.day
        self.assertEquals(instance.parameters.timestep, 0.5 | units.day)
        instance.commit_parameters()
        
        stars = datamodel.Particles(2)
        
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = [0.0,0.0,0.0] | units.m
        sun.velocity = [0.0,0.0,0.0] | units.ms
        sun.radius = units.RSun(1.0)
        
        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = [149.5e6, 0.0, 0.0] | units.km
        earth.velocity = [0.0, 29800, 0.0] | units.ms
        
        instance.particles.add_particles(stars)
        self.assertEquals(instance.parameters.timestep, 0.5 | units.day)
        
        postion_at_start = earth.position.x
        
        instance.evolve_model(365.0 | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        
        postion_after_full_rotation = earth.position.x
        
        self.assertAlmostRelativeEqual(postion_at_start, postion_after_full_rotation, 4)
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        
        postion_after_half_a_rotation = earth.position.x
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        
        postion_after_half_a_rotation = earth.position.x
        self.assertAlmostRelativeEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        
        postion_after_half_a_rotation = earth.position.y
        
        self.assertAlmostRelativeEqual(-postion_at_start, postion_after_half_a_rotation, 4)
        
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Test 2: testing Fi data directory"
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Fi(convert_nbody)
        instance.initialize_code()
        self.assertTrue('data/fi/input/' in instance.legacy_interface.get_fi_data_directory()['fi_data_directory'])
        self.assertEquals(instance.legacy_interface.get_fi_data_directory()['fi_data_directory'], 
            instance.get_data_directory()+'/')
        
        self.assertEquals(False, instance.parameters.radiation_flag)
        instance.parameters.radiation_flag = True
        self.assertEquals(True, instance.parameters.radiation_flag)
        
        self.assertEquals(False, instance.parameters.star_formation_flag)
        instance.parameters.star_formation_flag = True
        self.assertEquals(True, instance.parameters.star_formation_flag)
        
        instance.commit_parameters()
        stars = datamodel.Particles(2)
        stars.mass = [1.0, 3.0e-6] | units.MSun
        stars.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.AU
        stars.velocity = [[0.0,0.0,0.0], [0.0,29.8,0.0]] | units.km / units.s
        stars.radius = [1.0, 0.01] | units.RSun
        
        instance.particles.add_particles(stars)
        instance.evolve_model(1.0 | nbody_system.time)
        
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):
        print "Test 3: testing Fi boolean parameters"
        instance = Fi(nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        instance.initialize_code()
        
        for bool_par in ['radiation_flag','star_formation_flag',
            'square_root_timestep_flag','freeform_timestep_flag','quadrupole_moments_flag',
            'direct_sum_flag','fixed_halo_flag','adaptive_smoothing_flag','smooth_input_flag',
            'isothermal_flag','balsara_flag']:
            self.assertEquals(False, getattr(instance.parameters, bool_par))
            setattr(instance.parameters, bool_par, True)
            self.assertEquals(True, getattr(instance.parameters, bool_par))
        
        for bool_par in ['acc_timestep_flag','self_gravity_flag','gadget_cell_opening_flag',
            'use_hydro_flag','conservative_sph_flag','sph_dens_init_flag',
            'integrate_entropy_flag','eps_is_h_flag']:
            self.assertEquals(True, getattr(instance.parameters, bool_par))
            setattr(instance.parameters, bool_par, False)
            self.assertEquals(False, getattr(instance.parameters, bool_par))
        
        instance.cleanup_code()
        instance.stop()
    
    def test4(self):
        print "Test 4: testing Fi integer parameters"
        instance = Fi(nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        instance.initialize_code()
        
        for int_par, value in [('first_snapshot',0),('output_interval',5),('log_interval',5),
            ('maximum_time_bin',4096),('minimum_part_per_bin',1),('targetnn',32),
            ('verbosity',0),('n_smooth',64)]:
            self.assertEquals(value, getattr(instance.parameters,int_par))
            setattr(instance.parameters, int_par, 1)
            self.assertEquals(1, getattr(instance.parameters,int_par))
        
        instance.cleanup_code()
        instance.stop()
    
    def test5(self):
        print "Test 5: testing Fi double precision parameters"
        instance = Fi(nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        instance.initialize_code()
        
        par_names=['epsilon_squared','timestep','periodic_box_size','code_mass_unit','code_length_unit',
            'sqrt_timestep_crit_constant','acc_timestep_crit_constant','free_timestep_crit_constant_v',
            'free_timestep_crit_constant_a','free_timestep_crit_constant_vexp',
            'free_timestep_crit_constant_aexp','opening_angle','gadget_cell_opening_constant',
            'nn_tol','gas_epsilon','gamma','artificial_viscosity_alpha','beta','sph_artificial_viscosity_eps','courant',
            'min_gas_part_mass','sph_h_const','n_smooth_tol','grain_heat_eff','zeta_cr_ion_rate',
            'heat_par1','heat_par2','cool_par','optical_depth','star_form_delay_fac','star_form_mass_crit',
            'star_form_eff','supernova_duration','supernova_eff','t_supernova_start','max_density']
        defaults=[0.0 | nbody_system.length * nbody_system.length, 1.0 | nbody_system.time, 10000.0 | nbody_system.length, 
            1.0e9 | units.MSun, 1.0 | units.kpc, 1.0, 0.25, 0.5, 0.35, 
            0.0, -1.0, 0.5, 0.01, 0.1, 0.005 | nbody_system.length, 1.6666667, 0.5, 1.0, 0.01, 0.3, 
            0.25, 0.2 | nbody_system.length, 0.1, 0.05, 3.6 | 1.8e-17 * units.s**-1, 0.0, 0.0, 1.0, 
            0.0, 1.0, 1.0e5 | units.MSun, 0.25, 3.0e7 | units.Myr, 0.0, 3.e6 | units.Myr, 100.0]
        defaults = [instance.unit_converter.to_si(val) if not isinstance(val,float) and  nbody_system.is_nbody_unit(val.unit) else val for val in defaults]
        for double_par, value in zip(par_names, defaults):
            self.assertAlmostRelativeEquals(getattr(instance.parameters,double_par), value, 7)
            setattr(instance.parameters,double_par,2 * value)
            self.assertAlmostRelativeEquals(getattr(instance.parameters,double_par), 2*value, 7)
        instance.cleanup_code()
        instance.stop()
    
    def test6(self):
        print "Test 6: testing Fi string parameters"
        instance = Fi(nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        instance.initialize_code()
        
        par_names=['halofile','feedback','star_formation_mode','h_update_method','sph_viscosity','fi_data_directory']
        defaults=['none','fuv','gerritsen','mass','sph',instance.get_data_directory()+'/']
        new_values=['bct_02_10.halo','pres','nieuw','test','sphv','test']
        for string_par, value, new_value in zip(par_names, defaults, new_values):
            print instance.get_halofile(), getattr(instance.parameters, string_par), 
            self.assertEquals(getattr(instance.parameters, string_par), value)
            setattr(instance.parameters, string_par, new_value)
            self.assertEquals(getattr(instance.parameters, string_par), new_value)
        instance.cleanup_code()
        instance.stop()
    
    def test7(self):
        print "Test 7: testing Fi SPH particles"
        target_number_of_particles = 1000
        gas = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed = 1234)
        gas.h_smooth = 0.0 | nbody_system.length
        
        convert_nbody = nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        instance = Fi(convert_nbody,redirection="none")
        instance.initialize_code()
        instance.parameters.timestep = 0.025 | nbody_system.time
        instance.parameters.integrate_entropy_flag = True
        instance.gas_particles.add_particles(gas)
        
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.kinetic_energy),    0.00 | nbody_system.energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.potential_energy), -0.494 | nbody_system.energy, 2)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.05 | nbody_system.energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.total_energy),     -0.444 | nbody_system.energy, 2)
        
        instance.evolve_model(10.437 | units.Myr)
        instance.synchronize_model()
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.kinetic_energy),    0.0566 | nbody_system.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.potential_energy), -0.5661 | nbody_system.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.0635 | nbody_system.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.total_energy),     -0.444 | nbody_system.energy, 2)
        
        instance.cleanup_code()
        instance.stop()
    
    def test8(self):
        print "Test 8: testing Fi dark matter + SPH particles"
        target_number_of_particles = 100
        gas = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed = 1234)
        gas.h_smooth = 0.0 | nbody_system.length
        
        dark = datamodel.Particles(2)
        dark.mass = [1.0, 3.0e-6] | units.MSun
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.AU
        dark.velocity = [[0.0,0.0,0.0], [0.0,29.8,0.0]] | units.km / units.s
        dark.radius = [10.0, 0.1] | units.RSun
        
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Fi(convert_nbody)
        instance.parameters.eps_is_h_flag=False
        instance.parameters.gas_epsilon=0.4 | units.AU
        instance.parameters.timestep=0.001 | units.yr
#        instance.parameters.gadget_cell_opening_constant=0.001
        instance.dm_particles.add_particles(dark)
        instance.gas_particles.add_particles(gas)
        
        print convert_nbody.to_nbody(instance.kinetic_energy)
        print convert_nbody.to_nbody(instance.potential_energy)
        print convert_nbody.to_nbody(instance.thermal_energy)
        print convert_nbody.to_nbody(instance.total_energy)
        
        radius=instance.dm_particles.radius
        print instance.gas_particles.radius[0:10].in_(units.AU)
        self.assertEqual(dark[0].radius,radius[0])
        self.assertEqual(dark[1].radius,radius[1])
        
        
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.kinetic_energy),  
          convert_nbody.to_nbody(dark.kinetic_energy()), 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.potential_energy), -1.8656 | nbody_system.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.0500 | nbody_system.energy,3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.total_energy),     -1.8156 | nbody_system.energy, 3)
        
        instance.evolve_model(1.| units.yr)
        print instance.gas_particles.radius[0:10].in_(units.AU)

        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.kinetic_energy),    0.02507 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.potential_energy), -2.9630 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.thermal_energy),    1.122 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.total_energy),     -1.8156 | nbody_system.energy, 3)

        print convert_nbody.to_nbody(instance.total_energy)
        
        instance.cleanup_code()
        instance.stop()
    
    def test9(self):
        print "Test 9: testing Fi dark matter + SPH + star particles"
        target_number_of_particles = 100
        gas = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed = 1234)
        gas.h_smooth = 0.0 | nbody_system.length
        
        dark = datamodel.Particles(2)
        dark.mass = [0.4, 0.4] | nbody_system.mass
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.kpc
        dark.velocity = [[100.0,100.0,100.0], [1.0,1.0,1.0]] | units.km / units.s
        dark.radius = [0.0, 0.0] | units.RSun
        
        star = datamodel.Particles(2)
        star.mass = [0.02, 0.02] | nbody_system.mass
        star.position = [[0.1,0.2,0.3], [0.4,0.5,0.6]] | units.kpc
        star.velocity = [[-300.0,-200.0,-100.0], [-6.0,-5.0,-4.0]] | units.km / units.s
        star.radius = [0.0, 0.0] | units.RSun
        star.tform = [1000.0, 1000.0] | units.Myr
        
        convert_nbody = nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        instance = Fi(convert_nbody)
        instance.initialize_code()
        instance.parameters.verbosity = 0
        instance.parameters.timestep=0.05 | units.Myr
        instance.parameters.eps_is_h_flag=False
        instance.parameters.gas_epsilon=0.4 | units.kpc

        instance.dm_particles.add_particles(dark)
        instance.star_particles.add_particles(star)
        instance.gas_particles.add_particles(gas)
         
        print convert_nbody.to_nbody(instance.kinetic_energy)
        print convert_nbody.to_nbody(instance.potential_energy)
        print convert_nbody.to_nbody(instance.thermal_energy)
        print convert_nbody.to_nbody(instance.total_energy)
        print convert_nbody.to_nbody( (instance.dm_particles|instance.star_particles|instance.gas_particles).potential_energy(smoothing_length_squared = instance.parameters.epsilon_squared))
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.kinetic_energy),    1.7204 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.potential_energy), -1.3096 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.0500 | nbody_system.energy)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.total_energy),      0.461 | nbody_system.energy, 3)
        
        instance.evolve_model(100*instance.parameters.timestep)
        print convert_nbody.to_nbody(instance.kinetic_energy)
        print convert_nbody.to_nbody(instance.potential_energy)
        print convert_nbody.to_nbody(instance.thermal_energy)
        print convert_nbody.to_nbody(instance.total_energy)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.kinetic_energy),    1.6825 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.potential_energy), -1.277 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.05526 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(convert_nbody.to_nbody(instance.total_energy),      0.461 | nbody_system.energy, 3)
        
        instance.cleanup_code()
        instance.stop()
    
    def test10(self):
        print "Test 10: testing Fi star particles"
        stars = datamodel.Particles(2)
        stars.mass = [1.0, 3.00e-6] | units.MSun
        stars.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.AU
        stars.velocity = [[0.0,0.0,0.0], [0.0,29.8,0.0]] | units.km / units.s
        stars.radius = [1.0, 1.0] | units.RSun # results are nonsense if not (r1 == r2)... ?
        stars.tform = [-10.0, -50.0] | units.Myr
        
        instance = Fi(nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU))#, debugger='xterm')
        instance.initialize_code()
        instance.parameters.timestep = 0.5 | units.day
        instance.commit_parameters()
        instance.star_particles.add_particles(stars)
        self.assertAlmostEqual(instance.star_particles.tform, [-10.0, -50.0] | units.Myr)
        instance.star_particles.tform = [-100.0, -500.0] | units.Myr
        
        instance.synchronize_model()
        self.assertAlmostEqual(instance.star_particles.tform, [-100.0, -500.0] | units.Myr)
        self.assertAlmostEqual(instance.kinetic_energy, 2.6493e+33|units.J, 3, in_units=1e+33*units.J)
        self.assertAlmostEqual(instance.star_particles.kinetic_energy(), 2.6493e+33|units.J, 3, in_units=1e+33*units.J)

        instance.evolve_model(0.5 | units.yr)
        self.assertAlmostEqual(instance.star_particles[1].x, -1.0 | units.AU,2)
        instance.evolve_model(1.0 | units.yr)
        self.assertAlmostEqual(instance.star_particles[1].x, 1.0 | units.AU,3)
        self.assertAlmostEqual(instance.kinetic_energy, 2.6493e+33|units.J, 3, in_units=1e+33*units.J)
        
        instance.cleanup_code()
        instance.stop()
    
    def test11(self):
        print "Test 11: testing Fi (dm+sph+star) particles Superset"
        target_number_of_particles = 100
        gas = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed = 1234)
        gas.h_smooth = 0.0 | nbody_system.length
        number_sph_particles = len(gas)
        
        dark = datamodel.Particles(2)
        dark.mass = [0.4, 0.4] | nbody_system.mass
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.kpc
        dark.velocity = [[100.0,100.0,100.0], [1.0,1.0,1.0]] | units.km / units.s
        dark.radius = [0.0, 0.0] | units.RSun
        
        star = datamodel.Particles(2)
        star.mass = [0.02, 0.02] | nbody_system.mass
        star.position = [[0.1,0.2,0.3], [0.4,0.5,0.6]] | units.kpc
        star.velocity = [[-300.0,-200.0,-100.0], [-6.0,-5.0,-4.0]] | units.km / units.s
        star.radius = [0.0, 0.0] | units.RSun
        star.tform = [1000.0, 1000.0] | units.Myr
        
        convert_nbody = nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        instance = Fi(convert_nbody)
        instance.initialize_code()

        self.assertEquals(0, len(instance.particles))
        instance.dm_particles.add_particles(dark)
        instance.star_particles.add_particles(star)
        instance.gas_particles.add_particles(gas)
        self.assertEquals(number_sph_particles+4, len(instance.particles))
        
        print "The 'particles' superset can still be used as standard GD 'particles' set."
        default = datamodel.Particles(2)
        default.mass = [0.4, 0.4] | nbody_system.mass
        default.position = [[0.5,-0.5,0.04], [1.5,0.5,0.08]] | units.kpc
        default.velocity = [[10.0,10.0,10.0], [10.0,10.0,10.0]] | units.km / units.s
        default.radius = [0.0, 0.0] | units.RSun
        instance.particles.add_particles(default)
        
        self.assertEquals(number_sph_particles+6, len(instance.particles))
        print "'>>> print instance.particles' only prints those particle attributes the subsets have in common."
        string_produced_by_print = instance.particles.__str__()
        self.assertTrue("mass" in string_produced_by_print)
        self.assertTrue("vx" in string_produced_by_print)
        self.assertTrue("radius" in string_produced_by_print)
        self.assertFalse("tform" in string_produced_by_print)
        
        instance.synchronize_model()
        self.assertAlmostRelativeEqual(instance.particles.kinetic_energy(), 
            instance.dm_particles.kinetic_energy() +
            instance.star_particles.kinetic_energy() +
            instance.gas_particles.kinetic_energy())
        self.assertAlmostRelativeEqual(instance.particles.kinetic_energy(), instance.kinetic_energy)
        instance.evolve_model(1.0 | units.yr)
        instance.synchronize_model()
        self.assertAlmostRelativeEqual(instance.particles.kinetic_energy(), instance.kinetic_energy, -1)
        
        instance.cleanup_code()
        instance.stop()
    
    def test12(self):
        print "Testing Fi states"
        target_number_of_particles = 100
        gas = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed = 1234)
        gas.h_smooth = 0.0 | nbody_system.length
        dark = datamodel.Particles(2)
        dark.mass = [0.4, 0.4] | nbody_system.mass
        dark.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | nbody_system.length
        dark.velocity = [[0.1, 0.1, 0.1], [1.0, 1.0, 1.0]] | nbody_system.speed
        dark.radius = 0.0 | nbody_system.length
        
        print "First do everything manually:"
        convert_nbody = nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        instance = Fi(convert_nbody)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.gas_particles.add_particles(gas)
        instance.commit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        mass = instance.gas_particles[0].mass
        instance.evolve_model(0.001 | nbody_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()
        
        print "initialize_code(), commit_parameters(), (re)commit_particles(), and cleanup_code() should be called " \
            "automatically before setting parameters, new_xx_particle(), get_xx(), and stop():"
        instance = Fi(convert_nbody)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.timestep = 0.001 | nbody_system.time
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        self.assertEquals(instance.parameters.timestep, convert_nbody.to_si(0.001 | nbody_system.time))
        instance.gas_particles.add_particles(gas)
        instance.dm_particles.add_particles(dark)
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        mass = instance.gas_particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.synchronize_model() # model was already synchronized, but fi doesn't seem to get that...
        instance.dm_particles.remove_particles(dark)
        self.assertEquals(instance.get_name_of_current_state(), 'UPDATE')
        mass = instance.gas_particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.evolve_model(0.002 | nbody_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'STOPPED')
    
    def test13(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,10.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass

        instance = Fi(redirection="none")
        instance.initialize_code()
        instance.parameters.stopping_conditions_number_of_steps = 2
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps, 2)
        instance.parameters.epsilon_squared = (0.01 | nbody_system.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(10 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        self.assertTrue(instance.model_time < 10 | nbody_system.time)

        instance.stop()

    def test13a(self):
        particles = datamodel.Particles(2)
        particles.x = [0.0,10.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.005 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass

        very_short_time_to_evolve = 1 | units.s
        very_long_time_to_evolve = 1e9 | nbody_system.time

        instance = Fi()
        instance.initialize_code()
        instance.parameters.stopping_conditions_timeout = very_short_time_to_evolve
        self.assertEquals(instance.parameters.stopping_conditions_timeout, very_short_time_to_evolve)
        instance.parameters.epsilon_squared = (0.01 | nbody_system.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.timeout_detection.enable()
        start = time.time()
        instance.evolve_model(very_long_time_to_evolve)
        end = time.time()
        self.assertTrue(instance.stopping_conditions.timeout_detection.is_set())
        self.assertTrue((end-start)<very_short_time_to_evolve.value_in(units.s) + 2)#2 = some overhead compensation

        instance.stop()

    def test13b(self):
        print "Test out_of_box_detection stopping condition"
        particles = datamodel.Particles(3)
        particles.position = [0.0, 0.0, 0.0] | nbody_system.length
        particles.velocity = [[-1.0, 0, 0], [1.0, 0, 0], [0.5, 0, 0]] | nbody_system.speed
        particles.radius = 0.001 | nbody_system.length
        particles.mass = 1.0 | nbody_system.mass
        
        gas = datamodel.Particles(2)
        gas.position = [[0.0, 0.5, 0.0], [0.0, -0.5, 0.0]] | nbody_system.length
        gas.velocity = [[0.0, 0.5, 0.0], [0.0, -0.25, 0.0]] | nbody_system.speed
        gas.radius = 0.001 | nbody_system.length
        gas.mass = 1.0 | nbody_system.mass
        gas.u = 0.01 | nbody_system.specific_energy
        
        instance = Fi()
        instance.initialize_code()
        instance.parameters.stopping_conditions_out_of_box_size = 1.0 | nbody_system.length
        instance.parameters.self_gravity_flag = False
        instance.parameters.timestep = 0.25 | nbody_system.time
        self.assertEquals(instance.parameters.stopping_conditions_out_of_box_size,  1.0 | nbody_system.length)
        instance.gas_particles.add_particles(gas)
        instance.dm_particles.add_particles(particles)
        out_of_box_detection = instance.stopping_conditions.out_of_box_detection
        out_of_box_detection.enable()
        instance.evolve_model(2 | nbody_system.time)
        self.assertTrue(out_of_box_detection.is_set())
        self.assertAlmostEqual(instance.model_time, 1.0 | nbody_system.time, 3)
        self.assertEquals(len(out_of_box_detection.particles()), 3)
        self.assertEquals(out_of_box_detection.particles().position.lengths() >= 
            1.0 | nbody_system.length, [True, True, True])
        self.assertEquals(len(out_of_box_detection.particles(1)), 0)
        instance.particles.remove_particles(out_of_box_detection.particles())
        
        stopping_conditions_info = str(instance.stopping_conditions).split("\n")
        self.assertEquals(stopping_conditions_info[0], "Stopping conditions of a 'Fi' object")
        for i, word in enumerate(["supported", "enabled", "set"]):
            self.assertTrue(word in stopping_conditions_info[i+1])
            self.assertTrue("out_of_box_detection" in stopping_conditions_info[i+1])
        
        instance.evolve_model(1.5 | nbody_system.time)
        self.assertFalse(out_of_box_detection.is_set())
        self.assertAlmostEqual(instance.model_time, 1.5 | nbody_system.time, 3)
        
        instance.evolve_model(3.0 | nbody_system.time)
        self.assertTrue(out_of_box_detection.is_set())
        self.assertAlmostEqual(instance.model_time, 2.0 | nbody_system.time, 3)
        self.assertEquals(len(out_of_box_detection.particles()), 2)
        instance.stop()

    
    def test14(self):
        print "Testing Fi get_hydro_state_at_point"
        number_sph_particles = 100
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.kpc, 1.0e10 | units.MSun)
        gas = new_evrard_gas_sphere(number_sph_particles, convert_nbody, seed = 1234)
        gas.h_smooth = 0.01 | nbody_system.length
        
        instance = Fi(convert_nbody)
        instance.parameters.n_smooth = 64
        instance.parameters.n_smooth_tol = 0.01
        instance.gas_particles.add_particles(gas)
        instance.synchronize_model()
        
        coords = [0.0 | units.kpc]*3
        speeds = [0.0 | units.m / units.s]*3
        hydro_state = instance.get_hydro_state_at_point(*(coords + speeds))
        expected = [ 3.5540e-19 | units.kg * units.m**-3, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                     8.445e-10 | units.kg * units.m**-1 * units.s**-2]
        for value, expect in zip(hydro_state, expected):
            self.assertAlmostRelativeEqual(value, expect, places=3)
        
        coords = [0.1 | units.kpc]*3
        speeds = [0.0 | units.m / units.s]*3
        hydro_state = instance.get_hydro_state_at_point(*(coords + speeds))
        expected = [ 4.1789e-19 | units.kg * units.m**-3, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                     1.0868e-9 | units.kg * units.m**-1 * units.s**-2]
        for value, expect in zip(hydro_state, expected):
            self.assertAlmostRelativeEqual(value, expect, places=3)
        instance.stop()
    
    def test15(self):
        print "Testing Fi get_hydro_state_at_point II: uniform sphere"
        number_sph_particles = 1000
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.kpc, 1.0e10 | units.MSun)
        gas = new_uniform_spherical_particle_distribution(number_sph_particles, 1.0 | units.kpc, 1.0e10 | units.MSun)
        gas.velocity = [10.0, 20.0, 30.0] | units.km / units.s
        gas.h_smooth = 0.01 | nbody_system.length
        gas.u = 0.05 | nbody_system.specific_energy
        density = (1.0e10 | units.MSun) / (4.0/3.0 * numpy.pi * (1.0 | units.kpc)**3)
        
        instance = Fi(convert_nbody)
        instance.parameters.n_smooth = 64
        instance.parameters.n_smooth_tol = 0.01
        instance.gas_particles.add_particles(gas)
        instance.synchronize_model()
        
        coords = [0.0 | units.kpc]*3
        speeds = [0.0 | units.m / units.s]*3
        
        rho, rhovx, rhovy, rhovz, rhoe = instance.get_hydro_state_at_point(*(coords + speeds))
        self.assertAlmostRelativeEqual(rho,   density, places=3)
        self.assertAlmostRelativeEqual(rho,   max(instance.gas_particles.rho),      places=2)
        self.assertIsOfOrder(          rho,   instance.gas_particles.rho)
        self.assertAlmostRelativeEqual(rhovx, density*instance.gas_particles[0].vx, places=3)
        self.assertAlmostRelativeEqual(rhovy, density*instance.gas_particles[0].vy, places=3)
        self.assertAlmostRelativeEqual(rhovz, density*instance.gas_particles[0].vz, places=3)
        self.assertAlmostRelativeEqual(rhoe,  density * (instance.gas_particles[0].u + 
            0.5 * instance.gas_particles[0].velocity.length_squared()),  places=2)
        
        coords = [0.1 | units.kpc]*3
        rho, rhovx, rhovy, rhovz, rhoe = instance.get_hydro_state_at_point(*(coords + speeds))
        self.assertAlmostRelativeEqual(rho,   density, places=3)
        self.assertAlmostRelativeEqual(rhovx, density*instance.gas_particles[0].vx, places=3)
        self.assertAlmostRelativeEqual(rhovy, density*instance.gas_particles[0].vy, places=3)
        self.assertAlmostRelativeEqual(rhovz, density*instance.gas_particles[0].vz, places=3)
        self.assertAlmostRelativeEqual(rhoe,  density * (instance.gas_particles[0].u + 
             0.5 * instance.gas_particles[0].velocity.length_squared()),  places=2)
        instance.stop()
        
    
    def test16(self):
        instance = self.new_instance_of_an_optional_code(Fi)
        try:
            particles = new_plummer_model(100)
            more_particles = new_plummer_model(50)
          
            particles.h_smooth = 0.1 | nbody_system.length 
            particles.u = 0.1 | nbody_system.speed**2
            more_particles.h_smooth = 0.1 | nbody_system.length 
            more_particles.u = 0.1 | nbody_system.speed**2
            
            particles.move_to_center()
            more_particles.move_to_center()
            
            instance.gas_particles.add_particles(particles)
            instance.commit_particles()
            self.assertEquals(len(instance.particles), 100)
            instance.synchronize_model()
            instance.gas_particles.add_particles(more_particles)
            instance.recommit_particles()
            self.assertEquals(len(instance.particles), 150)
            instance.synchronize_model()
            selected1 = instance.particles.select_array(lambda x: x<0 | nbody_system.length, ["x",])
            instance.particles.remove_particles(selected1)
            self.assertEquals(len(instance.particles), 150 - len(selected1))
            instance.recommit_particles()
            instance.synchronize_model()
            
            selected2 = instance.particles.select_array(lambda x: x>0 | nbody_system.length, ["x",])
            instance.particles.remove_particles(selected2)
            self.assertEquals(len(instance.particles), 150 - len(selected1) - len(selected2))
            instance.recommit_particles()
        finally:
            instance.stop()
    
    def test17(self):
        UnitLength = 3.085678e21 | units.cm     # ~ 1.0 kpc
        UnitMass = 1.989e43 | units.g           # 1.0e10 solar masses
        UnitVelocity = 1e5 | units.cm / units.s # 1 km/sec
        convert_nbody = nbody_system.nbody_to_si(UnitLength, UnitMass)
        instance = Fi(convert_nbody, mode = FiInterface.MODE_PERIODIC_BOUNDARIES)
        self.assertEqual(instance.parameters.periodic_boundaries_flag, True)
        instance.parameters.use_hydro_flag = False
        instance.parameters.self_gravity_flag = False
        instance.parameters.periodic_box_size = 2.0 | nbody_system.length
        instance.parameters.timestep = 0.1 | nbody_system.time
        self.assertAlmostEqual(instance.parameters.periodic_box_size, 2.0 | units.kpc, places=6)
        
        three_particles_IC = datamodel.Particles(3)
        three_particles_IC.position = [[0.5, 0.0, 0.0], [0.0,-0.5, 0.0], [0.0, 0.0, 0.5]] | nbody_system.length
        three_particles_IC.velocity =[[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0,-1.0]] | nbody_system.speed
        three_particles_IC.radius = 0.0 | units.RSun
        three_particles_IC.mass = 1.0e10 | units.MSun
        
        instance.dm_particles.add_particles(three_particles_IC)
        self.assertAlmostEqual(instance.dm_particles.x, [0.5, 0.0, 0.0] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.y, [0.0,-0.5, 0.0] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.z, [0.0, 0.0, 0.5] | units.kpc, places=6)
        
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertAlmostEqual(instance.dm_particles.x, [0.4, 0.0, 0.0] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.y, [0.0,-0.4, 0.0] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.z, [0.0, 0.0, 0.4] | units.kpc, places=6)
        
        instance.evolve_model(1.0 | nbody_system.time)
        self.assertAlmostEqual(instance.dm_particles.x, [-0.5, 0.0, 0.0] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.y, [ 0.0, 0.5, 0.0] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.z, [ 0.0, 0.0,-0.5] | units.kpc, places=6)
        instance.stop()
        

    def test18(self):
        particles = datamodel.Particles(10)
        particles.x = (numpy.array(range(10)) * 1.0) | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.radius = 0.001 | nbody_system.length
        particles.vx =  0 | nbody_system.speed
        particles.vy =  0 | nbody_system.speed
        particles.vz =  0 | nbody_system.speed
        particles.mass = 1.0 | nbody_system.mass
    
        instance = Fi()
        instance.initialize_code()
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles[0:2]), 2)
        self.assertTrue(str(instance.particles[0:2]).find('key') > 0)
        instance.stop()
        
    def test19(self):
        
        instance = Fi()
        instance.dm_particles.add_particle(datamodel.Particle(
            x = 1.0  | nbody_system.length,
            y = 2.0 | nbody_system.length,
            z = 3.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 2.0 | nbody_system.speed,
            vz = 3.0 | nbody_system.speed,
            mass = 1.0 | nbody_system.mass,
        )) 
        instance.gas_particles.add_particle(datamodel.Particle(
            x = 1.0  | nbody_system.length,
            y = 2.0 | nbody_system.length,
            z = 4.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 2.0 | nbody_system.speed,
            vz = 3.0 | nbody_system.speed,
            mass = 1.0 | nbody_system.mass,
            u = 1.0 | nbody_system.potential
        )) 
        instance.star_particles.add_particle(datamodel.Particle(
            x = 1.0  | nbody_system.length,
            y = 2.0 | nbody_system.length,
            z = 5.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 2.0 | nbody_system.speed,
            vz = 3.0 | nbody_system.speed,
            mass = 1.0 | nbody_system.mass,
        )) 
        instance.commit_particles()
        self.assertEquals(instance.dm_particles[0].radius, 8.0 | nbody_system.length)
        self.assertAlmostRelativeEquals(instance.gas_particles[0].h_smooth, 116.754483949 | nbody_system.length, 6)
        self.assertEquals(instance.star_particles[0].radius, 8.0 | nbody_system.length)
        
        instance.stop()
        instance = Fi()
        instance.parameters.adaptive_smoothing_flag = False
        instance.dm_particles.add_particle(datamodel.Particle(
            x = 1.0  | nbody_system.length,
            y = 2.0 | nbody_system.length,
            z = 3.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 2.0 | nbody_system.speed,
            vz = 3.0 | nbody_system.speed,
            mass = 1.0 | nbody_system.mass,
            radius = 0.1 | nbody_system.length,
        )) 
        instance.gas_particles.add_particle(datamodel.Particle(
            x = 1.0  | nbody_system.length,
            y = 2.0 | nbody_system.length,
            z = 4.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 2.0 | nbody_system.speed,
            vz = 3.0 | nbody_system.speed,
            mass = 1.0 | nbody_system.mass,
            u = 1.0 | nbody_system.potential,
            h_smooth = 100 | nbody_system.length,
            
        )) 
        instance.star_particles.add_particle(datamodel.Particle(
            x = 1.0  | nbody_system.length,
            y = 2.0 | nbody_system.length,
            z = 5.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 2.0 | nbody_system.speed,
            vz = 3.0 | nbody_system.speed,
            mass = 1.0 | nbody_system.mass,
            radius = 0.2 | nbody_system.length,
        )) 
        instance.commit_particles()
        self.assertEquals(instance.dm_particles[0].radius, 0.1 | nbody_system.length)
        self.assertAlmostRelativeEquals(instance.gas_particles[0].h_smooth, 116.754483949 | nbody_system.length, 6)
        self.assertEquals(instance.star_particles[0].radius,0.2| nbody_system.length)
        instance.stop()
        
    def test20(self):
        instance=Fi()
        instance.initialize_code()
        instance.parameters.use_hydro_flag = 0
        instance.parameters.self_gravity_flag = 0
        instance.parameters.periodic_box_size = 1 | nbody_system.length
        instance.parameters.timestep = 0.2 | nbody_system.time
        
        p1 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 0.5 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.1 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        p2 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 0.1 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 0.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.001 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        
        instance.evolve_model(0.1 |nbody_system.time)
        instance.update_particle_set()
        self.assertEquals(len(instance.gas_particles), 2)
        
        self.assertAlmostRelativeEquals(instance.gas_particles[0].x, 0.5 | nbody_system.length)
        instance.evolve_model(0.15 |nbody_system.time)
        instance.update_particle_set()
        self.assertEquals(len(instance.gas_particles), 1)
        self.assertEquals(instance.gas_particles.mass, [0.001] | nbody_system.mass)
        self.assertEquals(instance.gas_particles[0].mass, 0.001 | nbody_system.mass)
        self.assertAlmostRelativeEquals(instance.gas_particles[0].x, 0.1 | nbody_system.length)
        self.assertAlmostRelativeEquals(instance.gas_particles.x, [0.1] | nbody_system.length)
        

    def test21(self):
        instance=Fi()
        instance.parameters.use_hydro_flag = 0
        instance.parameters.self_gravity_flag = 0
        instance.parameters.periodic_box_size = 1 | nbody_system.length
        instance.parameters.timestep = 0.2 | nbody_system.time
        instance.commit_parameters()
        instance.parameters.use_hydro_flag =1
        particles = datamodel.Particles(2)
        particles.mass = 0.5 | nbody_system.mass
        particles[0].position = [1.0,0.0,0.0] | nbody_system.length
        particles[1].position = [-1.0,0.0,0.] | nbody_system.length
        particles.velocity = [0,0,0] | nbody_system.speed
        particles.u = 0.0 | nbody_system.potential
        print particles
        
        try:
            t = instance.model_time 
            self.assertEquals(t, 0 | nbody_system.time)
            channel = instance.particles.new_channel_to(particles)
            self.assertEquals(0, len(instance.particles))
            particles.synchronize_to(instance.particles)
            self.assertEquals(2, len(instance.particles))
            instance.commit_particles()
            
            self.assertEquals(2, len(particles))
            instance.particles.synchronize_to(particles)
            self.assertEquals(2, len(particles))
            channel.copy()
        except Exception as ex:
            print ex
            self.fail("requesting the model time should not trigger a commit_particles")
            
    def test22(self):
        instance=Fi()
        instance.initialize_code()
        instance.parameters.use_hydro_flag = 0
        instance.parameters.self_gravity_flag = 0
        instance.parameters.periodic_box_size = 1 | nbody_system.length
        instance.parameters.timestep = 0.2 | nbody_system.time
        
        p1 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 0.5 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.1 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        p2 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 0.1 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 0.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.001 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        p3 = instance.gas_particles.add_particle(datamodel.Particle(
            x = -0.9 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = -2.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.002 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        
        instance.evolve_model(0.1 |nbody_system.time)
        instance.update_particle_set()
        self.assertEquals(len(instance.gas_particles), 3)
        
        self.assertAlmostRelativeEquals(instance.gas_particles[0].x, 0.5 | nbody_system.length)
        instance.evolve_model(0.15 |nbody_system.time)
        instance.update_particle_set()
        self.assertEquals(instance.gas_particles.mass, [0.001] | nbody_system.mass)
        self.assertEquals(instance.gas_particles[0].mass, 0.001 | nbody_system.mass)
        self.assertEquals(len(instance.gas_particles), 1)
        self.assertAlmostRelativeEquals(instance.gas_particles[0].mass, 0.001 | nbody_system.mass)
        self.assertAlmostRelativeEquals(instance.gas_particles[0].x, 0.1 | nbody_system.length)
        self.assertAlmostRelativeEquals(instance.gas_particles.x, [0.1] | nbody_system.length)
        

            
    def test23(self):
        instance=Fi()
        instance.initialize_code()
        instance.parameters.use_hydro_flag = 0
        instance.parameters.self_gravity_flag = 0
        instance.parameters.periodic_box_size = 1 | nbody_system.length
        instance.parameters.timestep = 0.2 | nbody_system.time
        
        p1 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 1.5 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.1 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        p2 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 0.1 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 0.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.001 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        p3 = instance.gas_particles.add_particle(datamodel.Particle(
            x = -1.5 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = -2.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.002 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        instance.commit_particles()
        self.assertEquals(len(instance.gas_particles), 3)
        self.assertEquals(instance.gas_particles.mass, [0.1, 0.001, 0.002] | nbody_system.mass)
        
        instance.evolve_model(0.15 |nbody_system.time)
        instance.update_particle_set()
        self.assertEquals(len(instance.gas_particles), 1)
        self.assertEquals(instance.gas_particles.mass, [0.001] | nbody_system.mass)
        self.assertEquals(instance.gas_particles[0].mass, 0.001 | nbody_system.mass)
        self.assertEquals(len(instance.gas_particles), 1)
        self.assertAlmostRelativeEquals(instance.gas_particles[0].mass, 0.001 | nbody_system.mass)
        self.assertAlmostRelativeEquals(instance.gas_particles[0].x, 0.1 | nbody_system.length)
        self.assertAlmostRelativeEquals(instance.gas_particles.x, [0.1] | nbody_system.length)
        
    def test24(self):
        instance=Fi()
        instance.initialize_code()
        instance.parameters.use_hydro_flag = 0
        instance.parameters.self_gravity_flag = 0
        instance.parameters.periodic_box_size = 1 | nbody_system.length
        instance.parameters.timestep = 0.2 | nbody_system.time
        
        p1 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 0.5 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.1 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        p2 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 0.1 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 0.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.001 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        #run more than one innerloop step in fi
        instance.evolve_model(1|nbody_system.time)
        instance.update_particle_set()
        self.assertEquals(len(instance.gas_particles), 1)
        self.assertEquals(instance.gas_particles.mass, [0.001] | nbody_system.mass)
        self.assertEquals(instance.gas_particles[0].mass, 0.001 | nbody_system.mass)
        self.assertAlmostRelativeEquals(instance.gas_particles[0].x, 0.1 | nbody_system.length)
        self.assertAlmostRelativeEquals(instance.gas_particles.x, [0.1] | nbody_system.length)
        
    def test25(self):
        print "Testing Fi density_limit_detection"
        number_gas_particles = 500
        UnitLength = 3.085678e21 | units.cm     # ~ 1.0 kpc
        UnitMass = 1.989e43 | units.g           # 1.0e10 solar masses
        convert_nbody = nbody_system.nbody_to_si(UnitLength, UnitMass)
        gas = new_evrard_gas_sphere(number_gas_particles, convert_nbody, do_scale=True, seed=12345)
        
        instance = Fi(convert_nbody)
        instance.initialize_code()
        instance.parameters.stopping_condition_maximum_density = 10.0 * UnitMass / UnitLength**3
        instance.parameters.timestep = 0.1 | units.Myr
        instance.gas_particles.add_particles(gas)
        stars = new_plummer_model(5, convert_nbody)
        stars.x += 1000 * UnitLength
        instance.dm_particles.add_particles(stars)
        self.assertIsOfOrder(max(instance.gas_particles.density), UnitMass / UnitLength**3)
        
        density_limit_detection = instance.stopping_conditions.density_limit_detection
        density_limit_detection.enable()
        
        instance.evolve_model(10.0 | units.Myr)
        print instance.model_time.as_quantity_in(units.Myr)
        print instance.stopping_conditions
        self.assertTrue(density_limit_detection.is_set())
        self.assertTrue(instance.model_time < 10.0 | units.Myr)
        self.assertEquals(len(density_limit_detection.particles()), 1)
        self.assertTrue((density_limit_detection.particles().density > 
                10 * UnitMass / UnitLength**3).all())
        
        instance.particles.remove_particles(density_limit_detection.particles())
        
        instance.evolve_model(10.0 | units.Myr)
        print instance.model_time.as_quantity_in(units.Myr)
        print instance.stopping_conditions
        self.assertTrue(density_limit_detection.is_set())
        self.assertTrue(instance.model_time < 10.0 | units.Myr)
        
        self.assertEquals(len(density_limit_detection.particles()), 1)
        self.assertTrue((density_limit_detection.particles().density > 
                10 * UnitMass / UnitLength**3).all())
        instance.stop()
    
    def test26(self):
        print "Testing Fi internal_energy_limit_detection"
        number_gas_particles = 500
        UnitLength = 3.085678e21 | units.cm     # ~ 1.0 kpc
        UnitMass = 1.989e43 | units.g           # 1.0e10 solar masses
        convert_nbody = nbody_system.nbody_to_si(UnitLength, UnitMass)
        gas = new_evrard_gas_sphere(number_gas_particles, convert_nbody, do_scale=True, seed=12345)
        initial_internal_energy = 0.05 * constants.G * UnitMass / UnitLength
        
        instance = Fi(convert_nbody)
        instance.initialize_code()
        instance.parameters.timestep = 0.1 | units.Myr
        instance.parameters.stopping_condition_maximum_internal_energy = 10 * initial_internal_energy
        instance.gas_particles.add_particles(gas)
        self.assertAlmostRelativeEquals(instance.gas_particles.u, initial_internal_energy, 8)
        
        internal_energy_limit_detection = instance.stopping_conditions.internal_energy_limit_detection
        internal_energy_limit_detection.enable()
        
        instance.evolve_model(10.0 | units.Myr)
        print instance.model_time.as_quantity_in(units.Myr)
        print instance.stopping_conditions
        self.assertTrue(internal_energy_limit_detection.is_set())
        self.assertTrue(instance.model_time < 10.0 | units.Myr)
        self.assertEquals(len(internal_energy_limit_detection.particles()), 8)
        self.assertTrue((internal_energy_limit_detection.particles().u > 
                10 * initial_internal_energy).all())
        
        instance.particles.remove_particles(internal_energy_limit_detection.particles())
        
        instance.evolve_model(10.0 | units.Myr)
        print instance.model_time.as_quantity_in(units.Myr)
        print instance.stopping_conditions
        self.assertTrue(internal_energy_limit_detection.is_set())
        self.assertTrue(instance.model_time < 10.0 | units.Myr)
        self.assertEquals(len(internal_energy_limit_detection.particles()), 5)
        self.assertTrue((internal_energy_limit_detection.particles().u > 
                10 * initial_internal_energy).all())
        instance.stop()
    
    def test27(self):
        print "Testing Fi SPH particle properties"
        number_sph_particles = 1000
        convert_nbody = nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        
        gas = new_evrard_gas_sphere(number_sph_particles, convert_nbody, seed = 1234)
        instance = Fi(convert_nbody)
        instance.parameters.timestep = 0.05 | nbody_system.time
        instance.parameters.integrate_entropy_flag = True
        instance.gas_particles.add_particles(gas)
        self.assertIsOfOrder(instance.gas_particles.h_smooth, 
            convert_nbody.to_si(0.5 | nbody_system.length) * 
            (instance.parameters.n_smooth*1.0/number_sph_particles)**(1.0/3))
        self.assertAlmostRelativeEqual(instance.gas_particles.u, 
            convert_nbody.to_si(0.05 | nbody_system.specific_energy))
        
        # the density of the cloud scales with 1/r:
        r_sort, rho_sort = instance.gas_particles.position.lengths().sorted_with(instance.gas_particles.rho)
        mean_density = convert_nbody.to_si(3.0/(4.0*numpy.pi) | nbody_system.density)
        select = slice(number_sph_particles/2) # select 50% particles closest to center to avoid boundaries
        self.assertIsOfOrder(rho_sort[select]/mean_density, r_sort.mean()/r_sort[select])
        
        self.assertAlmostEqual(instance.gas_particles.u * instance.gas_particles.rho, 
            1.5 * instance.gas_particles.pressure)
        
        self.assertAlmostEqual(instance.gas_particles.du_dt, 0 | units.m**2 * units.s**-3)
        u_0 = instance.gas_particles.u.sum()
        instance.evolve_model(0.1 | nbody_system.time)
        # Collapsing ==> heating:
        self.assertTrue(instance.gas_particles.du_dt.sum() >= 0 | units.m**2 * units.s**-3)
        self.assertIsOfOrder(instance.gas_particles.du_dt.sum() * convert_nbody.to_si(0.1 | nbody_system.time), 
            instance.gas_particles.u.sum() - u_0)
        instance.stop()

    def test28(self):
        print "Testing Fi u error"
        number_gas_particles = 500
        UnitLength = 3.085678e21 | units.cm     # ~ 1.0 kpc
        UnitMass = 1.989e43 | units.g           # 1.0e10 solar masses
        convert_nbody = nbody_system.nbody_to_si(UnitLength, UnitMass)
        gas = new_evrard_gas_sphere(number_gas_particles, convert_nbody, do_scale=True, seed=12345)
        
        instance = Fi(convert_nbody)
        instance.initialize_code()
        instance.parameters.stopping_condition_maximum_density = 10.0 * UnitMass / UnitLength**3
        instance.parameters.timestep = 0.1 | units.Myr
        instance.gas_particles.add_particles(gas)

        u=instance.gas_particles.u
        gas_particles=instance.gas_particles.copy()
        
        self.assertEqual(u,gas_particles.u)
        
        
    def test29(self):
        instance=Fi()
        instance.initialize_code()
        instance.parameters.use_hydro_flag = 0
        instance.parameters.self_gravity_flag = 1
        instance.parameters.epsilon_squared = 0.000000001 |  nbody_system.length**2
        instance.parameters.timestep = 0.2 | nbody_system.time
        
        x0 = 0.5 | nbody_system.length
        
        p1 = instance.star_particles.add_particle(datamodel.Particle(
            x = x0,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 0.0 | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.01 | nbody_system.mass,
            
        )) 
        p2 = instance.star_particles.add_particle(datamodel.Particle(
            x = 10| nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 0.0 | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 100 | nbody_system.mass,
            
        )) 
        #run more than one innerloop step in fi
        e0 = instance.kinetic_energy + instance.potential_energy
        instance.evolve_model(0.5|nbody_system.time)
        e1 = instance.kinetic_energy + instance.potential_energy
        self.assertAlmostRelativeEquals(e0, e1, 4)
        print instance.star_particles[0].position.x
        
        self.assertAlmostRelativeEquals(instance.star_particles[0].position.x, 0.588896329162| nbody_system.length, 12)
        
        for i in range(10):
            p3 = instance.star_particles.add_particle(datamodel.Particle(
                x = 100| nbody_system.length,
                y = 0.0 | nbody_system.length,
                z = 0.0 | nbody_system.length,
                vx = 0.0 | nbody_system.speed,
                vy = 0.0 | nbody_system.speed,
                vz = 0.0 | nbody_system.speed,
                mass = 100 | nbody_system.mass,
                
            )) 
            instance.recommit_particles()
            instance.star_particles.remove_particle(p3)
            instance.recommit_particles()
        
        e11 = instance.kinetic_energy + instance.potential_energy
        self.assertAlmostRelativeEquals(e1, e1, 8)
        instance.evolve_model(1.0|nbody_system.time)
        e2 = instance.kinetic_energy + instance.potential_energy
        self.assertAlmostRelativeEquals(e0, e1, 4)
        self.assertAlmostRelativeEquals(instance.star_particles[0].position.x, 1.06498184266 | nbody_system.length, 11)
       
    def test30(self):
        instance=Fi()
        instance.initialize_code()
        instance.parameters.use_hydro_flag = 0
        instance.parameters.self_gravity_flag = 1
        instance.parameters.epsilon_squared = 0.000000001 |  nbody_system.length**2
        instance.parameters.timestep = 0.2 | nbody_system.time
        instance.parameters.verbosity = 0
        
        x0 = 0.5 | nbody_system.length
        
        p1 = instance.star_particles.add_particle(datamodel.Particle(
            x = x0,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 0.0 | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.01 | nbody_system.mass,
            
        )) 
        p2 = instance.star_particles.add_particle(datamodel.Particle(
            x = 10| nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 0.0 | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 100 | nbody_system.mass,
            
        )) 
        instance.commit_particles()
        e0 = instance.kinetic_energy + instance.potential_energy
        print "e0:", e0
        for i in range(10):
            p3 = instance.star_particles.add_particle(datamodel.Particle(
                x = 100| nbody_system.length,
                y = 0.0 | nbody_system.length,
                z = 0.0 | nbody_system.length,
                vx = 0.0 | nbody_system.speed,
                vy = 0.0 | nbody_system.speed,
                vz = 0.0 | nbody_system.speed,
                mass = 100 | nbody_system.mass,
                
            )) 
            instance.recommit_particles()
            instance.star_particles.remove_particle(p3)
            instance.recommit_particles()
        #run more than one innerloop step in fi
        instance.evolve_model(0.5|nbody_system.time)
        
        e1 = instance.kinetic_energy + instance.potential_energy
        self.assertAlmostRelativeEquals(e0, e1, 4)
        print instance.star_particles[0].position.x
        
        self.assertAlmostRelativeEquals(instance.star_particles[0].position.x, 0.588896329162| nbody_system.length, 12)
        
        instance.evolve_model(1.0|nbody_system.time)
        e2 = instance.kinetic_energy + instance.potential_energy
        self.assertAlmostRelativeEquals(e0, e2,4)
        self.assertAlmostRelativeEquals(instance.star_particles[0].position.x, 1.06493331332 | nbody_system.length, 11)
       
    
    def test31(self):
        instance=Fi()
        instance.initialize_code()
        instance.parameters.use_hydro_flag = 0
        instance.parameters.self_gravity_flag = 0
        instance.parameters.periodic_box_size = 1 | nbody_system.length
        instance.parameters.timestep = 0.2 | nbody_system.time
        
        p1 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 0.5 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 1.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.1 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        p2 = instance.gas_particles.add_particle(datamodel.Particle(
            x = 0.1 | nbody_system.length,
            y = 0.0 | nbody_system.length,
            z = 0.0 | nbody_system.length,
            vx = 0.0  | nbody_system.speed,
            vy = 0.0 | nbody_system.speed,
            vz = 0.0 | nbody_system.speed,
            mass = 0.001 | nbody_system.mass,
            u = 0.0 | nbody_system.potential
            
        )) 
        
        instance.evolve_model(0.1 |nbody_system.time)
        instance.synchronize_model()
        self.assertEquals(len(instance.gas_particles), 2)
        
        self.assertAlmostRelativeEquals(instance.gas_particles[0].x, 0.5 | nbody_system.length)
        instance.evolve_model(0.15 |nbody_system.time)
        instance.synchronize_model()
        self.assertEquals(len(instance.gas_particles), 1)
        self.assertEquals(instance.gas_particles.mass, [0.001] | nbody_system.mass)
        self.assertEquals(instance.gas_particles[0].mass, 0.001 | nbody_system.mass)
        self.assertAlmostRelativeEquals(instance.gas_particles[0].x, 0.1 | nbody_system.length)
        self.assertAlmostRelativeEquals(instance.gas_particles.x, [0.1] | nbody_system.length)
        
        
        

        
        
        
