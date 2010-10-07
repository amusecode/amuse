import os
import sys
import numpy
import time
from amuse.test.amusetest import TestWithMPI

from amuse.legacy.fi.interface import FiInterface, Fi
from amuse.ext.evrard_test import new_evrard_gas_sphere
from amuse.ext.spherical_model import new_uniform_spherical_particle_distribution
from amuse.support.units import nbody_system as nbody
from amuse.support.units import units
from amuse.support.data import core
from amuse.support.legacy import channel

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
                    ('consph',1),('sphinit',1),('uentropy',0),('isotherm',0),
                    ('eps_is_h',1)]:
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
        
        instance.stop()
    
    def test3(self):
        instance=FiInterface()
        instance.initialize_code()
        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        instance.cleanup_code()
        instance.stop()
    
    def test4(self):
        instance=FiInterface()
        instance.initialize_code()
        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
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
        self.assertEqual(instance.get_number_of_particles().number_of_particles, 3)
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
        instance.commit_parameters()
        instance.new_sph_particle(mass, h, x, y, z, vx, vy, vz, u)
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

    def xtest7(self):
        instance=FiInterface()
        instance.initialize_code()
#        instance.set_periodic(1)
        instance.set_use_hydro(0)
        instance.set_selfgrav(0)
        instance.set_pboxsize(2.)
        instance.set_dtime(0.1)
        instance.commit_parameters()
        ids,err=instance.new_particle( 
           [1.0,1.0,1.0],
           [0.0,0.0,0.0],
           [0.5,0.0,0.0],
           [0.0,-0.5,0.0],
           [0.0,0.0,0.5],
           [1.0,0.0,0.0],
           [0.0,-1.0,0.0],
           [0.0,0.0,1.0])
        instance.commit_particles()
        
        instance.evolve(0.1)
        m,r,x,y,z,vx,vy,vz,err=instance.get_state(ids)
        self.assertAlmostEqual(x, [0.6,0.,0.], places=7)
        self.assertAlmostEqual(y, [0.,-0.6,0.], places=7)
        self.assertAlmostEqual(z, [0.,0.,0.6], places=7)
        instance.evolve(1.0)
        m,r,x,y,z,vx,vy,vz,err=instance.get_state(ids)
        self.assertAlmostEqual(x, [-0.5,0.,0.], places=7)
        self.assertAlmostEqual(y, [0.,0.5,0.], places=7)
        self.assertAlmostEqual(z, [0.,0.,-0.5], places=7)        
        instance.cleanup_code()
        instance.stop()


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
        
        nb.evolve(0.5)
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
        instance=Fi(nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
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
        convert_nbody = nbody.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = Fi(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00000001 | units.AU**2
        instance.parameters.timestep = 0.5 | units.day
        self.assertEquals(instance.parameters.timestep, 0.5 | units.day)
        instance.commit_parameters()
        
        stars = core.Particles(2)
        
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
        
        instance.particles.copy_values_of_state_attributes_to(stars)
        
        postion_after_full_rotation = earth.position.x
        
        self.assertAlmostRelativeEqual(postion_at_start, postion_after_full_rotation, 4)
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.particles.copy_values_of_state_attributes_to(stars)
        
        postion_after_half_a_rotation = earth.position.x
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.particles.copy_values_of_state_attributes_to(stars)
        
        postion_after_half_a_rotation = earth.position.x
        self.assertAlmostRelativeEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        instance.particles.copy_values_of_state_attributes_to(stars)
        
        postion_after_half_a_rotation = earth.position.y
        
        self.assertAlmostRelativeEqual(-postion_at_start, postion_after_half_a_rotation, 4)
        
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Test 2: testing Fi data directory"
        convert_nbody = nbody.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Fi(convert_nbody)
        error = instance.initialize_code()
        self.assertEquals(0, error)
        self.assertTrue('data/fi/input/' in instance.get_fi_data_directory().fi_data_directory)
        self.assertEquals(instance.get_fi_data_directory().fi_data_directory, 
            instance.get_data_directory()+'/')
        
        self.assertEquals(False, instance.parameters.radiation_flag)
        instance.parameters.radiation_flag = True
        self.assertEquals(True, instance.parameters.radiation_flag)
        
        self.assertEquals(False, instance.parameters.star_formation_flag)
        instance.parameters.star_formation_flag = True
        self.assertEquals(True, instance.parameters.star_formation_flag)
        
        stars = core.Particles(2)
        stars.mass = [1.0, 3.0e-6] | units.MSun
        stars.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.AU
        stars.velocity = [[0.0,0.0,0.0], [0.0,29.8,0.0]] | units.km / units.s
        stars.radius = [1.0, 0.01] | units.RSun
        
        instance.particles.add_particles(stars)
        instance.evolve_model(1.0 | nbody.time)
        
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):
        print "Test 3: testing Fi boolean parameters"
        instance = Fi(nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        self.assertEquals(0, instance.initialize_code())
        
        for bool_par in ['radiation_flag','star_formation_flag',
            'square_root_timestep_flag','freeform_timestep_flag','quadrupole_moments_flag',
            'direct_sum_flag','fixed_halo_flag','adaptive_smoothing_flag','smooth_input_flag',
            'integrate_entropy_flag','isothermal_flag']:
            self.assertEquals(False, eval("instance.parameters."+bool_par))
            exec("instance.parameters."+bool_par+" = True")
            self.assertEquals(True, eval("instance.parameters."+bool_par))
        
        for bool_par in ['acc_timestep_flag','self_gravity_flag','gadget_cell_opening_flag',
            'use_hydro_flag','conservative_sph_flag','sph_dens_init_flag','eps_is_h_flag']:
            self.assertEquals(True, eval("instance.parameters."+bool_par))
            exec("instance.parameters."+bool_par+" = False")
            self.assertEquals(False, eval("instance.parameters."+bool_par))
        
        instance.cleanup_code()
        instance.stop()
    
    def test4(self):
        print "Test 4: testing Fi integer parameters"
        instance = Fi(nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        self.assertEquals(0, instance.initialize_code())
        
        for int_par, value in [('first_snapshot',0),('output_interval',5),('log_interval',5),
            ('maximum_time_bin',4096),('minimum_part_per_bin',1),('targetnn',32),
            ('verbosity',0),('n_smooth',64)]:
            self.assertEquals(value | units.none, eval("instance.parameters."+int_par))
            exec("instance.parameters."+int_par+" = 1 | units.none")
            self.assertEquals(1 | units.none, eval("instance.parameters."+int_par))
        
        instance.cleanup_code()
        instance.stop()
    
    def test5(self):
        print "Test 5: testing Fi double precision parameters"
        instance = Fi(nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        self.assertEquals(0, instance.initialize_code())
        
        par_names=['epsilon_squared','timestep','pboxsize','code_mass_unit','code_length_unit',
            'sqrt_timestep_crit_constant','acc_timestep_crit_constant','free_timestep_crit_constant_v',
            'free_timestep_crit_constant_a','free_timestep_crit_constant_vexp',
            'free_timestep_crit_constant_aexp','opening_angle','gadget_cell_opening_constant',
            'nn_tol','gas_epsilon','gamma','artificial_viscosity_alpha','beta','sph_artificial_viscosity_eps','courant',
            'min_gas_part_mass','sph_h_const','n_smooth_tol','grain_heat_eff','zeta_cr_ion_rate',
            'heat_par1','heat_par2','cool_par','optical_depth','star_form_delay_fac','star_form_mass_crit',
            'star_form_eff','supernova_duration','supernova_eff','t_supernova_start','max_density']
        defaults=[0.0 | nbody.length * nbody.length, 1.0 | nbody.time, 300.0 | nbody.length, 
            1.0e9 | units.MSun, 1.0 | units.kpc, 1.0, 0.25, 0.5, 0.35, 
            0.0, -1.0, 0.5, 0.01, 0.1, 0.005 | nbody.length, 1.6666667, 0.5, 1.0, 0.01, 0.3, 
            0.25, 0.2 | nbody.length, 0.1, 0.05, 3.6 | 1.8e-17 * units.s**-1, 0.0, 0.0, 1.0, 
            0.0, 1.0, 1.0e5 | units.MSun, 0.25, 3.0e7 | units.Myr, 0.0, 3.e6 | units.Myr, 100.0]
        defaults = [val | units.none if isinstance(val,float) else val for val in defaults]
        defaults = [instance.unit_converter.to_si(val) if nbody.is_nbody_unit(val.unit) else val for val in defaults]
        for double_par, value in zip(par_names, defaults):
            self.assertAlmostRelativeEquals(eval("instance.parameters."+double_par), value, 7)
            exec("instance.parameters."+double_par+" = 2 * value")
            self.assertAlmostRelativeEquals(eval("instance.parameters."+double_par), 2*value, 7)
        instance.cleanup_code()
        instance.stop()
    
    def test6(self):
        print "Test 6: testing Fi string parameters"
        instance = Fi(nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        self.assertEquals(0, instance.initialize_code())
        
        par_names=['halofile','feedback','star_formation_mode','h_update_method','sph_viscosity','fi_data_directory']
        defaults=['none','fuv','gerritsen','mass','sph',instance.get_data_directory()+'/'] | units.string
        new_values=['bct_02_10.halo','pres','nieuw','test','sphv','test'] | units.string
        for string_par, value, new_value in zip(par_names, defaults, new_values):
            self.assertEquals(eval("instance.parameters."+string_par), value)
            exec("instance.parameters."+string_par+" = new_value")
            self.assertEquals(eval("instance.parameters."+string_par), new_value)
        instance.cleanup_code()
        instance.stop()
    
    def test7(self):
        print "Test 7: testing Fi SPH particles"
        target_number_of_particles = 1000
        gas = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed = 1234)
        gas.h_smooth = 0.0 | nbody.length
        
        convert_nbody = nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        instance = Fi(convert_nbody)
        self.assertEquals(0, instance.initialize_code())
        instance.parameters.timestep = 0.05 | nbody.time
        instance.parameters.integrate_entropy_flag = True
        instance.gas_particles.add_particles(gas)
        
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.kinetic_energy),    0.00 | nbody.energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.potential_energy), -0.50 | nbody.energy, 2)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.05 | nbody.energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.total_energy),     -0.45 | nbody.energy, 2)
        
        instance.evolve_model(10.437 | units.Myr)
        instance.synchronize_model()
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.kinetic_energy),    0.0566 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.potential_energy), -0.5661 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.0628 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.total_energy),     -0.45 | nbody.energy, 2)
        
        instance.cleanup_code()
        instance.stop()
    
    def test8(self):
        print "Test 8: testing Fi dark matter + SPH particles"
        target_number_of_particles = 100
        gas = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed = 1234)
        gas.h_smooth = 0.0 | nbody.length
        
        dark = core.Particles(2)
        dark.mass = [1.0, 3.0e-6] | units.MSun
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.AU
        dark.velocity = [[0.0,0.0,0.0], [0.0,29.8,0.0]] | units.km / units.s
        dark.radius = [1.0, 0.01] | units.RSun
        
        convert_nbody = nbody.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Fi(convert_nbody)
        self.assertEquals(0, instance.initialize_code())
        instance.dm_particles.add_particles(dark)
        instance.gas_particles.add_particles(gas)
        
        print convert_nbody.to_nbody(instance.kinetic_energy)
        print convert_nbody.to_nbody(instance.potential_energy)
        print convert_nbody.to_nbody(instance.thermal_energy)
        print convert_nbody.to_nbody(instance.total_energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.kinetic_energy),    0.0000 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.potential_energy), -1.3837 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.0500 | nbody.energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.total_energy),     -1.3337 | nbody.energy, 3)
        
        instance.evolve_model(1.113988| units.yr)
        instance.synchronize_model()
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.kinetic_energy),    0.00321 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.potential_energy), -1.8387 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.41868 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.total_energy),     -1.4167 | nbody.energy, 3)
        
        instance.cleanup_code()
        instance.stop()
    
    def test9(self):
        print "Test 9: testing Fi dark matter + SPH + star particles"
        target_number_of_particles = 100
        gas = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed = 1234)
        gas.h_smooth = 0.0 | nbody.length
        
        dark = core.Particles(2)
        dark.mass = [0.4, 0.4] | nbody.mass
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.kpc
        dark.velocity = [[100.0,100.0,100.0], [1.0,1.0,1.0]] | units.km / units.s
        dark.radius = [0.0, 0.0] | units.RSun
        
        star = core.Particles(2)
        star.mass = [0.02, 0.02] | nbody.mass
        star.position = [[0.1,0.2,0.3], [0.4,0.5,0.6]] | units.kpc
        star.velocity = [[-300.0,-200.0,-100.0], [-6.0,-5.0,-4.0]] | units.km / units.s
        star.radius = [0.0, 0.0] | units.RSun
        star.tform = [1000.0, 1000.0] | units.Myr
        
        convert_nbody = nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        instance = Fi(convert_nbody)
        self.assertEquals(0, instance.initialize_code())
        instance.parameters.verbosity = 1
        instance.dm_particles.add_particles(dark)
        instance.star_particles.add_particles(star)
        instance.gas_particles.add_particles(gas)
        
        print convert_nbody.to_nbody(instance.kinetic_energy)
        print convert_nbody.to_nbody(instance.potential_energy)
        print convert_nbody.to_nbody(instance.thermal_energy)
        print convert_nbody.to_nbody(instance.total_energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.kinetic_energy),    1.7204 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.potential_energy), -1.2935 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.0500 | nbody.energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.total_energy),      0.4768 | nbody.energy, 3)
        
        instance.evolve_model(instance.parameters.timestep)
        instance.synchronize_model()
        print convert_nbody.to_nbody(instance.kinetic_energy)
        print convert_nbody.to_nbody(instance.potential_energy)
        print convert_nbody.to_nbody(instance.thermal_energy)
        print convert_nbody.to_nbody(instance.total_energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.kinetic_energy),    1.5057 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.potential_energy), -1.281 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.thermal_energy),    0.147 | nbody.energy, 3)
        self.assertAlmostEqual(convert_nbody.to_nbody(instance.total_energy),      0.37178 | nbody.energy, 3)
        
        instance.cleanup_code()
        instance.stop()
    
    def test10(self):
        print "Test 10: testing Fi star particles"
        stars = core.Particles(2)
        stars.mass = [1.0, 3.00e-6] | units.MSun
        stars.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.AU
        stars.velocity = [[0.0,0.0,0.0], [0.0,29.8,0.0]] | units.km / units.s
        stars.radius = [1.0, 1.0] | units.RSun # results are nonsense if not (r1 == r2)... ?
        stars.tform = [-10.0, -50.0] | units.Myr
        
        instance = Fi(nbody.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU))#, debugger='xterm')
        self.assertEquals(0, instance.initialize_code())
        instance.parameters.timestep = 0.5 | units.day
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
    
    def xtest11(self):
        print "Test 11: testing Fi (dm+sph+star) particles Superset"
        target_number_of_particles = 100
        gas = new_evrard_gas_sphere(target_number_of_particles, do_scale=True, seed = 1234)
        gas.h_smooth = 0.0 | nbody.length
        number_sph_particles = len(gas)
        
        dark = core.Particles(2)
        dark.mass = [0.4, 0.4] | nbody.mass
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.kpc
        dark.velocity = [[100.0,100.0,100.0], [1.0,1.0,1.0]] | units.km / units.s
        dark.radius = [0.0, 0.0] | units.RSun
        
        star = core.Particles(2)
        star.mass = [0.02, 0.02] | nbody.mass
        star.position = [[0.1,0.2,0.3], [0.4,0.5,0.6]] | units.kpc
        star.velocity = [[-300.0,-200.0,-100.0], [-6.0,-5.0,-4.0]] | units.km / units.s
        star.radius = [0.0, 0.0] | units.RSun
        star.tform = [1000.0, 1000.0] | units.Myr
        
        convert_nbody = nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        instance = Fi(convert_nbody)
        self.assertEquals(0, instance.initialize_code())

        self.assertEquals(0, len(instance.particles))
        instance.dm_particles.add_particles(dark)
        instance.star_particles.add_particles(star)
        instance.gas_particles.add_particles(gas)
        self.assertEquals(number_sph_particles+4, len(instance.particles))
        
        print "The 'particles' superset can still be used as standard GD 'particles' set."
        default = core.Particles(2)
        default.mass = [0.4, 0.4] | nbody.mass
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
        gas.h_smooth = 0.0 | nbody.length
        dark = core.Particles(2)
        dark.mass = [0.4, 0.4] | nbody.mass
        dark.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | nbody.length
        dark.velocity = [[0.1, 0.1, 0.1], [1.0, 1.0, 1.0]] | nbody.speed
        dark.radius = 0.0 | nbody.length
        
        print "First do everything manually:"
        convert_nbody = nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        instance = Fi(convert_nbody)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.gas_particles.add_particles(gas)
        instance.commit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        mass = instance.gas_particles[0].mass
        instance.evolve_model(0.001 | nbody.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()
        
        print "initialize_code(), commit_parameters(), (re)commit_particles(), and cleanup_code() should be called " \
            "automatically before setting parameters, new_xx_particle(), get_xx(), and stop():"
        instance = Fi(convert_nbody)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.timestep = 0.001 | nbody.time
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        self.assertEquals(instance.parameters.timestep, convert_nbody.to_si(0.001 | nbody.time))
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
        instance.evolve_model(0.002 | nbody.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
    
    def test13(self):
        particles = core.Particles(2)
        particles.x = [0.0,10.0] | nbody.length
        particles.y = 0 | nbody.length
        particles.z = 0 | nbody.length
        particles.radius = 0.005 | nbody.length
        particles.vx =  0 | nbody.speed
        particles.vy =  0 | nbody.speed
        particles.vz =  0 | nbody.speed
        particles.mass = 1.0 | nbody.mass

        instance = Fi()
        instance.initialize_code()
        instance.parameters.stopping_conditions_number_of_steps = 2
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps, 2|units.none)
        instance.parameters.epsilon_squared = (0.01 | nbody.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(10 | nbody.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        self.assertTrue(instance.model_time < 10 | nbody.time)

        instance.stop()

    def test13a(self):
        particles = core.Particles(2)
        particles.x = [0.0,10.0] | nbody.length
        particles.y = 0 | nbody.length
        particles.z = 0 | nbody.length
        particles.radius = 0.005 | nbody.length
        particles.vx =  0 | nbody.speed
        particles.vy =  0 | nbody.speed
        particles.vz =  0 | nbody.speed
        particles.mass = 1.0 | nbody.mass

        very_short_time_to_evolve = 1 | units.s
        very_long_time_to_evolve = 1e9 | nbody.time

        instance = Fi()
        instance.initialize_code()
        instance.parameters.stopping_conditions_timeout = very_short_time_to_evolve
        self.assertEquals(instance.parameters.stopping_conditions_timeout, very_short_time_to_evolve)
        instance.parameters.epsilon_squared = (0.01 | nbody.length)**2
        instance.particles.add_particles(particles) 
        instance.stopping_conditions.timeout_detection.enable()
        start = time.time()
        instance.evolve_model(very_long_time_to_evolve)
        end = time.time()
        self.assertTrue(instance.stopping_conditions.timeout_detection.is_set())
        self.assertTrue((end-start)<very_short_time_to_evolve.value_in(units.s) + 2)#2 = some overhead compensation

        instance.stop()

    def test13b(self):
        particles = core.Particles(2)
        particles.x = [0.0,0.0] | nbody.length
        particles.y = 0 | nbody.length
        particles.z = 0 | nbody.length
        particles.radius = 0.001 | nbody.length
        particles.vx =  [-0.5,1.] | nbody.speed
        particles.vy =  0 | nbody.speed
        particles.vz =  0 | nbody.speed
        particles.mass = 1.0 | nbody.mass

        instance = Fi()
        instance.initialize_code()
        instance.parameters.stopping_conditions_out_of_box_size = 1. | nbody.length
        instance.parameters.self_gravity_flag=False
        instance.parameters.timestep=0.25 | nbody.time 
        self.assertEquals(instance.parameters.stopping_conditions_out_of_box_size,  1.| nbody.length)
        instance.parameters.epsilon_squared = (0.001 | nbody.length)**2
        instance.dm_particles.add_particles(particles) 
        instance.stopping_conditions.out_of_box_detection.enable()
        instance.evolve_model(2 | nbody.time)
        self.assertEquals(instance.stopping_conditions.out_of_box_detection.number_of_particles()[0],2)
        print instance.stopping_conditions.out_of_box_detection.particles(0, 'dm_particles')
        print instance.stopping_conditions.out_of_box_detection.particles(1, 'dm_particles')
        self.assertEquals(instance.stopping_conditions.out_of_box_detection.particles(2, 'dm_particles'),[])

        self.assertTrue(instance.stopping_conditions.out_of_box_detection.is_set())
        self.assertAlmostEqual(instance.model_time.number ,1.,3)
        
        instance.stop()

    
    def test14(self):
        print "Testing Fi get_hydro_state_at_point"
        number_sph_particles = 100
        convert_nbody = nbody.nbody_to_si(1.0 | units.kpc, 1.0e10 | units.MSun)
        gas = new_evrard_gas_sphere(number_sph_particles, convert_nbody, seed = 1234)
        gas.h_smooth = 0.01 | nbody.length
        
        instance = Fi(convert_nbody)
        instance.parameters.n_smooth     =   64 | units.none
        instance.parameters.n_smooth_tol = 0.01 | units.none
        instance.gas_particles.add_particles(gas)
        instance.synchronize_model()
        
        coords = [0.0 | units.kpc]*3
        speeds = [0.0 | units.m / units.s]*3
        hydro_state = instance.get_hydro_state_at_point(*(coords + speeds))
        expected = [ 3.5540e-19 | units.kg * units.m**-3, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                     7.6447e-10 | units.kg * units.m**-1 * units.s**-2]
        for value, expect in zip(hydro_state, expected):
            self.assertAlmostRelativeEqual(value, expect, places=3)
        
        coords = [0.1 | units.kpc]*3
        speeds = [0.0 | units.m / units.s]*3
        hydro_state = instance.get_hydro_state_at_point(*(coords + speeds))
        expected = [ 4.1789e-19 | units.kg * units.m**-3, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                     8.9889e-10 | units.kg * units.m**-1 * units.s**-2]
        for value, expect in zip(hydro_state, expected):
            self.assertAlmostRelativeEqual(value, expect, places=3)
        instance.stop()
    
    def test15(self):
        print "Testing Fi get_hydro_state_at_point II: uniform sphere"
        number_sph_particles = 1000
        convert_nbody = nbody.nbody_to_si(1.0 | units.kpc, 1.0e10 | units.MSun)
        gas = new_uniform_spherical_particle_distribution(number_sph_particles, 1.0 | units.kpc, 1.0e10 | units.MSun)
        gas.velocity = [10.0, 20.0, 30.0] | units.km / units.s
        gas.h_smooth = 0.01 | nbody.length
        gas.u = 0.05 | nbody.specific_energy
        density = (1.0e10 | units.MSun) / (4.0/3.0 * numpy.pi * (1.0 | units.kpc)**3)
        
        instance = Fi(convert_nbody)
        instance.parameters.n_smooth     =  64 | units.none
        instance.parameters.n_smooth_tol = 0.01 | units.none
        instance.gas_particles.add_particles(gas)
        instance.synchronize_model()
        
        coords = [0.0 | units.kpc]*3
        speeds = [0.0 | units.m / units.s]*3
        rho, rhovx, rhovy, rhovz, rhoe = instance.get_hydro_state_at_point(*(coords + speeds))
        self.assertAlmostRelativeEqual(rho,   density, places=3)
        self.assertAlmostRelativeEqual(rhovx, density*instance.gas_particles[0].vx, places=3)
        self.assertAlmostRelativeEqual(rhovy, density*instance.gas_particles[0].vy, places=3)
        self.assertAlmostRelativeEqual(rhovz, density*instance.gas_particles[0].vz, places=3)
        self.assertAlmostRelativeEqual(rhoe,  density * (instance.gas_particles[0].u + 
            instance.gas_particles[0].velocity.length_squared()),  places=3)
        
        coords = [0.1 | units.kpc]*3
        rho, rhovx, rhovy, rhovz, rhoe = instance.get_hydro_state_at_point(*(coords + speeds))
        self.assertAlmostRelativeEqual(rho,   density, places=3)
        self.assertAlmostRelativeEqual(rhovx, density*instance.gas_particles[0].vx, places=3)
        self.assertAlmostRelativeEqual(rhovy, density*instance.gas_particles[0].vy, places=3)
        self.assertAlmostRelativeEqual(rhovz, density*instance.gas_particles[0].vz, places=3)
        self.assertAlmostRelativeEqual(rhoe,  density * (instance.gas_particles[0].u + 
            instance.gas_particles[0].velocity.length_squared()),  places=3)
        instance.stop()
    
