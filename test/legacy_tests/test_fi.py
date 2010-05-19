import os
import sys
import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.legacy.fi.interface import fi, Fi
from amuse.ext.evrard_test import MakeEvrardTest

from amuse.support.units import nbody_system as nbody
from amuse.support.units import units
from amuse.support.data import core
from amuse.legacy.support import channel

class TestMPIInterface(TestWithMPI):

    def test1(self):
        instance=fi()  
        instance.setup_module()
        instance.stop()
    
    def test2(self):
        instance=fi()  
        instance.setup_module()
        
        for x,l in [('gravity_only',0),('radiate',0),('starform',0),('cosmo',1),
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
        instance=fi()  
        instance.setup_module()
        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        retrieved_state = instance.get_state(1)
        self.assertEquals(11.0,  retrieved_state['mass'])
        self.assertEquals(2.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 1)
        instance.cleanup_module()
        instance.stop()
    
    def test4(self):
        instance=fi()  
        instance.setup_module()
        instance.new_particle(11.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        instance.set_time_step(2.0)
        retrieved_state = instance.get_time_step()
        self.assertEquals(2.0, retrieved_state['time_step'])
        instance.cleanup_module()
        instance.stop()
    
    def test5(self):
        instance=fi()
        instance.setup_module()
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
        instance.initialize_particles(0.0)
        self.assertEqual(instance.get_number_of_particles().number_of_particles, 3)
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
        instance.stop()

class TestEvrard(TestWithMPI):

    def test0(self):
        evrard=MakeEvrardTest(1000)
        mass,x,y,z,vx,vy,vz,u=evrard.new_model()
        smooth=numpy.zeros_like(mass)
        nb = fi()
        nb.setup_module()

        nb.set_stepout(99999)
        nb.set_steplog(99999)
        nb.set_gravity_only(0)
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
        

class TestFiInterface(TestWithMPI):

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
        
        instance.update_particles(stars)
        
        postion_after_full_rotation = earth.position.x
        
        self.assertAlmostRelativeEqual(postion_at_start, postion_after_full_rotation, 4)
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
       
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.x
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.x
        self.assertAlmostRelativeEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
         
        instance.update_particles(stars)
        
        postion_after_half_a_rotation = earth.position.y
        
        self.assertAlmostRelativeEqual(-postion_at_start, postion_after_half_a_rotation, 4)
        
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Test 2: testing Fi data directory"
        convert_nbody = nbody.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = Fi(convert_nbody, redirection='none')
        error = instance.initialize_code()
        self.assertEquals(0, error)
        self.assertTrue('data/fi/input/' in instance.get_fi_data_directory().fi_data_directory)
        
        self.assertEquals(False, instance.parameters.radiation_flag)
        instance.parameters.radiation_flag = True
        self.assertEquals(True, instance.parameters.radiation_flag)
        
        self.assertEquals(False, instance.parameters.star_formation_flag)
        instance.parameters.star_formation_flag = True
        self.assertEquals(True, instance.parameters.star_formation_flag)
        instance.commit_parameters()
        
        stars = core.Particles(2)
        stars.mass = [1.0, 3.0e-6] | units.MSun
        stars.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.AU
        stars.velocity = [[0.0,0.0,0.0], [0.0,29.8,0.0]] | units.km / units.s
        stars.radius = [1.0, 0.01] | units.RSun
        
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model(365.0 | units.day)
        
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):
        print "Test 3: testing Fi boolean parameters"
        instance = Fi(nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc))
        self.assertEquals(0, instance.initialize_code())
        
        for bool_par in ['radiation_flag','star_formation_flag','gravity_only_flag',
            'square_root_timestep_flag','freeform_timestep_flag','quadrupole_moments_flag',
            'direct_sum_flag','fixed_halo_flag','adaptive_smoothing_flag','smooth_input_flag',
            'integrate_entropy_flag','isothermal_flag']:
            self.assertEquals(False, eval("instance.parameters."+bool_par))
            exec("instance.parameters."+bool_par+" = True")
            self.assertEquals(True, eval("instance.parameters."+bool_par))
        
        for bool_par in ['acc_timestep_flag','self_gravity_flag','gadget_cell_opening_flag',
            'conservative_sph_flag','sph_dens_init_flag','eps_is_h_flag']:
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
            ('verbosity',0),('nsmooth',64)]:
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
            'nn_tol','gas_epsilon','gamma','alpha','beta','sph_artificial_viscosity_eps','courant',
            'min_gas_part_mass','sph_h_const','n_neighbour_tol','grain_heat_eff','zeta_cr_ion_rate',
            'heat_par1','heat_par2','cool_par','optical_depth','star_form_delay_fac','star_form_mass_crit',
            'star_form_eff','supernova_duration','supernova_eff','t_supernova_start','max_density']
        defaults=[0.0 | nbody.length * nbody.length, 1.0 | nbody.time, 300.0 | nbody.length, 
            1.0e9 | units.MSun, 1.0 | units.kpc, 1.0, 0.25, 0.5, 0.35, 
            0.0, -1.0, 0.5, 0.01, 0.1, 0.005 | nbody.length, 1.6666667, 0.5, 1.0, 0.01, 0.3, 
            0.25, 0.2 | nbody.length, 0.1, 0.05, 3.6 | 1.8e-17 * units.s**-1, 0.0, 0.0, 1.0, 
            0.0, 1.0, 1.0e5 | units.MSun, 0.25, 3.0e7 | units.Myr, 0.0, 3.e6 | units.Myr, 100.0]
        defaults = [val | units.none if isinstance(val,float) else val for val in defaults]
        defaults = [instance.convert_nbody.to_si(val) if nbody.is_nbody_unit(val.unit) else val for val in defaults]
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
    
