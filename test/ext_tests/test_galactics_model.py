import os.path

from amuse.test.amusetest import TestWithMPI, get_path_to_results
from amuse.support.exceptions import AmuseException
from amuse.ext.galactics_model import new_galactics_model, new_galactics_gas_model
from amuse.units import nbody_system, generic_unit_converter, constants, units
from amuse.io import write_set_to_file

# testing the gas models is *very* slow

class NewGalactICsModelTests(TestWithMPI):
    
    def test1(self):
        halo_number_of_particles = 100
        particles = new_galactics_model(halo_number_of_particles, generate_bulge_flag = False, 
            generate_disk_flag = False, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = True)
        
        self.assertEqual(len(particles), halo_number_of_particles)
        self.assertAlmostEqual(particles.total_mass(), 1.0 | nbody_system.mass)
        self.assertAlmostEqual(particles.kinetic_energy(), 
            0.25 | nbody_system.energy)
    
    def test2(self):
        self.assertTrue("GalactICs documentation:" in new_galactics_model.__doc__)
        self.assertTrue("GalactICs allows to generate self-consistent disc-bulge-halo galaxy models." in new_galactics_model.__doc__)
        
        self.assertTrue("Parameters:" in new_galactics_model.__doc__)
        self.assertTrue("generate_disk_flag" in new_galactics_model.__doc__)
        self.assertTrue("halo_streaming_fraction" in new_galactics_model.__doc__)
        
        print(new_galactics_model.__doc__)
    
    def slowtest3(self):
        print("Generate a model for M31, using defaults (100k disk, 50k bulge, 200k halo) - Nbody units")
        halo_number_of_particles = 20000
        particles = new_galactics_model(halo_number_of_particles, do_scale = True, 
          bulge_number_of_particles=5000, disk_number_of_particles=10000)
        
        self.assertEqual(len(particles), 35000)
        self.assertAlmostEqual(particles.total_mass(), 1.0 | nbody_system.mass)
        self.assertAlmostEqual(particles.kinetic_energy(), 0.25 | nbody_system.energy)
        
        write_set_to_file(particles, os.path.join(get_path_to_results(), 'M31_galactICs.amuse'), 'amuse')
    
    def slowtest4(self):
        print("Generate a model for a disk galaxy (10k disk, no bulge, 20k halo) - SI units")
        halo_number_of_particles = 20000
        converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(constants.G, 1.0e12 | units.MSun, 50.0 | units.kpc)
        particles = new_galactics_model(halo_number_of_particles, disk_number_of_particles=10000, 
            generate_bulge_flag=False, do_scale=True, unit_system_converter=converter)
        
        self.assertEqual(len(particles), 30000)
        self.assertAlmostRelativeEquals(particles.total_mass(), 1.0e12 | units.MSun, 10)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), converter.to_si(0.25 | nbody_system.energy), 10)
        
        disk = particles[:10000]
        self.assertAlmostRelativeEquals(disk.total_mass(), 2.156e10 | units.MSun, 3)
        self.assertAlmostRelativeEquals(disk.position.lengths_squared().amax().sqrt().in_(units.kpc),
                                                 15.584 | units.kpc, 3)
                              
        self.assertAlmostRelativeEquals(disk.position.std(axis=0).in_(units.kpc), [3.5934, 3.6768, 0.17078] | units.kpc, 3)
        
        write_set_to_file(particles, os.path.join(get_path_to_results(), 'disk_galactICs.amuse'), 'amuse')
    
    def test5(self):
        halo_number_of_particles = 1000
        
        conv=nbody_system.nbody_to_si(1.e12 | units.MSun, 100. | units.kpc)
        
        particles1 = new_galactics_model(halo_number_of_particles, generate_bulge_flag = False, 
            generate_disk_flag = False, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = False)
        particles2 = new_galactics_model(halo_number_of_particles, conv, generate_bulge_flag = False, 
            generate_disk_flag = False, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = False)        
        particles3 = new_galactics_model(halo_number_of_particles, generate_bulge_flag = False, 
            generate_disk_flag = False, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = True)
        particles4 = new_galactics_model(halo_number_of_particles, conv,generate_bulge_flag = False, 
            generate_disk_flag = False, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = True)
            
                        
        self.assertEqual(len(particles1), halo_number_of_particles)
        self.assertEqual(len(particles2), halo_number_of_particles)
        self.assertEqual(len(particles3), halo_number_of_particles)
        self.assertEqual(len(particles4), halo_number_of_particles)
        
        self.assertAlmostEqual(conv.to_nbody(particles2.total_mass()), particles1.total_mass())
        self.assertAlmostEqual(conv.to_nbody(particles4.total_mass()), particles3.total_mass())
        self.assertAlmostEqual(particles3.total_mass(), 1. | nbody_system.mass)
        self.assertAlmostRelativeEquals(particles4.total_mass(), conv.to_si(1. | nbody_system.mass),12)

        r1=particles1.position.lengths().std()
        r2=particles2.position.lengths().std()
        r3=particles3.position.lengths().std()
        r4=particles4.position.lengths().std()
        self.assertAlmostEqual(conv.to_nbody(r2), r1)
        self.assertAlmostEqual(conv.to_nbody(r4), r3)
        self.assertTrue(r1/r3>100) # for the default parameters the scaling is quite drastic
        self.assertTrue(r2/r4>100)
        print(r1,r3)
        print(r2,r4)

    def test6(self):
        halo_number_of_particles = 1000
        
        conv=nbody_system.nbody_to_si(1.e12 | units.MSun, 100. | units.kpc)
        
        particles1 = new_galactics_gas_model(halo_number_of_particles, bulge_type_parameter=0, 
            disk_type_parameter=0, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = False, reuse_cached_model=False, verbose=True)
        particles2 = new_galactics_gas_model(halo_number_of_particles, conv, bulge_type_parameter=0, 
            disk_type_parameter=0, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = False, verbose=True)        
        particles3 = new_galactics_gas_model(halo_number_of_particles, bulge_type_parameter=0, 
            disk_type_parameter=0, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = True, verbose=True)
        particles4 = new_galactics_gas_model(halo_number_of_particles, conv, bulge_type_parameter=0, 
            disk_type_parameter=0, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = True, verbose=True)
            
                        
        self.assertEqual(len(particles1), halo_number_of_particles)
        self.assertEqual(len(particles2), halo_number_of_particles)
        self.assertEqual(len(particles3), halo_number_of_particles)
        self.assertEqual(len(particles4), halo_number_of_particles)
        
        self.assertAlmostEqual(conv.to_nbody(particles2.total_mass()), particles1.total_mass())
        self.assertAlmostEqual(conv.to_nbody(particles4.total_mass()), particles3.total_mass())
        self.assertAlmostEqual(particles3.total_mass(), 1. | nbody_system.mass)
        self.assertAlmostRelativeEquals(particles4.total_mass(), conv.to_si(1. | nbody_system.mass),12)

        r1=particles1.position.lengths().std()
        r2=particles2.position.lengths().std()
        r3=particles3.position.lengths().std()
        r4=particles4.position.lengths().std()
        self.assertAlmostEqual(conv.to_nbody(r2), r1)
        self.assertAlmostEqual(conv.to_nbody(r4), r3)
        self.assertTrue(r1/r3>30) # for the default parameters the scaling is quite drastic
        self.assertTrue(r2/r4>30)
        print(r1,r3)
        print(r2,r4)


        #~ self.assertAlmostEquals(particles.kinetic_energy(), 
            #~ 0.25 | nbody_system.energy)
