import os.path

from amuse.test.amusetest import TestWithMPI, get_path_to_results
from amuse.support.exceptions import AmuseException
from amuse.ext.galactics_model import new_galactics_model
from amuse.units import nbody_system, generic_unit_converter, constants, units
from amuse.io import write_set_to_file

class NewGalactICsModelTests(TestWithMPI):
    
    def test1(self):
        halo_number_of_particles = 100
        particles = new_galactics_model(halo_number_of_particles, generate_bulge_flag = False, 
            generate_disk_flag = False, order_of_multipole_expansion = 0, halo_random_seed = -1, 
            do_scale = True)
        
        self.assertEquals(len(particles), halo_number_of_particles)
        self.assertAlmostEquals(particles.total_mass(), 1.0 | nbody_system.mass)
        self.assertAlmostEquals(particles.kinetic_energy(), 
            0.25 | nbody_system.energy)
    
    def test2(self):
        self.assertTrue("GalactICs documentation:" in new_galactics_model.__doc__)
        self.assertTrue("GalactICs allows to generate self-consistent disc-bulge-halo galaxy models." in new_galactics_model.__doc__)
        
        self.assertTrue("Parameters:" in new_galactics_model.__doc__)
        self.assertTrue("generate_disk_flag" in new_galactics_model.__doc__)
        self.assertTrue("halo_streaming_fraction" in new_galactics_model.__doc__)
        
        print new_galactics_model.__doc__
    
    def slowtest3(self):
        print "Generate a model for M31, using defaults (100k disk, 50k bulge, 200k halo) - Nbody units"
        halo_number_of_particles = 200000
        particles = new_galactics_model(halo_number_of_particles, do_scale = True)
        
        self.assertEquals(len(particles), 350000)
        self.assertAlmostEquals(particles.total_mass(), 1.0 | nbody_system.mass)
        self.assertAlmostEquals(particles.kinetic_energy(), 0.25 | nbody_system.energy)
        
        write_set_to_file(particles, os.path.join(get_path_to_results(), 'M31_galactICs.amuse'), 'amuse')
    
    def slowtest4(self):
        print "Generate a model for a disk galaxy (10k disk, no bulge, 20k halo) - SI units"
        halo_number_of_particles = 20000
        converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(constants.G, 1.0e12 | units.MSun, 50.0 | units.kpc)
        particles = new_galactics_model(halo_number_of_particles, disk_number_of_particles=10000, 
            generate_bulge_flag=False, do_scale=True, unit_system_converter=converter)
        
        self.assertEquals(len(particles), 30000)
        self.assertAlmostRelativeEquals(particles.total_mass(), 1.0e12 | units.MSun, 10)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), converter.to_si(0.25 | nbody_system.energy), 10)
        
        disk = particles[:10000]
        self.assertAlmostRelativeEquals(disk.total_mass(), 2.156e10 | units.MSun, 3)
        self.assertAlmostRelativeEquals(disk.position.lengths_squared().amax().sqrt(), 15.5675526285 | units.kpc, 5)
        self.assertAlmostRelativeEquals(disk.position.std(axis=0), [3.64746422018, 3.62743124932, 0.171638441317] | units.kpc, 5)
        
        write_set_to_file(particles, os.path.join(get_path_to_results(), 'disk_galactICs.amuse'), 'amuse')
    

