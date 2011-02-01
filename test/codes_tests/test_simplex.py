import os.path
from amuse.test.amusetest import TestWithMPI
from amuse.community.simplex.interface import SimpleXInterface, SimpleX

# Change the default for some SimpleX(-Interface) keyword arguments:
default_options = dict(number_of_workers=2)
default_options = dict(number_of_workers=1, redirection="none")

class TestSimpleXInterface(TestWithMPI):

    def test1(self):
        instance = SimpleXInterface(**default_options)
        self.assertEqual(0, instance.set_output_directory(instance.output_directory))
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()
    
    def test2(self):
        instance = SimpleXInterface(**default_options)
        self.assertEqual(0, instance.set_output_directory(instance.output_directory))
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        
        input_file = os.path.join(instance.data_directory, 'vertices_10.txt')
        x, y, z, n_H, flux, X_ion = self.read_input_file(input_file)
        number_of_particles = len(x)
        inidices, errors = instance.new_particle(x, y, z, n_H, flux, X_ion)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertEqual(inidices, range(number_of_particles))
        self.assertEqual(0, instance.commit_particles())
        x, y, z, n_H, flux, X_ion, error = instance.get_state(0)
        for expected, received in zip([0.5, 0.5, 0.5, 1.0, 5.0, 0.0, 0], [x, y, z, n_H, flux, X_ion, error]):
            self.assertEqual(expected, received)
        self.assertEqual(0, instance.evolve(5.0, 1))
        print instance.get_state([0,1,2,3])
        x, y, z, n_H, flux, X_ion, error = instance.get_state(0)
        for expected, received in zip([0.5, 0.5, 0.5, 1.0, 5.0, 0], [x, y, z, n_H, flux, error]):
            self.assertAlmostEqual(expected, received, 5)
        self.assertTrue(0.0 < X_ion < 1.0)
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()
    
    def read_input_file(self, input_file):
        file = open(input_file, 'r')
        lines = file.readlines()
        lines.pop(0)
        x, y, z, nh, flux, xion = [], [], [], [], [], []
        for line in lines:
            l = line.strip().split()
            if len(l) >= 7:
                x.append(float(l[1]))
                y.append(float(l[2]))
                z.append(float(l[3]))
                nh.append(float(l[4]))
                flux.append(float(l[5]))
                xion.append(float(l[6]))
        return x, y, z, nh, flux, xion

