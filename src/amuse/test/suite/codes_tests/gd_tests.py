import numpy
from amuse.units import nbody_system
from amuse import datamodel
from amuse.ic.plummer import new_plummer_model


class _TestGravitationalDynamicsInterface:
    def gravity_code_interface(self):
        self.skip("abstract test")

    def reference_includes(self):
        self.skip("abstract test")

    def starting_particle_index(self):
        return 0

    def check_for_energy(self, *args, **kwargs):
        return self.assertEqual(*args, **kwargs)

    def almost_equal_precision(self):
        return 7

    def test_initialise(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()
        instance.stop()

    def test_literature_reference(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()
        reference_string = self.reference_includes()
        self.assertTrue(
            reference_string in instance.all_literature_references_string()
        )
        instance.stop()

    def test_add_and_retrieve_particle(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()

        index, error = instance.new_particle(
            11.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 2.0)
        self.assertEqual(error, 0)
        self.assertEqual(index, self.starting_particle_index())
        error = instance.commit_particles()
        self.assertEqual(error, 0)

        retrieved_state = instance.get_state(index)
        self.assertEqual(retrieved_state['__result'], 0)
        self.assertEqual(11.0, retrieved_state['mass'])
        self.assertEqual(2.0, retrieved_state['radius'])
        self.assertEqual(
            instance.get_number_of_particles()['number_of_particles'], 1
        )
        instance.stop()

    def test_add_and_retrieve_particles(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()

        instance.new_particle(
            [11.0, 12.0, 13.0, 14.0],
            [2.1, 3.1, 4.1, 5.1],
            [2.2, 3.2, 4.2, 5.2],
            [2.3, 3.3, 4.3, 5.3],
            [2.4, 3.4, 4.4, 5.4],
            [2.5, 3.5, 4.5, 5.5],
            [2.6, 3.6, 4.6, 5.6],
            [2.0, 3.0, 4.0, 5.0],
        ),
        error = instance.commit_particles()
        self.assertEqual(error, 0)
        retrieved_state = instance.get_state(self.starting_particle_index())
        self.assertEqual(11.0, retrieved_state['mass'])
        retrieved_state = instance.get_state(
            [
                self.starting_particle_index()+1,
                self.starting_particle_index()+2,
                self.starting_particle_index()+3,
            ]
        )
        self.assertEqual(12.0, retrieved_state['mass'][0])
        self.assertEqual(
            instance.get_number_of_particles()['number_of_particles'], 4
        )
        instance.stop()

    def test_add_and_retrieve_many_particles(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()
        n = 4000
        values = [1.0 * i for i in range(1, n)]
        instance.new_particle(
            values,
            values,
            values,
            values,
            values,
            values,
            values,
            values,
        )
        error = instance.commit_particles()
        self.assertEqual(error, 0)
        retrieved_state = instance.get_state(self.starting_particle_index())
        self.assertEqual(1.0, retrieved_state['mass'])
        retrieved_state = instance.get_state(
            self.starting_particle_index()+3998
        )
        instance.cleanup_code()
        instance.stop()
        self.assertEqual(3999.0,  retrieved_state['mass'])

    def test_epsilon_squared(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()
        for x in [0.101, 4.0]:
            error = instance.set_eps2(x)
            self.assertEqual(error, 0)
            value, error = instance.get_eps2()
            self.assertEqual(error, 0)
            self.assertEqual(x, value)
        instance.stop()

    def test_reversed_time_allowed(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()
        self.assertEqual(
            [0, 0], list(instance.get_is_time_reversed_allowed().values())
        )
        self.assertEqual(0, instance.set_is_time_reversed_allowed(1))
        self.assertEqual(
            [1, 0], list(instance.get_is_time_reversed_allowed().values())
        )
        instance.stop()

    def test_calculate_energies(self):
        interface = self.gravity_code_interface()
        instance = self.new_instance_of_an_optional_code(interface)
        instance.initialize_code()
        instance.set_eps2(0.0**2)
        instance.commit_parameters()

        instance.new_particle(
            [1.0, 1.0, 1.0],
            [1.0, 0.0, -1.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
        )
        instance.commit_particles()
        Ep = instance.get_potential_energy()['potential_energy']
        Ek = instance.get_kinetic_energy()['kinetic_energy']
        self.check_for_energy(Ek, 0.5)
        self.check_for_energy(Ep, -2.5)
        instance.delete_particle(self.starting_particle_index()+1)
        instance.recommit_particles()
        n = instance.get_number_of_particles()['number_of_particles']
        Ep = instance.get_potential_energy()['potential_energy']
        Ek = instance.get_kinetic_energy()['kinetic_energy']

        instance.cleanup_code()
        instance.stop()

        self.assertEqual(n, 2)
        self.check_for_energy(Ek, 0.)
        self.check_for_energy(Ep, -0.5)


class _TestGravityCodes:
    length_unit = nbody_system.length
    speed_unit = nbody_system.speed
    mass_unit = nbody_system.mass
    time_unit = nbody_system.time

    @property
    def nbody_converter(self):
        return None

    def gravity_code_factory(self):
        self.skip("abstract test")

    def test1(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        particles = new_plummer_model(100, convert_nbody=self.nbody_converter)
        particles.radius = 0 | self.length_unit
        particles.move_to_center()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        self.assertEqual(len(instance.particles), 100)
        outer = particles.select(lambda r: r.length() > (
            0.5 | self.length_unit), ["position"])
        print(len(outer))
        self.assertTrue(len(outer) > 0)
        self.assertTrue(len(outer) < 100)
        instance.synchronize_model()
        instance.particles.remove_particles(outer)
        instance.recommit_particles()
        self.assertEqual(len(instance.particles), 100-len(outer))
        number_of_particles_in_module = instance.legacy_interface.get_number_of_particles()[
            'number_of_particles']
        self.assertEqual(number_of_particles_in_module, 100-len(outer))
        instance.stop()

    def getset_attribute(self, attributename, attributevalue1, attributevalue2):

        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        try:
            particles = new_plummer_model(100)
            particles.move_to_center()
            particles.radius = 0.0 | nbody_system.length
            instance.particles.add_particles(particles)
            instance.commit_particles()
            setattr(instance.particles, attributename, attributevalue1)
            actual = getattr(instance.particles, attributename)
            self.assertAlmostRelativeEqual(actual, attributevalue1)
            setattr(instance.particles, attributename, attributevalue2)
            actual = getattr(instance.particles, attributename)
            print(actual.as_quantity_in(attributevalue2.unit))
            self.assertAlmostRelativeEqual(actual, attributevalue2)
        finally:
            instance.stop()

    def test2(self):
        self.getset_attribute(
            "radius", 1.0 | self.length_unit, 2.0 | self.length_unit)

    def test3(self):
        self.getset_attribute("position", [1.0, 2.0, 3.0] | self.length_unit, [
                              5.0, 6.0, 7.0] | self.length_unit)

    def test4(self):
        self.getset_attribute("velocity", [1.0, 2.0, 3.0] | self.speed_unit, [
                              5.0, 6.0, 7.0] | self.speed_unit)

    def test5(self):
        self.getset_attribute(
            "mass", 1.0 | self.mass_unit, 2.0 | self.mass_unit)

    def test6(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        try:
            particles = new_plummer_model(
                100, convert_nbody=self.nbody_converter)
            more_particles = new_plummer_model(
                50, convert_nbody=self.nbody_converter)
            particles.radius = 0 | self.length_unit
            more_particles.radius = 1.0 | self.length_unit
            particles.move_to_center()
            more_particles.move_to_center()

            instance.particles.add_particles(particles)
            instance.commit_particles()
            self.assertEqual(len(instance.particles), 100)
            instance.synchronize_model()
            instance.particles.add_particles(more_particles)
            instance.recommit_particles()
            number_of_particles_in_module = instance.legacy_interface.get_number_of_particles()[
                'number_of_particles']
            self.assertEqual(len(instance.particles), 150)
            self.assertEqual(number_of_particles_in_module, 150)
            instance.synchronize_model()
            instance.particles.remove_particles(particles)
            self.assertEqual(len(instance.particles), 50)
            instance.recommit_particles()
            instance.synchronize_model()
            number_of_particles_in_module = instance.legacy_interface.get_number_of_particles()[
                'number_of_particles']
            self.assertEqual(len(instance.particles), 50)
            self.assertEqual(number_of_particles_in_module, 50)
        finally:
            instance.stop()

    def test7(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        try:
            particles = new_plummer_model(
                100, convert_nbody=self.nbody_converter)
            new_particles = new_plummer_model(
                50, convert_nbody=self.nbody_converter)
            instance.particles.add_particles(particles)
            instance.commit_particles()

            self.assertEqual(len(instance.particles), 100)
            instance.synchronize_model()
            instance.particles.remove_particle(particles[40])
            instance.particles.remove_particle(particles[70])
            instance.particles.add_particle(new_particles[0])
            print(new_particles[0].key)
            instance.recommit_particles()
            # test the get_mass, get_position and get_velocity functions
            # if they are implemented for the code, otherwise will call
            # get_state multiple times
            # todo, fi fails, need to check with inti
            #self.assertAlmostRelativeEqual(instance.particles[-1].mass, new_particles[0].mass)
            #self.assertAlmostRelativeEqual(instance.particles[-1].velocity, new_particles[0].velocity)
            #self.assertAlmostRelativeEqual(instance.particles[-1].position, new_particles[0].position)
            instance.particles.synchronize_to(particles)
            self.assertEqual(len(particles), 99)
            self.assertEqual(particles[-1], new_particles[0])
        finally:
            instance.stop()

    def test8(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        random = numpy.random.mtrand.RandomState(3456)
        particles = new_plummer_model(
            10, convert_nbody=self.nbody_converter, random=random)
        particles.radius = 0.2 | self.length_unit
        particles.move_to_center()
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 10)
        collision_detection = instance.stopping_conditions.collision_detection
        collision_detection.enable()
        instance.evolve_model(1 | self.time_unit)
        self.assertTrue(collision_detection.is_set())
        instance.stop()

    def test9(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        instance.parameters.stopping_conditions_out_of_box_size = 2 | self.length_unit
        particles = datamodel.Particles(3)
        particles.x = [0, 1, 2] | self.length_unit
        particles.y = 0 | self.length_unit
        particles.z = 0 | self.length_unit
        particles.vx = [0, 0, 2] | self.speed_unit
        particles.vy = [1, -1, 0] | self.speed_unit
        particles.vz = 0 | self.speed_unit
        particles.mass = 1 | self.mass_unit
        particles.radius = 0.0 | self.length_unit
        instance.particles.add_particles(particles)
        stopping_condition = instance.stopping_conditions.out_of_box_detection
        stopping_condition.enable()
        instance.evolve_model(1 | self.time_unit)
        print(instance.particles)
        print(instance.particles.center_of_mass())
        print((instance.particles.position -
               instance.particles.center_of_mass()).lengths())
        self.assertTrue(stopping_condition.is_set())
        instance.stop()

    def test10(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        particle = datamodel.Particle()
        particle.position = [0, 0, 0] | self.length_unit
        particle.velocity = [1, -2, 3.0] | self.speed_unit
        particle.mass = 1 | self.mass_unit
        particle.radius = 0.0 | self.length_unit

        instance.particles.add_particle(particle)
        instance.evolve_model(1 | self.time_unit)
        self.assertAlmostEqual(instance.model_time, 1 | self.time_unit)
        self.assertAlmostEqual(instance.kinetic_energy,
                               7.0 | self.mass_unit * self.speed_unit**2)
        self.assertAlmostEqual(instance.potential_energy,
                               0.0 | self.mass_unit * self.speed_unit**2)
        self.assertAlmostEqual(instance.particles[0].position, [
                               1.0, -2.0, 3.0] | self.length_unit)
        instance.stop()

    def new_gravity_code(self):
        self.gravity_code_factory()
