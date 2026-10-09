import numpy as np
from numpy.typing import NDArray
import pytest
from uclchem.model import get_species_names as _get_species_names

from amuse.datamodel import Particle, Particles
from amuse.support.testing.amusetest import TestWithMPI
from amuse.units import units as u
from amuse_uclchem.interface import UclchemInterface, Uclchem


class TestUclchemInterface(TestWithMPI):

    def test_getters_and_setters(self):
        instance = self.new_instance_of_an_optional_code(UclchemInterface)
        assert instance is not None

        instance.new_particle(1, 2, 3, 4)
        instance.new_particle(11, 12, 13, 14)

        result = instance.get_state(0)
        self.assertEquals(result['number_density'], 1)
        self.assertEquals(result['temperature'], 2)
        self.assertEquals(result['ionrate'], 3)
        self.assertEquals(result['radfield'], 4)

        result = instance.get_state(1)
        self.assertEquals(result['number_density'], 11)
        self.assertEquals(result['temperature'], 12)
        self.assertEquals(result['ionrate'], 13)
        self.assertEquals(result['radfield'], 14)

        instance.set_state(0, 5, 6, 7, 8)
        result = instance.get_state(0)
        self.assertEquals(result['number_density'], 5)
        self.assertEquals(result['temperature'], 6)
        self.assertEquals(result['ionrate'], 7)
        self.assertEquals(result['radfield'], 8)

        instance.set_number_density(0, 10)
        result = instance.get_number_density(0)
        self.assertEquals(result['number_density'], 10)

        instance.set_temperature(0, 20)
        result = instance.get_temperature(0)
        self.assertEquals(result['temperature'], 20)

        instance.set_ionrate(0, 30)
        result = instance.get_ionrate(0)
        self.assertEquals(result['ionrate'], 30)

        instance.set_radfield(0, 40)
        result = instance.get_radfield(0)
        self.assertEquals(result['radfield'], 40)

        instance.set_chemical_model('prestellarcore')
        result = instance.get_chemical_model()
        self.assertEquals(result['chem_model'], 'prestellarcore')

        H_index = _get_species_names().index('H')
        result = instance.get_species_name(H_index)
        self.assertEquals(result['name'], 'H')
        result = instance.get_species_index('H')
        self.assertEquals(result['species_index'], H_index)

        result = instance.get_time()
        self.assertEquals(result['time'], 0)

        instance.commit_particles()
        instance.set_abundance(0, 0, 5)

        result = instance.get_abundance(0, 0)
        self.assertEquals(result['abundance'], 5)

        instance.stop()

class TestUclchem(TestWithMPI):

    def generate_single_particle(self):
        p = Particle()
        p.number_density = 1e4 | u.cm**-3
        p.temperature = 10 | u.K
        p.ionrate = 1.3e-17 | u.s**-1
        p.radfield = 1 | u.habing

        return p

    def generate_two_particles(self):
        p = Particles(2)
        p[0].number_density = 1e4 | u.cm**-3
        p[0].temperature = 10 | u.K
        p[0].ionrate = 1.3e-17 | u.s**-1
        p[0].radfield = 1 | u.habing

        p[1].number_density = 1e5 | u.cm**-3
        p[1].temperature = 20 | u.K
        p[1].ionrate = 1.3e-17 | u.s**-1
        p[1].radfield = 1 | u.habing

        return p

    def _validate_particle_state(self, particle1, particle2):
        attributes = [
            'key', 'number_density', 'temperature', 'ionrate', 'radfield'
        ]
        for attr in attributes:
            self.assertEquals(getattr(particle1, attr), getattr(particle2, attr))

    def test_parameters(self):
        """Test parameters defined for Uclchem."""
        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        self.assertEquals(instance.parameters.chem_model, 'cloud')
        instance.parameters.chem_model = 'jshock'
        self.assertEquals(instance.parameters.chem_model, 'jshock')

        instance.stop

    def test_methods(self):
        """Test methods defined for Uclchem."""
        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        p = self.generate_two_particles()
        instance.particles.add_particles(p)

        new_state = [10 | u.cm**-3, 20 | u.K, 30 | u.s**-1, 40 | u.habing]
        instance.set_state(0, *new_state)
        self.assertEquals(instance.get_state(0), new_state)

        instance.set_number_density(1, 200 | u.cm**-3)
        self.assertEquals(instance.get_number_density(1), 200 | u.cm**-3)

        instance.set_temperature(1, 200 | u.K)
        self.assertEquals(instance.get_temperature(1), 200 | u.K)

        instance.set_ionrate(1, 200 | u.s**-1)
        self.assertEquals(instance.get_ionrate(1), 200 | u.s**-1)

        instance.set_radfield(1, 200 | u.habing)
        self.assertEquals(instance.get_radfield(1), 200 | u.habing)

        self.assertEquals(instance.get_number_of_particles(), 2)
        self.assertEquals(len(instance.particles), 2)

        instance.stop()

    def test_add_particle(self):
        """Test add single particle."""
        p = self.generate_single_particle()
        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        instance.commit_parameters()
        instance.particles.add_particle(p)

        self.assertEquals(instance.get_number_of_particles(), 1)
        self._validate_particle_state(instance.particles, p)

        instance.evolve_model(1e3 | u.yr)
        self.assertAlmostEquals(instance.model_time, 1e3| u.yr)

        instance.evolve_model(2e3 | u.yr)
        self.assertAlmostEquals(instance.model_time, 2e3| u.yr)

        instance.stop()

    def test_add_particles(self):
        """Test add 2 particles."""
        p = self.generate_two_particles()
        instance = self.new_instance_of_an_optional_code(Uclchem, redirection='none')
        assert instance is not None

        instance.commit_parameters()
        instance.particles.add_particles(p)
        instance.commit_particles()

        self.assertEquals(instance.get_number_of_particles(), 2)

        self._validate_particle_state(instance.particles, p)

        instance.stop()

    def test_add_and_remove_particle(self):
        """Add then delete a particle."""
        p = self.generate_single_particle()
        instance = self.new_instance_of_an_optional_code(Uclchem, redirection='none')
        assert instance is not None

        instance.commit_parameters()
        instance.particles.add_particle(p)

        self.assertEquals(instance.get_number_of_particles(), 1)

        instance.particles.remove_particle(instance.particles[0])
        self.assertEquals(instance.get_number_of_particles(), 0)
        self.assertEquals(instance.particles.is_empty(), True)

        instance.stop()

    def test_add_and_remove_particles(self):
        """Test add and delete multiple particles."""
        p1 = self.generate_two_particles()
        instance = self.new_instance_of_an_optional_code(Uclchem, redirection='none')
        assert instance is not None

        instance.commit_parameters()
        instance.particles.add_particles(p1)

        self._validate_particle_state(instance.particles[0], p1[0])
        self._validate_particle_state(instance.particles[1], p1[1])

        self.assertEquals(instance.get_number_of_particles(), 2)

        instance.particles.remove_particle(instance.particles[1])

        self.assertEquals(instance.get_number_of_particles(), 1)
        self._validate_particle_state(instance.particles[0], p1[0])

        p2 = self.generate_two_particles()
        instance.particles.add_particles(p2)

        self.assertEquals(instance.get_number_of_particles(), 3)

        instance.particles.remove_particles(instance.particles)
        self.assertEquals(instance.get_number_of_particles(), 0)
        self.assertEquals(instance.particles.is_empty(), True)

        instance.stop()

    def test_setting_abundance(self):
        """Test setting a single chemical abundance for a particle."""
        instance = self.new_instance_of_an_optional_code(Uclchem, redirection='none')
        assert instance is not None

        p = self.generate_two_particles()
        instance.particles.add_particles(p)
        instance.commit_particles()

        instance.set_abundance(0, 0, 5)
        assert instance.get_abundance(0, 0) == 5

        instance.set_abundance(1, 230, 0.45283)
        assert instance.get_abundance(1, 230) == 0.45283

        instance.stop()

    def test_setting_abundances(self):
        """Test setting all abundances of a particle from an array."""
        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        p = self.generate_two_particles()
        instance.particles.add_particles(p)
        instance.commit_particles()

        abundances = np.random.rand(335)

        instance.set_abundances(0, abundances)
        assert np.all(instance.particles.abundances[0,:] == abundances)
        assert np.all(instance.particles.abundances[1,:] == np.zeros(335))

        instance.stop()

    def test_evolve_abundances(self):
        """
        Test evolving a cloud model and getting the abundances.
        The expected abundances are from evolving a cloud model
        in UCLCHEM for 1e6 years.
        """
        # expected abundances
        abundances = self._cloud_abundances()
        p = self.generate_single_particle()

        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        instance.set_chemical_model = 'cloud'
        instance.commit_parameters()

        instance.particles.add_particle(p)
        instance.commit_particles()

        instance.evolve_model(1e6 | u.yr)

        self.assertAlmostRelativeEquals(
            instance.particles.abundances,
            abundances,
            places=3
        )

        instance.stop()

    def test_evolve_model(self):
        """
        Test evolve model. The expected abundances are from
        a standalone identical Uclchem run.
        """
        p = Particle()
        p.number_density = 10010.000467300415 | u.cm**-3
        p.temperature = 10.0 | u.K
        p.ionrate = 1.3e-17 | u.s**-1
        p.radfield = 1 | u.habing

        instance = self.new_instance_of_an_optional_code(Uclchem, redirection='none')
        assert instance is not None

        instance.parameters.chem_model = 'cloud'
        instance.commit_parameters()

        instance.particles.add_particle(p)
        instance.commit_particles()

        instance.evolve_model(1e6 | u.yr)
        abund = instance.get_abundances_by_name(0, ['H', 'H2', 'H2O', 'CO', 'CH3OH'])

        expected_abundances = [
            0.33620836113812969,
            0.33172125743196351,
            8.0494187005664807e-08,
            1.8889925479868117e-05,
            4.2327403848735972e-08
        ]

        for a, ea in zip(abund[0], expected_abundances):
            self.assertAlmostEqual(a, ea, places=3)

        instance.stop()

    def test_delete_particles_and_abundances(self):
        """Test deleting particles deletes abundances correctly."""
        instance = self.new_instance_of_an_optional_code(Uclchem, redirection='none')
        assert instance is not None

        N_particles = 5
        parts = Particles(N_particles)
        parts.number_density = np.random.rand(N_particles)*1.e5 | u.cm**-3
        parts.temperature = np.random.rand(N_particles)*50 | u.K
        parts.ionrate = np.random.rand(N_particles)*2.e-17 | u.s**-1
        parts.radfield = np.random.rand(N_particles)* 1 | u.habing

        instance.particles.add_particles(parts)
        instance.commit_particles()

        instance.evolve_model(10 | u.yr)

        N_species = len(instance.species)
        abundances = instance.particles.abundances
        self.assertEquals(abundances.shape, (N_particles, N_species))
        self.assertEquals(instance.get_number_of_particles(), N_particles)

        instance.particles.remove_particle(parts[0])
        instance.recommit_particles()
        self.assertEquals(instance.get_number_of_particles(), N_particles-1)
        self.assertEquals(instance.particles.abundances.shape, (N_particles-1, N_species))
        self.assertEquals(instance.particles.abundances.shape, abundances[1:, :].shape)
        self.assertEquals(instance.particles.abundances, abundances[1:, :])

        instance.particles.remove_particle(parts[1:3])
        instance.recommit_particles()
        self.assertEquals(instance.get_number_of_particles(), N_particles-3)
        self.assertEquals(instance.particles.abundances.shape, (N_particles-3, N_species))
        self.assertEquals(instance.particles.abundances, abundances[3:, :])

        instance.particles.remove_particles(parts[3:])
        instance.recommit_particles()
        self.assertEquals(instance.get_number_of_particles(), 0)

        instance.stop()

    def test_species_index(self):
        """Test querying species."""
        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        species = instance.species
        for name, index in species.items():
            self.assertEquals(
                instance.species[name], instance.get_species_index(name)
            )

            self.assertEquals(
                name, instance.get_species_name(index)
            )

        first, last = instance.get_firstlast_species_index()
        for i in (first, last):
            name = instance.get_species_name(i)
            idx = instance.get_species_index(name)
            assert idx == i, f"boundary mismatch at i={i}: name={name}, got idx={idx}"

        for i in range(first, last + 1):
            name = instance.get_species_name(i)
            idx = instance.get_species_index(name)
            assert idx == i, f"mismatch at i={i}: name={name!r}, idx={idx}"

        instance.cleanup_code()
        instance.stop()

    def test_species_index_vectorized(self):
        print("Test querying species with arrays.")
        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        species = instance.species
        chem_names = list(species.keys())
        indices = list(species.values())

        vec_indices = instance.get_species_index(chem_names)
        assert len(vec_indices) == len(chem_names), (
            f'length mismatch: {len(vec_indices)} indices for {len(chem_names)} names'
        )
        for name, vec_idx, scalar_idx in zip(chem_names, vec_indices, indices):
            assert vec_idx == scalar_idx, (
                f'vectorized get_species_index mismatch for {name!r}: '
                f'vector={vec_idx}, scalar={scalar_idx}'
            )

        reversed_names = chem_names[::-1]
        vec_reversed = instance.get_species_index(reversed_names)
        assert list(vec_reversed) == [species[n] for n in reversed_names], (
            'vectorized get_species_index does not preserve input order'
        )

        np_names = np.array(chem_names)
        vec_np = instance.get_species_index(np_names)
        assert list(vec_np) == indices, (
            'get_species_index fails or diverges on np.ndarray input vs list input'
        )
        single = instance.get_species_index([chem_names[0]])
        assert hasattr(single, '__len__') and len(single) == 1, (
            f'single-element array input returned non-array result: {single!r}'
        )
        assert single[0] == species[chem_names[0]]

        instance.cleanup_code()
        instance.stop()

    def test_get_abundance(self):
        """Test getting abundances"""

        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        cloud = self.generate_two_particles()
        cloud.index = range(len(cloud))
        instance.particles.add_particles(cloud)
        instance.commit_particles()

        reference_abundances = {
            'H': 0.76875095999999998,
            'H2': 0.0023031866999999998,
            'H+': 5.8058537999999997e-09,
            'HE': 0.23894584999999999,
        }

        for p in cloud.index:
            for species, value in reference_abundances.items():
                instance.set_abundance(p, instance.species[species], value)


        for species, expected in reference_abundances.items():
            result = instance.get_abundance(1, instance.species[species])
            self.assertAlmostRelativeEquals(result, expected, 7)

        instance.cleanup_code()
        instance.stop()

    def test_set_abundance(self):
        """Test set_abundances."""
        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        cloud = self.generate_two_particles()
        instance.particles.add_particles(cloud)
        instance.commit_particles()
        instance.set_abundance(1, instance.species['H2'], 0.0023031866999999998)

        result = instance.get_abundance(1, instance.species['H2'])
        self.assertEquals(result, 0.0023031866999999998)

        instance.cleanup_code()
        instance.stop()

    def test_set_abundances(self):
        """Test setting abundances from an array."""
        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        cloud = self.generate_two_particles()
        instance.particles.add_particles(cloud)
        instance.commit_particles()

        n_species = len(instance.species)
        abundances = np.zeros(n_species)
        abundances[instance.species['H']] = 0.76875095999999998
        abundances[instance.species['H2']] = 0.0023031866999999998
        abundances[instance.species['H+']] = 5.8058537999999997e-09
        abundances[instance.species['HE']] = 0.23894584999999999

        instance.set_abundances(0, abundances)
        instance.set_abundances(1, abundances*2)

        for species, expected in [
            ('H', 0.76875095999999998),
            ('H2', 0.0023031866999999998),
            ('H+', 5.8058537999999997e-09),
            ('HE', 0.23894584999999999),
        ]:
            result = instance.get_abundance(0, instance.species[species])
            self.assertEquals(result, expected)
            result = instance.get_abundance(1, instance.species[species])
            self.assertAlmostRelativeEquals(result, expected*2, 7)

        instance.cleanup_code()
        instance.stop()

    def test_get_abundances_by_name(self):
        """Test getting abundances by name."""
        instance = self.new_instance_of_an_optional_code(Uclchem)
        assert instance is not None

        N_particles = 4
        cloud = Particles(N_particles)
        for i in range(N_particles):
            cloud[i].number_density = 1e4 | u.cm**-3
            cloud[i].temperature = 10 | u.K
            cloud[i].ionrate = 1.3e-17 | u.s**-1
            cloud[i].radfield = 1 | u.habing

        cloud.index = range(len(cloud))
        instance.particles.add_particles(cloud)
        instance.commit_particles()

        start_abund = instance.particles.abundances.copy()
        reference_abundances = {
            'H': 0.76875095999999998,
            'H2': 0.0023031866999999998,
            'H+': 5.8058537999999997e-09,
            'HE': 0.23894584999999999,
        }

        # give each particle increasing abundances
        for i, p in enumerate(cloud.index):
            for species, value in reference_abundances.items():
                instance.set_abundance(p, instance.species[species], value+i)

        species_names = list(reference_abundances.keys())
        species_indices = instance.get_species_index(species_names)

        abundances = instance.particles.abundances.copy()
        self.assertEquals(
            abundances.shape, (N_particles, len(instance.species))
        )

        # check that only requested abundances were mutated
        non_target_mask = np.ones(len(instance.species), dtype=bool)
        non_target_mask[species_indices] = False
        self.assertAlmostRelativeEquals(
            abundances[:, non_target_mask],
            start_abund[:, non_target_mask],
            7
        )

        # ####################################
        # CASE 1: single index, single species
        # ####################################
        idx = 0
        result_single = instance.get_abundances_by_name(idx, 'H')

        self.assertIsInstance(result_single, np.ndarray)
        # 1 particle, 1 abundance
        self.assertEquals(result_single.shape, (1,1))

        # check that result matches reference abundance
        self.assertAlmostRelativeEquals(
            result_single[0], reference_abundances['H'], 7
        )

        # check that particle abundances match abundances
        self.assertAlmostRelativeEquals(
            instance.particles.abundances[idx, :],
            abundances[idx, :],
            7
        )
        abunds_1 = start_abund.copy()
        abunds_1[idx, species_indices] = [ref for ref in reference_abundances.values()]
        self.assertAlmostRelativeEquals(
            instance.particles.abundances[idx, :],
            abunds_1[idx, :],
            7
        )

        # ######################################
        # CASE 2: single index, multiple species
        # ######################################
        idx = 1
        # sequence-of-names input, order must be preserved
        result_seq = instance.get_abundances_by_name(idx, species_names)

        self.assertIsInstance(result_seq, np.ndarray)
        # 1 particle, multiple abundances
        self.assertEquals(result_seq.shape, (1, len(species_names)))

        # check that particle abundances match abundances
        for res, name in zip(result_seq[0], species_names):
            self.assertAlmostRelativeEquals(
                res, reference_abundances[name]+1, 7
            )

        self.assertAlmostRelativeEquals(
            instance.particles.abundances[idx, :],
            abundances[idx, :],
            7
        )

        # ########################################
        # CASE 3: multiple indices, single species
        # ########################################
        idxs = [0, 1, 2, 3]
        result_seq = instance.get_abundances_by_name(idxs, 'H')

        self.assertIsInstance(result_seq, np.ndarray)
        # 4 particles, 1 abundance
        self.assertEquals(result_seq.shape, (4,1))

        for i, res in enumerate(result_seq):
            self.assertAlmostRelativeEquals(
                res, reference_abundances['H']+i, 7
            )

        for idx in idxs:
            self.assertAlmostRelativeEquals(
                instance.particles.abundances[idx, :],
                abundances[idx, :],
                7
            )

        # ##########################################
        # CASE 4: multiple indices, multiple species
        # ##########################################
        idxs = np.asarray([0, 1, 2, 3])
        shuffled_names = [species_names[2], species_names[0], species_names[3], species_names[1]]
        result_seq = instance.get_abundances_by_name(idxs, shuffled_names)

        self.assertIsInstance(result_seq, np.ndarray)
        # 4 particles, 4 abundances
        self.assertEquals(result_seq.shape, (4,4))

        # check that particle abundances match abundances
        for i, res in enumerate(result_seq):
            expected_row = [reference_abundances[name] + i for name in shuffled_names]
            self.assertAlmostRelativeEquals(res, expected_row, 7)

        # ##################
        # CASE 5: Edge cases
        # ##################
        with pytest.raises(ValueError):
            instance.get_abundances_by_name([], 'H')
        with pytest.raises(ValueError):
            instance.get_abundances_by_name(0, [])

        instance.cleanup_code()
        instance.stop()

    def _cloud_abundances(self) -> NDArray:
        """Expected abundances for `test_evolve_abundances`."""
        return np.array([
        3.36208361e-01, 4.61369488e-09, 3.31721257e-01,
        4.38197269e-13, 8.15598279e-11, 9.99488396e-02,
        8.12645596e-10, 7.08904893e-15, 5.16888465e-07,
        1.08596360e-09, 1.42864545e-08, 2.87411021e-15,
        7.51917181e-14, 4.28870164e-15, 2.43828929e-06,
        7.64367332e-12, 2.32751032e-09, 2.11880518e-13,
        2.07864774e-09, 1.69586361e-15, 6.26769743e-08,
        2.14176267e-14, 1.95447035e-09, 3.62119077e-15,
        9.95313910e-06, 1.08827712e-15, 3.19971826e-16,
        1.41425281e-08, 1.45296428e-13, 5.97751984e-08,
        2.34661031e-12, 8.04941870e-08, 6.61902375e-13,
        1.38106082e-15, 1.04967399e-14, 1.69083916e-10,
        1.02816814e-16, 2.42403960e-07, 3.04048256e-10,
        1.59373896e-10, 1.10012344e-16, 2.34185638e-09,
        5.97915707e-15, 1.52678652e-08, 1.68430327e-15,
        1.00000000e-30, 1.42401375e-07, 1.92584172e-15,
        2.17230482e-10, 1.88899255e-05, 2.37826162e-15,
        4.12739559e-14, 4.71140877e-14, 5.51363025e-07,
        4.29484423e-16, 1.62878151e-07, 1.74484429e-09,
        6.55590999e-10, 4.83742510e-09, 5.07433258e-13,
        7.18248760e-14, 1.48392378e-11, 3.31173363e-13,
        9.40929842e-17, 1.00000000e-30, 2.85293144e-08,
        4.71694221e-14, 4.19197231e-08, 6.62870254e-14,
        6.67440901e-18, 6.81309212e-17, 9.05161233e-08,
        4.68145842e-14, 5.78310798e-10, 2.64124780e-16,
        3.13549379e-23, 1.83731822e-21, 3.69728355e-09,
        4.23274038e-08, 7.13897444e-17, 8.09475813e-08,
        3.93543661e-14, 8.67122161e-07, 6.45547050e-10,
        5.42757174e-24, 1.00000000e-30, 3.52056723e-16,
        6.15576448e-09, 9.59852066e-15, 2.89098995e-15,
        4.44510918e-18, 2.40100604e-28, 1.08090881e-09,
        2.31894307e-15, 4.65021399e-09, 3.36479803e-11,
        2.64818515e-17, 1.46659285e-17, 6.61889156e-10,
        4.16306160e-17, 1.09076083e-22, 9.60083451e-12,
        3.03973148e-16, 1.00000000e-30, 6.39989322e-19,
        2.59686442e-14, 1.17614589e-13, 1.14610689e-14,
        4.23244839e-20, 8.20593623e-12, 4.72827238e-29,
        7.65633899e-13, 1.93349106e-17, 3.50836869e-14,
        2.75426473e-17, 9.16135047e-08, 1.04951960e-08,
        8.45826983e-15, 3.71549951e-08, 4.16848094e-14,
        1.24134269e-12, 1.72535687e-16, 3.42776908e-15,
        2.06186354e-14, 3.47867744e-13, 1.23347409e-11,
        1.00000000e-30, 1.70761558e-30, 5.04999946e-10,
        1.97571198e-15, 2.07955038e-10, 7.27834444e-11,
        8.99151127e-10, 1.13268490e-14, 1.00000000e-30,
        5.23636252e-16, 7.38363975e-17, 3.90547699e-12,
        1.28497122e-14, 1.00000000e-30, 2.45410589e-17,
        9.93967711e-14, 3.43220021e-11, 2.63594772e-19,
        3.72559009e-11, 2.84731830e-20, 1.82023994e-17,
        1.67475377e-09, 7.08004715e-16, 2.25659956e-09,
        1.83422292e-16, 1.27789499e-17, 5.74988257e-16,
        4.78376156e-20, 1.86954507e-17, 1.00000000e-30,
        1.00000000e-30, 1.18344510e-09, 1.99915901e-15,
        2.49119252e-21, 2.82446160e-15, 6.24740166e-13,
        6.34475884e-17, 2.84464002e-22, 3.39152297e-18,
        6.71433206e-13, 8.49966441e-08, 2.76360972e-08,
        5.95745469e-07, 1.79935805e-08, 1.18679841e-08,
        2.19166767e-08, 7.41495683e-08, 2.10906273e-08,
        6.98784562e-08, 9.09456185e-08, 6.07303266e-08,
        2.85491979e-07, 3.95060440e-07, 2.77333339e-07,
        1.78960527e-06, 3.82798327e-11, 4.22886005e-08,
        8.86610783e-11, 2.62887086e-10, 5.70752305e-10,
        2.32801138e-08, 1.55266295e-11, 1.00000000e-30,
        7.20995402e-07, 8.08703927e-08, 4.10848978e-15,
        5.03905135e-11, 7.08319460e-09, 7.22438067e-07,
        1.00000000e-30, 6.09926753e-09, 5.40185899e-09,
        3.12683239e-07, 5.40110042e-09, 9.44190223e-08,
        7.77051809e-11, 7.92534953e-08, 4.70742600e-09,
        9.32059348e-09, 9.87977467e-07, 1.04796423e-08,
        1.76637860e-08, 2.90738188e-16, 8.15181613e-09,
        4.99451493e-08, 1.01805087e-10, 6.52538037e-10,
        1.00000000e-30, 1.00000000e-30, 3.16973729e-15,
        1.25303750e-14, 8.38216460e-13, 5.56708716e-30,
        1.13282940e-12, 1.06141268e-15, 5.07243371e-16,
        1.64761325e-09, 3.53772976e-10, 1.00756299e-08,
        8.50669669e-14, 3.01479565e-10, 5.58077771e-10,
        1.00000000e-30, 1.00000000e-30, 5.50830415e-10,
        6.71010205e-12, 8.12120301e-11, 6.25144379e-10,
        1.00000000e-30, 2.89429591e-13, 1.00000000e-30,
        1.00000000e-30, 1.55826125e-12, 3.41505470e-12,
        1.54530315e-18, 2.14583928e-10, 7.48361616e-11,
        1.59369142e-18, 1.00000000e-30, 2.19347288e-16,
        1.44832408e-11, 9.51156923e-12, 2.31687436e-11,
        8.15517852e-05, 1.68746788e-06, 5.05638365e-05,
        6.93419081e-05, 1.17684170e-06, 2.77788296e-06,
        3.49903872e-05, 3.09449791e-06, 1.22639785e-05,
        1.06388505e-06, 1.98452588e-06, 1.32063889e-04,
        6.15115450e-06, 4.69674911e-05, 4.36811141e-05,
        1.07550079e-08, 1.97100339e-06, 6.84028094e-09,
        2.66491950e-09, 4.27052528e-08, 1.74467671e-07,
        2.38486698e-11, 1.79092206e-30, 4.14562085e-05,
        8.08921123e-07, 1.61743262e-14, 7.65974568e-11,
        1.32272263e-06, 2.26880314e-05, 1.00000000e-30,
        1.42281179e-08, 2.24008062e-08, 4.29249914e-06,
        1.14509845e-07, 2.02489782e-06, 4.42298831e-10,
        2.49669031e-07, 5.69087142e-08, 3.74370934e-08,
        5.42766847e-06, 9.11148355e-09, 1.79417221e-06,
        6.40840009e-16, 7.58330458e-08, 6.29698010e-07,
        1.92383281e-08, 8.56157760e-09, 1.00000000e-30,
        8.59568297e-29, 3.67385099e-14, 1.21079285e-13,
        4.21315677e-12, 1.98330137e-28, 1.23674959e-11,
        1.44708575e-14, 3.66834429e-16, 1.91927063e-08,
        2.56955051e-08, 6.66885986e-08, 1.59981919e-13,
        1.14872472e-08, 1.07678807e-10, 3.16141397e-30,
        8.65813992e-29, 3.66219021e-10, 3.09393311e-11,
        2.79837692e-10, 5.26125295e-10, 8.53322603e-29,
        6.80960009e-13, 4.82137645e-30, 2.79294102e-30,
        1.91833676e-12, 1.30380503e-11, 2.22045128e-17,
        1.62078139e-09, 9.22369486e-11, 2.58972913e-17,
        3.06973321e-29, 6.28934028e-16, 2.52306481e-11,
        8.02107879e-12, 2.35467269e-11, 9.30394227e-06,
        5.72748519e-04, 7.02930465e-06,
    ])
