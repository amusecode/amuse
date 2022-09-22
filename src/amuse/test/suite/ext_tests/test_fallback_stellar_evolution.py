import os
import os.path
import numpy

from amuse.units import units
from amuse.datamodel import Particles,Particle
from amuse.support.exceptions import AmuseException
from amuse.test.amusetest import TestCase, get_path_to_results

from amuse.community.mesa.interface import MESA
from amuse.community.evtwin.interface import EVtwin
from amuse.community.sse.interface import SSE

from amuse.couple.fallback_stellar_evolution import FallbackStellarEvolution

class TestFallbackStellarEvolution(TestCase):

    def test1(self):
        print("Testing FallbackStellarEvolution")
        instance = FallbackStellarEvolution()
        instance.stop()

    def xtest2(self):
        print("Testing FallbackStellarEvolution: evolve tests")


# results of original code  (not really check the numbers, tests have been relaxed because
# different timestepping of evolve_model and evolve_one_step)
        results=dict()
        results["10.0 MSun"]=dict(sse_age=3644655.52487 | units.yr,
                                      sse_mass=9.99289747724 | units.MSun,
                                      sse_rad=4.24764389443 | units.RSun,
                                      sse_L=5993.72678228 | units.LSun,
                                      evtwin_age=2999999.99732 | units.yr,
                                      evtwin_mass=9.99832250931 | units.MSun,
                                      evtwin_rad=4.1727904434 | units.RSun,
                                      evtwin_L=6175.5214329 | units.LSun)
        results["20.0 MSun"]=dict(sse_age= 2609167.84395| units.yr,
                                      sse_mass= 19.8895071502| units.MSun,
                                      sse_rad=7.07713598137| units.RSun,
                                      sse_L=51232.436257| units.LSun,
                                      evtwin_age= 2999999.99792| units.yr,
                                      evtwin_mass=19.877338008| units.MSun,
                                      evtwin_rad=6.77369109826| units.RSun,
                                      evtwin_L= 54123.7589452| units.LSun)
        results["40.0 MSun"]=dict(sse_age= 2627093.59096| units.yr,
                                      sse_mass=38.3702723816| units.MSun,
                                      sse_rad=13.6137662257 | units.RSun,
                                      sse_L=303498.877464 | units.LSun,
                                      evtwin_age= 2999999.99923| units.yr,
                                      evtwin_mass=37.0862227121 | units.MSun,
                                      evtwin_rad= 13.5982316053| units.RSun,
                                      evtwin_L=320099.608846| units.LSun)


        for m in [40.]| units.MSun:#,20.,10.] | units.MSun:
          instance = FallbackStellarEvolution()
          instance._main_se.parameters.max_age_stop_condition=3.| units.Myr

          star = Particle(1)
          star.mass = m
          star = instance.particles.add_particle(star)
          instance.commit_particles()

          while instance.ActiveModel[star].__class__.__name__=="EVtwin":
              instance.evolve_model()

          print("%s\t%s\t%s\t%s\t%s\t%s" % (star.age, star.mass, star.radius,
           star.luminosity, star.stellar_type, instance.ActiveModel[star].__class__.__name__))

          self.assertAlmostRelativeEqual(results[str(m)]["sse_age"],instance._FBTimeseries[star].particles[0].SSEAgeAtSwitch,7)
          self.assertAlmostRelativeEqual(results[str(m)]["sse_mass"],star.mass,7)
          self.assertAlmostRelativeEqual(results[str(m)]["sse_rad"],star.radius,7)
          self.assertAlmostRelativeEqual(results[str(m)]["sse_L"],star.luminosity,7)


          star=instance._main_se.particles[0]

          self.assertAlmostRelativeEqual(results[str(m)]["evtwin_age"],star.age,7)
          self.assertAlmostRelativeEqual(results[str(m)]["evtwin_mass"],star.mass,2)
          self.assertAlmostRelativeEqual(results[str(m)]["evtwin_rad"],star.radius,2)
          self.assertAlmostRelativeEqual(results[str(m)]["evtwin_L"],star.luminosity,2)

          print("%s\t%s\t%s\t%s\t%s" % (star.age, star.mass, star.radius,
           star.luminosity, star.stellar_type))

          instance.stop()

    def slowtest3(self):
        print("Testing FallbackStellarEvolution: evolve 3 stars at the same time")


# results of original code  (not really check the numbers, tests have been relaxed because
# different timestepping of evolve_model and evolve_one_step)
        results=dict()
        results["10.0 MSun"]=dict(sse_age=3644655.52487 | units.yr,
                                      sse_mass=9.99289747724 | units.MSun,
                                      sse_rad=4.35949272485 | units.RSun,
                                      sse_L=6180.22077675| units.LSun,
                                      evtwin_age=2999999.99732 | units.yr,
                                      evtwin_mass=9.99832250931 | units.MSun,
                                      evtwin_rad=4.1727904434 | units.RSun,
                                      evtwin_L=6175.5214329 | units.LSun)
        results["20.0 MSun"]=dict(sse_age= 2609167.84395| units.yr,
                                      sse_mass= 19.8895071502| units.MSun,
                                      sse_rad=7.30432858615| units.RSun,
                                      sse_L=51232.436257| units.LSun,
                                      evtwin_age= 2999999.99792| units.yr,
                                      evtwin_mass=19.877338008| units.MSun,
                                      evtwin_rad=6.77369109826| units.RSun,
                                      evtwin_L= 54123.7589452| units.LSun)
        results["40.0 MSun"]=dict(sse_age= 2869735.64001| units.yr, # for some reason this is not 2617093.59096
                                      sse_mass=37.9172423183| units.MSun,
                                      sse_rad=15.1388511451| units.RSun,
                                      sse_L=322674.980598| units.LSun,
                                      evtwin_age= 2999999.99923| units.yr,
                                      evtwin_mass=37.0862227121 | units.MSun,
                                      evtwin_rad= 13.5982316053| units.RSun,
                                      evtwin_L=320099.608846| units.LSun)


        instance = FallbackStellarEvolution()
        instance._main_se.parameters.max_age_stop_condition=3.| units.Myr

        stars = Particles(3,mass=[10.,20.,40.] | units.MSun)
        stars = instance.particles.add_particles(stars)
        stars.initial_mass=stars.mass
        instance.commit_particles()

        while instance.model_time<.05| units.Myr:
              instance.evolve_model()

        for star in stars:
          self.assertEqual(instance.ActiveModel[star].__class__.__name__,"EVtwin")

        while instance.model_time<=3.| units.Myr:
              instance.evolve_model()

        for star in stars:
          self.assertEqual(instance.ActiveModel[star].__class__.__name__,"SSE")

        for star in stars:
          m=star.initial_mass
          self.assertAlmostRelativeEqual(results[str(m)]["sse_age"],instance._FBTimeseries[star].particles[0].SSEAgeAtSwitch,3)

          evstar=star.as_particle_in_set(instance._main_se.particles)

          self.assertAlmostRelativeEqual(results[str(m)]["evtwin_age"],evstar.age,7)
          self.assertAlmostRelativeEqual(results[str(m)]["evtwin_mass"],evstar.mass,2)
          self.assertAlmostRelativeEqual(results[str(m)]["evtwin_rad"],evstar.radius,2)
          self.assertAlmostRelativeEqual(results[str(m)]["evtwin_L"],evstar.luminosity,2)


        instance.stop()

    def slowtest4(self):
        print("Testing FallbackStellarEvolution: evolve same 3 stars at the same time")

        instance = FallbackStellarEvolution()
        instance._main_se.parameters.max_age_stop_condition=3.| units.Myr

        stars = Particles(3,mass=[40.,40.,40.] | units.MSun)
        stars = instance.particles.add_particles(stars)
        instance.commit_particles()

        while instance.model_time<=3.| units.Myr:
              instance.evolve_model()

        for star in stars:
          self.assertEqual(instance.ActiveModel[star].__class__.__name__,"SSE")

        age1=instance._FBTimeseries[stars[0]].particles[0].SSEAgeAtSwitch
        age2=instance._FBTimeseries[stars[1]].particles[0].SSEAgeAtSwitch
        age3=instance._FBTimeseries[stars[2]].particles[0].SSEAgeAtSwitch

        self.assertEqual(age1,age2)
        self.assertEqual(age2,age3)

        instance.stop()

    def slowtest5(self):
        print("Testing FallbackStellarEvolution: evolve with end time")

        instance = FallbackStellarEvolution()
        instance._main_se.parameters.max_age_stop_condition=0.1| units.Myr

        stars = Particles(3,mass=[10.,20.,30.] | units.MSun)
        stars = instance.particles.add_particles(stars)
        instance.commit_particles()

        instance.evolve_model(.101 | units.Myr)

        for star in stars:
          self.assertEqual(instance.ActiveModel[star].__class__.__name__,"SSE")

        self.assertTrue(.101| units.Myr<instance.model_time)

        instance.stop()

    def test6(self):
        print("Testing FallbackStellarEvolution: enforce monotonic mass evolution")

        instance = FallbackStellarEvolution(enforce_monotonic_mass_evolution=True)
        instance._main_se.parameters.max_age_stop_condition=3.| units.Myr

        stars = Particles(1,mass=[40.] | units.MSun)
        stars = instance.particles.add_particles(stars)
        instance.commit_particles()

        while instance.ActiveModel[stars[0]].__class__.__name__=="EVtwin":
          emass=stars[0].mass
          instance.evolve_model()

        self.assertTrue(stars[0].mass<=emass)

        instance.stop()

    def slowtest7(self):
        print("Testing FallbackStellarEvolution: test mesa")

        instance = FallbackStellarEvolution(MESA)
        instance._main_se.parameters.max_age_stop_condition=3.| units.Myr

        stars = Particles(3,mass=[40.,40.,40.] | units.MSun)
        stars = instance.particles.add_particles(stars)
        instance.commit_particles()

        instance.evolve_model(3.01 | units.Myr)

        for star in stars:
          self.assertEqual(instance.ActiveModel[star].__class__.__name__,"SSE")

        age1=instance._FBTimeseries[stars[0]].particles[0].SSEAgeAtSwitch
        age2=instance._FBTimeseries[stars[1]].particles[0].SSEAgeAtSwitch
        age3=instance._FBTimeseries[stars[2]].particles[0].SSEAgeAtSwitch

        self.assertEqual(age1,age2)
        self.assertEqual(age2,age3)

        instance.stop()


