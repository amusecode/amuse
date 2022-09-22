import matplotlib.pyplot as plt
import numpy as np

from amuse.units import units
from amuse.community.mesa.interface import MESA
from amuse import datamodel

###BOOKLISTSTART###
def evolve_with_different_stellar_model():
    times = [10, 100, 1000] | units.Myr
    stars = datamodel.Particles(mass=[1, 2, 4]|units.MSun)
    stellars = [MESA(version='15140')]
    channels = []
    for star, stellar in zip(stars, stellars):
        stellar.particles.add_particle(star)
        channels.append(stellar.particles.new_channel_to(stars))

    for time, channel, stellar in zip(times, channels, stellars):
        stellar.evolve_model(time)
        channel.copy()

    for time, star in zip(times, stars):
        print("Time=", time, "M=", star.mass)
        stellar.stop()
###BOOKLISTSTOP###

evolve_with_different_stellar_model()
