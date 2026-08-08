import random
import numpy as np
from matplotlib import pyplot as plt
from amuse.community.arepo import Arepo
# from amuse.units import nbody_system
# import unit seconds as s

random.seed(123)
N_PLOT_PARTICLES = 300

# Check code runs without errors
x = Arepo(redirection="none")
x.initialize_code()

n_particles_total = x.get_number_of_particles()
tracked_ids = random.sample(range(n_particles_total), k=N_PLOT_PARTICLES)

positions = {}

# Get start position of tracked particles
for id in tracked_ids:
    positions[id] = [x.get_position(id)]

print("Evolving")
x.evolve_model(0.00001)

# Update positions of tracked particles
for id in tracked_ids:
    positions[id].append(x.get_position(id))

print("Evolving another step")
x.evolve_model(0.00002)

# Update positions of tracked particles
for id in tracked_ids:
    positions[id].append(x.get_position(id))


def dist(p0, p1):
    return np.linalg.norm(np.array(p1) - np.array(p0))


# Print the paths of 5 particles
for id, _ in zip(positions, range(5)):
    print(id)
    print(positions[id])
    print(dist(positions[id][0], positions[id][-1]))
    print()

# Plot1
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
for id, pos in positions.items():
    ax.plot(*np.array(pos).T)
fig.savefig('positions_1.png')

# Plot2
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
for id, pos in positions.items():
    ax.plot(*np.array(pos).T)
ax.set_xlim([-1000, 1000])
ax.set_ylim([-1000, 1000])
ax.set_zlim([-1000, 1000])
fig.savefig('positions_2.png')

# Plot3
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
for id, pos in positions.items():
    ax.plot(*np.array(pos).T)
ax.set_xlim([-100, 100])
ax.set_ylim([-100, 100])
ax.set_zlim([-100, 100])
fig.savefig('positions_3.png')

# x.run_sim()
# END_TIME = 1.0 | s
# x.evolve_model(END_TIME)

# x.(evolve for a single timestep)

# amuse tests to check for diverging behaviour
