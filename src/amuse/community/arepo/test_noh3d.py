import random
from amuse.community.arepo import Arepo
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
# import unit seconds as s

N_PLOT_PARTICLES = 3000

# Check code runs without errors
x = Arepo(redirection="none")
x.initialize_code()

n_particles_total = x.get_number_of_particles()
print('AMUSE: number of particles: {}'.format(n_particles_total))
random.seed(123)
tracked_ids = random.sample(range(n_particles_total), k=N_PLOT_PARTICLES)


def dist(p0, p1):
    return np.linalg.norm(np.array(p1) - np.array(p0))


positions = {id: [x.get_position(id)] for id in tracked_ids}
print(x.get_position(30))

print("Evolving")
x.evolve_model(0.00001)

for id in tracked_ids:
    positions[id].append(x.get_position(id))

print("Evolving another step")
x.evolve_model(0.00002)

final_densities = {}

for id in tracked_ids:
    positions[id].append(x.get_position(id))
    final_densities[id] = x.get_density(id)

max_density = max(final_densities.values())
min_density = min(final_densities.values())

for id, _ in zip(positions, range(5)):
    print(id)
    print(positions[id])
    print(dist(positions[id][0], positions[id][-1]))
    print()


# Plot1
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
cmap = cm.plasma
for id, pos in positions.items():
    ax.plot(*np.array(pos).T, color=cmap((final_densities[id] - min_density)/(max_density - min_density)))
fig.suptitle("Densities: t = {}".format(x.get_time()))
fig.savefig("noh_positions_0.png")

# Second plot at t = 1
print("Evolving another step")
x.evolve_model(1.0 - 0.00005)

for id in tracked_ids:
    positions[id].append(x.get_position(id))

x.evolve_model(1.0)

final_densities = {}

for id in tracked_ids:
    positions[id].append(x.get_position(id))
    final_densities[id] = x.get_density(id)

max_density = max(final_densities.values())
min_density = min(final_densities.values())

for id, _ in zip(positions, range(5)):
    print(id)
    print(positions[id])
    print(dist(positions[id][0], positions[id][-1]))
    print()

# Plot2
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
cmap = cm.plasma
for id, pos in positions.items():
    ax.plot(*np.array(pos[-2:]).T, color=cmap((final_densities[id] - min_density)/(max_density - min_density)))
fig.suptitle("Densities: t = {}".format(x.get_time()))
fig.savefig("noh_positions_1.png")

# Third plot at t = 2
print("Evolving another step")
x.evolve_model(1.9 - 0.00002)

for id in tracked_ids:
    positions[id].append(x.get_position(id))

x.evolve_model(1.9)

final_densities = {}

for id in tracked_ids:
    positions[id].append(x.get_position(id))
    final_densities[id] = x.get_density(id)

max_density = max(final_densities.values())
min_density = min(final_densities.values())

for id, _ in zip(positions, range(5)):
    print(id)
    print(positions[id])
    print(dist(positions[id][0], positions[id][-1]))
    print()

# Plot2
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'))
cmap = cm.plasma
for id, pos in positions.items():
    ax.plot(*np.array(pos[-2:]).T, color=cmap((final_densities[id] - min_density)/(max_density - min_density)))
fig.suptitle("Densities: t = {}".format(x.get_time()))
fig.savefig("noh_positions_2.png")
