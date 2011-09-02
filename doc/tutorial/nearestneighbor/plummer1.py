from interface import NearestNeighbor
from amuse.units import units, nbody_system
from amuse.ext import plummer
from amuse.io import text

number_of_particles = 1000
mass_per_particle = 1 | units.MSun

convert_nbody = nbody_system.nbody_to_si(number_of_particles * mass_per_particle, 1.0 | units.parsec)
uc = plummer.MakePlummerModel(number_of_particles, convert_nbody)
particles = uc.result

nn = NearestNeighbor()
nn.particles.add_particles(particles)
print "number of particles:", len(nn.particles)

nn.find_nearest_neighbors()

local_particles = nn.particles.copy()

for p in local_particles:
    delta =  p.neighbor1.position - p.position
    p.distance_to_neighbor = delta.length()
    p.dx = delta.x
    p.dy = delta.y
    p.dz = delta.z
    


output = text.TableFormattedText("output.txt", set = local_particles)
output.attribute_names = ['x','y','z', 'dx', 'dy','dz']
output.attribute_types = [units.parsec] * 6
output.store()

