from interface import NearestNeighbor
from amuse.lab import *
from amuse.io import text

if __name__ == '__main__':
    number_of_particles = 1000
    particles = new_plummer_sphere(1000)

    code = NearestNeighbor()
    code.set_maximum_number_of_particles(5000)
    code.commit_parameters
    code.particles.add_particles(particles)

    code.run()

    local_particles = code.particles.copy()
    delta = local_particles.neighbor1.position - local_particles.position
   
    local_particles.dx = delta[...,0]
    local_particles.dy = delta[...,1]
    local_particles.dz = delta[...,2]

    output = text.TableFormattedText("output.txt", set = local_particles)
    output.attribute_names = ['x','y','z', 'dx', 'dy','dz']
    output.store()

