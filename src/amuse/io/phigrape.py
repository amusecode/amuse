import re

from amuse.units import units
from amuse.units import nbody_system
from amuse.support import data
"""
fileformat:
===========

diskstep                # diskstep is used to number out, should be 0 for input dat
N                       # number of particles
time                    # current time, usually starts with 0.0
index   mass    x y z   vx vy vz        # N data lines where
                                        # index is a running particle index starting at 0!
                                        # mass is the mass of a particle
                                        # x y z are the x-, y-, and z- positions
                                        # vx vy vz are velocities components
"""
LINEFORMAT16 = "{particle_index:8d}   {mass:0.16E}   {x: 0.16E} {y: 0.16E} {z: 0.16E}   {vx: 0.16E} {vy: 0.16E} {vz: 0.16E}\n"
LINEFORMAT8 = "{particle_index:8d}   {mass:0.8E}   {x: 0.8E} {y: 0.8E} {z: 0.8E}   {vx: 0.8E} {vy: 0.8E} {vz: 0.8E}\n"

HEADER = "0\n{no_particles:d}\n{current_time}:f\n"

class Particles2Inp(object):

    def __init__(self):
        self.string = ""

    def __str__(self):
        return self.string

    def convert_to_inp(self, Particles, filename):

        self.no_particles = len(Particles)
        self.string +=  HEADER.format(no_particles=self.no_particles,
                                      current_time=0.0)

        masses = Particles.mass.value_in(nbody_system.mass)
        velocities = Particles.velocity.value_in(nbody_system.speed)
        positions = Particles.position.value_in(nbody_system.length)

        for i in range(self.no_particles):
            self.string += LINEFORMAT16.format(particle_index = i,
                                             mass = masses[i],
                                             x = positions[i][0],
                                             y = positions[i][1],
                                             z = positions[i][2],
                                             vx = velocities[i][0],
                                             vy = velocities[i][1],
                                             vz = velocities[i][2])

        f = open(filename, 'w')
        f.write(self.string)
        f.close()

class Inp2Particles(object):

    def __init__(self):
        pass

    def parse_int_from_group(self):
        pass

    def get_values_from_phase_line(self, line):
        """
           The regexp is grouping, therefore we add strings in the following code..
        """
        values =  re.findall('([-+])?([0-9]+(\.[0-9]*)?|\.[0-9]+)([eE][-+]?[0-9]+)?',line)
        args = ()
        n = int(values[0][1])
        m = float(values[1][1]+values[1][3])
        x = float(values[2][0]+values[2][1]+values[2][3])
        y = float(values[3][0]+values[3][1]+values[3][3])
        z = float(values[4][0]+values[4][1]+values[4][3])
        vx = float(values[5][0]+values[5][1]+values[5][3])
        vy = float(values[6][0]+values[6][1]+values[6][3])
        vz = float(values[7][0]+values[7][1]+values[7][3])
        return n, m, x, y, z, vx, vy, vz

    def read_to_ram(self, inputfile):
        f = open(inputfile,'r')
        lines = f.readlines()
        f.close()

        N_particles = int(lines[1])

        Particles = data.Particles(N_particles)

        for i, line in enumerate(lines):
            if i<3:
                pass#handle header
            else:
                n, m, x, y, z, vx, vy, vz = self.get_values_from_phase_line(line)
                Particles[n].mass = m | nbody_system.mass
                Particles[n].x = x |nbody_system.length
                Particles[n].y = y |nbody_system.length
                Particles[n].z = z |nbody_system.length
                Particles[n].vx = vx |nbody_system.speed
                Particles[n].vy = vy |nbody_system.speed
                Particles[n].vz = vz |nbody_system.speed

        return Particles

    def convert_to_particles(self, inputfile, converter = None):

        self.Particles = self.read_to_ram(inputfile)
        if not converter==None:
            self.Particles=data.ParticlesWithUnitsConverted(self.Particles,
                                                            converter.as_converter_from_si_to_nbody())
