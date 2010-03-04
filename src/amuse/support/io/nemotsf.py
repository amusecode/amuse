from amuse.support.data import core
from amuse.support.units import units
from amuse.support.units import nbody_system
import re
import numpy as np

TEMPLATE = \
"""set SnapShot
  set Parameters
    int Nobj %(n_obj)i
    double Time %(time)d
  tes
  set Particles
    int CoordSystem %(crdsystem)i
    double Mass[%(no_particles)i] %(mass_paragraph)s
    double PhaseSpace[%(no_particles)i][%(no_phases)i][%(no_dimensions)i] %(phase_paragraph)s
  tes
tes
"""

class Particles2Tsf(object):
    
    def __init__(self):
        self.tsf_string = ""
        self.no_particles = 0
        self.masses = ""
        self.phases = ""

    def init_string(self):
       
        fields = {'n_obj':self.no_particles,'time':0.0,'crdsystem':66306,
                  'no_particles':self.no_particles,
                  'mass_paragraph':self.masses,
                  'phase_paragraph':self.phases,'no_particles':self.no_particles, 'no_phases':2,'no_dimensions':3}
        self.tsf_string = TEMPLATE % fields
        print self.tsf_string

    def convert_to_tsf(self, Particles):
        self.no_particles = len(Particles)
        masses = Particles.mass.value_in(units.MSun)
        velocities = Particles.velocity.value_in(units.AUd)
        positions = Particles.position.value_in(units.AU)
        #self.masses = ' '.join([str(i) for i in masses])
        print self.masses
        for i, mass in enumerate(masses):
            self.masses += ' '+str(mass)
            if (i % 5) == 0:
                self.masses += '\n       '
        for i, phase in enumerate(zip(positions,velocities)):
            self.phases += ' '.join([str(j) for j in phase[0]]) + ' ' + ' '.join([str(k) for k in phase[1]]) 
            if (i % 1) == 0:
                self.phases += '\n       '

class Tsf2Particles(object):

    def __init__(self):
        self.number_of_particles = 0
        self.masses = []
        self.phases = []

    def return_numbers_in_brackets(self, line):
        numbers = []
        indices = re.findall('\[[0-9]*\]',line)
        for i in indices:
            numbers.append((int)(i.strip('[').strip(']')))
        return numbers    

    def return_numbers_from_paragraph(self, lines, start, end):

        numlist =[]
        
        for i in range(start, end):
            linechunks = lines[i].split(' ')
            for j in linechunks:
                try:
                    numlist.append(float(j))
                except:
                    pass

        return numlist

    def read_to_ram(self, inputfile):
        f = open(inputfile,'r')

        lines = f.readlines()
        start_masses = [i for i, oneline in enumerate(lines) if 'double Mass' in oneline][0]
        start_phasespace = [i for i, oneline in enumerate(lines) if 'double PhaseSpace' in oneline][0]
        massline_numbers = self.return_numbers_in_brackets(lines[start_masses])
        phaseline_numbers = self.return_numbers_in_brackets(lines[start_phasespace])
        
        no_particles = phaseline_numbers[0]
        no_phases = phaseline_numbers[1]
        no_dimensions = phaseline_numbers[2]

        self.number_of_particles = no_particles
        self.masses = self.return_numbers_from_paragraph(lines, start_masses, start_phasespace)
        self.phases = np.reshape(self.return_numbers_from_paragraph(lines, start_phasespace, len(lines)),
                                 [no_particles,no_phases,no_dimensions])

    def convert_to_particles(self, inputfile, converter = None):
        self.read_to_ram(inputfile)
        self.Particles = core.Particles(self.number_of_particles)
        
        self.Particles.mass = self.masses|nbody_system.mass
        self.Particles.position = [i[0] for i in self.phases]|nbody_system.length
        self.Particles.velocity = [i[1] for i in self.phases]|nbody_system.speed
        if not converter==None:
            self.Particles=core.ParticlesWithUnitsConverted(self.Particles,
                                                            converter.as_converter_from_si_to_nbody())
        
        
        
