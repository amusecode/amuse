import re
import numpy

from amuse.io import base
from amuse.units import units
from amuse.units import nbody_system
from amuse import datamodel
TEMPLATE = \
"""set SnapShot
  set Parameters
    int Nobj {0.number_of_particles:d}
    double Time {0.time:e}
  tes
  set Particles
    int CoordSystem {0.coordinate_system_id:d}
    double Mass[{0.number_of_particles:d}] {0.mass_paragraph:s}
    double PhaseSpace[{0.number_of_particles:d}][{0.number_of_phases:d}][{0.number_of_dimensions:d}] {0.phase_paragraph:s}
  tes
tes
"""

class Particles2Tsf(object):
    
    def __init__(self):
        self.number_of_particles = 0
        self.mass_paragraph = ""
        self.phase_paragraph = ""
        self.number_of_phases = 2
        self.number_of_dimensions = 3
        self.coordinate_system_id = 66306

    def convert_to_string(self, particles, converter = None):
        if not converter is None:
            particles=datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
            )
            
        self.time = (particles.get_timestamp() or (0.0|nbody_system.time)).value_in(nbody_system.time)
        self.particles = particles
        self.number_of_particles = len(particles)
        
        masses = particles.mass.value_in(nbody_system.mass)
        velocities = particles.velocity.value_in(nbody_system.speed)
        positions = particles.position.value_in(nbody_system.length)
       
        for i, mass in enumerate(masses):
            self.mass_paragraph += ' '+str(mass)
            if (i % 5) == 0:
                self.mass_paragraph += '\n       '
        for i, phase in enumerate(zip(positions,velocities)):
            self.phase_paragraph += ' '.join([str(j) for j in phase[0]]) + ' ' + ' '.join([str(k) for k in phase[1]]) 
            if (i % 1) == 0:
                self.phase_paragraph += '\n       '
                
        return TEMPLATE.format(self)

class Tsf2Particles(object):

    def __init__(self):
        self.number_of_particles = 0
        self.masses = []
        self.phases = []
        self.timestamp = None

    def return_numbers_in_brackets(self, line):
        numbers = []
        indices = re.findall(r'\[[0-9]*\]',line)
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

    def read_to_ram(self, string):
        lines = string.splitlines()
        
        timestamp_index = [i for i, oneline in enumerate(lines) if 'double Time' in oneline][0]
        self.timestamp = float(lines[timestamp_index].strip().split(' ')[2]) | nbody_system.time
        
        start_masses = [i for i, oneline in enumerate(lines) if 'double Mass' in oneline][0]
        start_phasespace = [i for i, oneline in enumerate(lines) if 'double PhaseSpace' in oneline][0]
        massline_numbers = self.return_numbers_in_brackets(lines[start_masses])
        phaseline_numbers = self.return_numbers_in_brackets(lines[start_phasespace])
        
        no_particles = phaseline_numbers[0]
        no_phases = phaseline_numbers[1]
        no_dimensions = phaseline_numbers[2]

        self.number_of_particles = no_particles
        self.masses = self.return_numbers_from_paragraph(lines, start_masses, start_phasespace)
        self.phases = numpy.reshape(self.return_numbers_from_paragraph(lines, start_phasespace, len(lines)),
                                 [no_particles,no_phases,no_dimensions])

    def convert_to_particles(self, string, converter = None):
        self.read_to_ram(string)
        result = datamodel.Particles(self.number_of_particles)
        
        result.mass = self.masses|nbody_system.mass
        result.position = [i[0] for i in self.phases]|nbody_system.length
        result.velocity = [i[1] for i in self.phases]|nbody_system.speed
        
        if not self.timestamp is None:
            result.savepoint(self.timestamp)
        
        if not converter is None:
            result=datamodel.ParticlesWithUnitsConverted(
                result,
                converter.as_converter_from_si_to_generic()
            )
        
        return result
                                                            
class NemoFileFormatProcessor(base.FullTextFileFormatProcessor):
    """
    Process a NEMO binary structured file
    """
    
    provided_formats = ['nemo', 'tsf']
    
    def __init__(self, filename = None, stream = None, set = None, format = None):
        base.FileFormatProcessor.__init__(self, filename, set, format)
        
    
    def load_string(self, string):
        x = Tsf2Particles()
        return x.convert_to_particles(string, self.nbody_to_si_converter)
        
    def store_string(self):
        x = Particles2Tsf()
        return x.convert_to_string(self.set, self.nbody_to_si_converter)

    @base.format_option
    def nbody_to_si_converter(self):
        "NEMO datafiles store nbody data, provide a converter to store si data (None means no converter)"
        return None
        
