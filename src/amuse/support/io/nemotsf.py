from amuse.support.data import core
from amuse.support.units import units

class Tsf2Particles(object):
    def __init__(self,filename=None):
        self.mass = []
        self.phase = []

    def convert_to_particles(self, inputfile):
        
        f = open(inputfile,'r')

        lines = f.readlines()
        start_masses = [i for i, oneline in enumerate(lines) if 'double Mass[' in oneline][0]
        start_phasespace = [i for i, oneline in enumerate(lines) if 'double PhaseSpace[' in oneline][0]
        lines[start_masses]=lines[start_masses].strip('double Mass')
        lines[start_phasespace]=lines[start_phasespace].strip('double PhaseSpace')
        print lines
if __name__ == "__main__":
    I = Tsf2Particles()
    I.convert_to_particles('p10.txt')
    
