import struct
import numpy

from collections import namedtuple

from amuse.io import base
from amuse.units import units
from amuse.units import nbody_system
from amuse.support.core import late

from amuse.support import data
class GadgetFileFormatProcessor(base.FortranFileFormatProcessor):
    """
    Process (read) a Gadget binary data file (saving to Gadget files not supported)
    """
    
    provided_formats = ['gadget']
    GAS = 0
    HALO = 1
    DISK = 2
    BULGE = 3
    STARS = 4
    BNDRY = 5
    
    @late
    def header_format(self):
        collate = []
        collate.append(self.endianness)        
        for name, times, formatchar in self.header_struct_format:
            for t in range(times):
                collate.append(formatchar)
        return ''.join(collate)
        
    @base.format_option
    def header_struct_format(self):
        return (
            ('Npart', 6, 'I'),
            ('Massarr', 6, 'd'),
            ('Time', 1, 'd'),
            ('Redshift', 1, 'd'),
            ('FlagSfr', 1, 'i'),
            ('FlagFeedback', 1, 'i'),
            ('Nall', 6, 'i'),
            ('FlagCooling', 1, 'i'),
            ('NumFiles', 1, 'i'),
            ('BoxSize', 1, 'd'),
            ('Omega0', 1, 'd'),
            ('OmegaLambda', 1, 'd'),
            ('HubbleParam', 1, 'd'),
            ('FlagAge', 1, 'd'),
            ('FlagMetals', 1, 'd'),
            ('NallHW', 6, 'd'),
            ('flag_entr_ics', 1, 'i'),
        )
        
    @base.format_option
    def has_potential_energy(self):
        """Set to true if the file has a potential energy block"""
        return False
        
    @base.format_option
    def has_acceleration(self):
        """Set to true if the file has a acceleration block"""
        return False
    
    @base.format_option
    def has_rate_of_entropy_production(self):
        """Set to true if the file has a block with the
        rate of change of the entropic function of each gas particle
        """
        return False
        
    @base.format_option
    def has_timestep(self):
        """Set to true if the file has a block with the
        individual timestep for each particle.
        """
        return False
    
    @base.format_option
    def is_initial_conditions_format(self):
        """Set to true if the file contains
        initial conditions. An initial conditions
        gadget file contains less data.
        """
        return True
        
    @late
    def header_size(self):
        return struct.calcsize(self.header_format)
    
    @late
    def struct_class_name(self):
        return 'GadgetHeader'
        
    @late
    def struct_class(self):
        collate = []       
        for name, times, formatchar in self.header_struct_format:
            collate.append(name)
        attribute_names = ' '.join(collate)
        return namedtuple(self.struct_class_name, attribute_names)
    
    
    
    @late
    def total_number_of_particles(self):
        return sum(self.header_struct.Npart)
        
    @late
    def total_number_of_particles_with_variable_masses(self):
        result = 0
        for x, n in zip(self.header_struct.Massarr, self.header_struct.Npart):
            if x == 0.0:
                result += n
        return result
    @late
    def number_of_gas_particles(self):
        return self.header_struct.Npart[self.GAS]
        
    def collect_values(self, values):
        offset = 0
        result = []
        for name, times, formatchar in self.header_struct_format:
            if times > 1:
                result.append(values[offset: offset+times])
            else:
                result.append(values[offset])
            offset += times
        return result
        
    def load_header(self, file):
        header_bytes = self.read_fortran_block(file)
        values = struct.unpack(self.header_format, header_bytes[0:self.header_size])
        values = self.collect_values(values)
        self.header_struct = self.struct_class(*values)
    
    def load_body(self, file):
        self.positions = self.read_fortran_block_float_vectors(file)
        self.velocities = self.read_fortran_block_float_vectors(file)
        self.ids = self.read_fortran_block_uints(file)
    
        if self.total_number_of_particles_with_variable_masses > 0:
            self.masses = self.read_fortran_block_floats(file)
        else:
            self.masses = None
    
        if self.number_of_gas_particles > 0:
            self.u = self.read_fortran_block_floats(file)
        
        if self.is_initial_conditions_format:
            self.density = None
            self.u = None
            self.hsml = None
            self.pot = None
            self.acc = None
            self.da_dt = None
            self.dt = None
            return
        
        if self.number_of_gas_particles > 0:
            self.density = self.read_fortran_block_floats(file)
            self.hsml = self.read_fortran_block_floats(file)
        else:
            self.u = None
            self.density  = None
            self.hsml = None
            
        if self.has_potential_energy:
            self.pot = self.read_fortran_block_floats(file)
        else:
            self.pot = None
            
        if self.has_acceleration:
            self.acc = self.read_fortran_block_floats(file)
        else:
            self.acc = None
            
        if self.has_rate_of_entropy_production:
            self.da_dt = self.read_fortran_block_floats(file)
        else:
            self.da_dt = None
            
        if self.has_timestep:
            self.dt = self.read_fortran_block_floats(file)
        else:
            self.dt = None
    
    def new_sets_from_arrays(self):
        offset = 0
        ids_per_set = []
        for x in self.header_struct.Npart:
            ids_per_set.append(self.ids[offset:offset+x])
            offset += x
            
        sets = [data.Particles(len(x), keys=x) for x in ids_per_set]
        offset = 0
        for x in sets:
            length = len(x)
            if length == 0:
                continue
            print offset,offset+length, self.positions[offset:(offset+length)]
            x.position = nbody_system.length.new_quantity(self.positions[offset:offset+length])
            x.velocity = nbody_system.speed.new_quantity(self.velocities[offset:offset+length])
            if not self.pot is None:
                x.potential_energy = nbody_system.energy.new_quantity(self.pot[offset:offset+length])
            if not self.acc is None:
                x.acceleration = nbody_system.acceleration.new_quantity(self.acc[offset:offset+length])
            if not self.dt is None:
                x.timestep = nbody_system.time.new_quantity(self.dt[offset:offset+length])
            offset += length
        
        offset = 0
        for x, mass in zip(sets, self.header_struct.Massarr):
            length = len(x)
            if length == 0:
                continue
            if mass == 0.0:
                x.mass = nbody_system.mass.new_quantity(self.masses[offset:offset+length])
                offset += length
            else:
                x.mass = nbody_system.mass.new_quantity(mass)
        
        if self.number_of_gas_particles > 0:
            gas_set = sets[self.GAS]
            unit = nbody_system.length / nbody_system.time ** 2
            gas_set.internal_energy = unit.new_quantity(self.u)
            unit = nbody_system.mass / nbody_system.length ** 3
            if not self.density is None:
                gas_set.rho = unit.new_quantity(self.density)
            if not self.hsml is None:
                gas_set.smoothing_length = nbody_system.length.new_quantity(self.hsml)
            
        
        return sets
        
            
        
    def load_file(self, file):
        self.load_header(file)
        self.load_body(file)
        return self.new_sets_from_arrays()
        
    
        
        
        
