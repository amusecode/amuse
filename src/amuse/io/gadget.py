import struct
import numpy

from collections import namedtuple

from amuse.io import base
from amuse.units import units
from amuse.units import nbody_system
from amuse.support.core import late

from amuse import datamodel

class GadgetFileFormatProcessor(base.FortranFileFormatProcessor):
    """
    Process a Gadget binary data file
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
        """The format of the header structure of the gadget file."""
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
            ('NallHW', 6, 'i'),
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
        
    @base.format_option
    def ids_are_keys(self):
        """Set to True if the file contains
        correct keys. Set to False to generate
        the keys in amuse and prove an id attribute for 
        the id's in the gadget file"""
        return True
        
    
    @base.format_option
    def equal_mass_array(self):
        """If filled with an array with masses > 0.0
        assume equal mass for the corresponding set"""
        return ([0.0] * 6) | nbody_system.mass
    
    @base.format_option
    def ids_are_long(self):
        """Set to true the ids will be written
        as longs in the gadget file"""
        return True
    
    @base.format_option
    def return_header(self):
        """Set to True to return both the particles and the header from the gadget file"""
        return False
    
    @base.format_option
    def write_header_from(self):
        """Pass a namedtuple to store non-default values in the header of the gadget file"""
        return None
    
    @base.format_option
    def convert_gadget_w_to_velocity(self):
        """Set to True to convert the w=sqrt(a)*dx/dt to (comoving) velocity in comoving integrations"""
        return False
    
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
        
        id_bytes = self.read_fortran_block(file)
        if len(id_bytes) == 4*self.total_number_of_particles:
            self.ids = numpy.frombuffer(id_bytes, dtype=self.uint_type)
        else:
            self.ids = numpy.frombuffer(id_bytes, dtype=self.ulong_type)
        
        
        if self.total_number_of_particles_with_variable_masses > 0:
            self.masses = self.read_fortran_block_floats(file)
        else:
            self.masses = None

        if self.number_of_gas_particles > 0:
            self.u = self.read_fortran_block_floats(file)[:self.number_of_gas_particles]
            #print self.u, self.number_of_gas_particles, len(self.u)
        else:
            self.u = None
        
        if self.is_initial_conditions_format:
            self.density = None
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
        
        if self.ids_are_keys:
            sets = [datamodel.Particles(len(x), keys=x) for x in ids_per_set]
        else:
            sets = [datamodel.Particles(len(x)) for x in ids_per_set]
            for set, x in zip(sets, ids_per_set):
                if len(set) > 0:
                    set.id = x
                
        offset = 0
        for x in sets:
            length = len(x)
            if length == 0:
                continue
                
            x.position = nbody_system.length.new_quantity(self.positions[offset:offset+length])
            x.velocity = nbody_system.speed.new_quantity(self.velocities[offset:offset+length])
            if self.convert_gadget_w_to_velocity:
                x.velocity *= numpy.sqrt(1.0 + self.header_struct.Redshift)
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
            unit = (nbody_system.length / nbody_system.time) ** 2
            gas_set.u = unit.new_quantity(self.u)
            unit = nbody_system.mass / nbody_system.length ** 3
            if not self.density is None:
                gas_set.rho = unit.new_quantity(self.density)
            if not self.hsml is None:
                gas_set.h_smooth = nbody_system.length.new_quantity(self.hsml)
            
        return sets
        
            
        
    def load_file(self, file):
        self.load_header(file)
        self.load_body(file)
        
        attribute_names = ["gas","halo","disk","bulge","stars","bndry"]
        values = self.new_sets_from_arrays()
        
        if self.return_header:
            attribute_names += [name for name, times, formatchar in self.header_struct_format]
            values += list(self.header_struct)
        return namedtuple("GadgetData", attribute_names)(*values)
    
    
    @late
    def sets_to_save(self):
        if isinstance(self.set, (tuple, list)):
            sets_to_save = self.set
        elif hasattr(self.set, 'key'):
            sets_to_save = (self.set, )
        else:
            raise Exception("The Gadget binary file writer can write a particle set or a list of sets but nothing else")
        
        if len(sets_to_save) < 6:
            sets_to_save = list(sets_to_save)
            sets_to_save.extend([()] * (6 - len(sets_to_save)))
        
        return sets_to_save
        
    def store_file(self, file):
        self.store_header(file)
        self.store_body(file)
        
    def header_field_value(self, name):
        if hasattr(self.write_header_from, name):
            return getattr(self.write_header_from, name)
        else:
            return self.default_header_field_value(name)
    
    def default_header_field_value(self, name):
        return self.default_header[name]
    
    def store_header(self, file):
        self.default_header = dict(
                Massarr = list(self.equal_mass_array.value_in(nbody_system.mass)),
                Time = 0.0,
                Redshift = 0.0,
                FlagSfr = 0,
                FlagFeedback = 0,
                Nall =  [len(x) for x in self.sets_to_save],
                FlagCooling = 0,
                NumFiles = 1,
                BoxSize = 0,
                Omega0 = 0,
                OmegaLambda = 0,
                HubbleParam = 0,
                FlagAge = 0.0,
                FlagMetals = 0.0,
                NallHW = [0]*6,
                flag_entr_ics = 0
        )
        if self.write_header_from is None:
            self.header_field_value = self.default_header_field_value
        
        self.header_struct = self.struct_class(
            Npart = [len(x) for x in self.sets_to_save],
            Massarr = list(self.header_field_value("Massarr")),
            Time = self.header_field_value("Time"),
            Redshift = self.header_field_value("Redshift"),
            FlagSfr = self.header_field_value("FlagSfr"),
            FlagFeedback = self.header_field_value("FlagFeedback"),
            Nall = list(self.header_field_value("Nall")),
            FlagCooling = self.header_field_value("FlagCooling"),
            NumFiles = self.header_field_value("NumFiles"),
            BoxSize = self.header_field_value("BoxSize"),
            Omega0 = self.header_field_value("Omega0"),
            OmegaLambda = self.header_field_value("OmegaLambda"),
            HubbleParam = self.header_field_value("HubbleParam"),
            FlagAge = self.header_field_value("FlagAge"),
            FlagMetals = self.header_field_value("FlagMetals"),
            NallHW = list(self.header_field_value("NallHW")),
            flag_entr_ics = self.header_field_value("flag_entr_ics")
        )
        
        bytes = bytearray(256)
        parts = list()
        for x in list(self.header_struct):
            if isinstance(x, list):
                parts.extend(x)
            else:
                parts.append(x)
                
        struct.pack_into(self.header_format, bytes, 0,  *parts)
        self.write_fortran_block(file, bytes)
    
    
    def _pack_sets(self, attributename, unit = None):
        result = []
        for x in self.sets_to_save:
            if len(x) > 0:
                if unit is None:
                    result.extend(getattr(x,attributename))
                else:
                    result.extend(getattr(x,attributename).value_in(unit))
        return result
    
    
    def store_body(self, file):
        self.positions = self._pack_sets('position', nbody_system.length)
        self.velocities = self._pack_sets('velocity', nbody_system.speed)
        
        if self.ids_are_keys:
            self.ids = self._pack_sets('key')
        else:
            self.ids = self._pack_sets('id')
            
        self.masses = []
        for equal_mass, x in zip(self.equal_mass_array, self.sets_to_save):
            if len(x) > 0 and not equal_mass > (0.0 | nbody_system.mass):
                self.masses.extend(x.mass.value_in(nbody_system.mass)) 
        if len(self.masses) == 0:
            self.masses = None
            
        
        number_of_gas_particles = len(self.sets_to_save[self.GAS])
        if number_of_gas_particles > 0:
            self.u = self.sets_to_save[0].u.value_in(nbody_system.potential)
        else:
            self.u = None
        
            
        self.write_fortran_block_float_vectors(file, self.positions)
        self.write_fortran_block_float_vectors(file, self.velocities)
        
        if self.ids_are_long:
            self.write_fortran_block_ulongs(file, self.ids)
        else:
            self.write_fortran_block_uints(file, self.ids)
        
        if not self.masses is None:
            self.write_fortran_block_floats(file, self.masses)
        
        if not self.u is None:
            self.write_fortran_block_floats(file, self.u)
            
        if self.is_initial_conditions_format:
            return
        
        if number_of_gas_particles > 0:
            self.density =  self.sets_to_save[0].rho.value_in(nbody_system.density)
            self.hsml = self.sets_to_save[0].h_smooth.value_in(nbody_system.length)
        else:
            self.density  = None
            self.hsml = None
        
        if self.has_potential_energy:
            self.pot = self._pack_sets('potential_energy', nbody_system.energy)
        else:
            self.pot = None
            
        if self.has_acceleration:
            self.acc = self._pack_sets('acceleration', nbody_system.acceleration)
        else:
            self.acc = None
            
        if self.has_rate_of_entropy_production:
            self.da_dt = None
        else:
            self.da_dt = None
            
        if self.has_timestep:
            self.dt = self._pack_sets('timestep', nbody_system.time)
        else:
            self.dt = None
        
        
        if not self.density is None:
            self.write_fortran_block_floats(file, self.density)
        
        if not self.hsml is None:
            self.write_fortran_block_floats(file, self.hsml)
            
        if not self.pot is None:
            self.write_fortran_block_floats(file, self.pot)
            
        if not self.acc is None:
            self.write_fortran_block_floats(file, self.acc)
            
        if not self.da_dt is None:
            self.write_fortran_block_floats(file, self.da_dt)
            
        if not self.dt is None:
            self.write_fortran_block_floats(file, self.dt)
            
