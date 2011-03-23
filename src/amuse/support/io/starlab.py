from amuse.support.data import core
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.units import generic_unit_converter
from amuse.support.io import base

import sys
import xml.dom.minidom 
import pdb #use with pdb.set_trace()
import re

starlab_stellar_types_to_amuse_stellar_type = {
    "planet":units.stellar_type("Unknown stellar type"),
    "proto_star":units.stellar_type("Unknown stellar type"),
    "brown_dwarf":units.stellar_type("Unknown stellar type"),
    "main_sequence":units.stellar_type("Main Sequence star"),
    "hyper_giant":units.stellar_type("Second Asymptotic Giant Branch"),
    "hertzsprung_gap":units.stellar_type("Hertzsprung Gap"),
    "sub_giant":units.stellar_type("First Giant Branch"),
    "horizontal_branch":units.stellar_type("Core Helium Burning"),
    "super_giant":units.stellar_type("First Asymptotic Giant Branch"),
    "thorne_zytkow":units.stellar_type("Unknown stellar type"),
    "carbon_star":units.stellar_type("Unknown stellar type"),
    "helium_star":units.stellar_type("Main Sequence Naked Helium star"),
    "helium_giant":units.stellar_type("Giant Branch Naked Helium star"),
    "helium_dwarf":units.stellar_type("Helium White Dwarf"),
    "carbon_dwarf":units.stellar_type("Carbon/Oxygen White Dwarf"),
    "oxygen_dwarf":units.stellar_type("Oxygen/Neon White Dwarf"),
    "xray_pulsar":units.stellar_type("Neutron Star"),
    "radio_pulsar":units.stellar_type("Neutron Star"),
    "neutron_star":units.stellar_type("Neutron Star"),
    "black_hole":units.stellar_type("Black Hole"),
    "Disintegrated":units.stellar_type("Massless Supernova"),
    "SPZDCH_star":units.stellar_type("Unknown stellar type"),
    "static_star":units.stellar_type("Unknown stellar type"),
    "star_cluster":units.stellar_type("Unknown stellar type")
}

class Dyn2Xml(object):
    
    def convert_startlab_string_to_xml_string(self, string):
        xml_lines = ['<?xml version="1.0"?>\n<System>']
        
        openpar  = 0
        closepar = 0
        
        for line in string.splitlines():
            line = line.rstrip()
                   
            if line.startswith("("):
                openpar += 1
                line = line.lstrip("(")
                newline = "<"+line+">"
            elif line.startswith(")"):
                closepar +=1
                line = line.lstrip(")")
                newline = "</"+line+">"
            else:
                newline = self.convert_storyline(line)
          
            xml_lines.append(newline)
    
        xml_lines.append("</System>")
    
       
        if closepar!=openpar:
            raise base.IoException("\nConversion failure, parenthesis mismatch. Return: no output string written\n")
    
        return '\n'.join(xml_lines)
        
    def convert_to_xml(self,inputfile):
        with open(inputfile,'r') as f:
            string = f.read()
        
        return self.convert_startlab_string_to_xml_string(string)
            
    def convert_storyline(self, line):
        if "===>" in line or line.startswith("       "):
            return line 
        
        parts = line.split("=",1)
        return "<pm "+ parts[0].lstrip()+'= "'+parts[1].lstrip()+'" />'
    
class Xml2Particles(object):
    dynamics_mass_units = nbody_system.mass
    dynamics_time_units = nbody_system.time
    dynamics_length_units = nbody_system.length
    
    def __init__(self):
        self.xmls = ""
        self.system = core.Particles()
        self.translator = {
            #'N':('number', lambda x : int(x) ),
            'm':('mass', lambda x : float(x)|self.dynamics_mass_units) ,
            't':('time', lambda x : float(x)|self.dynamics_time_units) ,
            'r':('position', lambda x : self.convert2vec(x)|self.dynamics_length_units),
            'v':('velocity', lambda x : self.convert2vec(x)|self.dynamics_length_units / self.dynamics_time_units),
            'a':('acceleration', lambda x : self.convert2vec(x)|self.dynamics_length_units / (self.dynamics_time_units ** 2)),
            'system_time':('timestamp', lambda x: float(x)|self.dynamics_time_units),
            'M_env': ('envelope_mass',  lambda x: float(x)|units.MSun),
            'M_rel': ('relative_mass',  lambda x: float(x)|units.MSun),
            'M_core': ('core_mass',  lambda x: float(x)|units.MSun),
            'T_eff' : ('effective_temperature', lambda x: float(x)|units.K),
            'T_cur' : ('age', lambda x: float(x)|units.Myr),
            'L_eff' : ('effective_luminocity', lambda x: float(x)|units.LSun),
            'Type'  : ('stellar_type', self.convert_starlab_stellar_type_to_amuse)
        }
        self.timestamp = None
        self.mass_scale = None
        self.size_scale = None
        self.time_scale = None
        
    def convert_starlab_stellar_type_to_amuse(self, string):
        if string in starlab_stellar_types_to_amuse_stellar_type:
            return starlab_stellar_types_to_amuse_stellar_type[string]
        else:
            return units.stellar_type("Unknown stellar type")
            
    def add_particle_with_parameters(self, subnode, parent):
        added_particle = self.system.add_particle(core.Particle())  
           
        self._recursive_parse_node_into_particles(
            subnode,
            added_particle, 
            parent = added_particle
        )
        
        if not parent is None:
            parent.add_child(added_particle)

    def _recursive_parse_node_into_particles(self, xml_node, particle_or_particleset, parent = None):
        node_list = xml_node.childNodes
        for subnode in node_list:
            if subnode.nodeType == subnode.ELEMENT_NODE:
                if subnode.tagName == "System": #overslaan
                    self._recursive_parse_node_into_particles(subnode, particle_or_particleset)
                elif subnode.tagName in ["Log", "Dynamics", "Star", "Hydro"]:  #overslaan
                    self._recursive_parse_node_into_particles(subnode, particle_or_particleset) 
                elif subnode.tagName == "Particle":
                    self.add_particle_with_parameters(subnode, parent)
                elif subnode.tagName == u"pm":                                
                    key = subnode.attributes.keys()[0]
                    value = subnode.getAttribute(key)
                    self.copy_starlab_parameter_to_star(key, value, particle_or_particleset)
                
    
    def copy_starlab_parameter_to_star(self, key, value, particle):
        if key == 'mass_scale':
            self.mass_scale = float(value)
        elif key == 'size_scale':
            self.size_scale = float(value)
        elif key == 'time_scale':
            self.time_scale = float(value)
        elif key in self.translator.keys():                                     
            amuse_attribute_name, conversion_function = self.translator[key]   
            setattr(particle,  amuse_attribute_name, conversion_function(value)) 
        
    def convert2vec(self, attribute):
        
        maybe_vec = attribute.split()

        if len(maybe_vec)>1:
            barout = [float(f) for f in maybe_vec]
        else:
            barout = float(maybe_vec[0])
        return barout

    def walk_me(self, SCL):

        self.xmls += self.elementstr(SCL,"open") +"\n"
            
        for attrlistitem in SCL.attriblist:

            attribute = attrlistitem[0]
            modified = attrlistitem[2]
            
            if attribute == "Sister" :
                if SCL.nsisters > 0:
                    for i in range(SCL.nsisters):
                        self.walk_me(SCL.Sister[i])
                        
            elif attribute in ["Log","Dynamics","Star","Hydro"]:
                self.walk_me(SCL.__getattribute__(attribute))
                
            elif not attribute in ["nsisters"]:
                foostr = str(SCL.__getattribute__(attribute))
                if not foostr == "" and modified:
                    #dont bother if empty
                    foostr = self.convert2propstr(foostr)
                    if attribute == "text":
                        #text is not a <pm />
                        self.xmls += foostr + "\n"
                    else:
                        self.xmls += '  <pm '+attribute+' = "'+ foostr +'" />\n'

        self.xmls += self.elementstr(SCL,"close")+"\n"
        
        
    def elementstr(self,SCL,oc):
        """
        write an element node,
        If particle is root then it is called System
        """
        if oc == "open":
            prestr = "<"
        elif oc == "close":
            prestr = "</"

        if SCL.__class__.__name__ == "Particle":
            if SCL.isroot:
                barstr =  prestr+"System>"
            else:
                barstr = prestr+SCL.__class__.__name__+">"        
        else:
            barstr = prestr+SCL.__class__.__name__.lstrip("_")+">"

        return barstr
    
    def convert2propstr(self,foostr):
        """
            if foostr is a vector we dont want the brackets and the ,'s
        """
        return foostr.replace("[","").replace("]","").replace(","," ")

    def makexml(self):
        """
            clear and init xmls
        """
        self.xmls = '<?xml version="1.0"?>\n'
        self.walk_me(self.System)
    
    def dumpxml(self):
        """
            Convert System to xml string and dump it on the screen
        """
        self.makexml()
        print self.xmls

    def savexml(self, file):
        """
            Convert System to xml string and put it in file with name file
        """
        f = open(file,'w')
        self.makexml()
        f.writelines(self.xmls)
        f.close()
       
    def loadxml(self,file):
        """
            Parse file into dom doc and then create particles accordingly 
        """
        f = open(file,'r')
        doc = xml.dom.minidom.parse(f)
        f.close()

        #copy the doc to a system object
        self._recursive_parse_node_into_particles(doc, self.system)

    def parse_xml(self, xmlstring):
        """
            Parse xml string into dom doc and then create particles accordingly 
        """
        doc = xml.dom.minidom.parseString(xmlstring)
        self._recursive_parse_node_into_particles(doc, self.system)
        if not self.timestamp is None:
            self.system.savepoint(self.timestamp)

class ParticlesFromDyn(object):
    
    def __init__(self, dyn_filename=None, convert_nbody=None):

        dyn2xml = Dyn2Xml()
        xml_string = dyn2xml.convert_to_xml(dyn_filename)
        
        xml2particles = Xml2Particles()
        err1 = xml2particles.parse_xml(xml_string)

        if convert_nbody is None:
            self.convert_nbody = None
            self.Particles = xml2particles.system
            return
        else:
            self.convert_nbody = convert_nbody
            self.Particles = core.ParticlesWithUnitsConverted(
                xml2particles.system,
                self.convert_nbody.as_converter_from_si_to_nbody()
            )
        
        
    def number_of_particles(self):
        return len(self.Particles)

class Particles2Dyn(object):
    
    def convert_to_string(self, particles):
        lines = []
        prefix = "  "
        lines.append("(Particle")
        lines.append(prefix + "N = " + str(len(particles)))
        lines.append("(Log")
        lines.append(")Log")
        #timestamp = particles.get_timestamp() # Timestamp is only saved if not None
        timestamp = particles.get_timestamp() or (0|nbody_system.time) # Timestamp is always saved
        if not timestamp is None:
            lines.append("(Dynamics")
            lines.append(prefix + "system_time  =  " + str(timestamp.value_in(nbody_system.time)))
            lines.append(")Dynamics")
        for index, x in enumerate(particles):
            lines.append("(Particle")
            lines.append(prefix + "i = " + str(index))
            lines.append(prefix + "N = " + str(1))
            
            float_format = '{0}'.format
            lines.append("(Dynamics")
            lines.append(prefix + "m = " + float_format(x.mass.value_in(nbody_system.mass)))
            r = x.position.value_in(nbody_system.length)
            lines.append(prefix + "r = " + " ".join(map(float_format, r)))
            v = x.velocity.value_in(nbody_system.speed)
            lines.append(prefix + "v = " + " ".join(map(float_format, v)))
            
            lines.append(")Dynamics")
            
            lines.append(")Particle")
            
        
        lines.append(")Particle")
        return '\n'.join(lines)

class StarlabFileFormatProcessor(base.FullTextFileFormatProcessor):
    """
    Process a Starlab binary structured file
    """
    
    provided_formats = ['starlab', 'dyn']
    
    def __init__(self, filename = None, stream = None, set = None, format = None):
        base.FileFormatProcessor.__init__(self, filename, set, format)
        
    
    def _is_valid_scaling_factor(self, factor):
        return not factor is None and not factor == -1.0
        
    def load_string(self, string):
        x = Dyn2Xml()
        xml_string = x.convert_startlab_string_to_xml_string(string)
        xml2particles = Xml2Particles()
        xml2particles.dynamics_mass_units = self.dynamics_mass_units
        xml2particles.dynamics_time_units = self.dynamics_time_units
        xml2particles.dynamics_length_units = self.dynamics_length_units
        xml2particles.parse_xml(xml_string)
        unit_converter = None
        if not self.nbody_to_si_converter is None:
            unit_converter = self.nbody_to_si_converter.as_converter_from_si_to_nbody()
        elif self.must_scale:
            if not self._is_valid_scaling_factor(xml2particles.mass_scale):
                unit_converter = None
            elif not self._is_valid_scaling_factor(xml2particles.time_scale):
                unit_converter = nbody_system.nbody_to_si(
                    (1.0 / xml2particles.mass_scale) | units.MSun,
                    (1.0 / xml2particles.size_scale) | units.RSun,
                ).as_converter_from_si_to_nbody()
            else:
                unit_converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
                    (1.0 / xml2particles.mass_scale) | units.MSun,
                    (1.0 / xml2particles.size_scale) | units.RSun,
                    (1.0 / xml2particles.time_scale) | units.Myr,
                ).as_converter_from_si_to_generic()
                
        if unit_converter is None:
            result = xml2particles.system
        else:
            result = core.ParticlesWithUnitsConverted(
                xml2particles.system,
                unit_converter
            )
        
        if self.return_children:
            return result[0].children()
        else:
            return result[0]
        
    def store_string(self):
        if not self.nbody_to_si_converter is None:
            particles = core.ParticlesWithUnitsConverted(
                self.set,
                self.nbody_to_si_converter.as_converter_from_nbody_to_si()
            )
        else:
            particles = self.set
        
        x = Particles2Dyn()
        return x.convert_to_string(particles)
       
    
    @base.format_option
    def return_children(self):
        """If True returns the children of the root node, if False returns the root node (defaults to True)"""
        return True
        

    @base.format_option
    def nbody_to_si_converter(self):
        """Starlab datafiles store stellar dynamics properties in scaled nbody values, 
        provide a converter to store si data (defaults to None). 
        Value None means no converter, or use scaling values provided in the file"""
        return None
        
    
    @base.format_option
    def dynamics_mass_units(self):
        """The m field in the dynamics section of a starlab file can be in MSun or in scaled units , defaults to scaled units (nbody_system.mass).
        When set to scaled units, AMUSE will convert the units if scaling parameters are also given in the file. See the `must_scale` option
        to turn this scaling off"""
        return nbody_system.mass
        
    @base.format_option
    def dynamics_time_units(self):
        """The time fields in the dynamics section of a starlab file can be in Myr or in scaled units , defaults to scaled units (nbody_system.time).
        When set to scaled units, AMUSE will convert the units if scaling parameters are also given in the file. See the `must_scale` option
        to turn this scaling off.
        """
        return nbody_system.time
        
    @base.format_option
    def dynamics_length_units(self):
        """The length fields in the dynamics section of a starlab file can be in parsec or in scaled units , defaults to scaled units (nbody_system.length)
        When set to scaled units, AMUSE will convert the units if scaling parameters are also given in the file. See the `must_scale` option
        to turn this scaling off.
        """
        return nbody_system.length
        
    @base.format_option
    def must_scale(self):
        """If True use the scaling values from the file, if False do not scale the stellar dynamics properties.
        Only used when no nbody to si converter has been set.
        """
        return True
        
