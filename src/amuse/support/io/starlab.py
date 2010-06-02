from amuse.support.data import core
from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.io import base

import sys
import xml.dom.minidom 
import pdb #use with pdb.set_trace()
import re

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
            raise Exception("\nConversion failure, parenthesis mismatch. Return: no output string written\n")

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

    def __init__(self):
        self.xmls = ""
        self.system = core.Particles()
        self.translator = {'N':'number','m':'mass','r':'position','v':'velocity','system_time':'timestamp'}
        self.timestamp = None

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
        if key in self.translator.keys():                                     
            amuse_key = self.translator[key]                                  
            if amuse_key == 'mass':                                           
                particle.mass = float(value)|nbody_system.mass                              
            if amuse_key == "position":                                       
                particle.position = self.convert2vec(value)|nbody_system.length                
            if amuse_key == "velocity":                                       
                particle.velocity = self.convert2vec(value)|nbody_system.speed
            if amuse_key == 'timestamp':
                self.timestamp = float(value)|nbody_system.time
        
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
        
    
    def load_string(self, string):
        x = Dyn2Xml()
        xml_string = x.convert_startlab_string_to_xml_string(string)
        xml2particles = Xml2Particles()
        xml2particles.parse_xml(xml_string)
        if not self.nbody_to_si_converter is None:
            result = core.ParticlesWithUnitsConverted(
                xml2particles.system,
                self.nbody_to_si_converter.as_converter_from_si_to_nbody()
            )
        else:
            result = xml2particles.system
        
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
        "if True returns the children of the root node, if False returns the root node"
        return True
        

    @base.format_option
    def nbody_to_si_converter(self):
        "starlab datafiles store nbody data, provide a converter to store si data (None means no converter)"
        return None
        
