from amuse.support.data import core
from amuse.support.units import units
from amuse.support.units import nbody_system

import sys
import xml.dom.minidom 
import pdb #use with pdb.set_trace()
import re

class Dyn2Xml(object):
    
    def __init__(self):

        self._xml = ""
    
    def convert_to_xml(self,inputfile):
        """
            Converts legacy starlab fileformat .dyn to xml string
       
        """

        f = open(inputfile,'r')

        fooline = f.readline()
            
        #lines either start with (,) followed by name of stories or Particle
        #Log, Dynamics, Hydro, Star
        #or with parameters
        #m = 0.5 etc...
        #Here we do not care about structure, whether the particle is a sister
        #or daughter etc...
        #..Thinking about reimplementing this using regularexpressions...

        #The xml string, init and add the root        
        #We need the root, an xml file with only sisters is not well fromed as it is
        #interpreted as multiple roots then. xml.dom yields an error
        xmls = '<?xml version="1.0"?>\n<System>\n'
        
        #counting opening and closing of parenthesis to check match at the end.
        #we do not want bad xml files.
        
        openpar  = 0
        closepar = 0
        
        while fooline:

            #State machine...
            
            fooline = fooline.rstrip()
                   
            if fooline.startswith("("):
                #Either Particle or stories
                #don't care about structure but do count parenthesis though
                openpar += 1
                #strip the (
                fooline = fooline.lstrip("(")
                newline = "<"+fooline+">"
                
            elif fooline.startswith(")"):
                closepar +=1
                #strip the )
                fooline = fooline.lstrip(")")
                newline = "</"+fooline+">"
            else:
                newline = self.convert_storyline(fooline)
          
            
            xmls += newline + "\n"                
            #Get next line from file...                
            fooline = f.readline()

        xmls += "</System>\n"    

        f.close()

        if closepar!=openpar:
            raise("\nConversion failure, parenthesis mismatch. Return: no output string written\n")

        self._xmls = xmls
            
    def convert_storyline(self, fooline):
        
        if "===>" in fooline or fooline.startswith("       "):
            newline = fooline 
        else:
            splithere = fooline.find("=")
            newline = "  <pm "+fooline[:splithere].lstrip()+'= "'+fooline[splithere+1:].lstrip()+'" />'
        
        return newline
    
class Xml2Particles(object):

    def __init__(self):
        self.xmls = ""
        self.system = core.Particles()
        self.translator = {'N':'number','m':'mass','r':'position','v':'velocity'}

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

class ParticlesFromDyn(object):
    
    def __init__(self, dyn_filename=None, convert_nbody=None):

        dyn2xml = Dyn2Xml()
        dyn2xml.convert_to_xml(dyn_filename)
        xml_string = dyn2xml._xmls
        
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
