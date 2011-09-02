import urllib2, urllib
import difflib
import os.path
import re

from hashlib import md5

from amuse.units import si
from amuse.units import derivedsi

NIST_URL = "http://132.229.222.6:9000/nistdata"
MESSAGE = \
"""
#This is an auto generated file, do not change manually. Instead if you want to add constants
#or change them, change the nist.txt file and run nist.py

import numpy
from amuse.units.units import *

"""
ADDITIONAL_DERIVED_CONSTANTS = \
"""
pi = numpy.pi | none
hbar=h/(2.0*numpy.pi)
four_pi_stefan_boltzmann = 4.0*numpy.pi*Stefan_hyphen_Boltzmann_constant
mu0=4*numpy.pi*1.e-7 | N/A**2
eps0=mu0**-1*c**-2
#machine constants
eps =numpy.finfo(numpy.double).eps
precision = int(numpy.log10(2/eps))
"""
class GetConstantsFromFiles(object):
    
    def __init__(self):
        self.nist_table = ""
        self.local_table = ""
        self.translator_table = []
        self.directory = os.path.dirname(__file__)

    def get_table_from_url(self):
        f = urllib2.urlopen(NIST_URL)
        self.nist_table = f.read()
        f.close()

    def save_table_as(self, filename):
        f = open(os.path.join(self.directory, 'nist.txt'), 'w')
        f.write(self.nist_table)
        f.close()
     
    def get_table_from_file(self):
        f = open(os.path.join(self.directory, 'nist.txt'), 'r')     
        self.nist_table = f.read()                             
        f.close()                                    

    def check_current_file_with_table(self):
        md5sum_local = md5()                         
        md5sum_local.update(self.local_table)              
        md5sum_local_val = md5sum_local.hexdigest()  
        md5sum_wgot = md5()                          
        md5sum_wgot.update(self.nist_table)          
        md5sum_wgot_val = md5sum_wgot.hexdigest()    
        return md5sum_local_val == md5sum_wgot_val   

    def compare_char_by_char(self):

        self.nist_table.lstrip('\n')
        mydiff = difflib.unified_diff(self.nist_table.splitlines(1), self.local_table.splitlines(1))
        for i in list(mydiff):
            print i
            
    def get_translator(self):
        f = open(os.path.join(self.directory, 'translator.txt'), 'r')     
        lines = f.readlines()

        for i, s in enumerate(lines):
            cols = s.split(',')
            self.translator_table.append(cols)

        f.close()         

class Constants(object):
    def __init__(self):
        self.I = GetConstantsFromFiles()
        #I.get_table_from_url()
        self.I.get_table_from_file()
        self.table = self.I.nist_table
        self.I.get_translator()
        self.translator = self.I.translator_table
        self.nistfile = MESSAGE
        
        self.nisttable = []
        self.nisttablederivedunits = []
        self.nisttablenoneunits = []
        self.nisttablebaseunits = []
        self.nisttabledependingunits = []

        self. siunits = dir(si)+dir(derivedsi)
        
    def test_regexp(self, regexp):
        lines =self.table.splitlines(1)
        for i,line in enumerate(lines):
            if i>80:
                break
            print re.findall(regexp, line)

    def translate(self, to_translate):
        list = [s[1] for s in self.translator if to_translate == s[0]]
        if list == []:
            return to_translate.lstrip(' ')
        else:
            return list[0].strip('\n')

    def list_constants(self):
        error =[]
        value = []
        name = []
        unit = []

        lines =self.table.splitlines(1)
        for tel in range(9): lines.pop(0)
        for i in lines:
            namestr = i[0:54]

            marker1 = 76
            marker2 = 98

            while 1:
                if i[marker1]=='\x20':
                    break
                else:
                    marker1+=1
            while 1:
                if i[marker2]=='\x20':
                    break
                else:
                    marker2+=1

            nrs=[]
            nrs.append(i[54:marker1])
            nrs.append(i[marker1:marker2])
            unitstr = i[marker2:-1]
                
            unitstr = unitstr.lstrip(' ').replace(' ','*').replace('^','**')

            new_name = self.translate(namestr.rstrip(' ').replace(' ','_').replace('.','').replace('{','X').replace('}','X').replace('(','X').replace(')','X').replace('-','_hyphen_').replace(',','_and_').replace('/','_div_'))
            error.append(nrs[1].replace(' ',''))
            if len(unitstr)==1:
                this_unit = "none\n"
            else:
                this_unit = unitstr

            self.nisttable.append([new_name, float(i[54:marker1].replace(' ','').replace('...','')), unitstr])

    def sort_units(self):

        for i in self.nisttable:
            if i[2].replace('\n','') in self.siunits:
                self.nisttablebaseunits.append(i)
            elif i[2] == '':
                self.nisttablenoneunits.append(i)
            else:
                derived_unit = i[2].replace('\n','')
                units_in_derived_unit = re.split('[*/^]',derived_unit)
                l = list(set(units_in_derived_unit)-set(self.siunits))
                for ind, j in enumerate(l): l[ind] = re.sub('-?[0-9.]*','',j)
                no_empty_l = []
                for ind, j in enumerate(l): 
                    if not j=='': no_empty_l.append(j)
                if len(no_empty_l)==0:
                    self.nisttablederivedunits.append(i)
                else:
                    self.nisttabledependingunits.append(i)
    
    def print_list_of_units(self, unitlist):
        for i in unitlist:
            i[2]=i[2].strip('\n')
            if i[2]=='':
                i[2] = 'none'
            self.nistfile += ("{0} = {1} | {2}\n".format(i[0], i[1], i[2]))

    def generate_constants(self):
        self.list_constants()
        self.sort_units()
        self.nistfile +=  "#BASE UNITS***********************************************\n"
        self.print_list_of_units(self.nisttablebaseunits)
        self.nistfile += "#DERIVED UNITS***********************************************\n"
        self.print_list_of_units(self.nisttablederivedunits)
        self.nistfile += "#RATIOS ***********************************************\n"
        self.print_list_of_units(self.nisttablenoneunits)
        self.nistfile += "#DERIVED CONSTANTS***********************************************"
        self.nistfile += ADDITIONAL_DERIVED_CONSTANTS
        self.nistfile += '#DROPPED UNITS***********************************************\n"""'
        self.print_list_of_units(self.nisttabledependingunits)
        self.nistfile +='"""\n'


        f = open(os.path.join(self.I.directory, 'constants.py'), 'w')
        f.write(self.nistfile)
        f.close()

if __name__  == "__main__":
    print "Generating constants.py...",
    Constants().generate_constants()
    print " done!"
