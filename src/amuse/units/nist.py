import urllib.request, urllib.error, urllib.parse, urllib.request, urllib.parse, urllib.error
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
from amuse.units.si import *
from amuse.units.derivedsi import *

"""
ADDITIONAL_DERIVED_CONSTANTS = \
"""
pi = numpy.pi
hbar = h / (2.0 * numpy.pi)
four_pi_stefan_boltzmann = 4.0 * numpy.pi * Stefan_hyphen_Boltzmann_constant
mu0 = 4 * numpy.pi * 1.e-7 | N/A**2
eps0 = mu0**-1 * c**-2
sidereal_day = 86164.100352 | s
#machine constants
eps = numpy.finfo(numpy.double).eps
precision = int(numpy.log10(2/eps))
"""
class GetConstantsFromFiles(object):
    
    def __init__(self):
        self.nist_table = ""
        self.local_table = ""
        self.translator_table = []
        self.directory = os.path.dirname(__file__)

    def get_table_from_url(self):
        f = urllib.request.urlopen(NIST_URL)
        self.nist_table = f.read()
        f.close()

    def save_table_as(self, filename):
        f = open(os.path.join(self.directory, 'nist.txt'), 'w')
        f.write(self.nist_table)
        f.close()
     
    def get_table_from_file(self):
        f = open(os.path.join(self.directory, 'nist.txt'), 'r') # CODATA 2006, for CODATA 2010 use 'nist2010.txt'
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
            print(i)
            
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

        self.siunits = dir(si)+dir(derivedsi)
        
    def test_regexp(self, regexp):
        lines =self.table.splitlines(1)
        for i,line in enumerate(lines):
            if i>80:
                break
            print(re.findall(regexp, line))

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
        for n, line in enumerate(lines):
            if "----------------------" in line:
                number_of_header_lines = n + 1
                break
        firstline = lines[number_of_header_lines]
        namestr_length = len(firstline) - len(firstline[firstline.find("   "):].lstrip())
        column_index_of_uncertainty = len(firstline) - len(firstline[namestr_length+21+firstline[namestr_length+21:].find(" "):].lstrip())
        column_index_of_unit = len(firstline) - len(firstline[column_index_of_uncertainty+21+firstline[column_index_of_uncertainty+21:].find(" "):].lstrip())
        for i in lines[number_of_header_lines:]:
            namestr = i[0:namestr_length]

            marker1 = column_index_of_uncertainty
            marker2 = column_index_of_unit

            while 1:
                if i[marker1-1]=='\x20':
                    break
                else:
                    marker1+=1
            while 1:
                if i[marker2-1]=='\x20':
                    break
                else:
                    marker2+=1

            nrs=[]
            nrs.append(i[namestr_length:marker1])
            nrs.append(i[marker1:marker2])
            unitstr = i[marker2:]
                
            unitstr = unitstr.strip().replace(' ','*').replace('^','**')

            new_name = self.translate(namestr.rstrip(' ').replace(' ','_').replace('.','').replace('{','X').replace('}','X').replace('(','X').replace(')','X').replace('-','_hyphen_').replace(',','_and_').replace('/','_div_'))
            error.append(nrs[1].replace(' ',''))
            if len(unitstr)==1:
                this_unit = "none\n"
            else:
                this_unit = unitstr

            self.nisttable.append([new_name, float(i[namestr_length:marker1].replace(' ','').replace('...','')), unitstr])

    def sort_units(self):
        for entry in self.nisttable:
            if entry[2] in self.siunits:
                self.nisttablebaseunits.append(entry)
            elif entry[2] == '':
                self.nisttablenoneunits.append(entry)
            elif set(re.split('[*/^]',re.sub('\*\*-?[0-9.]*','',entry[2]))).issubset(set(self.siunits)):
                self.nisttablederivedunits.append(entry)
            else:
                self.nisttabledependingunits.append(entry)
    
    def print_list_of_units(self, unitlist):
        for name, value, unit in unitlist:
            self.nistfile += ("{0} = {1} | {2}\n".format(name, value, unit or "none"))

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
    print("Generating constants.py...", end=' ')
    Constants().generate_constants()
    print(" done!")
