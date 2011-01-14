#!/usr/bin/env python

from optparse import OptionParser
import os.path
import sys
import re

def get_substitution_rules(
        value_of_periodic_flag = False, 
    ):
    result = []
    
    replace_periodic_re = re.compile(r'PARAMETER \(periodic=\.([A-Z]+)\.\)')
    fortran_boolean_value = "TRUE" if value_of_periodic_flag else "FALSE"
    sub_periodic = r'PARAMETER (periodic=.{0}.)'.format(fortran_boolean_value)
    
    result.append((replace_periodic_re, sub_periodic,))
    return result
    
def patch_global_h(
        value_of_periodic_flag = False, 
        is_dry_run = True,
        filename = "src/globals.h",
        name_of_backupfile = "src/globals.h.bck",
    ):
    
    rules = get_substitution_rules(value_of_periodic_flag)
    
    with open(filename, "r") as input:
        lines = input.readlines()
    
    if not is_dry_run:
        with open(name_of_backupfile, "w") as output:
            output.writelines(lines)
    
    output_lines = []
    for line in lines:
        for expression, substitution in rules:
            line =expression.sub(substitution, line)
        output_lines.append(line)
    
    
    if is_dry_run:
        sys.stdout.writelines(output_lines)
    else:
        with open(filename, "w") as output:
            output.writelines(output_lines)
    
    
    
def new_option_parser():
    result = OptionParser()
    result.add_option(
        "-p",
        "--periodic", 
        dest="value_of_periodic_flag",
        action="store_true",
        default = False,
        help="patch the globals.h file, set the periodic flag to true if set (false otherwise)"
    )
    result.add_option(
        "-n",
        "--try", 
        dest="is_dry_run",
        action="store_true",
        default = False,
        help="don't patch the file only output on standard input"
    )
    return result
    
if __name__ == "__main__":
    options, arguments = new_option_parser().parse_args()
    patch_global_h(**options.__dict__)

