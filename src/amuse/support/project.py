"""
Some general functions useful for AMUSE science projects
"""
import os
import os.path
import shutil
import pickle

def new_working_directory(script_filename=None, sub_directories=[]):
    """
    Call this function from your script to create a new directory and move 
    into it, for storing all your simulation output. Invoke it with:
        new_working_directory(__file__)
    to copy the current version of your script to this new directory for 
    book-keeping purposes.
    """
    i = 0
    while os.path.exists("run_{0:=03}".format(i)):
        i += 1
    new_directory = "run_{0:=03}".format(i)
    os.mkdir(new_directory)
    print "Created new directory for output:", new_directory
    for sub_directory in sub_directories:
        os.mkdir(os.path.join(new_directory, sub_directory))
    if not script_filename is None:
        shutil.copy(script_filename, new_directory)
    os.chdir(new_directory)

def store_results_in_file(results, datafile):
    with open(datafile, 'wb') as outfile:
        pickle.dump(results, outfile)

def load_results_from_file(datafile):
    with open(datafile, 'rb') as infile:
        results = pickle.load(infile)
    return results
