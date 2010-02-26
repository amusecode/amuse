import sys
import os

def get_path_to_test_results():
    dir = os.path.dirname(__file__)
    amuse_root_dir = os.path.dirname(dir)
    test_results_dir = os.path.join(amuse_root_dir, 'test_results')
    if os.path.exists(test_results_dir):
        return test_results_dir
    else:
        return './'
