from subprocess import call
from os import path, listdir
from amuse.community import get_amuse_root_dir
from amuse.test.amusetest import TestCase

amuse_sh_path = path.join(get_amuse_root_dir(), 'amuse.sh')
simple_examples_dir = path.dirname(path.abspath(__file__))
simple_scripts = [one_file for one_file in listdir(simple_examples_dir) if one_file[-3:] == '.py']

def path_to_script(filename):
    return path.join(simple_examples_dir, filename)

for script in simple_scripts:
    exec("def slowtest_" + script[:-3] + "(): " +
        "script_path = path_to_script('" + script + "') ; " +
        "call([amuse_sh_path, script_path])")


class TestTestingAllSimple(TestCase):
    
    def test_the_tests(self):
        for script in simple_scripts:
            self.assertTrue("slowtest_"+script[:-3] in globals())
        
