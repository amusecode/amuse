import unittest
import os.path
import shutil

from  path_to_test_results import get_path_to_test_results
from amuse.legacy.support import create_dir

class CreateADirectoryAndPopulateItWithFilesForALegacyCodeTest(unittest.TestCase):
    
    def test1(self):
        instance = create_dir.CreateADirectoryAndPopulateItWithFilesForALegacyCode()
        instance.name_of_the_code_interface_class = 'TestCode'
        self.assertEquals(instance.name_of_the_legacy_code , 'testcode')
        self.assertTrue(instance.path_of_the_legacy_code.endswith('testcode'))
    
    def test2(self):
        root = get_path_to_test_results()
        working_dir = os.path.join(root, 'testcode')
        
        print working_dir
        
        if os.path.exists(working_dir):
            shutil.rmtree(working_dir)
                        
        instance = create_dir.CreateADirectoryAndPopulateItWithFilesForALegacyCode()
        instance.name_of_the_code_interface_class = 'TestCode'
        instance.path_of_the_root_directory = root
        
        instance.start()
        
        print os.path.exists(working_dir)
        
        self.assertTrue(working_dir)
