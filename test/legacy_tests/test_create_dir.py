from amuse.test import amusetest
import os.path
import shutil
import subprocess

from amuse.legacy.support import create_dir

class CreateADirectoryAndPopulateItWithFilesForALegacyCodeTest(amusetest.TestCase):
    
    def test1(self):
        instance = create_dir.CreateADirectoryAndPopulateItWithFilesForALegacyCode()
        instance.name_of_the_code_interface_class = 'TestCode'
        self.assertEquals(instance.name_of_the_legacy_code , 'testcode')
        self.assertTrue(instance.path_of_the_legacy_code.endswith('testcode'))
    
    def test2(self):
        root = self.get_path_to_results()
        working_dir = os.path.join(root, 'testcode')
        
        print working_dir
        
        if os.path.exists(working_dir):
            shutil.rmtree(working_dir)
                        
        instance = create_dir.CreateADirectoryAndPopulateItWithFilesForALegacyCode()
        instance.name_of_the_code_interface_class = 'TestCode'
        instance.path_of_the_root_directory = root
        self.assertEquals(instance.reference_to_amuse_path ,'../..')
        
        instance.start()
        
        self.assertTrue(os.path.exists(working_dir))
        self.assertTrue(os.path.exists(os.path.join(working_dir,'interface.py')))
        self.assertTrue(os.path.exists(os.path.join(working_dir,'Makefile')))
        
        call = subprocess.Popen(
            'make',
            cwd=working_dir, 
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        outputstring, errorstring = call.communicate()
        
        print errorstring
    
        self.assertEquals(call.returncode, 0)
        
        self.assertTrue(os.path.exists(os.path.join(working_dir,'worker_code')))
        
        
        call = subprocess.Popen(
            ['make', 'clean'],
            cwd=working_dir, 
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )
        outputstring, errorstring = call.communicate()
        
        print errorstring
    
        self.assertEquals(call.returncode, 0)
        
        
            
        
