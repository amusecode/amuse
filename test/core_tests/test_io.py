from amuse.support.io import base
from amuse.test import amusetest


class TestFileFormatProcessor(base.FileFormatProcessor):
    """
    Save files in a test format
    
    long description    
    """
    provided_formats = ['test','123']
    instance = None
    
    def __init__(self, filename = None, set=None, format=None):
        base.FileFormatProcessor.__init__(self, filename, set, format)
        TestFileFormatProcessor.instance = self
        self.stored = False
        
    @base.format_option
    def add_comma(self):
        """if True will add a comma between each value"""
        return True
    
    @base.format_option
    def save_fast(self):
        """if True will save faster but less accurate"""
        return False

    def store(self):
        self.stored = True
        
    def load(self):
        return 10
        
class FrameworkTests(amusetest.TestCase):
    
    def test1(self):
        options = TestFileFormatProcessor.get_options()
        self.assertTrue('add_comma' in options)
        self.assertTrue('save_fast' in options)
        TestFileFormatProcessor.register()
        base.write_set_to_file("test.txt", None, format="123")
        self.assertTrue(TestFileFormatProcessor.instance.stored)
        self.assertTrue(TestFileFormatProcessor.instance.add_comma)
        
    def test2(self):
        TestFileFormatProcessor.register()
        base.write_set_to_file("test.txt", None, format="123", add_comma = False)
        self.assertFalse(TestFileFormatProcessor.instance.add_comma)
        
    def test3(self):
        TestFileFormatProcessor.register()
        documentation =  base.write_set_to_file.__doc__
        print documentation
        self.assertTrue("**123**,\n      Save files in a test format" in documentation)
    
    def test4(self):
        options = base.get_options_for_format('123')
        name, description, default = options[0]
        self.assertEquals(name, 'add_comma')
        self.assertEquals(description, 'if True will add a comma between each value')
        self.assertEquals(default, True)
        name, description, default = options[1]
        self.assertEquals(name, 'save_fast')
        self.assertEquals(description, 'if True will save faster but less accurate')
        self.assertEquals(default, False)
        
        
