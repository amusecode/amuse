from amuse.support import options
from amuse.test import amusetest
import StringIO
import textwrap

class OptionsTestsClass(options.OptionalAttributes):
    option_sections = ('amuse',)
    
    def string_option(self):
        return "default string"
    
    def boolean_option(self):
        return False
        
    def int_option(self):
        return 2
        
    def float_option(self):
        return 2.5
    
class OptionsTests(amusetest.TestCase):
    ini_contents = textwrap.dedent("""
    [amuse]
    string_option=a string
    boolean_option=1
    int_option=1
    float_option=1.5
    [bhtree]
    mode=gpu
    int_option=3
    
    """)
    
    
    def test1(self):
        global_options = options.GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        instance = OptionsTestsClass()
        option = options.option(instance.string_option, global_options=global_options)
        self.assertEquals(option.__get__(instance, type(instance)), "a string")
        
    def test2(self):
        global_options = options.GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        instance = OptionsTestsClass()
        option = options.option(instance.int_option, type="int", global_options=global_options)
        self.assertEquals(option.__get__(instance, type(instance)), 1)
        
        option = options.option(instance.boolean_option, type="boolean", global_options=global_options)
        self.assertEquals(option.__get__(instance, type(instance)), True)
        
        option = options.option(instance.float_option, type="float", global_options=global_options)
        
        self.assertEquals(option.__get__(instance, type(instance)), 1.5) 
        
    def test3(self):
        global_options = options.GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        instance = OptionsTestsClass()
        instance.option_sections =('unknown')
        
        option = options.option(OptionsTestsClass.string_option, global_options=global_options)
        self.assertEquals(option.__get__(instance, type(instance)), "default string")
        
        option = options.option(OptionsTestsClass.int_option, type="int", global_options=global_options)
        self.assertEquals(option.__get__(instance, type(instance)), 2)
        
        option = options.option(OptionsTestsClass.boolean_option, type="boolean", global_options=global_options)
        self.assertEquals(option.__get__(instance, type(instance)), False)
        
        option = options.option(OptionsTestsClass.float_option, type="float", global_options=global_options)
        self.assertEquals(option.__get__(instance, type(instance)), 2.5)
        
    
    def test4(self):
        global_options = options.GlobalOptions()
        global_options.config.readfp(StringIO.StringIO(self.ini_contents))
        instance = OptionsTestsClass()
        
        option = options.option(OptionsTestsClass.int_option, type="int", sections=("bhtree",), global_options=global_options)
        self.assertEquals(option.__get__(instance, type(instance)), 3)
        
    
