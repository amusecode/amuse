import ConfigParser
import os.path
import os

from amuse.support.core import late
from amuse.support.exception import CoreException

class GlobalOptions(object):
    INSTANCE = None
    
    def __init__(self):
        self.config=ConfigParser.RawConfigParser()
    
    def load(self):
        if not self.rcfilepath is None:
            self.config.read(self.rcfilepath)   
        
        
    @late
    def amuse_rootdirectory(self):
        result = os.path.abspath(__file__)
        while not os.path.exists(os.path.join(result,'build.py')):
            result = os.path.dirname(result)
        return result
    
    @late
    def rcfilepath(self):
        result = os.path.join(os.getcwd(), self.rcfilename)
        if os.path.exists(result): 
            return result
        
        result = os.path.join(self.homedirectory, '.' + self.rcfilename)
        if os.path.exists(result): 
            return result
        
        result = os.path.join(self.amuse_rootdirectory, self.rcfilename)
        if os.path.exists(result): 
            return result
        else:
            return None
    
    @late
    def rcfilename(self):
        return 'amuserc'
        
    @late
    def homedirectory(self):
        path=''
        try:
            path=os.path.expanduser("~")
        except:
            pass
        if not os.path.isdir(path):
            for evar in ('HOME', 'USERPROFILE', 'TMP'):
                try:
                    path = os.environ[evar]
                    if os.path.isdir(path):
                        break
                except: pass
        if path:
            return path
        else:
            raise RuntimeError('please define environment variable $HOME')
            
    @classmethod
    def instance(cls):
        if cls.INSTANCE is None:
            cls.INSTANCE = cls()
            cls.INSTANCE.load()
        return cls.INSTANCE
        
    def get_value_for_option(self, option, instance):
        for x in option.get_sections(instance):
            if self.config.has_option(x, option.name):
                return option.get_value(instance, x, self.config)
        return option.get_defaultvalue(instance)
        
class option(object):
    
    def __init__(self, function = None, type = "string", name = None, sections = (), global_options = None):
        self.specification_method = function
            
        if not name is None:
            self.name = name
        
        self.sections = sections
        
        if hasattr(self, type.upper()):
            self.valuetype = getattr(self, type.upper())
        else:
            raise CoreException("'{0}' is not a valid type for option".format(type))
            
        self.validator = self.default_validator
        self.choices = ()
        
        if global_options is None:
            self.global_options = GlobalOptions.instance()
        else:
            self.global_options = global_options
    
    @late
    def name(self):
        return self.specification_method.__name__
        
    def __call__(self, function):
        self.specification_method = function
        self.name = self.specification_method.__name__
        return self
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        
        if self.name in instance._local_options:
            return instance._local_options[self.name]
        else:
            return self.global_options.get_value_for_option(self, instance)
            
    def __set__(self, instance, value):
        instance._local_options[self.name] = value
        
    def get_value(self, instance, section, options):
        return self.validator(self.valuetype(section, options))
        
    def get_defaultvalue(self, instance):
        return self.specification_method(instance)
        
    def INT(self, section, options):
        return options.getint(section, self.name)
        
    def FLOAT(self, section, options):
        return options.getfloat(section, self.name)
        
    def BOOLEAN(self, section, options):
        return options.getboolean(section, self.name)
        
    def STRING(self, section, options):
        return options.get(section, self.name)
        
    def default_validator(self, value):
        return value
    
    def get_sections(self, instance):
        result = []
        result.extend(self.sections)
        result.extend(instance.option_sections)
        return result
        
    
        
        
class OptionalAttributes(object):
    option_sections = ()
    
    def __init__(self, **optional_keyword_arguments):
        for key, value in optional_keyword_arguments.iteritems():
            if self.hasoption(key):
                setattr(self, key, value)
        
        
    def hasoption(self, name):
        the_type = type(self)
        return hasattr(the_type, name)
        
    @late
    def _local_options(self):
        return {}
    
