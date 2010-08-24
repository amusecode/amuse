import ConfigParser
import os.path
import os

from amuse.support.core import late
from amuse.support import exceptions

class GlobalOptions(object):
    INSTANCE = None
    
    def __init__(self):
        self.config=ConfigParser.RawConfigParser()
    
    def load(self):
        if not self.rcfilepath is None:
            self.config.read(self.rcfilepath)   
        
        
    @late
    def amuse_rootdirectory(self):
        if 'AMUSE_ROOT_DIR' in os.environ:
            return os.environ['AMUSE_ROOT_DIR']
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
    """Decorator to define an option
    
    :argument type: Type of the value, used when reading from the configuration file. 
       Can be "string", "int", "float" or "boolean". Defaults to "string"
    :argument sections: Sections in the configuration file to search for
        the option value, must be an array of strings
    :argument choices: When given will check if the value of the option
        is in the array (must be a list or set of objects)
    :argument name: By default the name of the option in the configuration file
        is the same as the name of the function, use this argument to 
        use a different name (not recommended)
    
    Options can only be defined on subclasses of :class:`OptionalAttributes`
    
    
    
    """
    def __init__(self, function = None, type = "string", name = None, sections = (), choices = (), global_options = None):
        self.specification_method = function
            
        if not name is None:
            self.name = name
        
        if not self.specification_method is None:
            self.__doc__ = self.specification_method.__doc__
            
        self.sections = sections
        
        if hasattr(self, type.upper()):
            self.valuetype = getattr(self, type.upper())
        else:
            raise exceptions.CoreException("'{0}' is not a valid type for option".format(type))
            
        self.validator = self.default_validator
        self.choices = set(choices)
        if self.choices:
            self.validator = self.choice_validator
        
        if global_options is None:
            self.global_options = GlobalOptions.instance()
        else:
            self.global_options = global_options
    
    @late
    def name(self):
        return self.specification_method.__name__
        
    def __call__(self, function):
        self.specification_method = function
        self.__doc__ = self.specification_method.__doc__
        return self
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        
        if self.name in instance._local_options:
            return instance._local_options[self.name]
        else:
            return self.global_options.get_value_for_option(self, instance)
            
    def __set__(self, instance, value):
        instance._local_options[self.name] = self.validator(value)
        
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
        
    def choice_validator(self, value):
        if not value in self.choices:
            raise exceptions.CoreException("{0} is not a valid choice for option '{1}', valid values are: {2}".format(value, self.name, sorted(self.choices)))
        return value
    
    def get_sections(self, instance):
        result = []
        result.extend(instance.option_sections)
        result.extend(self.sections)
        return result
        
    
        
        
class OptionalAttributes(object):
    """
    Abstract superclass for all classes supporting optional
    attributes.
    
    To support optional attributes a class must inherit (directly
    or indirectly) from this class.
    
    To support setting the attributes when an object is created
    the class must define a *catch-all* keyword argument
    in the **__init__** function and send this argument to the
    __init__ of the superclass.
    
    The values of options are first searched for in the sections
    given in the **option_sections** attribute of the class (empty
    by default). Next the sections of the option are searched.
    
    For example::
        
        class MyInterface(OptionalAttributes):
            option_sections = ('mysection',)
            
            def __init__(self, **options):
                OptionalAttributes.__init__(self, **options)
            
            @option(type="int", choices=(5,10,15), sections=('try',))
            def number_of_tries(self):
                "number of times to try to connect"
                return 5
    
    To code will first search for the value of the option in the 
    *mysection*, if no value is found the *try* section is searched. 
    For the following configuration file the **number_of_tries** 
    attribute will be 10 as the *mysection* section is searched first.
    
    :: ini
    
        [mysection]
        number_of_tries = 10
        [try]
        number_of_tries = 5
        
    The value of the option can be overriden by specifying it when
    creating an object of the class.
    
    :: python
        
        x = MyInterface(number_of_tries = 15)
        print x.number_of_tries
        15
    
    """
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
    
