import configparser
import os.path
import os
import platform
from io import StringIO

from amuse.support.core import late
from amuse.support import exceptions

try:
    import pkg_resources
except ImportError:
    pkg_resources = None


class GlobalOptions(object):
    INSTANCE = None

    def __init__(self):
        self.config = configparser.ConfigParser()
        self.overriden_options = {}

    def load(self, preloadfp=None):
        if pkg_resources is not None:
            if pkg_resources.resource_exists('amuse', 'amuserc'):
                resourcerc = pkg_resources.resource_filename('amuse', 'amuserc')
                self.config.read(resourcerc)

        rootrc = os.path.join(self.amuse_rootdirectory, self.rcfilename)
        datarc = os.path.join(self.amuse_data_location, self.rcfilename)

        homedirrc = os.path.join(self.homedirectory, '.' + self.rcfilename)

        self.config.read(rootrc)
        self.config.read(datarc)

        if preloadfp is not None:
            self.config.read_file(preloadfp, "<amuserc>")

        self.config.read(homedirrc)
        self.config.read(os.path.join(self.homedirectory, '.' + platform.node() + '_' + self.rcfilename))
        if 'AMUSERC' in os.environ:
            self.config.read(os.environ['AMUSERC'])
        self.config.read(self.rcfilepath)

    @late
    def amuse_data_location(self):

        this = os.path.dirname(os.path.abspath(__file__))

        # installed
        result = os.path.abspath(os.path.join(this, "..", "..", "..", "..", "..", "share", "amuse"))
        if os.path.exists(os.path.join(result, 'build.py')):
            return result

        # for some virtualenv setups
        result = os.path.abspath(os.path.join(this, "..", "..", "..", "..", "..", "..", "share", "amuse"))
        if os.path.exists(os.path.join(result, 'build.py')):
            return result

        # in-place
        result = os.path.abspath(os.path.join(this, "..", "..", ".."))
        if os.path.exists(os.path.join(result, 'build.py')):
            return result

        raise exceptions.AmuseException("Could not locate AMUSE root directory! set the AMUSE_DIR variable")

    @late
    def amuse_rootdirectory(self):
        return os.path.dirname(os.path.dirname(__file__))

    @late
    def rcfilepath(self):
        return os.path.join(os.getcwd(), self.rcfilename)

    @late
    def rcfilename(self):
        return 'amuserc'

    @late
    def homedirectory(self):
        path = ''
        try:
            path = os.path.expanduser("~")
        except:
            pass
        if not os.path.isdir(path):
            for evar in ('HOME', 'USERPROFILE', 'TMP'):
                try:
                    path = os.environ[evar]
                    if os.path.isdir(path):
                        break
                except:
                    pass
        if path:
            return path
        else:
            raise RuntimeError('please define environment variable $HOME')

    @classmethod
    def instance(cls, preloadfp=None):
        if cls.INSTANCE is None:
            cls.INSTANCE = cls()
            cls.INSTANCE.load(preloadfp)
        return cls.INSTANCE

    def get_value_for_option(self, option, instance):
        if option.name in self.overriden_options:
            return self.overriden_options[option.name]

        for x in option.get_sections(instance):
            if self.config.has_option(x, option.name):
                return option.get_value(instance, x, self.config)
        return option.get_defaultvalue(instance)

    def to_ini_string(self):
        file = StringIO()
        self.config.write(file)
        return file.getvalue()

    def override_value_for_option(self, name, value):
        self.overriden_options[name] = value

    def read_from_ini_string(self, string):
        file = StringIO(string)
        self.config.read_file(file)


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

    def __init__(self, function=None, type="string", name=None, sections=(), choices=(), global_options=None):
        self.specification_method = function

        if name is not None:
            self.name = name

        if self.specification_method is not None:
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
            value = instance._local_options[self.name]
        else:
            value = self.global_options.get_value_for_option(self, instance)

        setattr(instance, self.name, value)

        return value

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
        if value not in self.choices:
            raise exceptions.CoreException("{0} is not a valid choice for option '{1}', valid values are: {2}".format(value, self.name, sorted(self.choices)))
        return value

    def get_sections(self, instance):
        result = []
        result.extend(instance.option_sections)
        result.append(instance.__class__.__name__)
        lastname = instance.__class__.__name__.split('.')[-1]
        if not lastname == instance.__class__.__name__:
            result.append(lastname)
            result.append(lastname.lower())
        else:
            result.append(instance.__class__.__name__.lower())
        result.extend(self.sections)
        return result

    def __setx__(self, instance, value):
        instance._local_options[self.name] = self.validator(value)
        try:
            delattr(instance, self.name)
        except AttributeError as ex:
            pass


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
        for key, value in optional_keyword_arguments.items():
            if self.hasoption(key):
                setattr(self, key, value)

    def hasoption(self, name):
        the_type = type(self)
        return hasattr(the_type, name)

    @late
    def _local_options(self):
        return {}

    def iter_options(self):
        cls = type(self)
        for x in dir(cls):
            if x.startswith('_'):
                continue
            value = getattr(cls, x)
            if isinstance(value, option):
                yield value
