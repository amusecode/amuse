"""

"""
import types 

class late(object):    
    """
    An attribute that is set at first access. 
    
    The value of the attribute will be determined from the *initializer*
    method. The name of the attribute is the same as the name of
    the *initializer* method.
    
    A late attribute is comparible with attributes set in the
    *__init__* method. Except the value of the late attribute is
    determined when first accessed and not when the class is
    instantiated.
    """ 
    def __init__(self, initializer):
        """
        Called when used as a decorator.
        
        Typical use to define a managed attribute x:
        
        >>> class C(object):
        ...    @late
        ...    def x(self):
        ...        return "i'm late!"
        ...
        >>> c = C()
        >>> print c.x
        i'm late!
        >>> c.x = "overridden"
        >>> print c.x
        overridden
        
        :argument initializer: function to determine the initial value of the property
        :returns: a descriptor to determine and set the value on first access
        """
        self.initializer = initializer
        self.__doc__ = self.initializer.__doc__
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        value = self.initializer(instance)
        setattr(instance,self.initializer.__name__, value)
        return value
        
        
class print_out(object):
    def __init__(self):
        self.parts = []
        self._indent = 0
        self.number_of_characters_on_current_line = 0
    
    def __add__(self, x):
        if self.isstring(x):
            self.parts.append(x)
            self.number_of_characters_on_current_line += len(x)
        elif self.isnumber(x):
            self.parts.append(str(x))
            self.number_of_characters_on_current_line += len(str(x))
        return self
        
    def n(self):
        if not self.parts:
            return self
        if self.parts[-1] == '\n':
            return self
        self.lf()
        return self
        
    def indent(self):
        self._indent += 1
        return self
        
    def dedent(self):
        self._indent -= 1
        return self
        
    def lf(self):
        self.parts.append('\n')
        self.number_of_characters_on_current_line = 0
        self.do_indent()
        return self

    def do_indent(self):
        for ignore in range(self._indent):
            self.parts.append(self.indent_characters())
            self.number_of_characters_on_current_line += len(self.indent_characters())

    def indent_characters(self):
        return '  '
        
    def __str__(self):
        return ''.join(self.parts)
        
    @property
    def string(self):
        return str(self)
        
    def isstring(self, x):
        return isinstance(x,types.StringType)
        
    def isnumber(self, x):
        return isinstance(x,types.IntType) or isinstance(x,types.FloatType) 
        
        
class OrderedDictionary(object):
    def __init__(self):
        self.mapping = {}
        self.orderedKeys = []
        
    def __setitem__(self, key, value):
        if key in self.mapping:
            raise Exception("key " + key + " already in the dictionary")
        self.orderedKeys.append(key)
        self.mapping[key] = value
        
    def __getitem__(self, key):
        return self.mapping[key]
        
    def  __contains__(self, key):
        return key in self.mapping
        
    def  __iter__(self):
        return self.values()
        
    def  __len__(self):
        return len(self.orderedKeys)
        
    def __str__(self):
        result = 'OrderedDictionary({'
        elements = []
        for x in self.keys():
            elements.append(str(x) + ':' + str(self[x]))
        result += ','.join(elements)
        result += '})'
        return result

    def __getattr__(self, key):
        return self.mapping[key]
        
    def keys(self):
        return iter(self.orderedKeys)
        
    def values(self):
        for x in self.orderedKeys:
            yield self.mapping[x]
