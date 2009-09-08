"""
usefull objects
"""
import types 

class late(object):    
    """  A late initialized property.
    The property is initialized at first access. 
    """ 
    def __init__(self, initializer):
        """late(initializer) -> late attribute
        
        initializer is a function to be used for getting the initial attribute value
        
        Typical use to define a managed attribute x:
        class C(object):
            @late
            def x(self):
                return "i'm late!"
        """
        self.initializer = initializer
    def __get__(self, instance, owner):
        """Initialize the value of the property using the *initializer* function
        """
        if instance is None:
            return self
        value = self.initializer(instance)
        setattr(instance,self.initializer.__name__,value)
        return value
        
        
class print_out(object):
    def __init__(self):
        self.parts = []
        self._indent = 0
    
    def __add__(self, x):
        if self.isstring(x):
            self.parts.append(x)
        elif self.isnumber(x):
            self.parts.append(str(x))
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
        self.do_indent()
        return self

    def do_indent(self):
        for i in range(self._indent):
            self.parts.append(self.indent_characters())

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
        
    def keys(self):
        return iter(self.orderedKeys)
        
    def values(self):
        for x in self.orderedKeys:
            yield self.mapping[x]
