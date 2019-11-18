"""

"""
import types 
import collections.abc
import re

def compare_version_strings(version1, version2):
    def normalize(v):
        return [int(x if x.isdigit() else 0) for x in re.sub(r'(\.0+)*$','', v).split(".")]
    version1 = normalize(version1)
    version2 = normalize(version2)
    
    return (version1 > version2) - (version1 < version2)

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
    def __init__(self, initializer):
        self.initializer = initializer
        self.__doc__ = self.initializer.__doc__
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        
        try:
            value = self.initializer(instance)
        except Exception as ex:
            raise AttributeError(ex)
            
        setattr(instance,self.initializer.__name__, value)
        
        return value
    
        
        
class print_out(object):
    """
    Efficient way to contruct large strings.
    
    Strings are build up out of parts. Objects of this
    class store these parts while building the string.
    Only on request the parts are concatenated into
    a large string.
    
    Strings and numbers can be added to the print_out.
    For other objects str(object) is called before
    adding it to the print_out.
    
    >>> p = print_out()
    >>> p + "number of counts : " + 10    #doctest: +ELLIPSIS
    <amuse.support.core.print_out object at 0x...>
    >>> print p.string
    number of counts : 10
    
    All methods return the print_out instance, so that
    calls can be chained.
    """
    
    def __init__(self):
        self.parts = []
        self._indent = 0
        self.number_of_characters_on_current_line = 0
    
    def __add__(self, x):
        """Add a new part to the print_out.
        """    
    
        if self.isstring(x):
            self.parts.append(x)
            self.number_of_characters_on_current_line += len(x)
        elif self.isnumber(x):
            self.parts.append(str(x))
            self.number_of_characters_on_current_line += len(str(x))
        elif isinstance(x, print_out):
            self.parts.extend(x.parts)
        else:
            part = str(x)
            self.parts.append(part)
            self.number_of_characters_on_current_line += len(part) 
        return self
        
    def n(self):
        """Start a new-line, if the current line is not-empty.
        
        >>> p = print_out()
        >>> for i in range(3):
        ...     p.n() + i #doctest: +ELLIPSIS
        ...
        <amuse.support.core.print_out object at 0x...>
        <amuse.support.core.print_out object at 0x...>
        <amuse.support.core.print_out object at 0x...>
        >>> print p.string
        0
        1
        2
        """
        if not self.parts:
            return self
        if self.parts[-1] == '\n':
            return self
        self.lf()
        return self
        
    def indent(self):
        """Increase the indent. The next and following lines
        will start indented.
        
        >>> p = print_out()
        >>> p + "01" #doctest: +ELLIPSIS
        <amuse.support.core.print_out object at 0x...>
        >>> p.indent().lf() + "2" #doctest: +ELLIPSIS
        <amuse.support.core.print_out object at 0x...>
        >>> p.lf() + "3" #doctest: +ELLIPSIS
        <amuse.support.core.print_out object at 0x...>
        >>> print p.string
        01
          2
          3
        """
        self._indent += 1
        return self
        
    def dedent(self):
        """Decrease the indent. The next line will start dedented.
        
        >>> p = print_out()
        >>> p + "01" #doctest: +ELLIPSIS
        <amuse.support.core.print_out object at 0x...>
        >>> p.indent().lf() + "2" #doctest: +ELLIPSIS
        <amuse.support.core.print_out object at 0x...>
        >>> p.dedent().lf() + "01" #doctest: +ELLIPSIS
        <amuse.support.core.print_out object at 0x...>
        >>> print p.string
        01
          2
        01
        """
        self._indent -= 1
        return self
        
    def lf(self):
        """Start a new-line"""
        self.parts.append('\n')
        self.number_of_characters_on_current_line = 0
        self.do_indent()
        return self
        
    def lf_noindent(self):
        """Start a new-line"""
        self.parts.append('\n')
        self.number_of_characters_on_current_line = 0
        return self

    def do_indent(self):
        for ignore in range(self._indent):
            self.parts.append(self.indent_characters())
            self.number_of_characters_on_current_line += len(self.indent_characters())

    def indent_characters(self):
        """ The indent characters, by default 2 spaces.
        
        Override this method to change the indent characters.
        """
        return '  '
        
    def __str__(self):
        return ''.join(self.parts)
        
    @property
    def string(self):
        """String version of the print_out.
        """
        return str(self)
        
    def isstring(self, x):
        return isinstance(x,bytes)
        
    def isnumber(self, x):
        return isinstance(x,int) or isinstance(x,float) 
        
        
class OrderedDictionary(object):
    """A dictionary that keeps the keys in the dictionary in order.
    
    Ordered dictionaries are just like regular dictionaries but they remember the
    order that items were inserted.  When iterating over an ordered dictionary,
    the values are returned in the order their keys were first added.

    >>> d = OrderedDictionary()
    >>> d["first"] = 0
    >>> d["second"] = 1
    >>> d["third"] = 2
    >>> [x for x in d]
    [0, 1, 2]
    """
    def __init__(self):
        self.mapping = {}
        self.orderedKeys = []
        
    def __setitem__(self, key, value):
        if key in self.mapping:
            self.mapping[key] = value
            return
        self.orderedKeys.append(key)
        self.mapping[key] = value
        
    def __getitem__(self, key):
        return self.mapping[key]
        
    def  __contains__(self, key):
        return key in self.mapping
        
    def  __iter__(self):
        return iter(self.values())
        
    def  __len__(self):
        return len(self.orderedKeys)
        
    def __str__(self):
        result = 'OrderedDictionary({'
        elements = []
        for x in self.keys():
            elements.append(repr(x) + ':' + repr(self[x]))
        result += ', '.join(elements)
        result += '})'
        return result
        
    def __repr__(self):
        return str(self)
        
    def iterkeys(self):
        return iter(self.orderedKeys)
        
    def itervalues(self):
        for x in iter(self.orderedKeys):
            yield self.mapping[x]
    
    def iteritems(self):
        for x in self.orderedKeys:
            yield x, self.mapping[x]
    
    def keys(self):
        return list(self.orderedKeys)
        
    def pop(self, key):
        index = self.orderedKeys.index(key)
        del self.orderedKeys[index]
        return self.mapping.pop(key)
    
    def values(self):
        return [self.mapping[x] for x in self.orderedKeys]
        
    def items(self):
        return [(x,self.mapping[x]) for x in self.orderedKeys]
        
    def copy(self):
        result = OrderedDictionary()
        result.mapping = self.mapping.copy()
        result.orderedKeys = list(self.orderedKeys)
        return result
            

        
class OrderedMultiDictionary(object):
    """A dictionary that keeps the keys in the dictionary in order and can store
    multiple items per key
    
    Ordered multi dictionaries remember the order that items were inserted
    and can store multple values per key.  When iterating over an ordered dictionary,
    the values are returned in the order their keys were first added.

    >>> d = OrderedMultiDictionary()
    >>> d["first"] = 0
    >>> d["second"] = 1
    >>> d["first"] = 2
    >>> [x for x in d]
    [0, 1, 2]
    >>> print d["first"]
    [0, 2]
    >>> print d["second"]
    [1]
    """
    
    def __init__(self):
        self.mapping = {}
        self.orderedKeys = []
        
    def __setitem__(self, key, value):
        if not key in self.mapping:
            self.mapping[key] = []
        self.mapping[key].append(value)
        self.orderedKeys.append((key, len(self.mapping[key]) - 1,))
        
    def __getitem__(self, key):
        return self.mapping[key]
        
    def  __contains__(self, key):
        return key in self.mapping
        
    def  __iter__(self):
        return list(self.values())
        
    def  __len__(self):
        return len(self.orderedKeys)
        
    def __str__(self):
        result = 'OrderedDictionary({'
        elements = []
        for x, index in self.orderedKeys:
            elements.append(repr(x) + ':' + repr(self[x][index]))
        result += ', '.join(elements)
        result += '})'
        return result
        
    def __repr__(self):
        return str(self)
        
    def __getattr__(self, key):
        return self.mapping[key]
        
    def keys(self):
        return [x for x,index in self.orderedKeys ]
                
    def values(self):
        for x, index in self.orderedKeys:
            yield self.mapping[x][index]
            
class CompositeDictionary(object):
    """A dictionary that defers to other dictionaries when an item is
    not found.
    
    Composite dictionaries are just like regular dictionaries but they
    get items from their parent dictionarkies when they do not contain
    the items.  
    
    >>> p = {'a':1, 'b':2}
    >>> d = CompositeDictionary(p)
    >>> d['a']
    1
    >>> p['a'] = 3
    >>> d['a']
    3
    >>> d['c'] = 2
    >>> 'c' in p
    False
    >>> 'b' in d
    True
    
    """
    def __init__(self, *parents):
        self.parents = parents
        self.mapping = {}
        
    def __setitem__(self, key, value):
        self.mapping[key] = value
        
    def __getitem__(self, key):
        if key in self.mapping:
            return self.mapping[key]
        for parent in self.parents:
            if key in parent:
                return parent[key]
        raise KeyError(key)
        
    def  __contains__(self, key):
        if key in self.mapping:
            return True
        for parent in self.parents:
            if key in parent:
                return True
        return False
        
    def  __iter__(self):
        return list(self.keys())
        
    def  __len__(self):
        return len(list(self.keys()))
        
    def __str__(self):
        result = 'CompositeDictionary({'
        elements = []
        for x in list(self.keys()):
            elements.append(str(x) + ':' + str(self[x]) )
        result += ','.join(elements)
        result += '})'
        return result

        
    def keys(self):
        keys = set(self.mapping.keys())
        
        for parent in self.parents:
            keys |= set(parent.keys())
        
        return iter(keys)
        
    def values(self):
        for x in list(self.keys()):
            yield self[x]
            
    def copy(self):
        result = type(self)(*self.parents)
        result.mapping = self.mapping.copy()
        return result
            
class OrderedSet(collections.abc.MutableSet):
    class Node(object):
        __slots__ = ['key', 'next', 'previous']
        
        def __init__(self, key, next = None, previous = None):
            self.key = key
            if next is None:
                next = self
            if previous is None:
                previous = self
                
            self.next = next
            self.previous = previous

            self.link()
            
        def link(self):
            self.next.previous = self
            self.previous.next = self
            
        def discard(self):
            self.previous.next = self.__next__
            self.next.previous = self.previous
            
    def __init__(self, iterable=None):
        self.end = self.Node(None, None, None)
        self.end.previous = self.end
        self.end.next = self.end
        self.map = {}
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def __iter__(self):
        end = self.end
        current = end.__next__
        while current is not end:
            yield current.key
            current = current.__next__

    def __reversed__(self):
        end = self.end
        current = end.previous
        while current is not end:
            yield current.key
            current = current.previous

    def __repr__(self):
        if not self:
            return '%s()' % (self.__class__.__name__,)
        return '%s(%r)' % (self.__class__.__name__, list(self))

    def __eq__(self, other):
        if isinstance(other, OrderedSet):
            return len(self) == len(other) and list(self) == list(other)
        return not self.isdisjoint(other)

    def __del__(self):
        self.clear()
    def add(self, key):
        if key not in self.map:
            self.map[key] = self.Node(key, self.end.previous, self.end.__next__)

    def discard(self, key):
        if key in self.map:        
            current = self.map.pop(key)
            current.discard()

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = next(reversed(self)) if last else next(iter(self))
        self.discard(key)
        return key

def memoize(f):
    def memof(*arg):
        try:
            return memof.d[arg]
        except:
            result=f(*arg)
            if len(memof.d)>5000:
                return result
            memof.d[arg]=result
            return result
    memof.d={}
    return memof


class MultitonMetaClass(type):
    def __new__(mcs, name, bases, dict):
        dict['__INSTANCES__'] = {}
        return type.__new__(mcs, name, bases, dict)
        
    def __call__(mcls, *arguments):
        if arguments in mcls.__INSTANCES__:
            return  mcls.__INSTANCES__[arguments]
        else:
            instance = type.__call__(mcls, *arguments)
            mcls.__INSTANCES__[arguments] = instance
            return instance

        
