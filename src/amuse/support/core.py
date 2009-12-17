"""

"""
import types 
import collections

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
        value = self.initializer(instance)
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
        return isinstance(x,types.StringType)
        
    def isnumber(self, x):
        return isinstance(x,types.IntType) or isinstance(x,types.FloatType) 
        
        
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
            
class OrderedSet(collections.MutableSet):
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
            self.previous.next = self.next
            self.next.previous = self.previous
            
    def __init__(self, iterable=None):
        self.end = Node(None, None, None)
        self.end.prev = end
        self.end.next = end
        self.map = {}
        if iterable is not None:
            self |= iterable

    def __len__(self):
        return len(self.map)

    def __contains__(self, key):
        return key in self.map

    def add(self, key):
        if key not in self.map:
            self.map[key] = self.Node(key, end.previous, end.next)

    def discard(self, key):
        if key in self.map:        
            current = self.map.pop(key)
            current.discard()

    def __iter__(self):
        end = self.end
        current = end.next
        while current is not end:
            yield current.key
            current = current.next

    def __reversed__(self):
        end = self.end
        current = end.previous
        while current is not end:
            yield current.key
            current = current.previous

    def pop(self, last=True):
        if not self:
            raise KeyError('set is empty')
        key = next(reversed(self)) if last else next(iter(self))
        self.discard(key)
        return key

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
