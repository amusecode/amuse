try:
    from docutils import core
except ValueError:
    import os
    import locale
    os.environ['LC_CTYPE'] = 'C'
    os.environ['LANG'] = 'C'
    from docutils import core

import docutils.nodes as nodes
from collections import namedtuple
from amuse.support import exceptions
import warnings

import atexit
import sys
import traceback

ClassWithLiteratureReferences = namedtuple(\
    "ClassWithLiteratureReferences", 
    "name_of_class_with_refs literature_references_of_class"
)
LiteratureReference = namedtuple(
    "LiteratureReference", 
    "id footnote"
)

class TrackLiteratureReferences(object):
    """
        .. [#] ** Portegies Zwart, S. et al., 2013, Multi-physics Simulations Using a Hierarchical Interchangeable Software Interface, Computer Physics Communications 183, 456-468 [2013CoPhC.183..456P]
        .. [#] ** Pelupessy, F. I. et al., 2013, The Astrophysical Multipurpose Software Environment, Astronomy and Astrophysics 557, 84 [2013A&A...557A..84P]
        .. [#] Portegies Zwart, S. et al., 2009, A multiphysics and multiscale software environment for modeling astrophysical systems, *New Astronomy*, **Volume 14**, **Issue 4**, 369-378 [2009NewA...14..369P]

    """
    INSTANCE = None
    
    def __init__(self):
        self.registered_classes = set([])
        self.must_show_literature_references_atexit = True
        self.original_excepthook = None
    
    @classmethod
    def default(cls):
        if cls.INSTANCE is None:
            cls.INSTANCE = cls()
            cls.INSTANCE.register()
        return cls.INSTANCE
    
    def register(self):
        self.original_excepthook = sys.excepthook
        sys.excepthook = self.exception_hook
        atexit.register(self.atexit_hook)
        
    @classmethod
    def suppress_output(cls):
        cls.default().must_show_literature_references_atexit = False
        
    def register_class(self, cls):
        self.registered_classes.add(cls)
    
    def exception_hook(self, *arguments):
        #print "exception", arguments, self.original_excepthook
        self.must_show_literature_references_atexit = False
        lines = traceback.format_exception(*arguments)
        #print ''.join(lines)
   
        self.original_excepthook(*arguments)
        
    def atexit_hook(self):
        if not self.original_excepthook is None:
            sys.excepthook = self.original_excepthook
        self.original_excepthook = None
        
        if self.must_show_literature_references_atexit:
            string = self.all_literature_references_string()
            if string:
                prefix = "\n\nYou have used the following codes, which contain literature references:\n"
                warnings.warn(prefix + self.all_literature_references_string(), exceptions.AmuseWarning)
        
    
    def get_literature_list_of_class(self, cls):
        """filter the refs form the docstring, if no refs there is no append"""

        result = []
        for current_class in cls.__mro__:
            docstring_in = current_class.__doc__
            if docstring_in:
                objectname = current_class.__name__
                doctree  = core.publish_doctree(source = docstring_in)
                ref_keys = doctree.ids.keys()
                natsort(ref_keys)
                ref_values = [doctree.ids[key] for key in ref_keys]
                literature_references_of_class = []
                for ikey, ival in zip(ref_keys, ref_values):
                    if isinstance(ival,nodes.footnote):
                        literature_references_of_class.append(
                            LiteratureReference(ikey, ival.rawsource)
                        )
                filled = len(literature_references_of_class) > 0
                if filled:
                    result.append(
                        ClassWithLiteratureReferences(
                            objectname, 
                            literature_references_of_class
                        )
                    )
        return result
    
    def get_literature_list(self):
        result = []
        for x in self.registered_classes:
            result.extend(self.get_literature_list_of_class(x))
        return result
        
   
    def all_literature_references_string(self):
        lines = []
        for s in self.get_literature_list():
            lines.append('\n\t"%s"' % s.name_of_class_with_refs)
            for literature_reference_of_class_item in s.literature_references_of_class:
                lines.append('\t\t%s' % (literature_reference_of_class_item.footnote))
        
        lines.append('\n\t"AMUSE"')
        amuse_list = self.get_literature_list_of_class(type(self))
        for x in amuse_list:
            for literature_reference_of_class_item in x.literature_references_of_class:
                lines.append('\t\t%s' % (literature_reference_of_class_item.footnote))
            
        return "\n".join(lines)
        
    def names_of_classes_with_references(self):
        return [x.name_of_class_with_refs for x in self.get_literature_list()]


def literature_references():
    return TrackLiteratureReferences.default().all_literature_references_string()
    
class LiteratureReferencesMixIn(object):


    def __init__(self):
        self.register_use()
 
    @classmethod
    def print_literature_references(cls):
        print "You are currently using the following codes, which contain literature references"
        print TrackLiteratureReferences.default().all_literature_references_string()
 
    @classmethod
    def export2html(cls):
        pass

    @classmethod
    def export2bibtex(cls):
        pass
   
    @classmethod
    def names_of_classes_with_references(cls):
        return TrackLiteratureReferences.default().names_of_classes_with_references()
 
    @classmethod
    def all_literature_references_string(cls):
        return TrackLiteratureReferences.default().all_literature_references_string()
        
    @classmethod
    def register_use(cls):
        TrackLiteratureReferences.default().register_class(cls)
        
   
# ------------------------------------------------------------------------------
# from natsort.py: Natural string sorting by Seo Sanghyeon and Connelly Barnes.
# ------------------------------------------------------------------------------

def try_int(s):
    "Convert to integer if possible."
    try: return int(s)
    except: return s

def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(a), natsort_key(b))

def natsort(seq, cmp=natcmp):
    "In-place natural string sort."
    seq.sort(cmp)
