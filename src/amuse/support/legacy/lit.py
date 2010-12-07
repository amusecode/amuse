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

class TracLiteratureReferences(object):
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
        
    def register_class(self, cls):
        self.registered_classes.add(cls)
    
    def exception_hook(self, *arguments):
        #print "exception", arguments, self.original_excepthook
        self.must_show_literature_references_atexit = False
        lines = traceback.format_exception(type, value, tb)
        print string.join(lines)
   
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
                Mydoctree  = core.publish_doctree(source = docstring_in)
                ref_keys = Mydoctree.ids.keys()
                natsort(ref_keys)
                ref_values = [Mydoctree.ids[key] for key in ref_keys]
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
        return "\n".join(lines)
        
    def names_of_classes_with_references(self):
        return [x.name_of_class_with_refs for x in self.get_literature_list()]

class LiteratureRefs(object):


    def __init__(self):
        self.register_use()
 
    @classmethod
    def print_refs(cls):
        print "You are currently using the following codes, which contain literature references"
        print TracLiteratureReferences.default().all_literature_references_string()
 
    @classmethod
    def export2html(cls):
        pass

    @classmethod
    def export2bibtex(cls):
        pass
   
    @classmethod
    def names_of_classes_with_references(cls):
        return TracLiteratureReferences.default().names_of_classes_with_references()
 
    @classmethod
    def all_literature_references_string(cls):
        return TracLiteratureReferences.default().all_literature_references_string()
        
    @classmethod
    def register_use(cls):
        TracLiteratureReferences.default().register_class(cls)
        
   
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
