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
try:
    from amuse.version import version as amuse_version
except ImportError:
    amuse_version = "unknown version"

import atexit
import sys
import traceback
import importlib

ClassWithLiteratureReferences = namedtuple(\
    "ClassWithLiteratureReferences", 
    "name_of_class_with_refs literature_references_of_class"
)
LiteratureReference = namedtuple(
    "LiteratureReference",
    "id footnote"
)

class TrackLiteratureReferences:
    """
        .. [#] [2018araa.book.....P] Portegies Zwart, S. & McMillan, S.L.W., 2018
        .. [#] [2013CoPhC.183..456P] ** Portegies Zwart, S. et al., 2013
        .. [#] [2013A&A...557A..84P] ** Pelupessy, F. I. et al., 2013
        .. [#] [2009NewA...14..369P] Portegies Zwart, S. et al., 2009
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
                prefix = """

Thank you for using AMUSE!
In this session you have used the modules below. Please cite any relevant articles:

"""
                print(prefix + self.all_literature_references_string())
        
    
    def get_literature_list_of_class(self, cls):
        """
        filter the refs from the docstring, if there are no refs nothing is appended
        """

        result = []
        for current_class in cls.__mro__:
            docstring_in = current_class.__doc__
            if docstring_in:
                if hasattr(current_class, "version"):
                    version = current_class.version()
                else:
                    version = amuse_version
                name = current_class.__name__
                if name.endswith("Interface"):
                    name = "AMUSE-" + name[:-9]
                objectname = "{name} ({version})".format(
                    name=name,
                    version=version,
                )
                doctree = core.publish_doctree(source = docstring_in)
                ref_keys = list(doctree.ids.keys())
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
        
        lines.append('\n\t"AMUSE (%s)"' % amuse_version)
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
    def version(cls):
        try:
            version = importlib.import_module(
                '..version',
                cls.__module__
            ).version
        except (ImportError, ValueError):
            try:
                from amuse.version import version
            except ImportError:
                version = "unknown"
        return version
 
    @classmethod
    def print_literature_references(cls):
        print("You are currently using the following codes, which contain literature references")
        print(TrackLiteratureReferences.default().all_literature_references_string())
 
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
    return list(map(try_int, re.findall(r'(\d+|\D+)', s)))

def natsort(seq):
    "In-place natural string sort."
    seq.sort(key=natsort_key)
