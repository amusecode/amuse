from docutils import core
import docutils.nodes as nodes
from collections import namedtuple

class LiteratureRefs(object):

    literature_list = []#Global containing all found lit refs (footnotes)

    def __init__(self):
        self.fill_literature_list()
 
    @classmethod
    def print_refs(cls):
        print "You are currently using the following codes, which contain literature references"
        print cls.all_literature_references_string()
 
    @classmethod
    def export2html(cls):
        pass

    @classmethod
    def export2bibtex(cls):
        pass
   
    @classmethod
    def names_of_classes_with_references(cls):
        return [x[0] for x in cls.literature_list]
 
    @classmethod
    def all_literature_references_string(cls):
        lines = []
        for s in cls.literature_list:
            lines.append('\n\t"%s"' % s.name_of_class_with_refs)
            for literature_reference_of_class_item in s.literature_references_of_class:
                lines.append('\t\t%s' % (literature_reference_of_class_item.footnote))
        return "\n".join(lines)
        
    @classmethod
    def fill_literature_list(cls):
        """filter the refs form the docstring, if no refs there is no append"""
        classes_with_references_item = namedtuple("cwrl", 
                                                  "name_of_class_with_refs literature_references_of_class"
                                                 )
        literature_references_of_class_item = namedtuple("lrocl", 
                                                         "id footnote"
                                                         )

        docstring_in = cls.__doc__
        objectname = cls.__name__
        Mydoctree  = core.publish_doctree(source = docstring_in)

        literature_references_of_class = []
        
        for ikey, ival in Mydoctree.ids.iteritems():
            if isinstance(ival,nodes.footnote):
                literature_references_of_class.append(literature_references_of_class_item(ikey, ival.rawsource))
        
        filled = bool(literature_references_of_class)
        
        if filled:
            cls.literature_list.append(classes_with_references_item(objectname, 
                                                                    literature_references_of_class
                                                                    )
                                       )
