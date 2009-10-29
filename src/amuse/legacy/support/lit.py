from docutils import core
import docutils.nodes as nodes

class LiteratureRefs(object):

    literature_list = []#Global containing all found lit refs (footnotes)

    def __init__(self):
        self.fill_literature_list()
 
    @classmethod
    def print_refs(cls, modus):
        print "You are currently using the following codes, which contain literature references"
        print cls.all_literature_references_string()
 
    @classmethod
    def names_of_classes_with_references(cls):
        return [x[0] for x in cls.literature_list]

    @classmethod
    def all_literature_references_string(cls):
        lines = []
        for s in cls.literature_list:
            lines.append('\n\t"%s"' % s[0])
            for key_and_note in s[1]:
                lines.append('\t\t%s' % (key_and_note[1]))
        return "\n".join(lines)
        
    @classmethod
    def fill_literature_list(cls):
        """filter the refs form the docstring, if no refs there is no append"""
        docstringin = cls.__doc__
        objectname = cls.__name__
        Mydoctree  = core.publish_doctree(source = docstringin)
        footnote_list = []
        for ikey, ival in Mydoctree.ids.iteritems():
            if isinstance(ival,nodes.footnote):
                footnote_list.append([ikey, ival.rawsource])
        
        filled = bool(footnote_list)
        
        if filled:
            cls.literature_list.append([objectname, footnote_list])
