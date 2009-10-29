from docutils import core
import docutils.nodes as nodes

class LiteratureRefs(object):

    literature_list = []#Global containing all found lit refs (footnotes)

    def __init__(self):
        self.fill_literature_list()
 
    @classmethod
    def print_refs(self, modus):
        print "Your are currently using the following codes, which contain literature references"
        for s in self.literature_list:
            print '\n\t"%s"' % s[0]
            for key_and_note in s[1]:
                print '\t\t%s: %s' % (key_and_note[0],key_and_note[1])
    
    @classmethod
    def fill_literature_list(self):
        """filter the refs form the docstring, if no refs there is no append"""
        docstringin = self.__doc__
        objectname = self.__name__
        Mydoctree  = core.publish_doctree(source = docstringin)
        footnote_list = []
        for ikey, ival in Mydoctree.ids.iteritems():
            if isinstance(ival,nodes.footnote):
                footnote_list.append([ikey, ival.rawsource])
        
        filled = bool(footnote_list)
        
        if filled:
            self.literature_list.append([objectname, footnote_list])
