from docutils import nodes
from docutils.parsers.rst import directives
from docutils.parsers.rst import Directive

import textwrap
import sys

from sphinx import addnodes

from amuse.rfi.core import is_mpd_running
from sphinx.ext.autodoc import AttributeDocumenter, ModuleLevelDocumenter
from sphinx.util import force_decode
from sphinx.util.docstrings import prepare_docstring

class ParametersAttributeDocumenter(AttributeDocumenter):
    """
    Specialized Documenter subclass for parameters attribute
    of interfaces
    """
    objtype = 'parametersattribute'
    directivetype = 'attribute'
    member_order = 60
    

    # must be higher than AttributeDocumenter
    priority = 11

    @classmethod
    def can_document_member(cls, member, membername, isattr, parent):
        return False

    def add_content(self, more_content, no_docstring=False):
        if not is_mpd_running():
            return
            
        try:
            
            cls = self.object
            instance = cls(must_start_worker = False, must_handle_state = False)
            try:
                #instance.initialize_code()
                parameter_documentation = self.get_sphinx_doc_for_parameters(instance.parameters)
            finally:
                instance.stop()
                
        except Exception as ex:
            print ex
            return
            
        if self.analyzer:
            # prevent encoding errors when the file name is non-ASCII
            filename = unicode(self.analyzer.srcname,
                               sys.getfilesystemencoding(), 'replace')
            sourcename = u'%s:docstring of %s' % (filename, self.fullname)
        else:
            sourcename = u'docstring of %s' % self.fullname
            
        encoding = self.analyzer and self.analyzer.encoding
        lines = prepare_docstring(force_decode(parameter_documentation, encoding))
        
        for i, line in enumerate(self.process_doc([lines,])):
                self.add_line(line, sourcename, i)
    
    def get_sphinx_doc_for_parameters(self, parameters):
        lines = []
        
        for parameter_definition in parameters._definitions:
            lines.append('.. py:attribute:: '+ self.objpath[-1] +'.' + parameter_definition.name)
            lines.append('')
            dedented = textwrap.dedent(parameter_definition.description)
            for x in dedented.splitlines():
                lines.append('    ' + x)
            try:
                lines.append('    ' + "(default value:" + str(parameters.get_default_value_for(parameter_definition.name)) + ")")
            except Exception as ex:
                lines.append('    ' + "(no default value)")
                
            lines.append('')
        
        lines.append('')
        
        return '\n'.join(lines)
        
    def import_object(self):
        """
        Import the object given by *self.modname* and *self.objpath* and sets
        it as *self.object*.

        Returns True if successful, False if an error occurred.
        """
        self._datadescriptor = False
        try:
            __import__(self.modname)
            parent = None
            obj = self.module = sys.modules[self.modname]
            for part in self.objpath[:-1]:
                parent = obj
                obj = self.get_attr(obj, part)
                self.object_name = part
            self.parent = parent
            self.object = obj
            return True
        # this used to only catch SyntaxError, ImportError and AttributeError,
        # but importing modules with side effects can raise all kinds of errors
        except Exception, err:
            if self.env.app and not self.env.app.quiet:
                self.env.app.info(traceback.format_exc().rstrip())
            self.directive.warn(
                'autodoc can\'t import/find %s %r, it reported error: '
                '"%s", please check your spelling and sys.path' %
                (self.objtype, str(self.fullname), err))
            self.env.note_reread()
            return False

def setup(app):
    app.add_autodocumenter(ParametersAttributeDocumenter)
