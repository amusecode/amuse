from docutils import nodes
from docutils.parsers.rst import directives
from docutils.parsers.rst import Directive

from amuse.support import io
import textwrap
from sphinx import addnodes


class IoOptions(Directive):

    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = False
    option_spec = {}
    has_content = False

    def run(self):
        options = io.get_options_for_format(self.arguments[0])
        field_list_node = nodes.definition_list()
        for name, description, value in options:
            item = nodes.definition_list_item()
            item.append(nodes.term(name + ' ',name+ ' '))
            item.append(nodes.definition('', nodes.paragraph('', description)))
            field_list_node.append(item)
        return [field_list_node]
        
def setup(app):
    directives.register_directive("iooptions", IoOptions)

