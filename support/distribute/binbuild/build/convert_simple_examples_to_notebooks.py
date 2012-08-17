import glob
import io
import os.path
import ast
from IPython.nbformat import current as notebook_format

def iter_python_files(path):
    globpath = os.path.join(path, '[a-zA-Z]*.py')
    return glob.iglob(globpath)

def convert_file(path, outputdir):
    filename = os.path.basename(path)
    notebookname, _ = os.path.splitext(filename)
    with open(path, 'r') as stream:
        sourcestring = stream.read()
        
    notebook = new_notebook_from_string(notebookname, filename, sourcestring)
    
    with open(os.path.join(outputdir, notebookname+'.ipynb'),'w') as stream:
        notebook_format.write(notebook, stream, u'json')

def new_notebook_from_string(notebookname, filename, sourcestring):
    root = ast.parse(sourcestring, filename=filename, mode='exec')
    x = DetermineBlocks()
    
    for child in ast.iter_child_nodes(root):
        print child.lineno, child
        x.visit(child)
    x.end()
    sourcelines = sourcestring.splitlines()
    cells = []
    for block in x.blocks:
        print block
        blocklines = sourcelines[block[1]-1:block[2]]
        blocksrc = '\n'.join(blocklines)
        if len(blocksrc) > 0:
            cell = notebook_format.new_code_cell(input=blocksrc)
            cells.append(cell)
            
    ws = notebook_format.new_worksheet(cells=cells)
    result = notebook_format.new_notebook(worksheets=[ws])
    result.metadata.name = notebookname
    return result
        
class DetermineBlocks(ast.NodeVisitor):
    def __init__(self):
        ast.NodeVisitor.__init__(self)
        self.kind = 'Start'
        self.start_lineno = 1
        self.blocks = []
        
    def generic_visit(self, node):
        pass
        
    def visit_Import(self,node):
        if not self.kind == 'Import':
            self.blocks.append( (self.kind, self.start_lineno, node.lineno - 1) )
            self.kind = 'Import'
            self.start_lineno = node.lineno
        
    def visit_ImportFrom(self,node):
        if not self.kind == 'Import':
            self.blocks.append( (self.kind, self.start_lineno, node.lineno - 1) )
            self.kind = 'Import'
            self.start_lineno = node.lineno
    
    def visit_FunctionDef(self,node):
        self.blocks.append( (self.kind, self.start_lineno, node.lineno - 1) )
        self.kind = 'FunctionDef'
        self.start_lineno = node.lineno
        
    def visit_If(self,node):
        self.blocks.append( (self.kind, self.start_lineno, node.lineno - 1) )
        self.kind = 'If'
        self.start_lineno = node.lineno
        
    def visit_Assign(self,node):
        if not self.kind == 'Assign':
            self.blocks.append( (self.kind, self.start_lineno, node.lineno - 1) )
            self.kind = 'Assign'
            self.start_lineno = node.lineno
    
    def visit_DDExpr(self,node):
        print node.col_offset
        if not self.kind == 'Expr':
            self.blocks.append( (self.kind, self.start_lineno, node.lineno - 1) )
            self.kind = 'Expr'
            self.start_lineno = node.lineno
            
    def end(self):
        self.blocks.append( (self.kind, self.start_lineno, None) )
        
if __name__ == '__main__':
    for x in iter_python_files('examples/simple'):
        convert_file(x, 'output')
