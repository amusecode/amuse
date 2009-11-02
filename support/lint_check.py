import subprocess
import xml.dom
import xml.dom.minidom

class InterfaceToPyLint(object):
    NAME_OF_THE_COMMAND = 'pylint' 
    
    def is_available(self):
        process = subprocess.Popen(
            ['which', self.NAME_OF_THE_COMMAND],
            stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE
        )
        process.communicate()
        return process.returncode == 0
        
    def run_onfile(self,path):
        process = subprocess.Popen(
            [self.NAME_OF_THE_COMMAND, '-f', 'parseable', path], 
            stdout = subprocess.PIPE, 
            stderr = subprocess.PIPE
        )
        stdout_string, stderr_string  = process.communicate()
        return stdout_string
        
class MakeHTMLDomFromFile(object):
    
    def __init__(self, path):
        self.path = path
        self.filecontents = self.get_filecontents()
        
    def get_filecontents(self):
        with open(self.path, "r") as file:
            return file.read()
            
    def start(self):
        self.document = self.new_document()
        
        self.head = self.document.createElement("head")
        self.document.documentElement.appendChild(self.head)
        
        table = self.new_table(1,2)
        
        td1 = table.firstChild.firstChild
        td2 = table.firstChild.lastChild
        td1.setAttribute("class", "numbers")
        td2.setAttribute("class", "lines")
        
        lines = self.filecontents.splitlines()
        counts = [str(x+1) for x in range(len(lines))]
        pre = self.document.createElement("pre")
        pre.appendChild(self.document.createTextNode('\n'.join(counts)))
        td1.appendChild(pre)
        
        pre = self.document.createElement("pre")
        file_data = self.document.createCDATASection(self.filecontents)
        pre.appendChild(file_data)
        self.document.documentElement.setAttribute("xmlns",xml.dom.XHTML_NAMESPACE)
        td2.appendChild(pre)
        self.document.documentElement.appendChild(table)
        
        self.add_stylesheet()
        
        print self.document.toxml()
        self.document.unlink()
        
    def new_table(self, number_of_rows, number_of_columns):
        table = self.document.createElement("table")
        for row in range(number_of_rows):
            tr = self.document.createElement("tr")
            table.appendChild(tr)
            for column in range(number_of_columns):
                td = self.document.createElement("td")
                tr.appendChild(td)
        return table
        
    def new_document(self):
        implementation = xml.dom.minidom.getDOMImplementation()
        return implementation.createDocument(xml.dom.XHTML_NAMESPACE, "xhtml", None)

    
    def add_stylesheet(self):
        self.style = self.document.createElement("style")
        self.style.setAttribute("type", "text/css")
        self.head.appendChild(self.style)
        stylesheet_string = ".lines {background-color: red;}"
        stylesheet_cdata = self.document.createCDATASection(stylesheet_string)
        self.style.appendChild(stylesheet_cdata)
        
    
if __name__ == '__main__':
    filename = "support/lint_check.py"
    i = InterfaceToPyLint()
    if i.is_available():
        pass #print i.run_onfile(filename)
        x = MakeHTMLDomFromFile(filename)
        x.start()

