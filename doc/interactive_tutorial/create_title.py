import sys
import os.path

def create_title(name):
    filename, ext = os.path.splitext(name)
    filename = os.path.basename(filename)
    filename = filename.replace('-',' - ')
    filename = filename.replace('_',' ')
    lines = []
    lines.append('='*len(filename))
    lines.append(filename)
    lines.append('='*len(filename))
    lines.append('')
    return '\n'.join(lines) 

if __name__ == '__main__':
    print create_title(sys.argv[1])
