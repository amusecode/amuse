from subprocess import call

import sys
import os.path
import time
import signal

def translate_filename_for_os(filename):
    if sys.platform == 'win32':
        if filename == '/dev/null':
            return 'nul'
        else:
            return filename
    else:
        return filename

if __name__ == '__main__':

    if sys.argv[1] == 'none':
        stdout = None
    else:
        fname=translate_filename_for_os(sys.argv[1])
        stdout = open(fname,'w')
        
    if sys.argv[2] == 'none':
        stderr = None
    else:
        fname=translate_filename_for_os(sys.argv[2])
        stderr = open(fname,'w')
    
    
    stdin = open(translate_filename_for_os('/dev/null'),'r')
    
    returncode = call(
        sys.argv[3:], 
        stdout = stdout,
        stderr = stderr,    
        stdin = stdin,
    #   close_fds = True
    )

    stdin.close()
    
    if not stdout is None:
        stdout.close()
        
    if not stderr is None:
        stderr.close()
    
    sys.exit(returncode)
