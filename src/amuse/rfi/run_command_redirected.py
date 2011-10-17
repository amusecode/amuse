from subprocess import call

import sys
import os.path
import time
import signal

if __name__ == '__main__':

    if sys.argv[1] == 'none':
        stdout = None
    else:
        fname=sys.argv[1]
        stdout = open(fname,'w')
        
    if sys.argv[2] == 'none':
        stderr = None
    else:
        fname=sys.argv[2]
        stderr = open(fname,'w')
    
    
    stdin = open('/dev/null','r')
    
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
