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
    stdoutfname = None
    if sys.argv[1] == 'none':
        stdout = None
    else:
        stdoutfname=translate_filename_for_os(sys.argv[1])
        stdout = open(stdoutfname,'w')
        
    if sys.argv[2] == 'none':
        stderr = None
    else:
        stderrfname=translate_filename_for_os(sys.argv[2])
        if sys.argv[2] != '/dev/null' and stdoutfname == stderrfname:
            stderr = open(stderrfname,'a')
        else:
            stderr = open(stderrfname,'w')
    
    
    stdin = open(translate_filename_for_os('/dev/null'),'r')
    
    returncode = call(
        sys.argv[3:], 
        stdout = stdout,
        stderr = stderr,    
        stdin = stdin,
        close_fds = False
    )

    stdin.close()
    
    if not stdout is None:
        stdout.close()
        
    if not stderr is None:
        stderr.close()
    
    sys.exit(returncode)
