import sys
import os
import os.path
import textwrap
from nose import main as nose_main
from optparse import OptionParser

_commands = []

def command(f):
    _commands.append(f)

@command
def test():
    """
    Run amuse-admin test, to perform units tests
    
    Without arguments this command will run the core units tests of
    amuse. All unit tests are run using nose and you can give all
    nose arguments here. Examples:
    
    amuse-admin test -v amuse.test.suite.codes_tests
    amuse-admin test -v amuse.test.suite.ext_tests
    
    """
        
    if len(sys.argv) == 1:
        arguments = list(sys.argv)
        arguments.append('-v')
        arguments.append('amuse.test.suite.core_tests')
        sys.argv = arguments
    
    if not os.path.exists('amuse-tests'):
        os.mkdir('amuse-tests')
        
    os.chdir('amuse-tests')
    try:
        nose_main()
    finally:
        os.chdir('..')
        
@command
def help():
    """
    Usage: amuse-admin help subcommand
    
    Type 'amuse-admin help <subcommand>' for help on a specific subcommand.
    """
    if len(sys.argv) > 1:
        commandname = sys.argv[1]
        for x in _commands:
            if x.__name__ == commandname:
                print textwrap.dedent(x.__doc__).strip()
    else:
        usage()

def usage():
    """
    Usage: amuse-admin command [options] [args]
    
    Type 'amuse-admin help <subcommand>' for help on a specific subcommand.
    """
    print textwrap.dedent(usage.__doc__).strip()
    print
    show_available_subcommands()
    
def show_available_subcommands():
    print "Available subcommands:"
    for x in _commands:
        print "    {0}".format(x.__name__)
    

def main():
    
    if len(sys.argv) == 1:
        usage()
    else:
        commandname = sys.argv[1]
        arguments = list(sys.argv)
        del arguments[1]
        sys.argv = arguments
        
        for x in _commands:
            if x.__name__ == commandname:
                x()
        



if __name__ == "__main__":
    main()
