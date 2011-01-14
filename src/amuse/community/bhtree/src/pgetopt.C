/*
 *  pgetopt.C: portable and simplified version of getopt() in UNIX, system V
 *.............................................................................
 *    version 1:  Nov 1989   Piet Hut               email: piet@iassns.bitnet
 *                           Institute for Advanced Study, Princeton, NJ, USA
 *    version 2:  Dec 1992   Piet Hut  --  adopted to the new C++-based starlab
 *    version 2.0A:  Dec 1998 Jun Makino  --  pskipopt() function added.
 *.............................................................................
 *  non-local function: 
 *    pgetopt
 *    pskipopt
 *.............................................................................
 *     Command line argument passing is done in System V UNIX style. Howver,
 *  instead of the function  getopt() , we use a private, portable version
 *  pgetopt(), given here.  The reason for providing our own function is that
 *  there are non-system V UNIX versions, such as some Berkeley UNICES, which
 *  do not provide  getopt().  By always using our own  pgetopt() we guarantee
 *  uniform behavior, independent of the UNIX version used.
 *     Restrictions: the length of an option is one character, and the argument
 *                   has to follow the corresponding option, but separated by
 *                   a space (or tab).
 *                   Options have to start with a minus sign, but more than one
 *                   option can be combined after the same minus sign.
 *  Examples: The following command lines all give the same effect:
 *  		mumble -a -b 10 -c
 *		mumble  -c -a -b 10
 *		mumble -b 10 -ca
 *		mumble -acb 10
 *  but the following versions are illegal, and will give error messages:      
 *		mumble -a -b10 -c
 *		mumble -ab10c
 * 		mumble -a -c -b
 *.............................................................................
 */

#include <stdlib.h>
#include <string.h>
//#include  <stdiostream.h>
#include <iostream>

using namespace std;

/*-----------------------------------------------------------------------------
 *  streq  --  a macro which returns 1 if two strings are equal, 0 otherwise
 *-----------------------------------------------------------------------------
 */
#define  streq(x,y)  (strcmp((x), (y)) == 0)

#define VERSION_OPTION_A   "-version"
#define VERSION_OPTION_B   "--version"

char *poptarg;

/*-----------------------------------------------------------------------------
 *  pgetopt  --  each call to  pgetopt()  returns the next option encountered
 *               on the command line, starting from the beginning. If the 
 *               option occurs in the  optstr  string as well, the option 
 *               itself is returned (as the int value of the character; 
 *               options are limited to one char). Otherwise, the character '?'
 *               is returned. If the end of the string is reached, the value 
 *               -1 is returned.  If an option is followed by the character ':'
 *               in  optstr , then a command line argument is expected, 
 *               separated by spaces. If such an argument if found, the pointer
 *                poptarg  is pointed at that string, so that the calling 
 *               function can access the argument.  If such an argument is not
 *               found, an error message is given.
 *
 *               NOTE: The option "-version" or "--version" to any Starlab
 *                     program is legal, and will result in the program
 *                     printing the current  version number on cerr and
 *                     terminating.
 *
 *----------------------------------------------------------------------------
 */
static int  argv_counter = 1;       /* skip argv[0], the command name */
static int  argv_offset = 0;
int  pgetopt(int argc, char ** argv,  char * optstr)
    {
    int  optstr_counter;
    char  option_char;
    
    if (argv_counter >= argc)
	return(-1);                  /* signal that we've run out of options */

    if (argv_offset == 0)
    	if (argv[argv_counter][argv_offset] != '-')
	    {
	    cerr << "pgetopt: command line argument does not begin with -\n";
	    exit(1);
	    }
	else
	    argv_offset++;

#ifdef IN_STARLAB    
//  We have a legal switch.  First check to see if all we want to
//  know is the STARLAB version number.

    if (streq(argv[argv_counter], VERSION_OPTION_A)
	|| streq(argv[argv_counter], VERSION_OPTION_B))
	{
	cerr << "Starlab version " << STARLAB_VERSION << endl;
	exit(0);
	}
#endif
    option_char = argv[argv_counter][argv_offset];

    optstr_counter = 0;
    while (optstr[optstr_counter] != option_char)
        if (optstr[optstr_counter] == '\0')
	    {
/*
 * prepare for the next call to pgetopt():
 */	
	    if (argv[argv_counter][++argv_offset] == '\0')
		{
		argv_counter++;
		argv_offset = 0;
		}
	    return('?');                                /* unexpected option */
	    }
	else
	    optstr_counter++;

    if (optstr[++optstr_counter] == ':')
        {
	if (argv[argv_counter][argv_offset + 1] != '\0')
	    {
	    cerr << "pgetopt: option -" << option_char
		 << " not followed by a space and argument\n";
            exit(1);
            }
	if (++argv_counter >= argc)
            {
            cerr << "pgetopt: option -" << option_char
                 << " requires a space-separated argument\n";
	    exit(1);
	    }
	poptarg = argv[argv_counter];
/*
 * prepare for the next call to pgetopt():
 */	
	argv_counter++;
	argv_offset = 0;
	}
    else
        {
/*
 * prepare for the next call to pgetopt():
 */	
	if (argv[argv_counter][++argv_offset] == '\0')
	    {
	    argv_counter++;
	    argv_offset = 0;
	    }
	}
    
    return(option_char);
    }



/*-----------------------------------------------------------------------------
 *  pskipopt  --  each call to  pskipopt()  inclements the internal pointer
 *                used by pgetopt by one, to skip next argument. This function
 *                is useful when an option takes more than one parameter.
 *----------------------------------------------------------------------------
 */

void pskipopt()
{
    argv_counter ++;
}

/* endof: pgetopt.C */





