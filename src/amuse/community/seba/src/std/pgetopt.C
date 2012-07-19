
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Portable and simplified version of getopt() in UNIX, system V.
//// This is a simple driver/test program for the library function.
//// The length of an option is one character, and the argument
//// has to follow the corresponding option, separated by a space (or tab).
//// Options must start with a minus sign, but more than one option can
//// be combined after the same minus sign.  Thus, the next four command
//// lines all have the same effect, while the following three are illegal
//// and will give error messages.
////
////              good -a -b 10 -c
////
////	          good -c -a -b 10
////
////	          good -b 10 -ca
////
////	          good -acb 10
////
////	          bad -a -b10 -c
////
////	          bad -ab10c 
////
//// 	          bad -a -c -b
////
//// Usage:   pgetopt
////
//// Options:
////          -a        test: no argument expected
////          -b        test: one argument expected
////          -c        test: two arguments expected
////          -d        test: three arguments expected
////          -e        test: gour arguments expected
////          -f        test: no argument expected
////          -g        test: one optional argument
////          -h        test: no argument expected
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

//.............................................................................
//    version 1:  Nov 1989   Piet Hut               email: piet@iassns.bitnet
//                           Institute for Advanced Study, Princeton, NJ, USA
//    version 2:  Dec 1992   Piet Hut  --  adopted to the new C++-based starlab
//    version 3:  Aug 1999   Steve McMillan
//				allow up to 16 arguments per option
//				allow single optional arguments
//.............................................................................
//  non-local function: 
//    pgetopt
//.............................................................................
//
//  Command line argument passing is done in System V UNIX style. However,
//  instead of the function  getopt(), we use a private, portable version
//  pgetopt(), given here.  The reason for providing our own function is that
//  there are non-system V UNIX versions, such as some Berkeley UNICES, which
//  do not provide  getopt().  By always using our own  pgetopt() we guarantee
//  uniform behavior, independent of the UNIX version used.
//
//  Restrictions: the length of an option is one character, and the argument
//                has to follow the corresponding option, but separated by
//                a space (or tab).
//                Options have to start with a minus sign, but more than one
//                option can be combined after the same minus sign.
//
//  Examples: The following command lines all give the same effect:
//                mumble -a -b 10 -c
//	          mumble  -c -a -b 10
//	          mumble -b 10 -ca
//	          mumble -acb 10
//  but the following versions are illegal, and will give error messages:      
//	          mumble -a -b10 -c
//	          mumble -ab10c
// 	          mumble -a -c -b
//.............................................................................

#include "stdinc.h"

#ifndef TOOLBOX

#define VERSION_OPTION_A   "-version"
#define VERSION_OPTION_B   "--version"

char *poptarg;				// global variable

// New version (Steve, 8/99): allow up to 16 numeric arguments to be
// associated with an option.  Store the pointers in the global array
// poparr[N_POP_ARG].  Retain poptarg for compatibility, but note that
// poparr[0] is the same thing.

#define N_POP_ARG	16

char *poparr[N_POP_ARG];		// global array

local inline char *get_version(const char *cvs_id)
{
    // Extract a version string from the CVS id.
    // CVS id format is "$Revision: 1.13 $".

    if (!cvs_id) return NULL;

    const char *start = strstr(cvs_id, "Revision:");
    if (!start) return NULL;

    start += 9;
    while (*start > 0 && *start <= ' ') start++;

    if (start - cvs_id > strlen(cvs_id)) return NULL;

    const char *end = start;
    while (*end > ' ') end++;

    if (*end == ' ') {
	int n = end-start+1;
	char *version = new char[n];
	strncpy(version, start, n-1);
	version[n-1] = 0;
	return version;
    } else {

	// Didn't find another space.  Assume no version was found.

	char *version = new char[4];
	strcpy(version, "0.0");
	return version;
    }
}

local inline char *get_name(const char *source)
{
    // Extract a program name from a source file specification.
    // Source format is /a/b/c/program.[cCfF].

    if (!source) return NULL;

    const char *start = source + strlen(source);
    while (start > source && *start != '/') start--;
    if (*start == '/') start++;

    const char *end = start;
    while (*end > 0 && *end != '.') end++;
    if (*end == '.') end--;

    if (end < start) return NULL;

    int n = end-start+2;
    char *name = new char[n];
    strncpy(name, start, n-1);
    name[n-1] = 0;

    return name;
}

//-----------------------------------------------------------------------------
//  pgetopt  --  each call to  pgetopt()  returns the next option encountered
//               on the command line, starting from the beginning. If the 
//               option occurs in the  optstr  string as well, the option 
//               itself is returned (as the int value of the character; 
//               options are limited to one char). Otherwise, the character '?'
//               is returned. If the end of the string is reached, the value 
//               -1 is returned.  If an option is followed by the character ':'
//               in  optstr , then a command line argument is expected, 
//               separated by spaces. If such an argument if found, the pointer
//                poptarg  is pointed at that string, so that the calling 
//               function can access the argument.  If such an argument is not
//               found, an error message is given.
//
//		 ..............................................................
//
//		 See above note for extension to this scheme, as of 8/99.
//		 Also added optional arguments, indicated by "." in optstr.
//
//		 Optional arguments *may not* start with "-", as that character
//		 is used to identify the start of the next option string...
//
//		 ..............................................................
//
//		 This function is quite unforgiving of errors in the format
//		 of the command line.  In typical use, there is no check in
//		 the calling program that poptarg (etc.) are properly set, or
//		 that the number of arguments is correct.  No default values
//		 are set, and errors generally cause the program to stop
//		  -- caveat emptor!
//
//               NOTE: The option "-version" or "--version" to any Starlab
//                     program is legal, and will result in the program
//                     printing the current Starlab and CVS versions on
//                     cout and terminating.
//
//----------------------------------------------------------------------------

int  pgetopt(int argc, char ** argv,
	     const char *optstr,
	     const char *cvs_id,	// default = NULL
	     const char *source)	// default = NULL
{
    static int argv_counter = 1;	// argument counter
					// skip argv[0], the command name
    static int argv_offset  = 0;	// character counter within argument
					// (multiple switches per argument
					//  are allowed)

    static char minus_chars[2] = {'-', 0};
    static char *minus = minus_chars;

    if (argv_counter >= argc)
	return -1;			// signal that we've run out of options

    if (argv_offset == 0) {
    	if (argv[argv_counter][argv_offset] != '-') {
	    cerr << "pgetopt: warning: command line argument \""
		 << argv[argv_counter]
		 << "\" does not begin with \"-\"\n";

	    // exit(1);		// too severe...

	    argv_counter++;
	    return '?';

	} else
	    argv_offset++;
    }

    char *version = get_version(cvs_id);
    char *name = get_name(source);

    // We have a legal switch.  First check to see if all we want to
    // know is the STARLAB version number.  New format and output to
    // cout is for use with help2man.

    if (streq(argv[argv_counter], VERSION_OPTION_A)
	|| streq(argv[argv_counter], VERSION_OPTION_B)) {

	// cerr << "Starlab version " << VERSION << endl;

	if (!version || !name)
	    cout << VERSION << endl;
	else {
	    cout << name << " (Starlab " << VERSION << "/CVS) " << version
		 << endl << endl;

	    // GNU boilerplate:

	    cout << "Copyright (C) 1994-2004, the Starlab development group."
		 << endl
		 << "This is free software; see the source for copying "
		 << "conditions.  There is NO"
		 << endl
		 << "warranty; not even for MERCHANTABILITY or FITNESS FOR "
		 << "A PARTICULAR PURPOSE."
		 << endl;
	}

	exit(0);
    }

    if (version) delete [] version;
    if (name) delete [] name;

    char option_char = argv[argv_counter][argv_offset];

    // Locate the next character in the argument list in optstr.

    int optstr_counter = 0;
    while (optstr[optstr_counter] != option_char) {

        if (optstr[optstr_counter] == '\0') {

	    // Prepare for the next call to pgetopt.

	    if (argv[argv_counter][++argv_offset] == '\0') {
		argv_counter++;
		argv_offset = 0;
	    }
	    return '?';			// couldn't find the specified
					// command-line option in optstr

	} else

	    optstr_counter++;
    }

    // Current command-line option is option_char, and it exists in
    // optstr, at position optstr_counter.

    poptarg = poparr[0] = NULL;		// default: no/optional argument

    // Look for arguments, and set up pointers to them.

    char opt;
    int narg = 0;

    while ((opt = optstr[++optstr_counter]) == ':' || opt == '.') {

	// Options must always be immediately followed by their arguments.

	if (narg == 0) {
	    if (argv[argv_counter][argv_offset + 1] != '\0') {
		cerr << "pgetopt: option \"-" << option_char
		     << "\" not followed by a space";
		if (opt == ':') cerr << " and argument";
		cerr << endl;
		exit(1);
	    }
	} else if (narg >= N_POP_ARG) {

	    // Shouldn't happen...

	    cerr << "pgetopt: too many arguments requested for option \"-"
		 << option_char << "\"" << endl;
            exit(1);
	}

	// Move to the next argument and check that it exists.

	if (++argv_counter >= argc && opt == ':') {
            cerr << "pgetopt: option \"-" << option_char
                 << "\" requires space-separated argument(s)\n";
	    exit(1);
	}

	// Check for a single optional argument.

	if (opt == '.'
	    && (argv_counter >= argc || argv[argv_counter][0] == '-')) {
	    argv_counter--;
	    break;
	}

	// Comptibility:

	if (narg == 0) poptarg = argv[argv_counter];

	poparr[narg++] = argv[argv_counter];
	poparr[narg] = minus;			// in case next arg is optional

    }

    // Prepare for the next call to pgetopt.

    if (poptarg || argv[argv_counter][++argv_offset] == '\0') {

	argv_counter++;
	argv_offset = 0;

    }
    
    return option_char;
}

#else

int main(int argc, char ** argv)
{
    extern char *poptarg;
    extern char *poparr[];				// new (8/99)
    int c;
    const char *param_string = "ab:c::d:::e::::fg.h";

    check_help();

    while ((c = pgetopt(argc, argv, param_string,
			"$Revision: 1.13 $", _SRC_)) != -1) {
	switch (c) {

	    case 'a':	cerr << "option a: no arguments"
			     << endl;
			break;

	    case 'b':	cerr << "option b: argument: "
			     << poptarg << " = " << poparr[0]
			     << endl;
			break;

	    case 'c':	cerr << "option c: arguments: "
			     << poparr[0] << "  " << poparr[1]
			     << endl;
			break;

	    case 'd':	cerr << "option d: arguments: "
			     << poparr[0] << "  " << poparr[1] << "  "
			     << poparr[2]
			     << endl;
			break;

	    case 'e':	cerr << "option e: arguments: "
			     << poparr[0] << "  " << poparr[1] << "  "
			     << poparr[2] << "  " << poparr[3]
			     << endl;
			break;

	    case 'f':	cerr << "option f: no arguments"
			     << endl;
			break;

	    case 'g':	if (poptarg)
			    cerr << "option g: optional argument: "
				 << poptarg << " = " << poparr[0]
				 << endl;
			else
			    cerr << "option g: no optional argument"
				 << endl;
			break;

	    case 'h':	cerr << "option h: no arguments"
			     << endl;
			break;

	}
    }
    return 0;
}

#endif
