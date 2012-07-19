/*
 *  gethist.c: for remembering the command line, and the date of execution
 *.............................................................................
 *    version 1:  May 1989   Piet Hut               email: piet@iassns.bitnet
 *                           Institute for Advanced Study, Princeton, NJ, USA
 *    version 2:  Dec 1992   Piet Hut  --  adopted to the new C++-based starlab
 *.............................................................................
 *  non-local function: 
 *    gethist
 *.............................................................................
 *     A string is returned which contains the date at which a command is
 *  invoked followed with the full command line, including all arguments.
 *.............................................................................
 */

//// Check Starlab function for saving command line and date of execution.
////
//// Options:
//// None.
////
//// Written by Piet Hut.
////
//// Report bugs to starlab@sns.ias.edu.


#include <time.h>
#include "stdinc.h"

#ifndef TOOLBOX

/*-----------------------------------------------------------------------------
 *  contains_white  --  returns TRUE of the string contains one or more
 *                      white spaces; returns FALSE otherwise.
 *-----------------------------------------------------------------------------
 */
local bool  contains_white(char * s)
    {
    while (*s != '\0')
        if (*s++ == ' ')
	    return(TRUE);

    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  eff_length  --  effective length, needed to stored a command line argument:
 *                  one for each character, and two extra for the surrounding
 *                  "" quotes, in case the string contains one or more 
 *                  internal white spaces.
 *          example:
 *                  "hi"        --> hi
 *                  "hi there"  --> "hi there"
 *                
 *                  eff_length("hi") = 2
 *                  eff_length("hi there") = 10
 *-----------------------------------------------------------------------------
 */
local int  eff_length(char * s)
    {
    if (contains_white(s))
        return(strlen(s) + 2);
    else
        return(strlen(s));
    }

/*-----------------------------------------------------------------------------
 *  DATE_TOT_LENGTH  --  total length, in characters, of the date plus the
 *                       date separator; this is the length reserved in total
 *                       for the date information by gethist().
 *  DATE_STR_LENGTH  --  The fixed standard length of the date containing 
 *                       string, as returned by the UNIX operating system,
 *                       apart from the last NULL character;
 *  DATE_SEPARATOR  --  This string can be changed whenever desired:
 *                      when chosing a new DATE_SEPARATOR string, all date
 *                      length macros will be automatically adjusted.
 *-----------------------------------------------------------------------------
 */
#define  INTRO_PART        " ===>  "
#define  INTRO_LENGTH      strlen(INTRO_PART)
#define  STARLAB_PART      "Starlab "
#define  STARLAB_LENGTH    strlen(STARLAB_PART)
#define  VERSION_PART      VERSION
#define  VERSION_LENGTH    strlen(VERSION_PART)
#define  DATE_STR_LENGTH   24
#define  DATE_SEPARATOR    "\n       "
#define  DATE_SEP_LENGTH   strlen(DATE_SEPARATOR)
#define  VERSION_SEPARATOR    " : "
#define  VERSION_SEP_LENGTH   strlen(VERSION_SEPARATOR)
#define  DATE_TOT_LENGTH   INTRO_LENGTH + STARLAB_LENGTH + VERSION_LENGTH \
                       + DATE_STR_LENGTH + DATE_SEP_LENGTH + VERSION_SEP_LENGTH

#define USER_HEAD " (user "

/*-----------------------------------------------------------------------------
 *  gethist  --  returns a newly allocated string containing the date at which
 *               the command of the calling program was given, and the full
 *               command line.
 *          note:
 *               simply printing successive command line arguments does not 
 *               provide the full information, since the following two commands
 *                 gobble this thing   and   gobble "this thing"
 *               would both be printed as   gobble this thing  .
 *               Therefore, all command line arguments are first checked for
 *               the presence of internal white space, in which case they are
 *               provided with the "" quotes.
 *-----------------------------------------------------------------------------
 */
char *gethist(int argc, char ** argv)
    {
    int  i, j, k;
    int  n;
    long int  clock;              /* will contain time in seconds since 1970 */
    char *hist_string;            /* will point to newly allocated string    */
/*
 * determine length of new character array, and allocate memory for it:
 */
    n = DATE_TOT_LENGTH;
    for (i = 0; i < argc; i++)    /* total length, apart from the white      */
        n += eff_length(argv[i]); /* spaces separating the commands          */
    n += argc;                    /* one for each white space delimiter,     */
				  /* and a final NULL to end the hist_string */

    // Added by Steve (7/01):

    char *user_part = getenv("USER");
    int user_length = 0;

    if (user_part) user_length = strlen(user_part) + strlen(USER_HEAD) + 1;
    //									')'

    n += user_length;
    hist_string = new char[n];
    
    if (hist_string == NULL)
	{
	cerr << "gethist: not enough memory left for hist_string\n";
	exit (1);
	}
/*
 * write the date:
 */
    strcpy(hist_string, INTRO_PART);
    clock = (long int)time(0);                 /* see time(2), UNIX manual   */
    strcpy(hist_string + INTRO_LENGTH, ctime((time_t*)&clock));
					       /* see ctime(3C), UNIX manual */
    strcpy(hist_string + INTRO_LENGTH + DATE_STR_LENGTH, DATE_SEPARATOR);
    strcpy(hist_string + INTRO_LENGTH + DATE_STR_LENGTH + DATE_SEP_LENGTH,
	   STARLAB_PART);
    strcpy(hist_string + INTRO_LENGTH + DATE_STR_LENGTH + DATE_SEP_LENGTH
	   + STARLAB_LENGTH, VERSION_PART);

    if (user_part) {

	// Add in the name of the user running the program:

	char tmp[1024];
	strcpy(tmp, USER_HEAD);
	strcat(tmp, user_part);
	strcat(tmp, ")");

	strcpy(hist_string + INTRO_LENGTH + DATE_STR_LENGTH + DATE_SEP_LENGTH
	   + STARLAB_LENGTH + VERSION_LENGTH, tmp);
    }

    strcpy(hist_string + INTRO_LENGTH + DATE_STR_LENGTH + DATE_SEP_LENGTH
	   + STARLAB_LENGTH + VERSION_LENGTH + user_length, VERSION_SEPARATOR);
/*
 * write the command line:
 */
    j = DATE_TOT_LENGTH + user_length;
    for (i = 0; i < argc; i++)
	{
	k = 0;
	if (contains_white(argv[i]))
	    hist_string[j++] = '"';
	while (argv[i][k] != '\0')
	    hist_string[j++] = argv[i][k++];
	if (contains_white(argv[i]))
	    hist_string[j++] = '"';
	hist_string[j++] = ' ';                     /* white space delimiter */
	}
    hist_string[--j] = '\0';              /* providing the string terminator */
				          /* by overwriting the last  ' '    */
/*
 * done:
 */
    return(hist_string);
    }

/*===========================================================================*/

#else

//-----------------------------------------------------------------------------
//  main  --  driver to test gethist() .
//            for example:
//                gethist these are words "but this is a single string"
//-----------------------------------------------------------------------------

main(int argc, char ** argv)
{
    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.5 $", _SRC_);

    printf("%s\n", gethist(argc, argv));
}

#endif

/*===========================================================================*/

/* endof: gethist.c */
