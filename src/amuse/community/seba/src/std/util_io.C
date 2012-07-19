
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// util_io.C:  Functions for formatted I/O.
//	       Note that xreal functions now reside in xreal.C.

#include "starlab_vector.h"
#include "util_io.h"
#include "story.h"
#include <ctype.h>

#undef isalnum		// hacks for Irix 6.5 <ctype.h> backward compatibility
#undef isspace

int get_line(istream & s, char * line)
{
    s.get(line, MAX_INPUT_LINE_LENGTH, '\n');	// don't read the '\n'

    if (s.fail())

	// Failed to read anything -- probably a blank line.
	// Just reset and continue.

	s.clear();

    // Read and discard all characters up to and including the next newline.

    char c;
    while (s.get(c))
	if (c == '\n') break;

    return strlen(line);
}

int check_input_line(istream &s, const char* reference_string)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    get_line(s,input_line);
    if(s.eof()){
	cerr << "check_input_line : unexpected EOF, expected '";
	cerr << reference_string <<"'\n";
	exit(1);
    }
    return matchbracket(reference_string, input_line);
}    

int check_and_skip_input_line(istream &s, const char* reference_string)
{
//  char dummy;
    char input_line[MAX_INPUT_LINE_LENGTH];

    while (! get_line(s,input_line) && !s.eof()) {
	int shhh = 1;           // to make the SG compiler shut up
	if (shhh)               // to make the SG compiler shut up
	    return 0;
    }

#if 0
    cerr << "line = {" << input_line << "}" << endl;
    for (int j = 0; j < strlen(input_line); j++)
	cerr << (int)input_line[j] << endl;
#endif

    if (!matchbracket(reference_string, input_line)) {

	if (s.eof()) {
	    return 0;
	} else {
	    cerr << "Input line must be '" << reference_string;
	    cerr << "', I got '" << input_line << "'" << endl;
	    exit(1);
	}
    }
    return 1;
}    

int get_data_line(istream & s,char * input_line)
{
    get_line(s,input_line);
    return strcmp(input_line,")");
}

// Returns 1 if either token matches line, or first two chars of token
// matches line.  So matchbracket("(Particle", line) matches "(Particle"
// or "(P".

int matchbracket(const char *token, const char *line)
{
    // Skip whitespace (actually, now anything <= ' '):

    // while(*line == ' ' || *line == '\t') line++;
    while(*line <= ' ') line++;

    if(token[0] != line[0] || token[1] != line[1])
	return 0;

    // Uncomment the next line for a stricter check: accept either "(P"
    // or "(Particle", but not "(Part".  Maybe too strict (Steve, 8/04).

    // return (line[2] == '\0') || (0 == strcmp(token+2, line+2));
    return 1;
}

const char *getequals(const char *input_line, char *keyword)
{
    const char *cp = input_line;

    // Grab first token from line, like sscanf %s.

    while (isspace(*cp)) cp++;
    int i;
    for (i = 0; isalnum(*cp) || *cp == '_'; )
	keyword[i++] = *cp++;
    keyword[i] = '\0';

    cp = strchr(cp, '=');
    if (cp == NULL) {

	// Maybe the warning is unnecessary?

//	cerr << "getequals: expected 'keyword = value', but got '"
//	     << input_line << "'" << endl;

	return cp;
    }
    cp++;
    while (isspace(*cp)) cp++;
    return cp;
}

void set_vector_from_input_line(vec & v, char * input_line)
{
    char *eq = strchr(input_line, '=');
    if(eq)
	set_vector_from_string( v, eq+1 );
}

void set_vector_from_string(vec & v, char *val)
{
    real component[3];
    char *cp, *ep;
    component[0] = strtod(val, &cp);
    component[1] = strtod(cp, &cp);
    component[2] = strtod(cp, &ep);
    if(cp == ep) {
	cerr << "Expected three reals, got: " << val << endl;
	exit(1);
    }
    v = vec(component[0],component[1],component[2]);
}

//
//---------------------------------------------------------------------
//
// Horrible kludges:
// ----------------
//

static bool short_story_keywords = false;

bool use_short_story_keywords( bool useshort ) {
  bool was = short_story_keywords;
  short_story_keywords = useshort;
  return was;
}


void put_story_header(ostream & s, const char * id)
{
#ifdef BAD_GNU_IO
    if (s == cout) {
	if(short_story_keywords) {
	    fprintf(stdout, "(%c\n", id[0]);
	} else {
	    fprintf(stdout, "(%s\n", id);
	}
    }
    else
#endif
    if(short_story_keywords) {
	s << '(' << id[0] << endl;
    } else {
	s << '(' << id << endl;
    }
}

void put_story_footer(ostream & s, const char * id)
{
#ifdef BAD_GNU_IO
    if (s == cout) {
	if(short_story_keywords) {
	    fprintf(stdout, ")%c\n", id[0]);
	} else {
	    fprintf(stdout, ")%s\n", id);
	}
    }
    else
#endif
    if(short_story_keywords) {
	s << ')' << id[0] << endl;
    } else {
	s << ')' << id << endl;
    }
}

// NOTE use of precision here.  If the STARLAB_PRECISION environment
// variable is set, the first call to set_starlab_precision will use
// its value (whatever it may be).  Subsequent calls will return the
// same value originally read from the environment.  If no environment
// variable is set, a default value (currently 18 -- see precision.C)
// is used.

// See also xreal version in xreal.C.

void put_real_number(ostream & s, const char * label, real x)
{
    // Note from Steve (7/04).  The precision of a stream isn't
    // exactly what we need here.  It simply determines the number
    // of digits printed after the decimal point, *including*
    // non-significant leading zeroes.  In many cases, we want
    // precision to set the number of significant digits printed
    // (especially in put_node).  The C printf "g" format seems
    // to do this properly, while it is unclear if there is an
    // equivalent state function in C++.  For now, use sprinf to
    // write a string which is then written to stream s.  The code
    // is essentially that used with the obsolete BAD_GNU_IO macro,
    // but this has *nothing* to do with that ancient G++ bug.

    int old_precision = set_starlab_precision(s);

#ifdef BAD_GNU_IO

    // Weird g++ bug means that we have to use fprintf instead of <<.
    // The problem is that we have no way to associate a steam with
    // a FILE pointer.  Deal with cout and cerr and simply use the
    // untrustworthy code for other streams.

    static FILE *f = NULL;
    if (s == cout)
	f = stdout;
    else
	f = stderr;

    if (f) {
	static int local_precision = -1;
	int precision = get_starlab_precision();
	static char format[64];
	if (local_precision != precision) {
	    local_precision = precision;
	    int p = precision;
	    if (p < 0) p = 5;
	    sprintf(format, "%%s%%.%dg\n", p);
	}
	fprintf(f, format, label, x);
    } else
	s << label << x << endl;

#else

# if 1

    // Even if the I/O is OK, we still have to deal with the precision
    // problem described above.  Code repeats much of the "BAD_GNU" code.

    static char format[64], *outstring = NULL;
    static int local_precision = -1, p = 0, nout = 0;

    int precision = get_starlab_precision();
    if (local_precision != precision) {

	// Recreate the format string.

	local_precision = precision;
	p = precision;
	if (p < 0) p = 5;
	sprintf(format,"%%.%dg", p);
    }

    // See if we need to resize.

    if (outstring && nout < p+10) delete [] outstring;

    // (Re)create outstring if necessary.

    if (!outstring) {
	nout = p+10;
	outstring = new char[nout];
    }

    // Finally, create the string...

    sprintf(outstring, format, x);

    // ...and print it along with the label.

    s << label << outstring << endl;
    
# else
    s << label << x << endl;
# endif

#endif

    // Restore the current precision.

    if (get_starlab_precision() != old_precision)
	s.precision(old_precision);
}

void put_real_vector(ostream & s, const char * label, vec v)
{
    // See various notes in put_real_number above...

    int old_precision = set_starlab_precision(s);

#ifdef BAD_GNU_IO

    static FILE *f = NULL;
    if (s == cout)
	f = stdout;
    else
	f = stderr;

    if (f){
	static int local_precision = -1;
	int precision = get_starlab_precision();
	static char format[64];
	if (local_precision != precision) {
	    local_precision = precision;
	    int p = precision;
	    if (p < 0) p = 5;
	    sprintf(format,"%%s%%.%dg %%.%dg %%.%dg\n", p, p, p);
	}
	fprintf(stdout, format,	label, v[0], v[1], v[2]);
    } else 
	s << label << v << endl;

#else

# if 1

    static char format[64], *outstring = NULL;
    static int local_precision = -1, p = 0, nout = 0;

    int precision = get_starlab_precision();
    if (local_precision != precision) {

	// Recreate the format string.

	local_precision = precision;
	p = precision;
	if (p < 0) p = 5;
	sprintf(format,"%%.%dg %%.%dg %%.%dg", p, p, p);
    }

    // See if we need to resize.

    if (outstring && nout < 3*p+30) delete [] outstring;

    // (Re)create outstring if necessary.

    if (!outstring) {
	nout = 3*p+30;
	outstring = new char[nout];
    }

    // Finally, create the string...

    sprintf(outstring, format, v[0], v[1], v[2]);

    // ...and print it along with the label.

    s << label << outstring << endl;
    
# else
    s << label << v << endl;
# endif

#endif

    // Restore the current precision.

    if (get_starlab_precision() != old_precision)
	s.precision(old_precision);
}

void put_integer(ostream & s, const char * label, int i)
{
#ifndef BAD_GNU_IO
    s << label << i << endl;
#else
    if (s == cout)
	fprintf(stdout, "%s%d\n", label, i);
    else
	s << label << i << endl;
#endif
}

void put_string(ostream & s, const char * label, const char * str)
{
#ifndef BAD_GNU_IO
    s << label << str << endl;
#else
    if (s == cout) 
	fprintf(stdout, "%s%s\n", label, str);
    else
	s << label << str << endl;
#endif
}
