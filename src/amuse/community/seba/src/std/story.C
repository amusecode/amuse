
//  story.C

//// Starlab story manipulation functions.
////
//// Options:
//// None.
////
//// Written by Piet Hut and Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

// The multiple instances of functions below are crying out
// to be templatized...  (Steve, 7/08)

#include "story.h"

#ifndef TOOLBOX

#define DEBUG false

// delete [] new_string lines added by Steve 8/5/97.

// New code to check if story exists before using member
// functions added by Steve 10/12/00.

story::~story()
{
    story * si = first_daughter_node;
    story * sn;
    while (si) {
	sn = si->next_story_node;
	delete si;
	si = sn;
    }
    if (DEBUG) cerr << "~story for " << this << endl << flush;
    if (text) {
	if (DEBUG) cerr << "~story: deleting text " << text << endl << flush;
	delete [] text;
    }
}

int  is_chapter(story * s)
{
    if (s)
	return s->get_chapter_flag();
    else
	return 0;
}

int  is_chapter_begin_line(const char * line)
{
    return (*line == chapter_begin_char);
}

int  is_chapter_end_line(const char * line)
{
    return (*line == chapter_end_char);
}

local int get_story_line(istream & str, char * line)
{
    str.get(line, MAX_STORY_LINE_LENGTH, '\n');

    if(str.eof())
        return 0;

    char c;
    if(str.get(c) && c!='\n') {
	cerr << "get_story_line : input line too long :'"<<line<<"'\n";
	exit(1);
    }

    return 1;
}

void add_daughter_story(story * s, story * d)
{
    if (!s || !d) return;

    story *fc = s->get_first_daughter_node();
    story *lc = s->get_last_daughter_node();

    // Redundant (?) check added by Steve 11/1/98.  Loses all
    // but the last daughter, but the best we can do...

    if (fc == NULL && lc != NULL) {
	warning("add_daughter_story: repairing inconsistent story list");
	s->set_first_daughter_node(lc);
    }

    if (lc) {
	lc->set_next_story_node(d);
	s->set_last_daughter_node(d);
    } else {
	s->set_first_daughter_node(d);
	s->set_last_daughter_node(d);
    }
}

void add_chapter(story * s, story * chap)
{
    if (!s || !chap) return;

    add_daughter_story(s, chap);
}

void rm_daughter_story(story * s, story * d)
{
    if (!s || !d) return;

    story * fd = s->get_first_daughter_node();
    story * ad;                                    // a daughter
    story * nd;                                    // next daughter

    if (fd == d) {
	s->set_first_daughter_node(d->get_next_story_node());

	// Added by Steve, 11/1/98, to deal with case of a single daughter:

	if (s->get_last_daughter_node() == d)
	    s->set_last_daughter_node(NULL);

    } else {

	ad = fd;
	nd = ad->get_next_story_node();
	while (nd != d)
	    {
	    ad = nd;
	    nd = ad->get_next_story_node();
	    }
	nd = d->get_next_story_node();
	ad->set_next_story_node(nd);
	if (nd == NULL)
	    s->set_last_daughter_node(ad);
    }

    delete d;    
}

story* mk_story_line()
{
    story* s = new story(0);
    return s;
}

story* mk_story_line(const char * line)
{
    story* s = new story(0);
    s->set_text(line);
    return s;
}

story* mk_story_chapter()       {story* s = new story(1); return s;}

story* mk_story_chapter(const char * title)
{
    story* s = new story(1);
    s->set_text(title);
    return s;
}

void add_story_line(story * s, const char * line)
{
    if (!s || !line) return;

    story * new_s = mk_story_line(line);
    add_daughter_story(s, new_s);
}

story* get_chapter(istream& str, const char* line)
{
    if (!is_chapter_begin_line(line)) {
	cerr << "get_chapter: first line not chapter_begin_line";
	exit(1);
    }

    story * chap = mk_story_chapter(++line);

    char new_line[MAX_STORY_LINE_LENGTH];

    while (get_story_line(str, new_line) && !is_chapter_end_line(new_line)) {
	if (is_chapter_begin_line(new_line))
	    add_chapter(chap, get_chapter(str, new_line));
	else
	    add_story_line(chap, new_line);
    }

    if (new_line == NULL) {
	cerr << "get_chapter: new_line == NULL before end of chapter\n";
	exit(1);
    }

    if (!streq(new_line+1, chap->get_text())) {
	cerr << "get_chapter: closing title ``" << new_line+1
	     << "'' differs from opening title ``" << chap->get_text()
	     << "''\n";
	exit(1);
    }

    return chap;
}

story* get_story(istream & str, const char *line)
{
    return get_chapter(str, line);
}

story* get_story(istream& str)
{ 
    char line[MAX_STORY_LINE_LENGTH];

    if (!get_story_line(str, line))
        return(NULL);

    if (!is_chapter_begin_line(line)) {
        cerr << "get_story: first line not chapter_begin_line\n";
	exit(1);
    }

    return get_chapter(str, line);
}

void put_headline(ostream& str, story& s)
{
#ifndef BAD_GNU_IO
    str << chapter_begin_char << s.get_text() << "\n";
#else
    fprintf(stdout, "%c%s\n", chapter_begin_char, s.get_text());
#endif
}

void put_tailline(ostream& str, story& s)
{
#ifndef BAD_GNU_IO
    str << chapter_end_char << s.get_text() << "\n";
#else
    fprintf(stdout, "%c%s\n", chapter_end_char,s. get_text());
#endif
}

void put_line_text(ostream& str, story& s)
{
#ifndef BAD_GNU_IO
    str << s.get_text() << "\n";
#else
    fprintf(stdout, "%s\n", s.get_text());
#endif
}

void put_chapter(ostream& str, story& s)
{
    if (!is_chapter(&s))
	{
	cerr << "put_chapter: not a story\n";
	exit(1);
	}
    put_headline(str, s);

    for (story * d = s.get_first_daughter_node(); d != NULL;
						  d = d->get_next_story_node())
	{
	if (is_chapter(d))
	    put_chapter(str, *d);
	else
	    put_line_text(str, *d);
	}
    str << flush;

    put_tailline(str, s);
}

void put_simple_story_contents(ostream& str, story& s,
			       const char *prefix)	// default = NULL
{
    // Assume s is a chapter and ignore any subchapters.
    // Intended for use with put_col() -- Log and Dyn stories,
    // without substructure.  Easy to extend if ever necessary.
    //						Steve, 5/03

    if (is_chapter(&s)) {
	for (story * d = s.get_first_daughter_node(); d != NULL;
	     d = d->get_next_story_node()) {

	    if (!is_chapter(d)) {
		if (prefix) str << prefix;
		put_line_text(str, *d);
	    }
	}
    }
}

void put_simple_story_contents(FILE *fp, story& s,
			       const char *prefix)	// default = NULL
{
    // Assume s is a chapter and ignore any subchapters.
    // Intended for use with put_col() -- Log and Dyn stories,
    // without substructure.  Easy to extend if ever necessary.
    //						Steve, 5/03

    if (is_chapter(&s)) {
	for (story * d = s.get_first_daughter_node(); d != NULL;
	     d = d->get_next_story_node()) {

	    if (!is_chapter(d)){
		if (prefix) 
		    fprintf(fp, "%s", prefix);
		fprintf(fp, "%s\n", d->get_text());
	    }
	}
    }
}

void put_story_contents(ostream& str, story& s,
			const char *prefix)	// default = NULL
{
    if (!is_chapter(&s))
	{
	cerr << "put_story_contents: not a story\n";
	exit(1);
	}

    for (story * d = s.get_first_daughter_node(); d != NULL;
	 d = d->get_next_story_node()) {

        if (is_chapter(d))
	    put_chapter(str, *d);
	else {
	    if (prefix) str << prefix;
	    put_line_text(str, *d);
	}
    }
}

void put_story(ostream& str, story& s)
{
    put_chapter(str, s);
}

// ostream & operator << (ostream & s, story * sptr)
//    {put_story(s, *sptr); return s;}

// istream & operator >> (istream & s, story * & sptr)
//    {sptr = get_story(s); return s;}

/*-----------------------------------------------------------------------------
 *     The following macros specify the buffer lengths for strings containing
 *  different types of quantities.
 *     The numbers given below are overestimates; if memory size is an 
 *  important issue, these numbers can be decreases somewhat, after careful
 *  analysis of how much space is actually needed in the worst case (beware of
 *  details such as minus signs, exponents, minus signs in exponents as well as
 *  in the mantissa, etc.).
 *-----------------------------------------------------------------------------
 */
#define  BYTE_LENGTH          8
#define  SAFE_INT_LENGTH     (5 + (BYTE_LENGTH*sizeof(int))/3)
#define  SAFE_REAL_LENGTH    (10 + (BYTE_LENGTH*sizeof(real))/3)
#define  SAFE_STRING_LENGTH   1    /* this should hold the string terminator */
#define  SAFE_VECTOR_LENGTH  (3 * (SAFE_REAL_LENGTH + 2))
#define  EXTRA_LENGTH  5           /* for initial "  " and " = " */

/*-----------------------------------------------------------------------------
 *  write_iq  --  write an integer quantity to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_iq(story * a_story_line, const char * name,
		     int value)
{
    if (!a_story_line || !name) return;

    char * new_string;
    int  new_string_length;
    
    new_string_length = EXTRA_LENGTH + strlen(name) + SAFE_INT_LENGTH;
    new_string  = new char[new_string_length];

    sprintf(new_string, "  %s = %d", name, value);

    a_story_line->set_text(new_string);
    delete [] new_string;
}

/*-----------------------------------------------------------------------------
 *  write_ulq  --  write an unsigned long integer quantity to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_ulq(story * a_story_line, const char * name,
		      unsigned long value)
{
    if (!a_story_line || !name) return;

    char * new_string;
    int  new_string_length;
    
    new_string_length = EXTRA_LENGTH + strlen(name) + SAFE_INT_LENGTH;
    new_string  = new char[new_string_length];

    sprintf(new_string, "  %s = %lu", name, value);

    a_story_line->set_text(new_string);
    delete [] new_string;
}

/*-----------------------------------------------------------------------------
 *  write_ullq  --  write an unsigned long long quantity to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_ullq(story * a_story_line, const char * name,
		       unsigned long long value)
{
    if (!a_story_line || !name) return;

    char * new_string;
    int  new_string_length;
    
    new_string_length = EXTRA_LENGTH + strlen(name) + SAFE_INT_LENGTH;
    new_string  = new char[new_string_length];

    sprintf(new_string, "  %s = %llu", name, value);

    a_story_line->set_text(new_string);
    delete [] new_string;
}

/*-----------------------------------------------------------------------------
 *  write_rq  --  write a real quantity to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_rq(story * a_story_line, const char * name, real value,
		     int p = STD_PRECISION)
{
    if (!a_story_line || !name) return;

    char * new_string;
    int  new_string_length;
    char format[128];
    
    new_string_length = EXTRA_LENGTH + strlen(name) + SAFE_REAL_LENGTH;
    new_string  = new char[new_string_length];

    // Allow variable precision (SLWM 6/99).

    sprintf(format, "  %%s = %%.%dg", p);

    // sprintf(new_string, "  %s = %lf", name, value);	// new format seems OK
    sprintf(new_string, format, name, value);		// here (SLWM 26/7/98)

    a_story_line->set_text(new_string);
    delete [] new_string;
}

/*-----------------------------------------------------------------------------
 *  write_sq  --  write a string quantity to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_sq(story * a_story_line, const char * name,
		     const char * value)
{
    if (!a_story_line || !name) return;

    char * new_string;
    int  new_string_length;
    
    new_string_length = EXTRA_LENGTH + strlen(name)
				     + strlen(value) + SAFE_STRING_LENGTH;
    new_string  = new char[new_string_length];

    sprintf(new_string, "  %s = %s", name, value);

    a_story_line->set_text(new_string);
    delete [] new_string;
}

//-----------------------------------------------------------------------------
//  write_vq  --  write a vector quantity to a line.
//-----------------------------------------------------------------------------

local void write_vq(story * a_story_line, const char * name, vec & value,
		    int p = STD_PRECISION)
{
    if (DEBUG) {
	cerr << "in write_vq" << endl;
	PRC(a_story_line); PRL(name);
	PRC(value); PRL(p);
    }

    if (!a_story_line || !name) return;

    char * new_string;
    int  new_string_length;
    char format[128];
    
    new_string_length = EXTRA_LENGTH + strlen(name) + SAFE_VECTOR_LENGTH;
    new_string  = new char[new_string_length];

    // Allow variable precision (SLWM 6/99).

    sprintf(format, "  %%s = %%.%dg %%.%dg %%.%dg",	// new format seems OK
	    p, p, p);
    // sprintf(new_string, "  %s = %.6g %.6g %.6g",	// old format...

    sprintf(new_string, format, name, value[0], value[1], value[2]);

    if (DEBUG) {
	PRL(new_string);
	PRL(1);
	char *c = a_story_line->get_text();
	PRL(10);
	//	unsigned int ic =(unsigned int) c;
	//	PRL(ic);
	if (c) {
	    PRL(*c);
	    PRL(c);
	    PRL(11);
	    PRL(a_story_line->get_text());
	    PRL(12);
	    put_line_text(cerr, *a_story_line);
	}
	PRL(2);
    }
    a_story_line->set_text(new_string);
    if (DEBUG) PRL(3);
    delete [] new_string;
    if (DEBUG) PRL(4);
}

/*-----------------------------------------------------------------------------
 *  write_ra  --  write a real array to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_ra(story * a_story_line, const char * name,
		     real * value, int n)
{
    if (!a_story_line || !name) return;

    char *new_string, *tmp;
    int  new_string_length;

    new_string_length = EXTRA_LENGTH + strlen(name)
      				+ n * (SAFE_REAL_LENGTH + 2);
    new_string  = new char[new_string_length];
    tmp = new char[SAFE_REAL_LENGTH + 2];

    sprintf(new_string, "  %s =", name);

    for (int i = 0; i < n; i++) {

//	sprintf(tmp, " %.6g", value[i]);  	// %g seems to cause problems
	sprintf(tmp, " %lf", value[i]);   	// here, for unknown reasons

	strcat(new_string, tmp);
    }

    a_story_line->set_text(new_string);
    delete [] new_string;
    delete [] tmp;
}

/*-----------------------------------------------------------------------------
 *  write_ia  --  write an integer array to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_ia(story * a_story_line, const char * name,
		     int * value, int n)
{
    if (!a_story_line || !name) return;

    char *new_string, *tmp;
    int  new_string_length;

    new_string_length = EXTRA_LENGTH + strlen(name)
      				+ n * (SAFE_INT_LENGTH + 2);
    new_string  = new char[new_string_length];
    tmp = new char[SAFE_INT_LENGTH + 2];

    sprintf(new_string, "  %s =", name);

    for (int i = 0; i < n; i++) {

	sprintf(tmp, " %d", value[i]);

	strcat(new_string, tmp);
    }

    a_story_line->set_text(new_string);
    delete [] new_string;
    delete [] tmp;
}

/*-----------------------------------------------------------------------------
 *  write_ia  --  write an unsigned long array to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_ia(story * a_story_line, const char * name,
		     unsigned long * value, int n)
{
    if (!a_story_line || !name) return;

    char *new_string, *tmp;
    int  new_string_length;

    new_string_length = EXTRA_LENGTH + strlen(name)
      				+ n * (SAFE_INT_LENGTH + 2);
    new_string  = new char[new_string_length];
    tmp = new char[SAFE_INT_LENGTH + 2];

    sprintf(new_string, "  %s =", name);

    for (int i = 0; i < n; i++) {

	sprintf(tmp, " %d", value[i]);

	strcat(new_string, tmp);
    }

    a_story_line->set_text(new_string);
    delete [] new_string;
    delete [] tmp;
}

/*-----------------------------------------------------------------------------
 *  write_ia  --  write an unsigned long long array to a line.
 *-----------------------------------------------------------------------------
 */
local void  write_ia(story * a_story_line, const char * name,
		     unsigned long long * value, int n)
{
    if (!a_story_line || !name) return;

    char *new_string, *tmp;
    int  new_string_length;

    new_string_length = EXTRA_LENGTH + strlen(name)
      				+ n * (SAFE_INT_LENGTH + 2);
    new_string  = new char[new_string_length];
    tmp = new char[SAFE_INT_LENGTH + 2];

    sprintf(new_string, "  %s =", name);

    for (int i = 0; i < n; i++) {

	sprintf(tmp, " %d", value[i]);

	strcat(new_string, tmp);
    }

    a_story_line->set_text(new_string);
    delete [] new_string;
    delete [] tmp;
}

/*-----------------------------------------------------------------------------
 *  qmatch  --  checks whether a line contains a physical quantity with a
 *              matching name; if so, returns TRUE, otherwise FALSE.
 *         note: 
 *              leading blanks are discarded, both in the name and in the text
 *              of a_line ; trailing blanks in the text up to the equal sign
 *              are also discarded; if the equal sign is missing (i.e. if the
 *              line is not a proper quantity line) FALSE is returned.
 *-----------------------------------------------------------------------------
 */
local bool qmatch(story * a_story_line, const char * name)
{
    if (!a_story_line || !name) return false;

    int  i, j;
    char *s;

    if ((s = a_story_line->get_text()) == NULL)
        return(FALSE);
    i = 0;
    while (s[i] == ' ')
        i++;

    j = 0;
    while (name[j] != '\0')
        {
	if (name[j] == ' ')
            j++;
	else 
	    break;
	}

    while (name[j] != '\0')
        {
	if (s[i] != name[j])
            return(FALSE);
	i++;
	j++;
	}

    while (s[i] == ' ')
        i++;

    if (s[i] != '=')
        return(FALSE);
    if (s[++i] != ' ')
        return(FALSE);

    return(TRUE);
}

/*-----------------------------------------------------------------------------
 *  find_qmatch  --  within a story, find a line which contains a physical
 *                   quantity with a matching name; or return NULL if such a
 *                   line is not found.
 *-----------------------------------------------------------------------------
 */
story * find_qmatch(story * a_story, const char * name)
{
    if (!a_story || !name) return NULL;

    story * a_chapter;

    a_chapter = a_story->get_first_daughter_node();

    while (a_chapter != NULL)
        {
	if (qmatch(a_chapter, name))
	    return(a_chapter);
	else
	    a_chapter = a_chapter->get_next_story_node();
	}

    return(NULL);
}

/*-----------------------------------------------------------------------------
 *  get_qstring  --  
 *               note:
 *                    no check for presence of '=' ; BEWARE
 *-----------------------------------------------------------------------------
 */
local char * get_qstring(story * a_story_line)
{
    if (!a_story_line) return NULL;

    char *c;
    
    c = a_story_line->get_text();
    
    while (*(++c) != '=')
	;
    while (*(++c) == ' ')
        ;

    return(c);
}

/*-----------------------------------------------------------------------------
 *  getiq  --  reads an integer quantity from a line in a story
 *-----------------------------------------------------------------------------
 */
int getiq(story * a_story, const char * name, bool verbose)
{
    if (!a_story || !name) return 0;

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	if (verbose) cerr << "getiq: no quantity found with name \""
	                  << name << "\"" << endl;
	return -VERY_LARGE_INTEGER;
	// exit(1);
	}

    return atoi(get_qstring(story_line));
}

/*-----------------------------------------------------------------------------
 *  getulq  --  reads an unsigned long integer quantity from a line in a story
 *-----------------------------------------------------------------------------
 */
unsigned long getulq(story *  a_story, const char * name, bool verbose)
{
    if (!a_story || !name) return 0;

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	if (verbose) cerr << "getulq: no quantity found with name \""
	                  << name << "\"" << endl;
	return 0;
	// exit(1);
	}

    return strtoul(get_qstring(story_line), (char**)NULL, 10);
}

/*-----------------------------------------------------------------------------
 *  getullq  --  reads an unsigned long long quantity from a line in a story
 *-----------------------------------------------------------------------------
 */
unsigned long long getullq(story *  a_story, const char * name, bool verbose)
{
    if (!a_story || !name) return 0;

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	if (verbose) cerr << "getulq: no quantity found with name \""
	                  << name << "\"" << endl;
	return 0;
	// exit(1);
	}

    return strtoull(get_qstring(story_line), (char**)NULL, 10);
}

/*-----------------------------------------------------------------------------
 *  getrq  --  reads a real quantity from a line in a story
 *-----------------------------------------------------------------------------
 */
real getrq(story * a_story, const char * name, bool verbose)
{
    if (!a_story || !name) return 0.0;

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL) {

	if (verbose) cerr << "getrq: no quantity found with name \""
			  << name << "\"" << endl;
	return -VERY_LARGE_NUMBER;
    }

    return atof(get_qstring(story_line));
}

/*-----------------------------------------------------------------------------
 *  getsq  --  reads a string quantity from a line in a story
 *             BEWARE: a pointer is returned to the original string.
 *                     this has two dangerous aspects:
 *                     1) the user may inadvertantly change the story line;
 *                     2) the story line may be changed before the user uses
 *                        the string.
 *             note: when needed, we can easily provide an alternative function
 *                   `getnsq' that returns a new string (as a copy).
 *-----------------------------------------------------------------------------
 */
char *getsq(story * a_story, const char * name, bool verbose)
{
    if (!a_story || !name) return NULL;

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	if (verbose) cerr << "getiq: no quantity found with name \""
	                  << name << "\"" << endl;
	return (char*)NULL;
	// exit(1);
	}

    return get_qstring(story_line);
}

/*-----------------------------------------------------------------------------
 *  getvq  --  reads a vector quantity from a line in a story
 *-----------------------------------------------------------------------------
 */
vec getvq(story *  a_story, const char * name, bool verbose)
{
    if (!a_story || !name) return vec(-VERY_LARGE_NUMBER,
					      -VERY_LARGE_NUMBER,
					      -VERY_LARGE_NUMBER);

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL) {
	if (verbose) cerr << "getiq: no quantity found with name \""
	                  << name << "\"" << endl;
	return vec(-VERY_LARGE_NUMBER,
		      -VERY_LARGE_NUMBER,
		      -VERY_LARGE_NUMBER);
	// exit(1);
    }

    vec v;
    sscanf(get_qstring(story_line), "%lf %lf %lf", &v[0], &v[1], &v[2]);

    return v;
}

#define GET_ARR_BUFSIZ 8192

static char s[GET_ARR_BUFSIZ];

local void get_array(char * sin, real * x, int n)
{
    // Extract real array x from (a copy of) string sin.

    // Note: SPZ says this code as it stands causes problems on HP/g++...
    //       Certainly, it isn't very elegant!

    // Work with a copy of the input string.

    // Hmmm... problems with ccmalloc if we allocate space for s
    // and delete it later (no problem with just g++).  For now, use
    // static space and hope it is big enough.

    // char *s = new char[strlen(sin)];	 // possible that this "new" is not
    					 // really static, so the string is
					 // automatically deleted on exit?

    if (n > GET_ARR_BUFSIZ)
	warning("get_array: truncating input string");

    strncpy(s, sin, GET_ARR_BUFSIZ);

    char *s1 = s, *s2 = s;
    int i = 0;

    while (*(s2++) > '\0') {
        if (*(s2-1) > ' ' && *s2 <= ' ') {
	    char save = *s2;

	    *s2 = '\0';
	    x[i++] = atof(s1);
	    *s2 = save;

	    if (i >= n) break;
	    s1 = s2;
	}
    }

    // delete [] s;
}

local void get_array(char * sin, int * x, int n)	// overloaded!
{
    // Copy of get_array, for int instead of real.  Previous comments apply.

    if (n > GET_ARR_BUFSIZ)
	warning("get_array: truncating input string");

    strncpy(s, sin, GET_ARR_BUFSIZ);

    char *s1 = s, *s2 = s;
    int i = 0;

    while (*(s2++) > '\0') {
        if (*(s2-1) > ' ' && *s2 <= ' ') {
	    char save = *s2;

	    *s2 = '\0';
	    x[i++] = atoi(s1);				// <-- only change!
	    *s2 = save;

	    if (i >= n) break;
	    s1 = s2;
	}
    }
}

local void get_array(char * sin, unsigned long * x, int n)  // overloaded!
{
    // Copy of get_array, for int instead of real.  Previous comments apply.

    if (n > GET_ARR_BUFSIZ)
	warning("get_array: truncating input string");

    strncpy(s, sin, GET_ARR_BUFSIZ);

    char *s1 = s, *s2 = s;
    int i = 0;

    while (*(s2++) > '\0') {
        if (*(s2-1) > ' ' && *s2 <= ' ') {
	    char save = *s2;

	    *s2 = '\0';
	    x[i++] = atoi(s1);				// <-- only change!
	    *s2 = save;

	    if (i >= n) break;
	    s1 = s2;
	}
    }
}

/*-----------------------------------------------------------------------------
 *  getra  --  reads a real array of length n from a line in a story
 *-----------------------------------------------------------------------------
 */
void getra(story *  a_story, const char * name, real * x, int n, bool verbose)
{
    if (!a_story || !name) {
	*x = -VERY_LARGE_NUMBER;
	return;
    }

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL) {

	if (verbose) cerr << "getra: no quantity found with name \""
	                  << name << "\"" << endl;

	*x = -VERY_LARGE_NUMBER;
	return;

    }

    get_array(get_qstring(story_line), x, n);
}

/*-----------------------------------------------------------------------------
 *  getia  --  reads an integer array of length n from a line in a story
 *-----------------------------------------------------------------------------
 */
void getia(story *  a_story, const char * name, int * x, int n, bool verbose)
{
    if (!a_story || !name) {
	*x = -VERY_LARGE_INTEGER;
	return;
    }

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL) {

	if (verbose) cerr << "getia: no quantity found with name \""
	                  << name << "\"" << endl;

	*x = -VERY_LARGE_INTEGER;
	return;

    }

    get_array(get_qstring(story_line), x, n);
}

/*-----------------------------------------------------------------------------
 *  getia  --  reads an unsigned long array of length n from a line in a story
 *-----------------------------------------------------------------------------
 */
void getia(story *  a_story, const char * name, unsigned long * x,
	   int n, bool verbose)
{
    if (!a_story || !name) {
	*x = 0;
	return;
    }

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL) {

	if (verbose) cerr << "getia: no quantity found with name \""
	                  << name << "\"" << endl;

	*x = 0;
	return;

    }

    get_array(get_qstring(story_line), x, n);
}

/*-----------------------------------------------------------------------------
 *  putiq  --  write an integer quantity at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putiq(story * a_story, const char * name, int value)
{
    if (!a_story || !name) return;

    story * story_line;
    
    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	story_line = new story;	
	add_daughter_story(a_story, story_line);
	}

    write_iq(story_line, name, value);
}

/*-----------------------------------------------------------------------------
 *  putiq  -- write an unsigned long integer quantity at the end of a
 *            story; or, if a previous quantity of the same name is found,
 *            overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putiq(story * a_story, const char * name, unsigned long value)
{
    if (!a_story || !name) return;

    story * story_line;
    
    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	story_line = new story;	
	add_daughter_story(a_story, story_line);
	}

    write_ulq(story_line, name, value);
}

/*-----------------------------------------------------------------------------
 *  putiq  -- write an unsigned long long quantity at the end of a
 *            story; or, if a previous quantity of the same name is found,
 *            overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putiq(story * a_story, const char * name, unsigned long long value)
{
    if (!a_story || !name) return;

    story * story_line;
    
    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	story_line = new story;	
	add_daughter_story(a_story, story_line);
	}

    write_ullq(story_line, name, value);
}


// Handy function for debugging (Steve, 10/98)...

void dump_story(story* s, int indent)
{
    if (!s) return;

    if (indent == 0) cerr << endl;

    PRI(indent); cerr << "story s = " << s << endl;
    PRI(indent); cerr << "flag = " << s->get_chapter_flag() << endl;
    PRI(indent); cerr << "text:  \"" << s->get_text() << "\"" << endl;
    PRI(indent); cerr << "next = " << s->get_next_story_node() << endl;
    PRI(indent); cerr << "first = " << s->get_first_daughter_node() << endl;
    PRI(indent); cerr << "last = " << s->get_last_daughter_node() << endl;

    if (indent >= 0) {
	if (s->get_first_daughter_node())
	    for (story* ss = s->get_first_daughter_node(); ss != NULL;
		 ss = ss->get_next_story_node())
		dump_story(ss, indent+4);
    }
}

/*-----------------------------------------------------------------------------
 *  putrq  --  write a real quantity at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putrq(story * a_story, const char * name, real value,
	   int precision)			// default = STD_PRECISION
{
    if (!a_story || !name) return;

    story * story_line;

    bool dbg = false;
    if (name[0] == '-') {
	name++;
	dbg = true;
    }

    if (dbg) {
	cerr << endl << "in putrq..." << endl;
	dump_story(a_story);
    }

    if ((story_line = find_qmatch(a_story, name)) == NULL) {
	if (dbg) cerr << "making new story" << endl;
	story_line = new story;	
	add_daughter_story(a_story, story_line);
    }


    if (dbg) {
	dump_story(a_story);
	PRL(story_line);
    }

    write_rq(story_line, name, value, precision);

    if (dbg) {
	cerr << "wrote " << name << " = " << value << " to story" << endl;
	put_story(cerr, *a_story);
	PRL(story_line->get_text());
	PRL(get_qstring(story_line));
	dump_story(a_story);
    }
}

/*-----------------------------------------------------------------------------
 *  putra  --  write a real array at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putra(story * a_story, const char * name, real * value, int n)
{
    if (!a_story || !name) return;

    story * story_line;
    
    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	story_line = new story;	
	add_daughter_story(a_story, story_line);
	}

    write_ra(story_line, name, value, n);
}

/*-----------------------------------------------------------------------------
 *  putia  --  write an integer array at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putia(story * a_story, const char * name, int * value, int n)
{
    if (!a_story || !name) return;

    story * story_line;
    
    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	story_line = new story;	
	add_daughter_story(a_story, story_line);
	}

    write_ia(story_line, name, value, n);
}

/*-----------------------------------------------------------------------------
 *  putia  --  write an unsigned long array at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putia(story * a_story, const char * name, unsigned long * value, int n)
{
    if (!a_story || !name) return;

    story * story_line;
    
    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	story_line = new story;	
	add_daughter_story(a_story, story_line);
	}

    write_ia(story_line, name, value, n);
}

/*-----------------------------------------------------------------------------
 *  putia  --  write an unsigned long long array at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putia(story * a_story, const char * name, unsigned long long * value, int n)
{
    if (!a_story || !name) return;

    story * story_line;
    
    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	story_line = new story;	
	add_daughter_story(a_story, story_line);
	}

    write_ia(story_line, name, value, n);
}

/*-----------------------------------------------------------------------------
 *  putsq  --  write a string quantity at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putsq(story * a_story, const char * name, const char * value)
{
    if (!a_story || !name) return;

    story * story_line;
    
    if ((story_line = find_qmatch(a_story, name)) == NULL)
	{
	story_line = new story;	
	add_daughter_story(a_story, story_line);
	}

    write_sq(story_line, name, value);
}

/*-----------------------------------------------------------------------------
 *  putvq  --  write a vector quantity at the end of a story;
 *             or, if a previous quantity of the same name is found,
 *             overwrite the line containing that quantity.
 *-----------------------------------------------------------------------------
 */
void putvq(story * a_story, const char * name, vec & value,
	   int precision)			// default = STD_PRECISION
{
    if (!a_story || !name) return;

    story * story_line;

    if ((story_line = find_qmatch(a_story, name)) == NULL) {
	story_line = new story;	
	add_daughter_story(a_story, story_line);
    }

    write_vq(story_line, name, value, precision);
}

/*-----------------------------------------------------------------------------
 *  rmq  --  removes a quantity line from a story.
 *           returns 1 when a line is actually removed,
 *           returns 0 when no line is found for quantity name `name'.
 *-----------------------------------------------------------------------------
 */
int  rmq(story *  a_story, const char * name)
{
    if (!a_story || !name) return 0;

    story * story_line;

    if (story_line = find_qmatch(a_story, name)) {
	rm_daughter_story(a_story, story_line);
	return 1;
    } else
	return 0;
}

/*-----------------------------------------------------------------------------
 *  is_quantity_name  --  checks whether there exists a quantity by the name
 *                        'name' in story 'a_story'.
 *                        returns 1 when that quantity is found,
 *                        returns 0 when that quantity is not found.
 *-----------------------------------------------------------------------------
 */
int  is_quantity_name(story * a_story, const char * name)
{
    if (!a_story || !name) return 0;

    if (find_qmatch(a_story, name))
	return 1;
    else
	return 0;
}


#else

main(int argc, char** argv)
{
    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.15 $", _SRC_);

    story * s;

//    while (cin >> s)
//	{
/*
	putiq(s, "test message", 41);
	putiq(s, "test message", 42);
	putiq(s, "yet another test message", 137);
	putrq(s, "pi", 3.14);
	putrq(s, "pi", 3.14, 2);
	putrq(s, "pi", 3.14, 10);
	putsq(s, "star cluster", "47 Tuc");
	vec * tmp = new vec(1,2,3);
	putvq(s, "1-2-3-test vector", *tmp);
	cout << "test message = " << getiq(s, "test message") << endl;
	cout << "yet another test message = "
	     << getiq(s, "yet another test message") << endl;
	cout << "pi = " << getrq(s, "pi") << endl;
	cout << "star cluster = " << getsq(s, "star cluster") << endl;
	cout << "1-2-3-test vector = " << getvq(s, "1-2-3-test vector")
	     << endl;
*/
//        cout << s;
//	}

//    cerr << "story.C : TOOLBOX : Normal exit\n";
}

#endif
