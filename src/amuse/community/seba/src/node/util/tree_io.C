
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// tree_io.C
//
//	The most important (pipeable) output routines have been
//	modified to work on HP/HPUX-8.x -- S. McMillan, 3/3/93
//
//	Modified to read/write compressed snapshots (cout only),
//	under control of environment variable STARLAB_ZIP.
//				-- S. McMillan, X. Zhuge, 6/22/98
//	Minor leak fix		-- PJT , 11/24/98

#include "starlab_vector.h"
#include "util_io.h"
#include "story.h"
#include "node.h"

void  node::log_comment(char * comment)
{
    if (log_story)
	add_story_line(log_story, comment);
}

void  node::log_history(int argc, char ** argv)	
{
    if (log_story) {
	char *hist = gethist(argc, argv);

	// Note from Steve, 05/30/03.  Gethist returns a string
	// containing a line break.  This doesn't coexist well
	// with the "col" output format, which is line oriented.
	// Rather than rewrite gethist, split the line here...

	int n = strlen(hist);
	if (n > 0) {
	    char *tmp = new char[n+1];
	    int j = 0;
	    for (int i = 0; i < strlen(hist); i++) {
		if (hist[i] == '\n') {
		    strncpy(tmp, hist+j, i-j);
		    tmp[i-j] = '\0';
		    add_story_line(log_story, tmp);
		    j = i + 1;
		}
	    }
	    if (j < n)
		add_story_line(log_story, hist+j);
	    delete tmp;
	}
	delete hist;
    }
}

ostream& node::print_log_story(ostream& s)
{
    story * the_log_story = node::get_log_story();
    if (the_log_story)
        put_story(s, *the_log_story);
    return s;
}

ostream& node::print_hydro_story(ostream& s)
{
    if (hbase)
	hbase->print_hydro_story(s);
    return s;
}

ostream& node::print_star_story(ostream& s,
				int short_output)	// default = 0
{
    if (sbase)
	sbase->print_star_story(s, short_output);
    return s;
}

istream & node::scan_log_story(istream& s, char * line)
{
    node::set_log_story(get_story(s, line));
    return s;
}

istream & node::scan_hydro_story(istream& s)
{
    if (hbase)
	hbase->scan_hydro_story(s);
    else {

	// No hydrobase is defined, but we have already read the
	// START_HYDRO string.  Must keep reading the input stream
	// until the matching END_HYDRO is found.

	char input_line[MAX_INPUT_LINE_LENGTH];
	while(get_line(s,input_line), !matchbracket(END_HYDRO, input_line))
	    ;
    }
    return s;
}

istream & node::scan_star_story(istream& s, int level)
{
    if (sbase)
	sbase->scan_star_story(s, level);
    else {

	// (See hydrobase note above.)

	char input_line[MAX_INPUT_LINE_LENGTH];
	while (get_line(s,input_line), !matchbracket(END_STAR, input_line));
    }
    return s;
}

#include <string.h>
#define BUF_SIZE 1024

static char format_string[BUF_SIZE];

// Note added by Steve 7/6/98.  Do NOT use format_label() to print out
// labels of more than one particle in the same print statement if the
// particles are identified by index rather than by label.  Because we
// use a single string here, only one format can be stored at any
// time.  Therefore, since all calls to format_label() will be executed
// before the string is sent to the output stream, a statement of the
// form
//
//	cerr << bi->format_label() << bj->format_label() << endl;
//
// will actually print out the label of bj twice!
//
// Better to use print_label() instead in those circumstances.

char* node::format_label() const
{
    if (is_valid()) {

	// Precedence:	print name string if defined
	//		otherwise, print index if non-negative
	//		otherwise, print '?'

	if (name != NULL) {
	    strncpy(format_string, name, BUF_SIZE-1);
	    format_string[BUF_SIZE-1] = '\0';
	} else if(index >= 0) {
	    sprintf(format_string, "%d", index);	// SLWM removed leading
							// "#", July 1998
	} else {
	    sprintf(format_string, "?");
	}
    } else
	sprintf(format_string, "(invalid)");
	
    return format_string;
}

bool node_contains(node * b, int i)		// overloaded
{
    if (b->is_parent()) {
        for_all_nodes(node, b, bi)
	    if (bi->get_index() == i)
		return true;
    } else
        if (b->get_index() == i)
	    return true;

    return false;
}

bool node_contains(node * b, const char* s)		// overloaded
{
    if (b->is_parent()) {
        for_all_nodes(node, b, bi)
	    if (bi->name_is(s))
		return true;
    } else
        if (b->name_is(s))
	    return true;

    return false;
}

bool clump_contains(node * b, int i)		// overloaded
{
    return node_contains(b->get_top_level_node(), i);
}

bool clump_contains(node * b, const char *s)		// overloaded
{
    return node_contains(b->get_top_level_node(), s);
}

bool node::name_is(const char* s) const
{
    return streq(format_label(), s);
}

void node::print_label(ostream & s) const
{
    s << format_label();
}

void node::pretty_print_node(ostream & s) const		// default = cerr
{
    print_label(s);					// no trailing endl
}

void node::pretty_print_tree(ostream & s,		// default = cerr
			     int level)			// default = 0
{
    int  k = level;
    while (k--)	s << "  ";

    pretty_print_node(s);

    if (mass != 1) s << "        m = " << mass;
    s << endl;

    if (is_parent())
	for_all_daughters(node, this, d)
	    d->pretty_print_tree(s, level + 1);
}

void pp(const node *b, ostream & s)
{
    s << "(";
    b->pretty_print_node(s);
    for_all_daughters(node, b, daughter)
	pp(daughter, s);	
    s << ")";
}

void pp2(const node *b, ostream & s,  int level)
{
    for (int i = 0; i<level*2; i++) {s << " ";}
    b->pretty_print_node(s);
    s << "\n";
    for_all_daughters(node, b, daughter)
	pp2(daughter, s, level + 1);	
}

//------------------------------------------------------------------------


local node *get_node_recursive(istream& s,
			       npfp the_npfp,
			       hbpfp the_hbpfp,
			       sbpfp the_sbpfp,
			       bool use_stories,
			       int level)
{
    // Functions the_hbpfp and the_sbpfp allow the possibility that the
    // hydrobase and starbase pointers point to something more complex
    // than the (default) simple base classes.
    //
    // As of 1/04, only the default is used (SLWM).

    node *b = (*the_npfp)(the_hbpfp, the_sbpfp, use_stories);
    char line[MAX_INPUT_LINE_LENGTH];

    get_line(s, line);

    node *elder_sister = (node *)42;	// to make some compilers happy

    // Would be highly desirable to have the code ignore unexpected data
    // when seeking e.g. a START_PARTICLE line.  Currently, the input is
    // very intolerant of even extra whitespace in the input stream.
    //
    // However, the logic below makes this a little tricky... (Steve, 6/01)

    while (!matchbracket(END_PARTICLE, line)) {
	if (matchbracket(START_DYNAMICS, line)) {

	    // Virtual function scan_dyn_story() allows the calling
	    // class to control handling of Dyn data format.

	    b->scan_dyn_story(s);

	} else if (matchbracket(START_HYDRO, line)) {

	    // Virtual function scan_hydro_story() allows the calling
	    // class to control handling of Hydro data.  As of 1/04,
	    // all classes inherit the node version defined above, and
	    // hence use the hydrobase version of scan_hydro_story(),
	    // which simply saves all data as story text.

	    b->scan_hydro_story(s);

	} else if (matchbracket(START_STAR, line)) {

	    // Virtual function scan_star_story() allows the calling
	    // class to control handling of Star data.  As of 1/04, all
	    // classes except pdyn (and tdyn) inherit the node version
	    // defined above, and hence use the starbase version of
	    // scan_star_story(), which simply saves all data as story
	    // text.
	    //
	    // The pdyn version reads Star data, but saves it as pdyn
	    // data, without using a star pointer.

	    b->scan_star_story(s, level);

	} else if (matchbracket(START_LOG, line)) {

	    // bug: every node gets a log story, but when you see
	    // one from input, you set this to be the new one
	    // and thus never dealloc the old one ???

	    b->scan_log_story(s, line);

	} else if (matchbracket(START_PARTICLE, line)) {
	    node *daughter =
		get_node_recursive(s, the_npfp, the_hbpfp, the_sbpfp,
				   use_stories, level+1);
	    if (b->get_oldest_daughter() == NULL) {
		b->set_oldest_daughter(daughter);
	    } else {
		daughter->set_elder_sister(elder_sister);
		elder_sister->set_younger_sister(daughter);
	    }
	    daughter->set_parent(b);
	    elder_sister = daughter;
	} else {

	    char keyword[MAX_INPUT_LINE_LENGTH];
	    const char *val = getequals(line, keyword);

	    if (val) {
		if (!strcmp("i",keyword)) {
		    int index = strtol(val, NULL, 10);
		    b->set_label(index);
		} else if (!strcmp("name",keyword)) {
		    char cptr[MAX_INPUT_LINE_LENGTH];
		    sscanf(val,"%s",cptr);
		    b->set_label(cptr);
		} else if (!strcmp("N",keyword)) {   // N is not read in here;
		    ;                                // instead N is recomputed
		}                                    // at output time.
		else {
		    cerr << line <<" unexpected\n";
		    exit(1);
		}
	    }
        }
        get_line(s, line);
    }
    return b;
}

local node *get_node_init(istream& s,
			  npfp the_npfp,
			  hbpfp the_hbpfp,
			  sbpfp the_sbpfp,
			  bool use_stories)
{

    if (!check_and_skip_input_line(s, START_PARTICLE))
	return NULL;

    node *root = get_node_recursive(s, the_npfp, the_hbpfp, the_sbpfp,
				    use_stories, 0);
    root->set_root(root);
    return root;
}

#ifdef HAS_GZIP
#include <pfstream.h>	// (May not exist on all systems...)
#endif

node *get_node(istream& s,		// default = cin
	       npfp the_npfp,		// default = new_node
	       hbpfp the_hbpfp,		//	     etc.
	       sbpfp the_sbpfp,	 	//
	       bool use_stories)	// default = true
{
    // If STARLAB_USE_GZIP is defined, the input will be decompressed
    // by gzip first.  This is accomplished by redefining the stream s.

    // Note that, for compression to be carried out, we must (1) compile
    // HAS_GZIP set, then (2) run with STARLAB_USE_GZIP defined.

    node *n = NULL;

#ifdef HAS_GZIP
    if (char * zip = getenv("STARLAB_USE_GZIP")) {
	ipfstream sz("|gzip -d -f");
	if (sz) {
	    n = get_node_init(sz, the_npfp, the_hbpfp, the_sbpfp,
			      use_stories);
	}
    }
#endif

    if (!n)
	n = get_node_init(s, the_npfp, the_hbpfp, the_sbpfp, use_stories);

    if (n)
	n->check_and_correct_node();		// action depends strongly on
						// which class is involved

    return n;
}

//----------------------------------------------------------------------


static bool first_log = true;

inline local void put_node_body(node *b,
				ostream &s = cout,
				bool print_xreal = true,
				int short_output = 0)
{
    // Now have one "short" option that basically is full output,
    // for restart purposes.  Check for it explicitly.

    bool short_short = (short_output && short_output != 4);

    if (short_short) {

	// Need i for scatter3; see if this breaks anything...

	// put_string(s, "  name = ", b->format_label());
	if (b->get_index() >= 0) put_integer(s, "  i = ", b->get_index());
	if (b->get_name() != NULL) put_string(s, "  name = ", b->get_name());

    } else {
	if (b->get_index() >= 0) put_integer(s, "  i = ", b->get_index());
	if (b->get_name() != NULL) put_string(s, "  name = ", b->get_name());
	put_integer(s, "  N = ", b->n_leaves());
    }

    if (!short_short || (b->is_root() && first_log)) {
	b->print_log_story(s);
	first_log = false;
    }

    // ------------------------------------------------------------
    // *** Changed the way output is done (Steve, 5/01) ***
    //
    // Virtual print_dyn_story() functions now just print out
    // "known" quantities in the form
    //
    //		keyword = value
    //
    // and each uses its base class function to avoid repetition.
    // We add enclosing "parens" and any additional story output HERE.

    put_story_header(s, DYNAMICS_ID);			// new

    b->print_dyn_story(s, print_xreal, short_output);

    if (!short_short && b->get_dyn_story())
        put_story_contents(s, *b->get_dyn_story());	// new

    put_story_footer(s, DYNAMICS_ID);			// new

    // ------------------------------------------------------------

    // For now, all short output is handled as part of Dyn.

   if (!short_short) {
       b->print_hydro_story(s);
       b->print_star_story(s, short_output);
   }
}

inline local void put_node_recursive(node *b,
				     ostream & s = cout, 
				     bool print_xreal = true,
				     int short_output = 0)
{
    put_story_header(s, PARTICLE_ID);

    put_node_body(b, s, print_xreal, short_output);

    for_all_daughters(node, b, daughter)
	put_node_recursive(daughter, s, print_xreal, short_output);

    put_story_footer(s, PARTICLE_ID);
}

void put_single_node(node *b,
		     ostream & s,		// default = cout
		     bool print_xreal,		// default = true
		     int short_output)		// default = 0
{
    // Same as put_node, but without recursion.

    put_story_header(s, PARTICLE_ID);
    put_node_body(b, s, print_xreal, short_output);
    put_story_footer(s, PARTICLE_ID);
}

// put_node: the function that does all the work of outputting nodes
//	     of all sorts...

void put_node(node *b,
	      ostream & s, 		// default = cout
	      bool print_xreal,		// default = true
	      int short_output)		// default = 0
{
    // If STARLAB_USE_GZIP is set, everything written to cout will be
    // compressed by gzip.

#ifdef HAS_GZIP
    if (&s == (ostream *) &cout) {
	if (char * zip = getenv("STARLAB_USE_GZIP")) {
	    opfstream sz("|gzip -c -f");
	    if (sz) {
		put_node_recursive(b, sz, print_xreal, short_output);
		return;
	    }
	}
    }
#endif
    put_node_recursive(b, s, print_xreal, short_output);
}

//----------------------------------------------------------------------


local void forget_node_recursive(istream& s)
{
    char line[MAX_INPUT_LINE_LENGTH];

    get_line(s, line);

    while(!matchbracket(END_PARTICLE, line)) {
	if(matchbracket(START_PARTICLE, line))                 
	    forget_node_recursive(s);
	get_line(s, line);
    }
}

// forget_node: reads in a complete node structure, just as get_node does, but 
//              without storing anything.  This is useful in the function
//              snapprune, which prunes a long list of snapshots.
//              Piet, 941125.

bool forget_node(istream& s)
{
    if (!check_and_skip_input_line(s, START_PARTICLE)) {
	return FALSE;
    }
    forget_node_recursive(s);
    return TRUE;
}
