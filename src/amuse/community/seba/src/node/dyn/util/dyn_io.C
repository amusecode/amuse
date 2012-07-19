
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Test Starlab dyn class I/O functions.  Define scan_dyn_story and
//// print_dyn_story for the dyn class.
////
//// Usage: dyn_io [OPTIONS] < input > output
////
//// Options:
//// None.
////
//// Written by the Starlab development group.
////
//// Report bugs to starlab@sns.ias.edu.

#include "dyn.h"
#include "util_io.h"

#ifndef TOOLBOX

// Initialize all static dyn data here...

xreal dyn::system_time          = 0.0;
real dyn::real_system_time      = 0.0;
bool dyn::use_sstar	        = false;

bool dyn::ignore_internal	= false;

unsigned int  dyn::external_field	= 0;

int  dyn::tidal_type	= 0;
real dyn::omega	= 0;
real dyn::omega_sq	= 0;
real dyn::alpha1	= 0;
real dyn::alpha3	= 0;
vec  dyn::tidal_center = vec(0,0,0);

real dyn::p_mass = 0;
real dyn::p_scale_sq = 0;
vec  dyn::p_center = vec(0,0,0);
bool dyn::p_friction = false;

real dyn::pl_coeff = 0;
real dyn::pl_scale = 0;
real dyn::pl_exponent = 0;
vec  dyn::pl_center = vec(0,0,0);

FILE* dyn::ifp = 0;
FILE* dyn::ofp = 0;

bool dyn::col_output = false;

void dyn::print_static(ostream& s)		// default = cerr
{
    node::print_static(s);

    s << "system_time = " << system_time << endl;
    s << "real_system_time = " << real_system_time << endl;

    s << "use_sstar = " << use_sstar << endl;
    s << "ignore_internal = " << ignore_internal << endl;

    s << "external_field = " << external_field << endl;
    s << "tidal_type = " << tidal_type << endl;

    s << "omega = " << omega << endl;
    s << "omega_sq = " << omega_sq << endl;
    s << "alpha1 = " << alpha1 << endl;
    s << "alpha3 = " << alpha3 << endl;
    s << "tidal_center = " << tidal_center << endl;

    s << "p_mass = " << p_mass << endl;
    s << "p_scale_sq = " << p_scale_sq << endl;
    s << "p_center = " << p_center << endl;

    s << "pl_coeff = " << pl_coeff << endl;
    s << "pl_scale = " << pl_scale << endl;
    s << "pl_exponent = " << pl_exponent << endl;
    s << "pl_center = " << pl_center << endl;
}

static bool read_xreal = false;

bool dyn::check_and_correct_node(bool verbose)	// default = false
{
    // cerr << "dyn::check_and_correct_node: "; PRL(this);

    bool ok = true;

    if (oldest_daughter) {		// for dyn, check that masses,
					// pos and vel are all consistent
	real m = 0;
	vec p = 0, v = 0;
	bool low = false;

	for_all_daughters(dyn, this, bb) {
	    if (bb->oldest_daughter)
		ok &= bb->check_and_correct_node(verbose);
	    real mm = bb->get_mass();
	    m += mm;
	    p += mm*bb->get_pos();
	    v += mm*bb->get_vel();
	}
	if (!ok) low = true;

	if (!twiddles(m, mass)) ok = false;
	if (!twiddles(abs(p), 0)) ok = false;
	if (!twiddles(abs(v), 0)) ok = false;

	mass = m;
	if (m > 0) {
	    p /= m;
	    v /= m;

	    // Force daughters to have zero CM quantities relative to parent.
	    // Impose an arbitrary limit on abs1(p) to prevent modification
	    // of systems that were (presumably) already intended to be in
	    // the CM frame.

//	    if (abs1(p) > 1.e-12) {
		for_all_daughters(dyn, this, bb) {
		    bb->inc_pos(-p);
		    bb->inc_vel(-v);
		}
	    }
// 	}

	if (!ok && verbose && parent == NULL) {
	    cerr << "check_and_correct_node: applied ";
	    if (low)
		cerr << "low-level";
	    else
		cerr << "top-level";
	    cerr << " mass/pos/vel correction" << endl;
	}
    }

    return ok;
}

istream & dyn::scan_dyn_story(istream& s)
{
    char input_line[MAX_INPUT_LINE_LENGTH];
    real last_real = false;

    while (get_line(s,input_line), !matchbracket(END_DYNAMICS, input_line)) {

	char keyword[MAX_INPUT_LINE_LENGTH];
	const char *val = getequals(input_line, keyword);

	if (val) {

	    if (!strcmp("real_system_time", keyword)) {

		read_xreal = true;
		last_real = true;

		// We don't know anything about parent nodes yet, so it is
		// not easy to know if we are the root node.  Rule: if we
		// find real_system_time, assume that we should read an
		// xreal as system_time.  Otherwise, read it as real.
		// The tortuous logic is to keep the determination of
		// which choice we should make completely local.
		//
		// Unfortunately, this logic must be duplicated in all
		// other *dyn::scan_dyn_story functions (see _dyn_io.C,
		// hdyn_io.C, sdyn3_io.C)...

	    } else if (!strcmp("system_time", keyword)) {

		// Check input format before reading.

		if (!last_real) read_xreal = false;

		if (read_xreal) {

		    //cerr << "dyn::scan_dyn_story: input "
		    //     << "time data type is xreal"
		    //     << endl;

		    // The following should set real_system_time too...

		    set_system_time(get_xreal_from_input_line(input_line));

		} else {

		    //cerr << "dyn::scan_dyn_story: input "
		    //     << "time data type is real"
		    //     << endl;

		    real_system_time = system_time = strtod(val, NULL);
		}
	    } else {

		last_real = false;

		if (!strcmp("m", keyword))
		    mass = strtod(val, NULL);
		else if (!strcmp("r", keyword))
		    set_vector_from_input_line(pos, input_line);
		else if (!strcmp("v", keyword))
		    set_vector_from_input_line(vel, input_line);
		else
		    add_story_line(dyn_story, input_line);
	    }
	}
    }

    return s;
}

ostream& dyn::print_dyn_story(ostream& s,
			      bool print_xreal,		// default = true
			      int short_output)		// default = 0
{
    // Modifications by Steve (5/01) to streamline output.

    // Print system time first (root node only).

    if (!parent) {

#ifdef USE_XREAL
	if (print_xreal) {

	    // Note (Steve 5/00): system_time is now xreal and hard to read.
	    // For convenience, also print out a real version of the time.
	    // By printing out real_system_time first, we set a flag that
	    // allows scan_dyn_story to read xreal, rather than real, input.

	    put_real_number(s, "  real_system_time  =  ", (real)system_time);
	    put_real_number(s, "  system_time  =  ", system_time);

	} else

	    put_real_number(s, "  system_time  =  ", (real)system_time);
#else

	put_real_number(s, "  system_time  =  ", system_time);

#endif
    }

    // Mass is now printed by node::print_dyn_story().

    node::print_dyn_story(s, print_xreal, short_output);

    put_real_vector(s, "  r  =  ", pos);
    put_real_vector(s, "  v  =  ", vel);

    return s;
}



local inline dyn* _fget_dyn(FILE* fp) {
#ifdef STARLAB_USE_FDYN
  // cerr << "get_dyn:  using fget_dyn." << endl;
  return fget_dyn(fp);
#else
  cerr << "get_dyn: fget_dyn not available." << endl;
  return NULL;
#endif
}

// Newest get_dyn() is a switch between several input formats.

dyn *get_dyn(istream & s,		// default = cin
	     hbpfp the_hbpfp,		// default = new_hydrobase,
	     sbpfp the_sbpfp,		// default = new_starbase,
	     bool use_stories)		// default = true
{

  // Currently the fast input fget_dyn() only uses stdin (and not cin) as
  // its input file pointer.  Get_col() uses get_ifp().  This is a mess,
  // which will be fixed soon...  (Steve and Ernie, 1/04)

  static enum { NUL, DYN, COL } format = NUL;
  static bool fast = false;
  FILE* ifp = dyn::get_ifp() ?: stdin;
  switch (format) {
  case DYN:
    return fast ? _fget_dyn(ifp) :
      (dyn*)get_node(s, new_dyn, the_hbpfp, the_sbpfp, use_stories);
  case COL: return get_col(s, new_dyn, the_hbpfp, the_sbpfp, use_stories);
  case NUL:
    fast = getenv("STARLAB_USE_FDYN");		// fast input off by default
    format = ungetc(getc(ifp), ifp) == '(' ? DYN : COL;
    return get_dyn(s, the_hbpfp, the_sbpfp, use_stories);
  }
}

// Called by get_dyn when input format is columns of numbers.

dyn* get_col(istream& s,
	     npfp the_npfp,		// note: return value is always node*
	     hbpfp the_hbpfp,
	     sbpfp the_sbpfp,
	     bool use_stories)
{
  dyn* root = (dyn*)the_npfp(the_hbpfp, the_sbpfp, use_stories);
  root->set_index(0);

  static unsigned lineno = 0;
  bool first_dyn = true;
  dyn* bo;

  bool first_data = true;

  while (true) {
    char line[256];
    if (dyn::get_ifp()) {
      if (!fgets(line, sizeof line, dyn::get_ifp())) break;
    } else {
      if (s.getline(line, sizeof line - 1)) strcat(line, "\n");
      else break;
    }

    // Only set col output here, once we know we aren't at EOF.

    if (first_data) dyn::set_col_output(true);

    switch (++lineno, line[0]) {
      case '\n': return root;
      case '#': continue;
      case ';':
	if (use_stories) {	      // deal with Log (;) and Dyn (;;) stories
	  line[strlen(line)-1] = '\0';
	  if (line[1] == ';') {
	    if (strlen(&line[2])) {
	      story *st = root->get_dyn_story();
	      if (st)
		add_story_line(st, line+2);
	    }
	  } else
	    if (strlen(&line[1])) root->log_comment(&line[1]);
        }
      continue;
    }

    // OK, looks like we have real data.  ASSUME that all comments and
    // non-particle data precede the data, and update the root node from
    // the Dyn story before proceeding.

    if (first_data) {

	// Upgrade root-node Dyn story data, if any, to dyn class data.
	// Delete story entries as we promote them.

	// Choices follow data in scan_dyn_story() above, except that we
	// never write xreal time, so omit real_system_time.  Also allow
	// "time" as a synomym for "system_time" (for non-Starlab users).

	story *st = root->get_dyn_story();

	if (find_qmatch(st, "time")) {
	    real system_time = getrq(st, "time");
	    root->set_system_time(system_time);
	    rmq(st, "time");
	}

	if (find_qmatch(st, "system_time")) {
	    real system_time = getrq(st, "system_time");
	    root->set_system_time(system_time);
	    rmq(st, "system_time");
	}

	if (find_qmatch(st, "r")) {
	    vec pos = getvq(st, "r");
	    root->set_pos(pos);
	    rmq(st, "r");
	}

	if (find_qmatch(st, "v")) {
	    vec vel = getvq(st, "v");
	    root->set_vel(vel);
	    rmq(st, "v");
	}

	first_data = false;
    }

    int id;
    double m, x[3], v[3];
    if (sscanf(line, "%i %lg %lg %lg %lg %lg %lg %lg\n", &id, &m, &x[0], &x[1],
               &x[2], &v[0], &v[1], &v[2]) != 8)
      cerr << "Malformed input on line #" << lineno << " at `" << line << "'\n",
	exit(1);
    dyn* const b = (dyn*)the_npfp(the_hbpfp, the_sbpfp, use_stories);
    b->set_index(id);
    b->set_mass(m);
    b->set_pos(vec(x[0], x[1], x[2]));
    b->set_vel(vec(v[0], v[1], v[2]));
    b->set_parent(root);
    if (first_dyn) root->set_oldest_daughter(b), first_dyn = false;
    else  b->set_elder_sister(bo), bo->set_younger_sister(b);
    bo = b;
  }

  if (!first_dyn) return root;
  else { delete root; return NULL; }
}

// This I/O problem is getting to be a real pain in the (my) ass...

void put_col(dyn* root, ostream& s, bool put_time) {

  int p = adjust_starlab_precision(-1);	  // odd name, but establishes
					  // the precision of the output

  story *st;		// note the elegant use of replicated code below...

  if (dyn::get_ofp()) {

    // Define format strings for use below.

    char fmt1[80], fmt2[80], fmt3[80];
    sprintf(fmt1, "%%s%%.%dg\n", p);
    sprintf(fmt2, "%%s%%.%dg %%.%dg %%.%dg\n", p, p, p);
    sprintf(fmt3, "%%d %%.%dg %%.%dg %%.%dg %%.%dg %%.%dg %%.%dg %%.%dg\n",
	    p, p, p, p, p, p, p);

    // Special treatment to preserve the root node member data
    // and Log and Dyn stories.

    st = root->get_log_story();
    if (st) put_simple_story_contents(dyn::get_ofp(), *st, ";");

    // Root node member data won't be properly written out by the loop below.
    // Save it as Dyn story data (see get_col above) -- less than elegant.
    // Don't bother with xreal output in this format, so don't print a
    // real_system_time line.

    fprintf(dyn::get_ofp(), fmt1, ";;  system_time  =  ",
	    (real)root->get_system_time());
    fprintf(dyn::get_ofp(), fmt2, ";;  r  =  ",
	    root->get_pos()[0], root->get_pos()[1], root->get_pos()[2]);
    fprintf(dyn::get_ofp(), fmt2, ";;  v  =  ",
	    root->get_vel()[0], root->get_vel()[1], root->get_vel()[2]);

    st = root->get_dyn_story();
    if (st) put_simple_story_contents(dyn::get_ofp(), *st, ";;");

    for_all_daughters(dyn, root, i)
      fprintf(dyn::get_ofp(), fmt3, i->get_index(),
	      i->get_mass(), i->get_pos()[0], i->get_pos()[1], i->get_pos()[2],
	      i->get_vel()[0], i->get_vel()[1], i->get_vel()[2]);
    putc('\n', dyn::get_ofp());

  } else {

    int oldp = s.precision(p);

    // Special treatment to preserve the root node member data
    // and Log and Dyn stories.

    st = root->get_log_story();
    if (st) put_simple_story_contents(s, *st, ";");

    // Root node member data won't be properly written out by the loop below.
    // Save it as Dyn story data (see get_col above) -- less than elegant.
    // Don't bother with xreal output in this format, so don't print a
    // real_system_time line.

    s << ";;  system_time  =  " <<  root->get_system_time() << endl;
    s << ";;  r  =  " << root->get_pos() << endl;
    s << ";;  v  =  " << root->get_vel() << endl;

    st = root->get_dyn_story();
    if (st) put_simple_story_contents(s, *st, ";;");

    for_all_daughters(dyn, root, i)
      s << i->get_index() << ' ' << i->get_mass() << ' ' << i->get_pos() << ' '
	<< i->get_vel() << '\n';
    s << '\n';

    s.precision(oldp);
  }

}



#else
main(int argc, char** argv)
{
    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.32 $", _SRC_);

    dyn *b;

    while (b = get_dyn()) {
	cout << "TESTING put_dyn:" << endl;
        put_dyn(b);
	cout << "TESTING pp2()   :" << endl;
	pp2(b);
	rmtree(b);
    }
    cerr << "Normal exit\n";
}
#endif
