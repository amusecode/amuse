
       //=======================================================//   _\|/_
      //  __  _____           ___                    ___       //     /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //         _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //           /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //          _\|/_
//=======================================================//            /|\ ~

// Analysis routines for use with smallN.  The only externally visible
// functions are
//
//	bool check_structure(hdyn *bin, real eps2)
//	dyn* get_tree(hdyn *bin, real eps2)
//	void my_sys_stats(dyn *b)
//
// NEED 1. better treatment of stability (Mardling)
//	2. "quarantining" of marginal cases that slip through
//	3. removal of escapers -- simplifies the quarantining logic
//	4. time and size limits
//	5. reporting and storage of details

#include  "smallN2.h"

dyn *flat_copy_tree(dyn *bin)
{
  dyn *b = new dyn(), *bo = NULL, *by = NULL;
  b->set_parent(NULL);
  b->set_name("root");
  b->set_mass(0);
  // b->set_root(b);			// avoid set/get_root() here!

  for_all_leaves(dyn, bin, bb) {
    bo = by;
    by = new dyn();
    if (!b->get_oldest_daughter())
      b->set_oldest_daughter(by);
    if (bo)
      bo->set_younger_sister(by);
    by->set_parent(b);
    by->set_elder_sister(bo);
    by->set_younger_sister(NULL);

    by->set_index(bb->get_index());	// worry about strings later
    by->set_mass(bb->get_mass());
    b->set_mass(b->get_mass()+bb->get_mass());
    by->set_pos(bb->get_pos());
    by->set_vel(bb->get_vel());
    dyn *pp = bb->get_parent();
    while (pp != bin) {
      by->inc_pos(pp->get_pos());
      by->inc_vel(pp->get_vel());
      pp = pp->get_parent();
    }
  }

  // pp2(b);
  // my_sys_stats(b);

  return b;
}

local vec relative_pos(dyn *b, dyn *i)
{
  // Return the position of i relative to its ancestor b.

  vec pos = i->get_pos();
  while (i != b) {
    dyn *p = i->get_parent();
    pos += p->get_pos();
    i = p;
  }
  return pos;
}

local void get_potentials(dyn *bi, dyn *bj, real eps2, real &phi, real &phi2)
{
  // Compute the potential and softened potential of node bi relative
  // to its sister node node bj.

  phi = phi2 = 0;
  dyn *p = bi->get_parent();	// should also be bj->get_parent()

  for_all_leaves(dyn, bi, ii) {
    vec posi = relative_pos(p, ii);
    for_all_leaves(dyn, bj, jj) {
      vec posj = relative_pos(p, jj);
      phi  -= ii->get_mass()*jj->get_mass()/abs(posi-posj);
      phi2 -= ii->get_mass()*jj->get_mass()/sqrt(square(posi-posj)+eps2);
    }
  }
}

local void decompose(dyn *b)
{
  // Construct a binary tree by recursively pairing the top-level
  // particles with the greatest mutual potential.  Use the correct
  // potential in the decomposition, not the 2-body approximation.

  while (b->n_daughters() > 2) {

    // cerr << endl << "start new decompose pass" << endl;
    // pp2(b);

    real potmin = 0;
    dyn *bimin = NULL, *bjmin = NULL;
    for_all_daughters(dyn, b, bi) {
      for (dyn *bj = bi->get_younger_sister();
	   bj != NULL;
	   bj = bj->get_younger_sister()) {
#if 0	
	real pot = -bi->get_mass()*bj->get_mass()
			/abs(bi->get_pos()-bj->get_pos());
#else
	real pot, phi2;
	get_potentials(bi, bj, 0, pot, phi2);
#endif

	// cerr << "potential " << bi->format_label();
	// cerr << " " << bj->format_label() << " " << pot << endl;

	if (pot < potmin) {
	  potmin = pot;
	  bimin = bi;
	  bjmin = bj;
	}
      }
    }
    if (bimin && bjmin) {

      // Make these particles into a binary.

#if 0
      cerr << "decompose: merging " << bimin->format_label();
      cerr << " and " << bjmin->format_label() << " "
	   << abs(bimin->get_pos()-bjmin->get_pos()) << endl;
#endif

      real mtot = bimin->get_mass() + bjmin->get_mass();
      vec cmpos = (bimin->get_mass()*bimin->get_pos()
		    + bjmin->get_mass()*bjmin->get_pos()) / mtot;
      vec cmvel = (bimin->get_mass()*bimin->get_vel()
		    + bjmin->get_mass()*bjmin->get_vel()) / mtot;
      dyn *cm = new dyn();
      cm->set_mass(mtot);
      cm->set_pos(cmpos);
      cm->set_vel(cmvel);

      // Create the binary (bimin,bjmin).

      detach_node_from_general_tree(bimin);
      bimin->set_younger_sister(NULL);
      bimin->set_elder_sister(NULL);

      bimin->inc_pos(-cm->get_pos());
      bjmin->inc_pos(-cm->get_pos());
      bimin->inc_vel(-cm->get_vel());
      bjmin->inc_vel(-cm->get_vel());

      insert_node_into_binary_tree(bimin, bjmin, cm);
      label_binary_node(cm);

    } else
      break;
  }

  // pp2(b);
  // my_sys_stats(b);
}

local void compute_energies(dyn *b, real eps2)
{
  for_all_nodes(dyn, b, bb)
    if (bb->is_parent()) {		// bb may be root, note
      dyn *bi = bb->get_oldest_daughter();
      dyn *bj = bi->get_younger_sister();

      // Note: 2-body energy will be too crude if (e.g.) the
      // components of one sister node are widely separated and the CM
      // happens to lie close to the other sister.  Better to compute
      // the relative potential energy exactly.  Note also that the
      // energies stored are per unit reduced mass.

      real E = 0.5*square(bi->get_vel()-bj->get_vel()), E2 = E;
#if 0
      E  -= bb->get_mass()/abs(bi->get_pos()-bj->get_pos());
      E2 -= bb->get_mass()/sqrt(square(bi->get_pos()-bj->get_pos())+eps2);
#else
      real mu_inv = bb->get_mass()/(bi->get_mass()*bi->get_mass());
      real phi, phi2;
      get_potentials(bi, bj, eps2, phi, phi2);
      E += mu_inv*phi;
      E2 += mu_inv*phi2;
#endif

      putrq(bb->get_dyn_story(), "E", E);
      putrq(bb->get_dyn_story(), "E2", E2);
      // PRC(bb->format_label()); PRC(E); PRL(E2);
    }
}

local inline real separation_squared(dyn *bi, dyn *bj)
{
  return square(bi->get_pos()-bj->get_pos());
}

local inline real size_squared(dyn *b)
{
  if (!b->is_parent()) return 0;
  dyn *od = b->get_oldest_daughter();
  dyn *yd = od->get_younger_sister();
  return separation_squared(od, yd);
 }

local void flag_escapers(dyn *b)
{
  // Components can be flagged as escaping only if all ancestors are
  // unbound.  Use the true energy E here.

  if (b->is_parent()) {
    real E = getrq(b->get_dyn_story(), "E");
    if (E >= 0) {

      // Components are formally unbound, before accepting them as
      // escaping, check that the separation is much greater than the
      // size of any bound component.

      bool esc = true;
      dyn *od = b->get_oldest_daughter();
      dyn *yd = od->get_younger_sister();

      if (od->is_parent() || yd->is_parent()) {
	real sep2 = separation_squared(od, yd);
	for_all_daughters(dyn, b, bb) {
	  if (bb->is_parent()) {
	    if (find_qmatch(bb->get_dyn_story(), "E")) {
	      real EE = getrq(bb->get_dyn_story(), "E");
	      if (EE < 0 && sep2 < 25*size_squared(bb)) {
		esc = false;
		break;
	      }
	    } else {
	      esc = false;
	      break;
	    }
	  }
	}
      }

      if (esc) {
	//	cerr << "flagging escaping components for "
	//	     << b->format_label() << endl;
	putiq(b->get_dyn_story(), "escaping components", 1);
      }

      flag_escapers(od);
      flag_escapers(yd);
    }
  }
}

const real STABLE_FAC = 2.8;

local bool widely_separated(real sma, real sep)
{
  return sep > STABLE_FAC*sma;			// to be improved
}

local real get_semi(dyn *b)
{
  // Return the semi-major-axis of the binary formed by the
  // (presumably bound) components of b.  By construction, we should
  // have an energy already stored in the dyn story.

  real semi = -1;

  if (b->is_parent()) {
    if (find_qmatch(b->get_dyn_story(), "E")) {
      real E = getrq(b->get_dyn_story(), "E");
      if (E < 0) semi = -0.5*b->get_mass()/E;
    } else {

      // Looks like we have to compute the energy from scratch.  We
      // really do use the 2-body expression here.

      dyn *od = b->get_oldest_daughter();
      dyn *yd = od->get_younger_sister();

      real E = 0.5*square(od->get_vel()-yd->get_vel())
		- b->get_mass()/abs(od->get_pos()-yd->get_pos());
      if (E < 0) semi = -0.5*b->get_mass()/E;
    }
  }

  return semi;
}

local bool check_stability(dyn *b)
{
  // Stablility criterion for a CM node is that both daughters are
  // stable and sufficiently widely separated.  We currently don't
  // include external perturbations in the decision (although we
  // probably should...).

  bool stable = false;

  if (b->is_leaf())				// single stars are stable

    stable = true;

  else {

    real E = getrq(b->get_dyn_story(), "E");    // NB E not E2

    if (E < 0) {
      dyn *od = b->get_oldest_daughter();
      dyn *yd = od->get_younger_sister();

      if (od->is_leaf() && yd->is_leaf())	// binaries are stable too
	stable = true;
      else {

	// Check both daughters (don't && these tests in the if
	// statement because this could skip the second check).

	bool od_stab = check_stability(od);
	bool yd_stab = check_stability(yd);
	if (od_stab && yd_stab) {

	  // Identify an "inner" binary.  Presumably the only
	  // multiples under b are stable.  However, both components
	  // could be multiples, so choose the larger one.  This
	  // stability criterion should be improved.

	  real a_in = max(get_semi(od), get_semi(yd));
	  real r_out = abs(od->get_pos() - yd->get_pos());
	  if (widely_separated(a_in, r_out)) {

	    // Apply the same test to the periastron of the relative orbit.

	    kepler k;
	    initialize_kepler_from_dyn_pair(k, od, yd);
	    stable = widely_separated(a_in, k.get_periastron());
	    putrq(b->get_dyn_story(), "separation ratio",
		  k.get_periastron()/a_in);
// 	    real a_out = k.get_semi_major_axis();
// 	    real e_out = k.get_eccentricity();
// 	    real peri = k.get_periastron();
	  }
	}
      }
    }

    // If we want to apply a perturbation condition, do it here...

  }
  if (stable && b->is_parent()) putiq(b->get_dyn_story(), "stable", 1);
  return stable;
}

local void print_recursive(dyn *b, int indent = 0)
{
  for (int i = 0; i < indent; i++) cerr << " ";
  cerr << b->format_label() << ":";
  int esc = 0;
  if (b->is_parent()) {
    if (b->n_leaves() > 2) cerr << " multiple";
    if ((esc = getiq(b->get_dyn_story(), "escaping components") == 1))
      cerr << " escaping components";
    real E;
    if ((E = getrq(b->get_dyn_story(), "E")) < 0) {
      cerr << " bound components " << E;
      real sep = getrq(b->get_dyn_story(), "separation ratio");
      if (sep > 0) cerr << " (" << sep << ")";
    } else if (!esc)
      cerr << " unbound components " << E;
    if (getiq(b->get_dyn_story(), "stable") == 1 && b->n_leaves() > 2)
      cerr << " stable multiple";
  } else {
    cerr << " single";
  }
  if (b->get_parent()
      && getiq(b->get_parent()->get_dyn_story(), "escaping components") == 1)
    cerr << " escaper";
  cerr << endl;

  for_all_daughters(dyn, b, bb) print_recursive(bb, indent+2);
}

local void set_string_node_name(dyn *b, char *s)
{
  if (b->get_name())
    strcpy(s, b->get_name());		// NB no size check...
  else if (b->get_index() >= 0)
    sprintf(s, "%d", b->get_index());
  else
    strcpy(s, "-");
}

local char *construct_binary_label(dyn *bi, dyn *bj, char b[])
{
  // Based on node/util/node_tt.C...
  // ASSUMES that indices are set (i.e. leaf names are all numbers).

  static char new_name[256];	// shared by all function calls
  char si[128], sj[128];

  set_string_node_name(bi, si);
  set_string_node_name(bj, sj);

  sprintf(new_name,"%c%s,%s%c", b[0], si, sj, b[1]);
  return new_name;
}

local void create_label(dyn *b)
{
  // Construct strings describing the dynamical state of the system.
  // We will use the usual notation for multiples, except that the
  // choice of symbols will reflect nore accurately the state of each
  // node:
  //
  //	(a,b)	a and b are bound, but this node is not necessarily stable
  //	[a,b]	a and b are bound, and this node is flagged as stable
  //	<a,b>	a and b are unbound
  //
  // Note that we use the unsoftened energy E in determining the
  // structure.

  for_all_daughters(dyn, b, bb)
    create_label(bb);

  if (b->is_parent()) {
    dyn *od = b->get_oldest_daughter();
    dyn *yd = od->get_younger_sister();
    if (getrq(b->get_dyn_story(), "E") >= 0)
      b->set_label(construct_binary_label(od, yd, "<>"));
    else {
      if (getiq(b->get_dyn_story(), "stable") == 1
	  || b->n_leaves() <= 2)
	b->set_label(construct_binary_label(od, yd, "[]"));
      else
	b->set_label(construct_binary_label(od, yd, "()"));
    }
  }
}

local void print_summary_string(dyn *b)
{
  create_label(b);
  cerr << "summary: " << b->get_name() << endl;
}

local bool is_over(dyn *b)
{
  // "Over" means that b consists only of unbound stable objects.
  // This means that, as we descend the tree, an object whose
  // components aren't escaping has to be stable.

  if (b->is_leaf()) return true;
  if (b->n_leaves() <= 2 && getrq(b->get_dyn_story(), "E") < 0) return true;
  if (getiq(b->get_dyn_story(), "stable") == 1) return true;
  if (getiq(b->get_dyn_story(), "escaping components") != 1) return false;

  // Not stable, but components are escaping.  Check the components.

  dyn *od = b->get_oldest_daughter();
  dyn *yd = od->get_younger_sister();

  return is_over(od) && is_over(yd);
}



// Externally visible functions:

bool check_structure(hdyn *bin,
		     real eps2,			// default = 0
		     bool verbose)		// default = true
{
  if (verbose)
    cerr << endl << "checking structure at time "
	 << bin->get_time() << endl;

  // Simple approach: start by copying the tree to a flat dyn tree;
  // then identify structure by proximity (actually potential) and
  // dynamics (bound or not).  Then classify bound systems by their
  // stability properties.  At present, we compute and store the true
  // and softened energies, but we only use the true energy in out
  // classification.  The softened energy has meaning only in relating
  // the multiple to the larger system in which it resides.

  dyn *b = flat_copy_tree(bin);
  decompose(b);

  // pp2(b);
  // my_sys_stats(b);

  // First determine the binding energy of each node.  Note that the
  // order of for_all_nodes() means that this function and the next
  // could probably be merged, but don't assume that here.

  compute_energies(b, eps2);

  // Then descend the tree to look for escapers (simply unbound, since
  // that is what the dynamical code will care about).  Note that we
  // use the softening from the external program in making this
  // decision.

  flag_escapers(b);

  // Then ascend the tree looking for widely separated binaries.  For
  // now, just use a simple distance criterion.  Eventually, we will
  // use the Mardling expression.

  check_stability(b);

  // Print the results.

  if (verbose) print_recursive(b);

  // Construct and print a summary string describing the state of the stsrem.

  if (verbose) print_summary_string(b);

  bool over = false;

  if (is_over(b)) {
    if (verbose) cerr << endl << "interaction is over" << endl;
    over = true;
  }

  // Clean up and exit.

  rmtree(b);
  return over;
}

dyn *get_tree(hdyn *bin,
	      real eps2,			// default = 0
	      bool verbose)			// default = false
{
  // Same as check structure, but return the annotated tree
  // and don't print anything.

  dyn *b = flat_copy_tree(bin);
  decompose(b);
  compute_energies(b, eps2);
  flag_escapers(b);
  check_stability(b);
  if (verbose) print_summary_string(b);
  return b;
}

void my_sys_stats(dyn *b,
		  int verbose)			// default = 1
{
  if (verbose > 1) verbose -= 1;		// so 1 and 2 both produce
						// 1 for sys_stats purposes
  sys_stats(b,
	    0.1,				// energy cutoff
	    verbose,				// verbosity
	    true,				// binaries
	    true,				// long binary output
	    2,					// which_lagr
	    false,				// print_time
	    true,				// compute_energy
	    true);				// allow_n_sq_ops
  cout << flush;
  cerr << flush;
}
