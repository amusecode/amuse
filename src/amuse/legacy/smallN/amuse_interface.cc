#include "amuse_interface.h"
#include <vector>
#include <algorithm>

// Externally visible interface functions:
//
//	report_multiples(int level)
//      set_rmax_system(int rm)
//	get_rmax_system()
//	set_name_ranges(int nstart, int delta_n)
//	add_binary(int id, int id1, int id2, real mass1, real mass2,
//		   real period, real eccentricity)
//	remove_multiple(int id)
//	get_energy(int id)
//	get_total_energy()
//	get_top_level_energy_error()
//	get_mass(int id)
//	get_radius(int id)
//	get_n(int id)
//	int integrate_multiple(real eps2 = 0);
//
// Functions to load/retrieve data on an individual interaction.
//
//	void clear_multiple()
//	void add_to_interaction(int i, real m, real x[3], real v[3])
//	int get_status(int i)
//	void get_particle(int k, int &i, int &isn, real &m,
//			  real x[3], real v[3])

// Internal data management:

static vector<int> ident;	// multiple IDs
static vector<dyn*> msys;	// multiple dyn trees
				// note: store dyns but integrate with hdyns

static int nstart=100, delta_n=100, n_next[4]={100,200,300,400};

static real rmax_system = VERY_LARGE_NUMBER;
static real top_level_energy_error = 0;

void set_rmax_system(real rm)
{
  rmax_system = rm;
}

real get_rmax_system()
{
  return rmax_system;
}

local real get_energy(dyn *b, real eps2 = 0)
{
  if (!b || !b->is_valid() || !b->is_parent()) return 0;

  dyn *d = flat_copy_tree(b);
  real kinetic = 0, potential = 0;
  for_all_daughters(dyn, d, di) {
    kinetic += di->get_mass()*square(di->get_vel());
    real dpot = 0;
    for (dyn *dj = di->get_younger_sister();
	 dj != NULL; dj = dj->get_younger_sister())
      dpot -= dj->get_mass()/sqrt(square(di->get_pos()-dj->get_pos())+eps2);
    potential += di->get_mass()*dpot;
  }
  kinetic /= 2;

  rmtree(d);
  return kinetic+potential;
}

local real get_energy_by_index(int i)
{
  if (i < (int)ident.size())
    return get_energy(msys[i]);
  else
    return 0;
}

void report_multiples(int level)
{
  // Return summary information about the multiples in the system.

  cout << "report_multiples: " << ident.size() << " multiple";
  if (ident.size() > 1) cout << "s";
  cout << " in the system" << endl;

  cout << "id(n): ";
  for (unsigned int i = 0; i < ident.size(); i++) {
    cout << ident[i] << "(";
    if (msys[i] && msys[i]->is_valid())
      cout << msys[i]->n_leaves();
    else
      cout << "invalid";
    cout << ") ";
  }
  cout << endl;

  if (level > 0) {
    for (unsigned int i = 0; i < ident.size(); i++) {
      cout << "multiple " << i << ", id = " << ident[i]
	   << "  E = " << get_energy_by_index(i);
      if (level < 2)
	cout << "  root name = " << msys[i]->format_label() << endl;
      else {
	cout << endl;
	pp2(msys[i]);
      }
    }
  }

  // Consistency checks: no duplication of IDs; leaves only appear
  // once.

  for (unsigned int i = 0; i < ident.size(); i++)
    for (unsigned int j = i+1; j < ident.size(); j++)
      if (ident[i] == ident[j])
	cout << "warning: duplicate index " << ident[i] << endl;

  for (unsigned int i = 0; i < ident.size(); i++) {
    vector<int>ileaves;
    for_all_leaves(dyn, msys[i], di) ileaves.push_back(di->get_index());
    for (unsigned int j = i+1; j < ident.size(); j++) {
      vector<int>jleaves;
      for_all_leaves(dyn, msys[j], dj) jleaves.push_back(dj->get_index());

      // Compare leaves of i and j.  There should be no overlap.

      for (int ki = 0; ki < (int)ileaves.size(); ki++)
	for (int kj = 0; kj < (int)ileaves.size(); kj++)
	  if (ileaves[ki] == jleaves[kj])
	    cout << "warning: duplicate leaf " << ileaves[ki]
		 << " in " << ident[i] << " and " << ident[j] << endl;
    }
  }
  cout << flush;
  cerr << flush;
}

// Need a naming scheme for new multiples, but the outside program
// must be able to influence it.  For now, adopt the following simple
// scheme: binaries start at nstart, triples at nstart + delta_n,
// quadruples at nstart+2*delta_n, etc.  We do *not* reuse names.

void set_name_ranges(int ns, int del)
{
  nstart = ns;
  delta_n = del;

  // Keep track of names to avoid having to search the list.
  // Categories are binaries, triples, quadruples, and "higher."

  for (int k = 0; k < 4; k++)
    n_next[k] = ns + k*delta_n;
}

local int get_new_name(int n)
{
  if (n < 2) return -1;
  if (n > 5) n = 5;
  return n_next[n-2]++;
}

// The only way to make anything other than a binary is to do it
// dynamically, i.e. within this module, so for now there seems to be
// no particular reason for anything more complex.

vec random_unit_vec()
{
  real costh = randinter(-1,1), phi = randinter(0, 2*M_PI);
  real sinth = sqrt(1-costh*costh);
  if (randinter(-1,1) > 0) sinth = -sinth;
  return vec(sinth*cos(phi), sinth*sin(phi), costh);
}

real semi_from_period(real m1, real m2, real period)
{
  return pow((m1+m2)*period*period/(4*M_PI*M_PI), 1./3);
}

int add_binary(int id1, int id2, real mass1, real mass2,
	       real period, real eccentricity)
{
  // Create a new dyn tree and store it.  Arguments follow those in
  // the stellar package.  Phase and orientation are chosen
  // randomly.  Time and center of mass position and velocity will
  // be derived from the dynamical module when this binary is
  // selected for integration.
  //
  // Return value is the new id of the binary (not the total number of
  // objects).

  real mtot = mass1+mass2;

  dyn *cm = new dyn();
  dyn *od = new dyn();
  dyn *yd = new dyn();

  od->set_parent(cm);
  od->set_younger_sister(yd);
  od->set_index(id1);
  od->set_mass(mass1);

  yd->set_parent(cm);
  yd->set_elder_sister(od);
  yd->set_index(id2);
  yd->set_mass(mass2);

  cm->set_oldest_daughter(od);
  cm->set_label(construct_binary_label(od, yd));
  cm->set_mass(mtot);

  kepler k;
  k.set_time(0);
  k.set_total_mass(mtot);

  real semi = semi_from_period(mass1, mass2, period);
  k.set_semi_major_axis(semi);
  k.set_eccentricity(eccentricity);
  k.set_mean_anomaly(randinter(-M_PI, M_PI));
  vec l = random_unit_vec();
  vec t = l^random_unit_vec();
  t /= abs(t);
  vec n = l^t;
  k.set_orientation(l, t, n);
  k.initialize_from_shape_and_phase();	// expects a, e [, q [, E]]

  od->set_pos(-mass2*k.get_rel_pos()/mtot);
  yd->set_pos(mass1*k.get_rel_pos()/mtot);
  od->set_vel(-mass2*k.get_rel_vel()/mtot);
  yd->set_vel(mass1*k.get_rel_vel()/mtot);

  int id = get_new_name(2);
  ident.push_back(id);
  msys.push_back(cm);

  return id;
}

local int id_to_index(int id)
{
  return (int)(find(ident.begin(), ident.end(), id) - ident.begin());
}

// Removal may be prompted by binary evolution.

local int remove_multiple_by_index(int i)
{
  // Remove the specified dyn tree.

  if (i < (int)ident.size()) {
    ident.erase(ident.begin()+i);
    if (msys[i]) rmtree(msys[i]);
    msys.erase(msys.begin()+i);
  }
  return (int)ident.size();
}

int remove_multiple(int id)
{
  // Locate and remove the specified dyn tree.

  int i = id_to_index(id);
  return remove_multiple_by_index(i);
}

real get_energy(int id)
{
  int i = id_to_index(id);
  return get_energy_by_index(i);
}

real get_total_energy()
{
  real total_energy = 0;
  for (int i = 0; i < (int)ident.size(); i++)
    total_energy += get_energy_by_index(i);
  return total_energy;
}

real get_top_level_energy_error()
{
  return top_level_energy_error;
}

real get_mass(int id)
{
  // Locate id, return its total mass.

  int i = id_to_index(id);

  if (i < (int)ident.size())
    return msys[i]->get_mass();
  else
    return -1;
}

#define FAC 3

local real get_effective_radius(dyn *b)
{
  if (!b->is_parent()) return 0;

  // Return an effective (interaction) radius for object b.  This is
  // FAC times the binary semi-major axis of the top-level pair.
  // Stored multiples are always saved in hierarchical form.

  dyn *od = b->get_oldest_daughter();
  dyn *yd = od->get_younger_sister();

  real mass = od->get_mass() + yd->get_mass();
  real energy = 0.5*square(od->get_vel()-yd->get_vel())
			- mass/abs(od->get_pos()-yd->get_pos());

  if (energy < 0)
    return -0.5*FAC*mass/energy;
  else
    return 0;
}

real get_radius(int id)
{
  // Locate id, return its effective radius.

  int i = id_to_index(id);

  if (i < (int)ident.size())
    return get_effective_radius(msys[i]);
  else
    return -1;
}

real get_n(int id)
{
  // Locate id, return the number of stars in it.

  int i = id_to_index(id);

  if (i < (int)ident.size())
    return msys[i]->n_leaves();
  else
    return -1;
}

bool is_multiple(int id)
{
  return (id_to_index(id) < (int)ident.size());
}


local void duplicate_daughters(dyn *d, hdyn *h)
{
  hdyn *hh = NULL, *prevhh = NULL;
  for_all_daughters(dyn, d, dd) {
    hh = new hdyn();
    if (!h->get_oldest_daughter())
      h->set_oldest_daughter(hh);
    hh->set_parent(h);
    hh->set_elder_sister(prevhh);
    if (prevhh)
      prevhh->set_younger_sister(hh);
    duplicate_daughters(dd, hh);
    prevhh = hh;
  }
}

local void copy_data(dyn *d, dyn *h)
{
  h->set_mass(d->get_mass());
  h->set_pos(d->get_pos());
  h->set_vel(d->get_vel());
  h->set_index(d->get_index());
  h->set_name(d->get_name());
}

local void copy_dyn_to_hdyn(dyn *d, hdyn *h)
{
  // Make an hdyn copy h of a dyn tree d.  On entry, h points to an
  // empty node.  Note similarities to functions in analyze.cc and
  // kira_smallN.C...

  // First create parallel tree structures.

  duplicate_daughters(d, h);

  // Then copy node data from d to h.

  hdyn *hh = h;
  for_all_nodes(dyn, d, dd) {
    copy_data(dd, hh);
    hh = (hdyn*)hh->next_node(h);
  }
}

local hdyn *create_multiple(vector<int> id, vector<real> mass,
			    vector<vec> pos, vector<vec> vel,
			    real &mtot, vec &cmpos, vec &cmvel)
{
  // Create a multiple system from the input data.

  hdyn *b = new hdyn();
  b->set_name("root");
  b->set_root(b);

  int nin = id.size();
  hdyn *last = NULL;

  for (int k = 0; k < nin; k++) {
    hdyn *bb = new hdyn();

    int i = id_to_index(id[k]);
    if (i < (int)ident.size()) copy_dyn_to_hdyn(msys[i], bb);

    bb->set_index(id[k]);	// always store the id in the node.
    bb->set_mass(mass[k]);	// should check that bb->get_mass() = mass[k]
    bb->set_pos(pos[k]);
    bb->set_vel(vel[k]);

    mtot += mass[k];
    cmpos += mass[k]*pos[k];
    cmvel += mass[k]*vel[k];

    // Add bb to the tree.

    bb->set_parent(b);
    bb->set_elder_sister(last);
    if (last)
      last->set_younger_sister(bb);
    else
      b->set_oldest_daughter(bb);
    last = bb;
  }

  b->set_mass(mtot);
  cmpos /= mtot;
  cmvel /= mtot;

  // Move to the center of mass frame.

  for_all_daughters(hdyn, b, bb) {
    bb->inc_pos(-cmpos);
    bb->inc_vel(-cmvel);
  }

  return b;
}

local void set_top_level_escapers(dyn *d, bool verbose = true)
{
  // Restructure d into a flat array of stable objects.  This function
  // should probably be merged into analyze.cc.

  // Currently d is a binary tree, but don't assume this here, as the
  // tree construction may already have assigned escapers to the top
  // level.  (We do assume binarity at lower levels.)

  // Every node in d should either be flagged as an escaper or be
  // stable.  Systematically check each top-level CM node, moving its
  // components to the top level (replacing the CM) if they are
  // flagged as escapers.

  // The top-level nodes should already be flagged as escaping.  If
  // they are not, then we must have created (been given?) a
  // hierarchically bound system.  This shouldn't happen, but we
  // should be able to cope with it...

  if (getiq(d->get_dyn_story(), "escaping components") != 1) {
    if (getiq(d->get_dyn_story(), "stable") == 1) {
      if (verbose)
	cout << "top-level node " << d->format_label()
	     << " is stable" << endl << flush;

      // This can only happen if the interaction was bound to begin
      // with (which is not impossible).  Install the entire system
      // under a new single top-level node of d.

      dyn *dnew = new dyn;
      dnew->set_index(99999);	// get_new_name(d->n_leaves()));
      dnew->set_parent(d);
      dnew->set_oldest_daughter(d->get_oldest_daughter());
      for_all_daughters(dyn, d, dd) dd->set_parent(dnew);
      d->set_oldest_daughter(dnew);
      dnew->set_name(d->get_name());
      dnew->set_mass(d->get_mass());
      dnew->set_pos(d->get_pos());
      dnew->set_vel(d->get_vel());

      return;

    } else if (verbose)
      cout << "warning: top-level nodes of " << d->format_label()
	   << " are not escaping" << endl << flush;
  }

  dyn *dd = d->get_oldest_daughter();
  while (dd) {
    dyn *next = dd->get_younger_sister();

    if (dd->is_parent()) {

      if (getiq(dd->get_dyn_story(), "escaping components") == 1) {

	// Components are escaping; promote them to the top level.
	// The first replaces dd; the second is inserted as dd's
	// younger sister.

	dyn *os = dd->get_elder_sister();
	dyn *ys = dd->get_younger_sister();
	dyn *od = dd->get_oldest_daughter();
	dyn *yd = od->get_younger_sister();

	od->inc_pos(dd->get_pos());
	od->inc_vel(dd->get_vel());
	yd->inc_pos(dd->get_pos());
	yd->inc_vel(dd->get_vel());

	if (!os)
	  d->set_oldest_daughter(od);
	else
	  os->set_younger_sister(od);

	od->set_elder_sister(os);
	od->set_parent(d);

	if (ys)
	  ys->set_elder_sister(yd);
	yd->set_younger_sister(ys);
	yd->set_parent(d);

	dd->set_oldest_daughter(NULL);
	delete dd;

	next = od;

      } else {

	// Components are bound.  Node should already have been
	// flagged as stable...

	if (dd->n_leaves() > 2
	    && getiq(dd->get_dyn_story(), "stable") != 1)
	  cout << "warning: top-level bound node " << dd->format_label()
	       << " is not stable" << endl;
      }
    }
	
    dd = next;
  }
  cout << flush;
  cerr << flush;
}

local void split_wide_binaries(dyn *d, real rlimit, bool verbose)
{
  // Recursively disassemble binaries whose effective sizes are
  // greater than the specified limit.

  dyn *dd = d->get_oldest_daughter();
  while (dd) {
    real r_eff = get_effective_radius(dd);
    if (r_eff > rlimit) {

      // Replace dd by its first component, and insert the other
      // component into the top-level tree.

      if (verbose) {
	cout << "splitting " << dd->format_label()
	     << " into components, " << flush;
	PRC(r_eff); PRL(rlimit);
      }

      dyn *od = dd->get_oldest_daughter();
      dyn *yd = od->get_younger_sister();

      od->set_parent(d);
      yd->set_parent(d);

      od->inc_pos(dd->get_pos());
      od->inc_vel(dd->get_vel());
      yd->inc_pos(dd->get_pos());
      yd->inc_vel(dd->get_vel());

      dyn *os = dd->get_elder_sister();
      od->set_elder_sister(os);
      if (os) os->set_younger_sister(od);

      dyn *ys = dd->get_younger_sister();
      yd->set_younger_sister(ys);
      if (ys) ys->set_elder_sister(yd);

      if (d->get_oldest_daughter() == dd) d->set_oldest_daughter(od);
      dd->set_oldest_daughter(NULL);
      delete dd;

      dd = od;

    } else
      dd = dd->get_younger_sister();
  }    
}

local bool have_same_components(dyn *d1, dyn *d2)
{
  // Return true iff d1 and d2 have the same (possibly permuted)
  // components.

  bool same = false;
  if (d1->n_leaves() == d2->n_leaves()) {
    vector<int> init, final;
    for_all_leaves(dyn, d1, dd) init.push_back(dd->get_index());
    for_all_leaves(dyn, d2, dd) final.push_back(dd->get_index());
    sort(init.begin(), init.end());
    sort(final.begin(), final.end());
    same = (init == final);
  }
  return same;
}

local void label_new_nodes(vector<int> idlist, dyn *d,
			   bool verbose = true)
{
  // Label "new" top-level multiples not on the original list idlist,
  // and make sure that surviving multiples have the proper label.  An
  // initial multiple is regarded as surviving (and its ID preserved),
  // if all of its components are found in a multiple of the same
  // order in the final system (even if the components have been
  // reorganized).  Also optionally print a list of new nodes as we
  // go.

  if (verbose) cout << "new nodes: ";
  int n_new = 0;

  for_all_daughters(dyn, d, dd) {
    if (dd->is_parent()) {

      // This is a multiple.  Find its first component.

      int first = 0;
      for_all_leaves(dyn, dd, ddd) {first = ddd->get_index(); break;}

      // Find the initial multiple containing first and compare its
      // content with that of dd.

      bool is_new = true;

      for (int k = 0; k < (int)idlist.size(); k++) {
	int i = id_to_index(idlist[k]);
	if (i < (int)ident.size())
	  if (node_contains(msys[i], first)) {

	    if (have_same_components(msys[i], dd)) {

	      // This is an old multiple.  Make sure the name is
	      // consistent.

	      dd->set_index(ident[i]);
	      is_new = false;
	    }
	  }
      }

      if (is_new) {

	// We have a new multiple.  Give it a name and save it.

	int id = get_new_name(dd->n_leaves());
	dd->set_index(id);
	if (verbose) cout << id << " ";
	n_new++;
      }

    } else {

      // This is a single star.  Is it a new top-level node?
      // Note that we only care if verbose = true.

      if (verbose) {
	bool new_node = true;
	for (int k = 0; k < (int)idlist.size(); k++)
	  if (idlist[k] == dd->get_index()) {
	    new_node = false;
	    break;
	  }

	if (new_node) {
	  cout << dd->get_index() << " ";
	  n_new++;
	}
      }
    }
  }  
  if (verbose) {
    if (n_new == 0) cout << "(none)";
    cout << endl;
  }
  cout << flush;
  cerr << flush;
}

local void flag_status(vector<int> idlist, dyn *d, vector<int> &status,
		       bool verbose = true)
{
  // Compare the old (idlist) and new (d) nodes in this interaction,
  // and construct a status array for the original nodes.  An object
  // on the old list will be flagged as surviving (status = 1) if it
  // is now a top-level node in d, as indicated by get_index().
  // Otherwise, it is flagged for removal from the external
  // integration (status = 0).

  // Code here is inefficient, but the arrays are very small.

  int nin = idlist.size();
  status.clear();
  for (int k = 0; k < nin; k++)
    status.push_back(0);

  for (int k = 0; k < nin; k++)
    for_all_daughters(dyn, d, dd)
      if (dd->get_index() == idlist[k]) {
	status[k] = 1;
	break;
      }

  if (verbose) {
    cout << "status: ";
    for (int k = 0; k < nin; k++) cout << idlist[k] << ":" << status[k] << " ";
    cout << endl;
  }
  cout << flush;
  cerr << flush;
}

local void remove_old_multiples(vector<int> idlist, vector<int> status)
{
  // Remove *all* of the original multiple nodes from storage (we will
  // replace the survivors in a moment).

  for (int k = 0; k < (int)idlist.size(); k++) {
    int i = id_to_index(idlist[k]);
    if (i < (int)ident.size()) {

      // This is a multiple.  Remove it completely if status = 0,
      // otherwise delete the tree and leave a stub.

      if (status[k] == 0)
	remove_multiple_by_index(i);
      else {
	rmtree(msys[i]);
	msys[i] = NULL;
      }
    }
  }
}

local real store_new_multiples(dyn *d)
{
  // Move all multiples from the new list into storage.  As each node
  // is processed, it is removed from the d tree and replaced by its
  // center of mass, so that, on return, d is a flat tree containing
  // only top-level information, and all nodes mave been flagged and
  // removed or saved.  Return value is the effective radius of the
  // largest node.

  real max_size = 0;

  for_all_daughters(dyn, d, dd) {
    if (dd->is_parent()) {
      dyn *cm = new dyn();

      // Sever dd's links and insert cm instead.

      cm->set_parent(d);
      dd->set_parent(NULL);
      if (d->get_oldest_daughter() == dd)
	d->set_oldest_daughter(cm);

      dyn *os = dd->get_elder_sister();
      dyn *ys = dd->get_younger_sister();
      cm->set_elder_sister(os);
      if (os) os->set_younger_sister(cm);
      cm->set_younger_sister(ys);
      if (ys) ys->set_elder_sister(cm);
      dd->set_elder_sister(NULL);
      dd->set_younger_sister(NULL);

      // Copy all data and clear dd's pos and vel.

      copy_data(dd, cm);
      dd->set_pos(0);
      dd->set_vel(0);

      // Store dd.

      int i;
      if ((i=id_to_index(dd->get_index())) < (int)ident.size()) {

	// The name is already on the list.  The msys pointer should be
	// NULL, but check just in case.

	if (msys[i]) rmtree(msys[i]);
	msys[i] = dd;

      } else {

	// Add the new multiple to the list.

	ident.push_back(dd->get_index());
	msys.push_back(dd);

      }

      real size = get_effective_radius(dd);
      if (size > max_size) max_size = size;

      dd = cm;
    }
  }

  return max_size;
}


// Avoid the need for a large number of variable-length arguments to
// integrate_multiple by providing direct access to a set of internal
// temporary arrays.

// Input arrays: the initial multiple system.  Used to create the
// interaction sent to integrate().

static vector<int> input_id;
static vector<real> input_mass;
static vector<vec> input_pos;
static vector<vec> input_vel;

static vector<int> status;	// the fate of the original CMs

// Output arrays: the final multiple system.  Unclear what the
// best/most convenient quantities are.  For now, the data are
// somewhat redundant.  To be revisited once a working driver
// application is in place.

static vector<int> new_id;
static vector<real> new_mass;
static vector<vec> new_pos;
static vector<vec> new_vel;
static vector<int> is_new;	// 1 iff this is a new CM

// Accessors:

void clear_multiple()
{
  input_id.clear();
  input_mass.clear();
  input_pos.clear();
  input_vel.clear();

  // Unnecessary, but...

  status.clear();
  new_id.clear();
  new_mass.clear();
  new_pos.clear();
  new_vel.clear();
  is_new.clear();
}

void add_to_interaction(int i, real m, real x[3], real v[3])
{
  input_id.push_back(i);
  input_mass.push_back(m);
  input_pos.push_back(vec(x[0],x[1],x[2]));
  input_vel.push_back(vec(v[0],v[1],v[2]));
}

void add_to_interaction(int i, real m, real x, real y, real z,
			real vx, real vy, real vz)
{
  input_id.push_back(i);
  input_mass.push_back(m);
  input_pos.push_back(vec(x,y,z));
  input_vel.push_back(vec(vx,vy,vz));
}

int get_status(int i)		// i is an id, not an index
{
  for (int k = 0; k < (int)input_id.size(); k++)
    if (input_id[k] == i) return status[k];
  return -1;
}

local void list_final(dyn *d)
{
  // Save the final top-level structure in accessible form.  Restore
  // the center of mass position and velocity in new_pos and new_vel.

  new_id.clear();
  new_mass.clear();
  new_pos.clear();
  new_vel.clear();
  is_new.clear();

  for_all_daughters(dyn, d, dd) {
    new_id.push_back(dd->get_index());
    new_mass.push_back(dd->get_mass());
    new_pos.push_back(dd->get_pos());
    new_vel.push_back(dd->get_vel());
    int isn = 1;
    for (int k = 0; k < (int)input_id.size(); k++)
      if (dd->get_index() == input_id[k]) {
	isn = 0;
	break;
      }
    is_new.push_back(isn);
  }
}

void get_particle_result(int k, int *id, real *mass, real *x, real *y, real *z, real *vx, real *vy, real *vz)
{
  *id = new_id[k];
  *mass = new_mass[k];
  *x = new_pos[k][0];
  *y = new_pos[k][1];
  *z = new_pos[k][2];
  *vx = new_vel[k][0];
  *vy = new_vel[k][1];
  *vz = new_vel[k][2];
}

void get_particle_original(int k, int *id, real *mass, real *x, real *y, real *z, real *vx, real *vy, real *vz)
{
  *id = input_id[k];
  *mass = input_mass[k];
  *x = input_pos[k][0];
  *y = input_pos[k][1];
  *z = input_pos[k][2];
  *vx = input_vel[k][0];
  *vy = input_vel[k][1];
  *vz = input_vel[k][2];
}



int is_new_particle(int k)
{
  return is_new[k];
}

local real local_potential(dyn *d, dyn *dd, real eps2 = 0)
{
  real pot = 0;
  for_all_daughters(dyn, d, ddd)
    if (ddd != dd)
      pot -= ddd->get_mass()
		/ sqrt(square(ddd->get_pos()-dd->get_pos())+eps2);
  return pot;
}

local void expand_nodes(dyn *d, real rscale,
			bool verbose = false,
			real eps2 = 0)
{
  // Attempt to expand the top-level node system so that the rms
  // separation is rscale, but preserving the top-level energy.  Do
  // nothing if the separation already exceeds rscale.  If the nodes
  // are bound, there may be a limit on how much expansion is
  // possible.  Should also check for (and correct) any accidental
  // close pairs created, but not done yet.

  // Determine the rms separation (no mass weighting).

  real sep2 = 0;
  int np = 0;
  for_all_daughters(dyn, d, dd)
    for (dyn *ddd = dd->get_younger_sister();
	 ddd != NULL; ddd = ddd->get_younger_sister()) {
      sep2 += square(dd->get_pos()-ddd->get_pos());
      np++;
    }

  real sep = sqrt(sep2/np);
  bool expand = false;

  // Proceed if there is something to expand.

  if (sep < rscale) {

    // Start by recomputing the center of mass pos and vel (should
    // both be zero).

    int n = 0;
    real mass = 0;
    vec cmpos = 0, cmvel = 0;

    for_all_daughters(dyn, d, dd) {
      mass += dd->get_mass();
      cmpos += dd->get_mass()*dd->get_pos();
      cmvel += dd->get_mass()*dd->get_vel();
      n++;
    }
    cmpos /= mass;
    cmvel /= mass;

    // Move to the center of mass frame, but save the original values,
    // just in case they aren't zero for some application.

    vec cmpos_init = cmpos;
    vec cmvel_init = cmvel;

    for_all_daughters(dyn, d, dd) {
      dd->inc_pos(-cmpos_init);
      dd->inc_vel(-cmvel_init);
    }

    // Compute the initial energies.

    real ke0 = 0, pe0 = 0;
    real pot[n], energy[n];
    n = 0;
    bool bound = false;
    for_all_daughters(dyn, d, dd) {
      pot[n] = local_potential(d, dd, eps2);
      real v2 = 0.5*square(dd->get_vel());
      energy[n] = v2 + pot[n];
      if (energy[n] < 0) bound = true;
      ke0 += dd->get_mass() * v2;
      pe0 += dd->get_mass() * pot[n];
      n++;
    }

    pe0 /= 2;
    real e0 = ke0 + pe0;

    if (bound) {

      // Select a smaller rscale based on the bound energies.

      n = 0;
      for_all_daughters(dyn, d, dd) {
	if (energy[n] < 0)
	  rscale = min(rscale, -0.5*mass/energy[n]);
	n++;
      }
    }

    if (sep < rscale) {

      // Scale the positions.

      expand = true;
      real rfac = rscale/sep;

      for_all_daughters(dyn, d, dd)
	dd->set_pos(rfac*dd->get_pos());

      // Position scaling will keep the center of mass at 0, but the
      // following velocity scaling in general won't.  Recompute and
      // correct the center of mass velocity after scaling the
      // velocities.  (DOESN'T PRESERVE ANGULAR MOMENTUM -- *** TO
      // COME *** USE THE KEPLER PACKAGE TO ACCOMPLISH THIS IN THE
      // 2-BODY CASE ***).

      cmvel = 0;
      n = 0;
      for_all_daughters(dyn, d, dd) {
	real p = local_potential(d, dd, eps2);
	real dpot = p - pot[n];
	pot[n] = p;
	real vfac = sqrt(1-dpot/square(dd->get_vel()));
	dd->set_vel(vfac*dd->get_vel());
	cmvel += dd->get_mass()*dd->get_vel();
	n++;
      }
      cmvel /= mass;

      // Reset velocities to the center of mass frame and recompute
      // the energies.

      real ke1 = 0, pe1 = 0;
      n = 0;
      for_all_daughters(dyn, d, dd) {
	dd->inc_vel(-cmvel);
	ke1 += dd->get_mass()*square(dd->get_vel());
	pe1 += dd->get_mass()*pot[n++];
      }
      ke1 /= 2;
      pe1 /= 2;
      real e1 = pe1;

      // Final velocity rescaling should restore the desired energy.

      real vfac = sqrt((e0-pe1)/ke1);
      for_all_daughters(dyn, d, dd) {
	dd->set_vel(vfac*dd->get_vel());
	e1 += 0.5*dd->get_mass()*square(dd->get_vel());
      }

      if (verbose) {
	cout << "expand_nodes: rfac = " << rfac << endl << flush;
	cerr << "              "; PRC(e0); PRC(e1); PRL(e1-e0);
      }
    }

    // Restore the initial center of mass pos and vel.

    for_all_daughters(dyn, d, dd) {
      dd->inc_pos(cmpos_init);
      dd->inc_vel(cmvel_init);
    }
  }

  if (!expand && verbose)
    cout << "expand_nodes: no expansion performed" << endl;

  cout << flush;
  cerr << flush;
}

local void contract_nodes(dyn *d, real rscale,
			 bool verbose = false,
			 real eps2 = 0)
{
  // Contract the top-level escapers onto a sphere of radius rscale.
  // Correct velocities to preserve the top-level energy.  Ultimately,
  // we will preserve angular momentum too.  *** TO COME. *** Should
  // also check for (and correct) any accidental close pairs created,
  // but not done yet.  In this case, we make no attempt to preserve
  // relative locations.  The goal is to create a compact system for
  // reinsertion into the external integration.

  // Determine the maximum separation (no mass weighting).

  real sep2 = 0;
  for_all_daughters(dyn, d, dd)
    for (dyn *ddd = dd->get_younger_sister();
	 ddd != NULL; ddd = ddd->get_younger_sister()) {
      sep2 = max(sep2, square(dd->get_pos()-ddd->get_pos()));
    }

  real sep = sqrt(sep2);
  bool contract = false;

  // Proceed if there is something to shrink.

  if (sep > rscale) {

    contract = true;

    // Start by recomputing the center of mass pos and vel (should
    // both be zero).

    int n = 0;
    real mass = 0;
    vec cmpos = 0, cmvel = 0;

    for_all_daughters(dyn, d, dd) {
      mass += dd->get_mass();
      cmpos += dd->get_mass()*dd->get_pos();
      cmvel += dd->get_mass()*dd->get_vel();
      n++;
    }
    cmpos /= mass;
    cmvel /= mass;

    // Move to the center of mass frame, but save the original values,
    // just in case they aren't zero for some application.

    vec cmpos_init = cmpos;
    vec cmvel_init = cmvel;

    for_all_daughters(dyn, d, dd) {
      dd->inc_pos(-cmpos_init);
      dd->inc_vel(-cmvel_init);
    }

    // Compute initial energies.

    real ke0 = 0, pe0 = 0;
    real pot[n], energy[n];
    n = 0;
    for_all_daughters(dyn, d, dd) {
      pot[n] = local_potential(d, dd, eps2);
      real v2 = 0.5*square(dd->get_vel());
      energy[n] = v2 + pot[n];
      ke0 += dd->get_mass() * v2;
      pe0 += dd->get_mass() * pot[n];
      n++;
    }

    pe0 /= 2;
    real e0 = ke0 + pe0;

    // Scale the positions.  This rescaling will in general preserve
    // neither the center of mass position or velocity, so recompute
    // both and correct.

    cmpos = 0;
    cmvel = 0;

    for_all_daughters(dyn, d, dd) {
      real rr = abs(dd->get_pos());
      if (rr > rscale)
	dd->set_pos((rscale/rr)*dd->get_pos());
      cmpos += dd->get_mass()*dd->get_pos();
    }

    // Scale the velocities.  (*** DOESN'T PRESERVE ANGULAR MOMENTUM
    // -- TO COME -- USE THE KEPLER PACKAGE IN THE 2-BODY CASE ***).

    n = 0;
    for_all_daughters(dyn, d, dd) {
      real p = local_potential(d, dd, eps2);
      real dpot = p - pot[n];
      pot[n] = p;
      if (dpot < 0) {
	real vfac = sqrt(1-dpot/square(dd->get_vel()));
	dd->set_vel(vfac*dd->get_vel());
      }
      cmvel += dd->get_mass()*dd->get_vel();
      n++;
    }

    cmpos /= mass;
    cmvel /= mass;

    // Reset velocities to the center of mass frame and recompute
    // the energies.

    real ke1 = 0, pe1 = 0;
    n = 0;
    for_all_daughters(dyn, d, dd) {
      dd->inc_pos(-cmpos);
      dd->inc_vel(-cmvel);
      ke1 += dd->get_mass()*square(dd->get_vel());
      pe1 += dd->get_mass()*pot[n++];
    }

    ke1 /= 2;
    pe1 /= 2;
    real e1 = pe1;

    // Final velocity rescaling should restore the desired energy.

    real vfac = sqrt((e0-pe1)/ke1);
    for_all_daughters(dyn, d, dd) {
      dd->set_vel(vfac*dd->get_vel());
      e1 += 0.5*dd->get_mass()*square(dd->get_vel());
    }

    if (verbose) {
      cerr << "contract_nodes: "; PRL(rscale);
      cerr << "                "; PRC(e0); PRC(e1); PRL(e1-e0);
    }

    // Restore the initial center of mass pos and vel.

    for_all_daughters(dyn, d, dd) {
      dd->inc_pos(cmpos_init);
      dd->inc_vel(cmvel_init);
    }
  }

  if (!contract && verbose)
    cout << "contract_nodes: no contraction performed" << endl;

  cout << flush;
  cerr << flush;
}

local void get_top_level_mult_energies(dyn *b,
				       real &kinetic,
				       real &potential, real &potential2,
				       real eps2 = 0)
{
  kinetic = potential = potential2 = 0;

  for_all_daughters(dyn, b, bi) {
    kinetic += 0.5*bi->get_mass()*square(bi->get_vel());
    real dpot = 0, dpot2 = 0;
    for (dyn *bj = bi->get_younger_sister(); bj != NULL;
	 bj = bj->get_younger_sister()) {
      dpot += bj->get_mass() / abs(bi->get_pos()-bj->get_pos());
      dpot2 += bj->get_mass()
		/ sqrt(square(bi->get_pos()-bj->get_pos())+eps2);
    }
    potential -= bi->get_mass()*dpot;
    potential2 -= bi->get_mass()*dpot2;
  }
}


//--------------------------------------------------------------------
//
// integrate_true_multiple() -- follow an unsoftened multiple to
// completion.  Let integrate_multiple() take care of the interface
// with MUSE.

local dyn *integrate_true_multiple(hdyn *b,
				   real &de_internal,
				   real rscale,
				   int verbose = 1,
				   real eps2 = 0)
{
  // Do the entire UNSOFTENED multiple calculation.  Start from the
  // hdyn specification of the multiple system.  Internally save new
  // low-level structure and return a dyn structure describing the
  // final top level syetem, scaled back to a compact sphere, but
  // still unsoftened.  Also compute the change in internal energy of
  // the system.  Don't do any rescaling.

  // Note that the input value of eps2 is ONLY used in the calls to
  // integrate() (for get_structure) and get_tree().  Maybe should
  // always use eps2 = 0 in analyze.cc, and drop from the argument
  // list.  *** TO BE DETERMINED. ***

  //--------------------------------------------------------------------
  //
  // The integration/check time scale should not be allowed to grow
  // indefinitely for very wide systems, as this may cause us to take
  // far too many steps in a wide system.
  //
  // Start by determining a system time scale.  Use the top-level
  // dynamical time, modified for cases far from virial equilibrium.
  // (This calculation is somewhat redundant, as we just did something
  // very similar in the calling function, but repeat it here for
  // modularity.)  Limit the time scale to the "transit time" of the
  // to-level system.

  real mtot = b->get_mass();
  real kinetic = 0, potential = 0, potential2 = 0;
  get_top_level_mult_energies(b, kinetic, potential, potential2, 0);
  real escale = -potential/2;
  if (kinetic > escale) escale = min(kinetic, 4*escale);

  real tscale = 0.125*sqrt(pow(mtot,5)/pow(escale,3));
  real tscale2 = 2*sqrt(pow(rscale,3)/mtot);
  if (tscale > tscale2) tscale = tscale2;

  real tcheck = tscale;
  int ncheck = 10;
  int ntime = 2000;

  real ttrans = 2*sqrt(pow(rmax_system,3)/mtot);
  if (tcheck > ttrans) {
    tcheck = ttrans;
    ncheck = 1;
    ntime = 100;		// *** these numbers may need adjustment
  }

  if (verbose > 0) {
    PRL(tscale); PRC(tcheck); PRC(ntime); PRL(ncheck);
  }

  //--------------------------------------------------------------------
  //
  // Compute the initial total internal energy of all multiples.

  real ie_internal = 0;
  for_all_daughters(hdyn, b, bb)
    if (bb->is_parent()) ie_internal += get_energy(bb);

  //--------------------------------------------------------------------
  //
  // Absorb the tidal error due to the multiples into the top-level
  // kinetic energy (if possible) -- this is not strictly necessary,
  // but we do it to define a conserved energy for the top-level
  // system.  We will force the total energy to equal the incoming
  // top-level + internal energy.

  real ext_energy = kinetic + potential + ie_internal;
  real true_pot, true_kin, true_energy;
  calculate_energies(b, 0, true_pot, true_kin, true_energy);
  real de_tidal = ext_energy - true_energy;
  if (verbose > 0)
    cerr << "tidal error 0 = " << de_tidal << endl << flush;

  // Note: this code is replicated below (second tidal correction).

  real vfac = 1 + de_tidal/kinetic;
  if (vfac > 0) {
    vfac = sqrt(vfac);
    if (verbose > 0) {
      cerr << "velocity correction: "; PRL(vfac);
    }
    for_all_daughters(hdyn, b, bb)
      bb->set_vel(vfac*bb->get_vel());
    kinetic += de_tidal;
    true_energy = ext_energy;
  } else if (verbose > 0)
      cerr << "...no correction applied" << endl << flush;

  //--------------------------------------------------------------------
  //
  // Integrate the system to completion (note that integrate() starts
  // by flattening b).

  kira_options ko;
  ko.perturber_criterion = 2;
  b->set_kira_options(&ko);

  b->set_system_time(0);
  for_all_nodes(hdyn, b, bi)
    bi->set_time(0);

  if (verbose > 0) my_sys_stats(b, verbose);

  // *** Note that we need better (consistent) treatment of
  // *** unperturbed multiples in kira_smallN in order to efficiently
  // *** handle hierarchical systems.

  // PRL(verbose);
  // bool over = 
  integrate(b, ntime, ncheck, tcheck, 0, verbose>1);
  // PRL(over);

  if (verbose > 0) {
    cout << "  End of interaction at " << flush;
    PRL(b->get_time());
    my_sys_stats(b, verbose);
  }

  if (verbose > 2) {
    cout << "hdyn tree:" << endl;
    pp2(b);
    for_all_daughters(hdyn, b, bb)
      cout << bb->get_index() << " " << bb->get_pos() << endl;
  }

  //--------------------------------------------------------------------
  //
  // Construct a list d of top-level stable, escaping nodes, then
  // restructure d into a flat array of stable objects.  Probably
  // should merge these two functions.

  if (verbose > 0) PRL(get_energy(b));
  dyn *d = get_tree(b, eps2, verbose>0);	// this is the only place eps2
						// is used in this function
						// -- set eps2 = 0 always?

  if (verbose > 2) {
    cout << "binary tree:" << endl;
    pp2(d);
    for_all_nodes(dyn, d, dd)
      cout << dd->format_label() << " " << dd->get_pos() << endl;
  }

  set_top_level_escapers(d, verbose>0);

  if (verbose > 2) {
    cout << "escaper tree:" << endl;
    pp2(d);
    for_all_daughters(dyn, d, dd)
      cout << dd->get_index() << " " << dd->get_pos() << endl;
  }

  //--------------------------------------------------------------------
  //
  // Modify the tree to exclude nodes that are too large.  Replace
  // them by their components at the top level, even if the components
  // are bound.

   split_wide_binaries(d, rmax_system, verbose>0);

  //--------------------------------------------------------------------
  //
  // Add appropriate labels (get_index) to all new top-level nodes.

  label_new_nodes(input_id, d, verbose>0);

  //--------------------------------------------------------------------
  //
  // Determine the status of all old nodes.

  flag_status(input_id, d, status, verbose>0);

  //--------------------------------------------------------------------
  //
  // Compute the final total internal energy of multiples.

  real fe_internal = 0;
  for_all_daughters(dyn, d, dd)
    if (d->is_parent()) fe_internal += get_energy(dd);

  de_internal = fe_internal - ie_internal;
  if (verbose > 0) {
    PRC(ie_internal); PRC(fe_internal); PRL(de_internal);
    PRL(get_energy(d));
  }

  // Calculate the total energy of the system before stripping the
  // multiples.

  d->set_root(d);
  calculate_energies(d, 0, true_pot, true_kin, true_energy);

  //--------------------------------------------------------------------
  //
  // Remove all old multiples from the storage arrays.

  remove_old_multiples(input_id, status);

  //--------------------------------------------------------------------
  //
  // Store all new multiples (effectively updates the survivors),
  // leaving d as a flat tree of top-level centers of mass.  On
  // return, only top-level nodes are left under d.  Note that the
  // logic here should also be able to cope with the possibility that
  // the system is now a single bound object.

  real max_node_size = store_new_multiples(d);

  if (verbose > 1) {
    cout << "multiple initial " << flush;
    PRC(rscale); PRL(max_node_size);
  }
  if (rscale < max_node_size) rscale = max_node_size;

  if (verbose > 2) {
    cout << "final top-level nodes:" << endl;
    pp2(d);
    for_all_daughters(dyn, d, dd)
      cout << dd->get_index() << " " << dd->get_pos() << endl;
  }

  //--------------------------------------------------------------------
  //
  // Again correct for the tidal error due to multiples.

  get_top_level_mult_energies(d, kinetic, potential, potential2, 0);
  ext_energy = kinetic + potential + fe_internal;
  de_tidal = ext_energy - true_energy;
  if (verbose > 0)
    cerr << "tidal error 1 = " << de_tidal << endl << flush;

  // Note: this code is a near-copy of previous (first tidal
  // correction).

  vfac = 1 - de_tidal/kinetic;
  if (vfac > 0) {
    vfac = sqrt(vfac);
    if (verbose > 0) {
      cerr << "velocity correction: "; PRL(vfac);
    }
    for_all_daughters(dyn, d, dd)
      dd->set_vel(vfac*dd->get_vel());
    kinetic += de_tidal;
    true_energy = ext_energy;
  } else if (verbose > 0)
      cerr << "...no correction applied" << endl << flush;

  //--------------------------------------------------------------------

  cout << flush;
  cerr << flush;
  return d;
}

local void print_total_energy(dyn *d, real de_internal = 0)
{
  // Compute and print the total top-level energy and the change in
  // the internal energy.

  real kinetic = 0, potential = 0, potential2 = 0;
  get_top_level_mult_energies(d, kinetic, potential, potential2);
  cout << "top-level energy = " << kinetic+potential
       << ", de_internal = " << de_internal
       << ", total = " << kinetic+potential+de_internal
       << endl << flush;
}


//--------------------------------------------------------------------
//
// integrate_multiple() -- the MUSE interface function.  Note that
// most of the actual work is done by integrate_true_multiple().

int integrate_multiple(int  verbose,	// default = 1
		       real eps2)	// default = 0
		       
{
  // This function (1) constructs an N-body system representing the
  //                   interaction,
  //		   (2) converts softened to unsoftened potentials and
  //		       rescales all velocities accordingly,
  //		   (3) rescales the system to a standard radius (in case
  //		       of errors in the calling function), preserving
  //		       unsoftened energy,
  //               (4) sends the system to integrate(),
  //               (5) recovers the system once the interaction is over,
  //	     	       with all internal bookkeeping updated and with
  //		       only top-level escaping nodes remaining,
  //		   (6) rescales the system to the standard radius, again
  //		       preserving unsoftened energy,
  //		   (7) converts unsoftened to softened potentials,
  //		       rescaling all velocities accordingly,
  //               (8) returns to the external calculation.
  //
  // Some of this functionality may change or migrate as the interface
  // with the external dynamical integrator evolves.

  //--------------------------------------------------------------------
  //
  // Construct the system to integrate.  Determine the center of mass
  // and work in the center of mass frame.  We will restore the center
  // of mass properties at the end.

  // *** Might be cleaner to deal with the top-level scaling, etc.,
  // *** before adding internal structure...

  real mtot = 0;
  vec cmpos, cmvel;
  hdyn *b = create_multiple(input_id, input_mass, input_pos, input_vel,
			    mtot, cmpos, cmvel);

  if (verbose > 1) {
    for_all_daughters(hdyn, b, bb) {
      PRC(bb->format_label()); PRL(get_effective_radius(bb));
      PRL(bb->get_pos()+cmpos);
      PRL(bb->get_vel()+cmvel);
      for (hdyn *bbb=bb->get_younger_sister();
	   bbb != NULL; bbb = bbb->get_younger_sister()) {
	PRC(bbb->format_label()); PRL(abs(bb->get_pos()-bbb->get_pos()));
      }
    }
  }

  b->set_mass(mtot);

  // (This is actually the last time the input mass, pos, and vel are
  // used.  However, don't delete the input arrays yet, as we will
  // want to use input_id later to create and access the status
  // array.)

  // At this stage, b is structured, with the top-level nodes
  // representing the input objects from the external system.  We are
  // in the center of mass frame, with mtot, cmpos, cmvel the total
  // mass and center of mass position and velocity.

  //--------------------------------------------------------------------
  //
  // Determine a suitable radius for the system.

  real rscale = 0;
  for_all_daughters(hdyn, b, bb) {
    real r2 = square(bb->get_pos());
    if (r2 > rscale) rscale = r2;
  }
  rscale = sqrt(rscale);

  // Include the radii of all top-level nodes (not quite the same
  // thing).

  for_all_daughters(hdyn, b, bb)
    if (bb->is_parent()) {
      real r_eff = get_effective_radius(bb);
      if (r_eff > rscale) rscale = r_eff;
    }

  // We don't want systems growing to arbitrary size.  Impose a hard
  // upper limit on the size of any resultant system, and don't
  // combine nodes whose effective size would become greater than
  // this.  Actually simpler to combine nodes as usual, then split
  // nodes that are too large to be practical (currently done in
  // integrate_true_multiple).

  if (rscale > rmax_system) rscale = rmax_system;

  if (verbose > 0) {
    PRC(rscale); PRL(rmax_system);
  }

  //--------------------------------------------------------------------
  //
  // Compute the top-level energies (to preserve the total energy of
  // the external system), with and without softening.  Note that
  // tidal effects are not yet included -- this is the responsibility
  // of the calling function.

  real ikinetic = 0, ipotential = 0, ipotential2 = 0;
  get_top_level_mult_energies(b, ikinetic, ipotential, ipotential2, eps2);
  real ikinetic2 = ikinetic;

  //--------------------------------------------------------------------
  //
  // --> A gory detail...  We must make a concession here to the basic
  // inconsistency between the external (softened) and the internal
  // (unsoftened) descriptions of the motion.  In a typical situation,
  // the top-level nodes will be unbound in the external system.
  // However, when we switch off softening, it may happen that they
  // become bound in the internal description.  In principle, we can
  // cope with a bound system, but this is usually not what we
  // intended.
  //
  // The underlying picture here is that the external integrator
  // delivers the smallN subsystem a series of encounters to follow.
  // If the external system really had been unsoftened, then the
  // relative top-level velocities would have been greater.  We mimic
  // this by scaling up the incoming velocities to follow the
  // unsoftened potential.  We will do something similar, in reverse,
  // at the end, and will take care to ensure that the
  // scaling/unscaling doesn't bias the overall energetics of the
  // interaction.
  //
  // Scale top-level velocities to preserve T/U as U goes from
  // softened to unsoftened.  This keeps unbound systems unbound, and
  // vice versa.  Use an unsoftened potential from here until the end
  // of the calculation.

  real vfac2 = sqrt(ipotential/ipotential2);
  real vfac = sqrt(vfac2);
  if (verbose > 0) {
    PRC(ipotential); PRL(ipotential2);
    cout << "scaling velocities to unsoftened potential: " << flush;
    PRL(vfac);
  }
  for_all_daughters(hdyn, b, bi)
    bi->set_vel(vfac*bi->get_vel());
  ikinetic *= vfac2;

  // Note: ikinetic2 is the original kinetic energy, before scaling to
  // the unsoftened potential; ikinetic refers to the scaled system.

  //--------------------------------------------------------------------
  //
  // Expand top-level nodes so that their mean separation is rscale,
  // if possible.  Correct the velocities to preserve the top-level
  // unsoftened energy.  DON'T rescale them all onto the same sphere
  // -- try to preserve the encounter.

  if (verbose > 0) print_total_energy(b);
  expand_nodes(b, rscale, verbose>1, 0);
  if (verbose > 0) print_total_energy(b);

  //--------------------------------------------------------------------
  //
  // Compute the multiple calculation, using an unsoftened potential,
  // classifying and saving all internal structure, but performing no
  // rescaling.  Return an unsoftened tree d of top-level escaping
  // nodes, probably widely separated.

  real de_internal = 0;
  dyn *d = integrate_true_multiple(b, de_internal, rscale, verbose, eps2);

  if (verbose > 0) {
    cout << "after integrate_true_multiple..." << endl << flush;
    pp2(d);
    cerr << flush;
  }

  // All rescaling to compact configurations and between softened and
  // unsoftened potentials should be handled at THIS level, not in
  // integrate_true_multiple().

  //--------------------------------------------------------------------
  //
  // The following scalings only make sense if there is more than one
  // top-level node, and if the top-level nodes are mutually unound.

  if (d->n_daughters() > 1) {

    //------------------------------------------------------------------
    //
    // Map escapers back to a sphere of radius rscale.  Correct the
    // velocities to preserve the top-level unsoftened energy.

    if (verbose > 0) print_total_energy(d, de_internal);
    contract_nodes(d, rscale, verbose>1, 0);
    if (verbose > 0) print_total_energy(d, de_internal);

    //------------------------------------------------------------------
    //
    // Restore the top-level softened potential.  Scale top-level
    // velocities to preserve T/U as U goes from unsoftened to softened.
    // This keeps unbound systems unbound, and vice versa.  (See the
    // earlier "gory" discussion.)

    real fkinetic = 0, fpotential = 0, fpotential2 = 0;
    get_top_level_mult_energies(d, fkinetic, fpotential, fpotential2, eps2);

    vfac = sqrt(fpotential2/fpotential);
    if (verbose > 0) {
      cout << "scaling velocities back to softened potential: " << flush;
      PRL(vfac);
    }
    for_all_daughters(dyn, d, di)
      di->set_vel(vfac*di->get_vel());

    fkinetic *= vfac*vfac;

    if (verbose > 1) {
      PRC(ipotential); PRL(ipotential2);
      PRC(ikinetic); PRL(ikinetic2);
      PRC(ikinetic+ipotential); PRL(ikinetic2+ipotential2);
      PRC(fpotential); PRC(fpotential2); PRL(fkinetic);
      PRC(fkinetic+fpotential); PRL(fkinetic+fpotential2);
      PRL(fkinetic+fpotential2-ikinetic-ipotential2);
      PRC(de_internal);
    }

#if 1
    //------------------------------------------------------------------
    //
    // The scaling to and from rscale and the internal integration
    // should have preserved total energy (top-level + internal)
    // exactly, but it is not guaranteed that when we restore
    // softening we will still have conservation of external softened
    // energy + internal energy.

    real de_top_level = fkinetic + fpotential2 - ikinetic2 - ipotential2;
    real dde = -de_top_level - de_internal;	// additional top-level
						// energy change needed
    if (verbose > 0) {
      PRC(de_top_level); PRL(dde);
    }

    // Rescale so that the change in the (softened) top-level energy
    // matches the change in internal energy.  Flag a problem if the
    // energy change required is greater than the top-level kinetic
    // energy fkinetic. ***

    // *** What is the proper treatment here?  TO COME... ***

    vfac = 1 + dde/fkinetic;
    if (vfac < 0) {
      if (verbose > 0) cout << "warning: dde < -fkinetic" << endl << flush;
      if (verbose == 1) PRL(fkinetic);
      vfac = 0.1;				// arbitrary

      top_level_energy_error -= dde + (1-vfac)*fkinetic;
    }

    vfac = sqrt(vfac);

    if (verbose > 0) {
      cout << "scaling velocities to reflect internal energy change: "
	   << flush;
      PRL(vfac);
    }
    for_all_daughters(dyn, d, dd)
      dd->set_vel(vfac*dd->get_vel());

    if (verbose > 1) {
      get_top_level_mult_energies(d, fkinetic, fpotential, fpotential2, eps2);
      PRC(fpotential); PRC(fpotential2); PRL(fkinetic);
      PRC(fkinetic+fpotential); PRL(fkinetic+fpotential2);
      PRL(fkinetic+fpotential2-ikinetic2-ipotential2);
    }
#endif

  }

  //--------------------------------------------------------------------
  //
  // Restore center of mass positions and velocities to the top-level
  // nodes.

  for_all_daughters(dyn, d, dd) {
    dd->inc_pos(cmpos);
    dd->inc_vel(cmvel);
  }

#if 0
  real fkinetic, fpotential, fpotential2;
  get_top_level_mult_energies(d, fkinetic, fpotential, fpotential2, eps2);
  PRC(fpotential); PRC(fpotential2); PRL(fkinetic);
  PRC(fkinetic+fpotential); PRL(fkinetic+fpotential2);
  PRL(fkinetic+fpotential2-ikinetic2-ipotential2);
#endif

  //--------------------------------------------------------------------
  //
  // Create a final list of all ids, masses, positions, and velocities
  // present at the end of the interaction.  In addition, set the
  // (somewhat redundant) is_new flag.

  list_final(d);

  // Clean up.

  rmtree(d);
  rmtree(b);

  cout << flush;
  cerr << flush;

  return new_id.size();
}
