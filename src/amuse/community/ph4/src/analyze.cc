
// Analysis routines for use with smallN.  Based on the Starlab
// version, with a custom class to make it self-contained.  The only
// externally visible functions are
//
//	bool  check_structure(hdyn *bin, bool verbose = true)
//	hdyn* get_tree(hdyn *bin)
//
// NOTES: The original code code used Starlab stories to track
// critical information.  Stories are handy, but maybe too much to
// implement here.  The relevant story quantities are "E",
// "escaping_cpt", "escaper", "periastron ratio", and "stable":
//
//	- in a CM node, E is the relative energy (per unit reduced
//	  mass) of the two (bound) components
//
//	- in a CM node, escaping_cpt is true if one or more daughter
//	  nodes is escaping; in practice, the CM node is the root node
//
//	- a (top-level) node is flagged with escaper = true if it is
//	  escaping
//
//	- in a multiple CM node, periastron_ratio is the ratio of outer
//	  periastron to inner semi-major axis (forms the basis for a
//	  simple stability criterion)
//
//	- in a multiple CM node, stable = true means that the multiple
//	  can be regarded as stable; a multiple with stable = false
//	  may be promoted to stable if its configuration remains
//	  unchanged for an extended quarantine period
//
// Implement these here as extra public variables in the derived hdyn2
// class.
//
// The code also uses labels in addition to indices to represent node
// names.  Again, these are implemented in class hdyn2, which is local
// to this file.

#include "hdyn.h"
#include <cstring>
#include "nstab.h"

#ifndef TOOLBOX

#define local static

class hdyn2 : public hdyn
{
    private:
	string name;
    public:
	real relative_energy;
	real periastron_ratio;
	bool escaping_cpt;
	bool escaper;
	bool stable;

    hdyn2() : hdyn() {
	relative_energy = periastron_ratio = 0;
	escaping_cpt = escaper = stable = false;
    }
    hdyn2 *get_parent()			const {return (hdyn2*)parent;}
    hdyn2 *get_oldest_daughter()	const {return (hdyn2*)oldest_daughter;}
    hdyn2 *get_younger_sister()		const {return (hdyn2*)younger_sister;}
    hdyn2 *get_older_sister()		const {return (hdyn2*)older_sister;}

    void set_name(char *id)		{name = id;}
    void set_name(string id)		{name = id;}
    const string get_name()		const {return name;}

    const char *format_label() const;
    int  n_leaves() {int n = 0; for_all_leaves(hdyn2, this, b) n++; return n;}
};

#define BUF_SIZE 1024
static char format_string[BUF_SIZE];
const char* hdyn2::format_label() const
{
    // Precedence:	print name string if defined
    //			otherwise, print index if non-negative
    //			otherwise, print '?'

    if (name.size() > 0) {
	strncpy(format_string, name.c_str(), BUF_SIZE-1);
	format_string[BUF_SIZE-1] = '\0';
    } else if(index >= 0)
	sprintf(format_string, "%d", index);
    else
	sprintf(format_string, "?");

    return format_string;
}

const char *string_index_of_node(hdyn2 *n)
{
    if (n->get_name().size() > 0) {
	return n->get_name().c_str();
    } else {
	static char integer_id[20];
	sprintf(integer_id, "%d", n->get_index());
	n->set_name(integer_id);
	return n->get_name().c_str();
    }
}

local const char *construct_binary_label(hdyn2 *ni, hdyn2 *nj)
{
    static char new_name[256];
    sprintf(new_name,"(%s,%s)", string_index_of_node(ni),
	    string_index_of_node(nj));
    return new_name;
}

local void label_binary_node(hdyn2 *n)
{
    if (n->is_root()) return;

    for_all_daughters(hdyn2, n, nn) label_binary_node(nn);

    if (!n->is_leaf()) {
        hdyn2* od = n->get_oldest_daughter();
        hdyn2* yd = od->get_younger_sister();
        n->set_name(construct_binary_label(od, yd));
    }
}



hdyn2 *flat_copy_tree(hdyn *bin)
{
    // Flatten the input tree into a new tree, and return a pointer to
    // the result.

    hdyn2 *b = new hdyn2(), *bo = NULL, *by = NULL;
    b->set_parent(NULL);
    b->set_mass(0);
    b->set_name((char*)"root");

    for_all_leaves(hdyn, bin, bb) {	// bb is hdyn, others are hdyn2
	bo = by;
	by = new hdyn2();
	if (!b->get_oldest_daughter())
	    b->set_oldest_daughter(by);
	if (bo)
	    bo->set_younger_sister(by);
	by->set_parent(b);
	by->set_older_sister(bo);
	by->set_younger_sister(NULL);

	by->set_index(bb->get_index());
	by->set_mass(bb->get_mass());
	b->set_mass(b->get_mass()+bb->get_mass());
	by->set_pos(bb->get_pos());
	by->set_vel(bb->get_vel());
	hdyn *pp = bb->get_parent();
	while (pp != bin) {
	    by->inc_pos(pp->get_pos());
	    by->inc_vel(pp->get_vel());
	    pp = pp->get_parent();
	}
    }

    b->set_system_time(bin->get_system_time());

    return b;
}

local vec relative_pos(hdyn2 *b, hdyn2 *bb)
{
    // Return the position of bb relative to its ancestor b.

    vec pos = bb->get_pos();
    while (bb != b) {
	hdyn2 *p = bb->get_parent();
	pos += p->get_pos();
	bb = p;
    }
    return pos;
}

local vec relative_vel(hdyn2 *b, hdyn2 *bb)
{
    // Return the velocity of bb relative to its ancestor b.

    vec vel = bb->get_vel();
    while (bb != b) {
	hdyn2 *p = bb->get_parent();
	vel += p->get_vel();
	bb = p;
    }
    return vel;
}

local inline real get_potential(hdyn2 *bi, hdyn2 *bj)
{
    // Compute the total potential energy of node bi relative to its
    // sister node node bj.  Compute the potential energy exactly, by
    // including all leaves, not using the center of mass
    // approximation.

    real phi = 0;
    hdyn2 *p = bi->get_parent();	// should also be bj->get_parent()

    for_all_leaves(hdyn2, bi, ii) {
	vec posi = relative_pos(p, ii);
	for_all_leaves(hdyn2, bj, jj) {
	    vec posj = relative_pos(p, jj);
	    phi -= ii->get_mass()*jj->get_mass()/abs(posi-posj);
	}
    }
    return phi;
}

local real get_relative_energy(hdyn2 *bi, hdyn2 *bj)
{
    // Return the relative energy (per unit reduced mass) of nodes bi
    // and bj.

    hdyn2 *p = bi->get_parent();	// p should be root

    vec xi = relative_pos(p, bi);
    vec xj = relative_pos(p, bj);
    vec vi = relative_vel(p, bi);
    vec vj = relative_vel(p, bj);
    real mi = bi->get_mass();
    real mj = bj->get_mass();
    real mu_inv = (mi+mj)/(mi*mj);

    return 0.5*square(vi-vj) + mu_inv*get_potential(bi, bj);
}

local void combine_nodes(hdyn2 *bi, hdyn2 *bj)
{
    // Make these particles into a binary node.

    // Create the CM node.

    real mtot = bi->get_mass() + bj->get_mass();
    vec cmpos = (bi->get_mass()*bi->get_pos()
		 + bj->get_mass()*bj->get_pos()) / mtot;
    vec cmvel = (bi->get_mass()*bi->get_vel()
		 + bj->get_mass()*bj->get_vel()) / mtot;
    hdyn2 *cm = new hdyn2();
    cm->set_mass(mtot);
    cm->set_pos(cmpos);
    cm->set_vel(cmvel);

    // Offset components to the center of mass frame.

    bi->inc_pos(-cm->get_pos());
    bj->inc_pos(-cm->get_pos());
    bi->inc_vel(-cm->get_vel());
    bj->inc_vel(-cm->get_vel());

    // Restructure the tree to create the binary (bi,bj).  The new CM
    // node is cm.

    create_binary_node(cm, bi, bj);
    label_binary_node(cm);
}

local inline real separation_squared(hdyn2 *bi, hdyn2 *bj)
{
    return square(bi->get_pos()-bj->get_pos());
}

local inline real size_squared(hdyn2 *b)
{
    if (!b->is_parent()) return 0;
    hdyn2 *od = b->get_oldest_daughter();
    hdyn2 *yd = od->get_younger_sister();
    return separation_squared(od, yd);
}

// The perturbation defined here as "light" effectively determines the
// threshold perturbation for unperturbed multiple motion.  Must make
// it consistent with the criteria applied in smallN.

const real MAX_PERT_SQ = 1.e-4;

local bool is_lightly_perturbed(hdyn2 *b, hdyn2 *bi, hdyn2 *bj)
{
    real mtot = bi->get_mass() + bj->get_mass();
    vec cm = (bi->get_mass()*bi->get_pos() + bj->get_mass()*bj->get_pos())
		/ mtot;
    real mr3_ij_sq = 0.25 * mtot * mtot
			  / pow(square(bi->get_pos()-bj->get_pos()), 3);

    real max_pert_sq = 0;
    hdyn2 *max_pert_bb = NULL;  
    for_all_daughters(hdyn2, b, bb)
	if (bb != bi && bb != bj) {
	    real r3_sq = square(bb->get_pos() - cm);
	    real pert_sq = pow(bb->get_mass(), 2) / pow(r3_sq, 3);
	    if (pert_sq > max_pert_sq) {
		max_pert_sq = pert_sq;
		max_pert_bb = bb;
	    }
	}

    return (max_pert_sq < MAX_PERT_SQ*mr3_ij_sq);
}

#define PRINT_DETAILS 0

local void decompose(hdyn2 *b)
{
    // Construct a binary tree by recursively pairing the top-level
    // particles with the greatest mutual potential.  Don't pair if
    // strongly perturbed.
  
    while (true) {

	// Look for the most bound pair among the current top-level nodes.

	if (PRINT_DETAILS) {
	    cout << "top level:";
	    for_all_daughters(hdyn2, b, bi)
		cout << " " << bi->format_label();
	    cout << endl;
	}

	real min_energy = 0;
	hdyn2 *bimin = NULL;
	hdyn2 *bjmin = NULL;

	// Slightly wasteful to recompute all pairwise energies every
	// time around (really only a problem for N >> 3)...

	for_all_daughters(hdyn2, b, bi)
	    for (hdyn2 *bj = bi->get_younger_sister(); bj != NULL;
		 bj = bj->get_younger_sister()) {
		real Eij = get_relative_energy(bi, bj);
		if (PRINT_DETAILS) {
		    cout << "    " << bi->format_label();
		    cout << " " << bj->format_label() << " ";
		    PRL(Eij);
		}
		if (Eij < min_energy && is_lightly_perturbed(b, bi, bj)) {
		    min_energy = Eij;
		    bimin = bi;
		    bjmin = bj;
		}
	    }

	if (bimin && bjmin) combine_nodes(bimin, bjmin);
	else break;
    }
}

local void compute_relative_energy(hdyn2 *b)
{
    for_all_nodes(hdyn2, b, bb)
	if (bb->get_parent() && bb->is_parent()) {    // bb may be root, note
	    hdyn2 *od = bb->get_oldest_daughter();
	    hdyn2 *yd = od->get_younger_sister();

	    // Note: 2-body energy will be too crude if (e.g.) the
	    // components of one sister node are widely separated and
	    // the CM happens to lie close to the other sister.
	    // Better to compute the relative potential energy
	    // exactly.  Note also that the energies stored are per
	    // unit reduced mass.

      	    real mu_inv = bb->get_mass()/(od->get_mass()*yd->get_mass());
	    real E = 0.5*square(od->get_vel()-yd->get_vel())
			+ mu_inv*get_potential(od, yd);
	    bb->relative_energy = E;				  // STORY
	}
}

local bool is_escaper(hdyn2 *b, hdyn2 *bi)
{
    // An escaper is a (top-level) particle having positive energy
    // relative to all other top-level nodes as well as their CM, and
    // moving away from the CM.
  
    bool esc1 = true;

    real mi = bi->get_mass();
    vec xi = bi->get_pos();
    vec vi = bi->get_vel();
    float mcm = 0;
    vec rcm = 0;
    vec vcm = 0;
    real phicm = 0;

    real Emin = 0;
    hdyn2 *bmin = NULL;
    for_all_daughters(hdyn2, b, bj)
	if (bj != bi) {
	    real mj = bj->get_mass();
	    vec vj = bj->get_vel();
	    mcm += mj;
	    rcm += mj*bj->get_pos();
	    vcm += mj*vj;
	    real mu_inv = (mi+mj)/(mi*mj);
	    real dphi = get_potential(bi, bj);
	    phicm += dphi;
	    real Eij = 0.5*square(vi-vj) + mu_inv*dphi;
	    esc1 &= (Eij > 0);
	    if (Eij < 0 && Eij < Emin) {
		Emin = Eij;
		bmin = bj;
	    }
	}

    bool esc = esc1;

    if (esc && mcm > 0) {
	rcm = rcm/mcm;
	vcm = vcm/mcm;
	esc &= ((vi-vcm)*(xi-rcm) > 0
		  && 0.5*mi*mcm*square(vi-vcm)/(mi+mcm) + phicm > 0);
    }

#if 0
    if (!esc) {
	if (!esc1) {
	    cout << "particle " << bi->format_label() << " is bound to ";
	    cout << bmin->format_label() << " (" << Emin << ")"
		 << endl << flush;
	} else {
	    cout << "particle " << bi->format_label()
		 << " is bound to the system CM"
		 << endl << flush;
	}
    }
#endif

    return esc;
}

local void flag_escapers(hdyn2 *b)
{
    for_all_daughters(hdyn2, b, bi)
	if (is_escaper(b, bi)) {
	    b->escaping_cpt = true;				  // STORY
	    bi->escaper = true;					  // STORY
	}
}



local real get_semi(hdyn2 *b)
{
    // Return the semi-major-axis of the binary formed by the
    // (presumably bound) components of b.  By construction, we should
    // have an energy already stored in the hdyn2 "story."

    real semi = -1;

    if (b->is_parent()) {
	if (b->relative_energy != 0) {
	    real E = b->relative_energy;			  // STORY
	    if (E < 0) semi = -0.5*b->get_mass()/E;

	} else {

	    // Looks like we have to compute the energy from scratch.
	    // Really do use the 2-body expression here.

	    hdyn2 *od = b->get_oldest_daughter();
	    hdyn2 *yd = od->get_younger_sister();

	    real E = 0.5*square(od->get_vel()-yd->get_vel())
			- b->get_mass()/abs(od->get_pos()-yd->get_pos());
	    if (E < 0) semi = -0.5*b->get_mass()/E;
	}
    }

    return semi;
}

// Code from Starlab:

local inline bool harrington_stable(real m12, real m3,
				    real peri_fac, real cos_i)
{
    // Use a simple distance criterion with an extra mass term that
    // goes to 1 for equal-mass systems.  The Harrington (1972)
    // criterion ranges from outer_peri/inner_semi ~ 3.5 for planar
    // prograde orbits (cos_i = 1) to ~2.75 for planar retrograde
    // (cos_i = -1).  The overall factor can be adjusted, but the
    // inclination variation preserves this ratio.

    const real HARRINGTON_FAC = 3.5;		// (NB Starlab uses 3.25)
    const real HARRINGTON_RATIO = 25.0/3;

    // Ratio quantifies "closeness to stability."

    real ratio = (peri_fac * pow(m12/m3, 1./3))
		  / (HARRINGTON_FAC * (HARRINGTON_RATIO+cos_i)
				    / (HARRINGTON_RATIO+1));
    return (ratio > 1);
}

local inline bool aarseth_stable(real m12, real m3,
				 real e_outer, real peri_fac, real cos_i)
{
    // The Mardling and Aarseth (1999) stability criterion is
    //
    //    outer_peri/inner_semi > AARSETH_STABLE_FAC
    //				* [ (1 + q_outer)
    //				    * (1 + e_outer)
    //				    / sqrt(1 - e_outer) ] ^(2/5).
    //
    // This criterion omits an inclination factor that should
    // effectively reduce the critical outer periastron by about 30%
    // for retrograde orbits relative to prograde orbits.  Aarseth
    // (2003) recommends an additional factor linear in the
    // inclination to correct this.

    const real AARSETH_FAC = 2.8;		// don't change!
    const real AARSETH_RATIO = 23.0/3;		// (basically the same as the
						//  Harrington ratio above)

    real q_outer = m3/m12;
    real ratio = (peri_fac * pow((1+q_outer)*(1+e_outer)/sqrt(1-e_outer), -0.4))
		  / (AARSETH_FAC * (AARSETH_RATIO+cos_i) / (AARSETH_RATIO+1));
    return (ratio > 1);
}

local inline bool mardling_stable(kepler& outerkep, kepler& innerkep,
				  real m1, real m2, real m3, real cos_i)
{
    // Apply the Mardling (2007) criterion.  No free parameters!

    // First apply a couple of extra checks not in the Mardling
    // function, but included in the Aarseth nbody6 version:

    // 1. Reject outer pericenter inside inner apocenter.

    if (outerkep.get_periastron() < innerkep.get_apastron()) return false;

    // 2. Reject period ratio < 1.

    real sigma = outerkep.get_period()/innerkep.get_period();
    if (sigma < 1) return false;

    // Create parameters for the Mardling function:
    //
    // sigma      =  period ratio (outer/inner) ** should be > 1 **
    // ei0        =  initial inner eccentricity
    // eo         =  outer eccentricity
    // relinc     =  relative inclination (radians)
    // m1, m2, m3 =  masses (any units; m3=outer body)

    real ei0 = innerkep.get_eccentricity();
    real eo = outerkep.get_eccentricity();
    real relinc = acos(cos_i);

    // A return value of 1 from nstab_ means that we can't treat this
    // object as stable.  No in between!

    return (nstab_(&sigma, &ei0, &eo, &relinc, &m1, &m2, &m3) != 1);
}

local bool is_stable(hdyn2 *b)
{
    // Determine the stability of the multiple with top-level center
    // of mass b.

    hdyn2 *od = b->get_oldest_daughter();
    if (!od) return false;
    hdyn2 *yd = od->get_younger_sister();
    if (!yd) return false;

    // Identify an "inner" binary.  Presumably any multiples under b
    // are stable.  However, both components could be multiples, so
    // choose the wider one.  This stability criterion should perhaps
    // be improved.

    hdyn2 *b_in = od;
    real a_in = get_semi(od);
    if (get_semi(yd) > a_in) {
	b_in = yd;
	a_in = get_semi(yd);
    }
    if (a_in < 0) return false;

    // Create a kepler structure describing the outer orbit...

    kepler outerkep;
    outerkep.set_time(0);
    outerkep.set_total_mass(b->get_mass());
    outerkep.set_rel_pos(yd->get_pos()-od->get_pos());
    outerkep.set_rel_vel(yd->get_vel()-od->get_vel());
    outerkep.initialize_from_pos_and_vel(true);
    b->periastron_ratio = outerkep.get_periastron() / a_in;	  // STORY

    // ...and another describing the inner orbit.

    hdyn2 *od_in = b_in->get_oldest_daughter();
    if (!od_in) return false;
    hdyn2 *yd_in = od_in->get_younger_sister();
    if (!yd_in) return false;

    kepler innerkep;
    innerkep.set_time(0);
    innerkep.set_total_mass(b_in->get_mass());
    innerkep.set_rel_pos(od_in->get_pos() - yd_in->get_pos());
    innerkep.set_rel_vel(od_in->get_vel() - yd_in->get_vel());
    innerkep.initialize_from_pos_and_vel();

    real cos_i = innerkep.get_normal_unit_vector()
			* outerkep.get_normal_unit_vector();

    // Harrington criterion.

    bool hstable = harrington_stable(od->get_mass(), yd->get_mass(),
				     b->periastron_ratio, cos_i);

    // Aarseth criterion.

    bool astable = aarseth_stable(od->get_mass(), yd->get_mass(),
				  outerkep.get_eccentricity(),
				  b->periastron_ratio, cos_i);

    // Mardling criterion.

    bool mstable = mardling_stable(outerkep, innerkep,
				   od_in->get_mass(), yd_in->get_mass(),
				   b_in->get_binary_sister()->get_mass(),
				   cos_i);

    if ((mstable && (!hstable || !astable))
	|| (!mstable && (hstable || astable))) {
	PRC(hstable); PRC(astable); PRL(mstable);
    }

    // Choose Mardling.

    return mstable;
}

local bool check_stability(hdyn2 *b)
{
    // Stablility criterion for a CM node is that both daughters are
    // stable and sufficiently widely separated.  We currently don't
    // include external perturbations explicitly in the decision, but
    // note that we have already applied a perturbation condition in
    // constructing the tree.

    bool stable = false;

    if (b->is_leaf())				// single stars are stable
	stable = true;

    else {
	hdyn2 *od = b->get_oldest_daughter();
	hdyn2 *yd = od->get_younger_sister();

	if (od->is_leaf() && yd->is_leaf()) {
	    if (b->relative_energy < 0)				  // STORY
		stable = true;

	} else {

	    // Check both daughters (don't && these tests in the if
	    // statement because this could skip the second check).

	    bool od_stab = check_stability(od);
	    bool yd_stab = check_stability(yd);

	    if (od_stab && yd_stab)
		if (b->relative_energy < 0)			 // STORY
		    stable = is_stable(b);
	}
    }

    if (stable) b->stable = true;				 // STORY
    return stable;
}

local void check_stability_master(hdyn2 *b)
{
    for_all_daughters(hdyn2, b, bd)
	if (bd->is_leaf()) {
	    if (is_escaper(b, bd)) bd->stable = true;		 // STORY
	} else
	    if (check_stability(bd)) bd->stable = true;		 // STORY
}



local void print_recursive(hdyn2 *b, bool verbose, int indent = 0)
{
    for (int i = 0; i < indent; i++) cout << " ";
    cout << b->format_label() << ":";

    // Notes: (1) relative_energy is stored in a binary CM and refers
    //		  to the components
    //	      (2) escaping_cpt means that one or more components of
    //		  a binary CM is escaping
    //	      (3) stable refers to the CM

    if (b->is_parent()) {

	if (b->n_leaves() > 2) cout << " multiple";
	else if (b->n_leaves() == 2) cout << " binary";

	if (b->escaping_cpt)					 // STORY
	    cout << " escaping_components";

	real E = b->relative_energy;				 // STORY
	if (E < 0) {
	    cout << " bound_components E = " << E;
	    real sep = b->periastron_ratio;			 // STORY
	    if (sep > 0) cout << " (" << sep << ")";
	} else if (E > 0)
	    cout << " unbound_components E = " << E;

	if (b->stable && b->n_leaves() > 2)			 // STORY
	    cout << " stable";

    } else
	cout << " single";

    if (b->escaper)						 // STORY
	cout << " escaper";

    if (verbose && b->is_parent() && !b->is_root()) {
	hdyn2 *od = b->get_oldest_daughter();
	hdyn2 *yd = od->get_younger_sister();
	if (od && yd) {
	    cout << endl; PRI(indent);
	    print_orbital_elements(od, yd, false);
	}
    } else
	cout << endl;

    for_all_daughters(hdyn2, b, bb) print_recursive(bb, verbose, indent+2);
}

local inline void set_string_node_name(hdyn2 *b, char *s)
{
    if (b->get_name().size() > 0)
	strcpy(s, b->get_name().c_str());	// NB no size check...
    else if (b->get_index() >= 0)
	sprintf(s, "%d", b->get_index());
    else
	strcpy(s, "-");
}

local const char *string_index_copy(hdyn2 *n)
{
    if (n->get_name().size() == 0) {
	char integer_id[20];
	sprintf(integer_id, "%d", n->get_index());
	n->set_name(integer_id);
    }
    return n->get_name().c_str();
}

local void create_summary_strings(hdyn2 *b,
				  string& summary, string& this_state)
{
    // Construct summary and state strings describing the system.

    for_all_daughters(hdyn2, b, bi) {
	if (!bi->is_parent()) {
	    string index = string_index_copy(bi);
	    this_state += "s";
	    this_state += index;
	    this_state += " ";
	    if (bi->relative_energy >= 0) {			  // STORY
		summary += "<";
		summary += index;
		summary += "> ";
	    } else if (bi->stable) {				  // STORY
		summary += "[";
		summary += index;
		summary += "] ";
	    } else {
		summary += "(";
		summary += index;
		summary += ") ";
	    }
	} else {
	    int mult = bi->n_leaves();
	    if (mult==2) this_state += "b";
	    if (mult==3) this_state += "t";
	    if (mult==4) this_state += "q";
	    if (mult==5) this_state += "v";
	    if (mult >5) this_state += "x";
	    int low = 10000;
	    int last = 0;
	    int track = 0;
	    
	    while (track == 0) {
		track = 1;
		string next_num;
		for_all_leaves(hdyn2, bi, l) {
		    if (l->get_index() < low && l->get_index() > last) {
			low = l->get_index();
			track = 0;
			next_num = string_index_copy(l);
		    }
		}
		last = low;
		low = 10000;
		this_state += next_num;
	    }
      
	    this_state += " ";

	    if (bi->relative_energy >= 0)			 // STORY
		summary += "<";
	    else if (bi->stable || bi->n_leaves() <=2)		 // STORY
		summary += "[";
	    else
		summary += "(";     

	    hdyn2 *od = bi->get_oldest_daughter();
	    hdyn2 *yd = od->get_younger_sister();
	    int index_od = od->get_index();
	    int index_yd = yd->get_index();
	    string name1 = string_index_copy(od);
	    string name2 = string_index_copy(yd);
	    if ( index_od > 0 && index_yd > 0) {
		summary += name1;
		summary += ",";
		summary += name2;
	    } else if ( index_od > 0) {;
		summary += name1;
		summary += ",";
		summary += yd->get_name();
	    } else if ( index_yd > 0) {
		summary += od->get_name();
		summary += ",";
		summary += name2;
	    } else {;
		summary += od->get_name();
		summary += ",";
		summary += yd->get_name();
	    }

	    if (bi->relative_energy >= 0)			 // STORY
		summary += "> ";
	    else if (bi->stable || bi->n_leaves() <= 2)		 // STORY
		summary += "] ";
	    else
		summary += ") ";
	}
    }
}



local hdyn2 *get_tree2(hdyn *bin)
{
    // Determine the current structure and return the annotated tree.

    // Start by copying the tree to a flat hdyn2 tree; then identify
    // structure by proximity (actually potential) and dynamics (bound
    // or not).  Then classify bound systems by their stability
    // properties.

    hdyn2 *b = flat_copy_tree(bin);
    decompose(b);

    // Determine the binding energy of each node.  Note that the order
    // of for_all_nodes() means that this function and the next could
    // probably be merged, but don't assume that here.  Function sets
    // the relative_energy "story" value.
  
    for_all_daughters(hdyn2, b, bi) compute_relative_energy(bi);

    // Descend the tree to look for escapers (simply unbound, since
    // that is what the dynamical code will care about).  Note that we
    // use the softening from the external program in making this
    // decision.  Function sets the escaping_cpt and escaper "story"
    // values.

    flag_escapers(b);

    // Ascend the tree looking for widely separated binaries.
    // Function sets the stable and periastron_ratio "story" values.

    check_stability_master(b);

    return b;
}

// Manage quarantining of quasi-stable systems.  If a multiple cannot
// be confirmed as stable, but the interaction is otherwise over,
// declare the multiple stable if the system configuration survives
// without changing state for some specified period of time.
// Currently, we simply count calls to check_structure, but probably
// we want to follow the system for some number of outer multiple
// periods.

static string last_state;
static int state_count = 0;
static real state_time = 0;
const int min_qstable_count = 15;

local inline bool is_quasi_stable(hdyn2 *b)
{
    if (state_count >= min_qstable_count) {
	real max_period = 0;
	for_all_daughters(hdyn2, b, bi) {
	    if (!bi->stable) {					 // STORY
		hdyn2 *od = bi->get_oldest_daughter();
		hdyn2 *yd = od->get_younger_sister();
		kepler k;
		k.set_time(0);
		k.set_total_mass(od->get_mass()+yd->get_mass());
		k.set_rel_pos(yd->get_pos()-od->get_pos());
		k.set_rel_vel(yd->get_vel()-od->get_vel());
		k.initialize_from_pos_and_vel(true);
		if (k.get_period() > max_period) max_period = k.get_period();
	    }
	}
	return (b->get_system_time()-state_time > 10*max_period);
    }
    return false;
}

local inline bool is_over(hdyn2 *b, bool verbose)
{
    bool over = true;
    for_all_daughters(hdyn2, b, bi)
	if (!is_escaper(b, bi)) over = false;

    bool stable = true;
    for_all_daughters(hdyn2, b, bi)
	if (!bi->stable) stable = false;			 // STORY

    if (over) {
	if (stable) {
	    if (verbose) cout << "normal termination: "
			      << last_state << endl << flush;
	} else {
	    over = is_quasi_stable(b);
	    if (over) cout << "quase-stable termination: "
			   << last_state << endl << flush;
	}
    }

    return over;
}



//----------------------------------------------------------------------
// Externally visible functions:

bool check_structure(hdyn *bin,			// input root node
                     int verbose)		// default = 1
{
    if (verbose)
	cout << endl << "checking structure at time "
	     << bin->get_system_time() << endl;

    // Decompose the system.

    hdyn2 *b = get_tree2(bin);

    string summary, curr_state;
    create_summary_strings(b, summary, curr_state);

    // Print the results.

    if (verbose) {

	// Overall structure.

	print_recursive(b, verbose > 1);

	// Summary strings.

	cout << "summary: " << summary << endl;
	cout << "state:   " << curr_state << endl;
    }

    // Count the latest run of unchanged state strings.

    if (curr_state == last_state)
	state_count += 1;
    else {
	state_count = 0;
	state_time = bin->get_system_time();
    }
    last_state = curr_state;

    // See if the interaction is over.

    bool over = is_over(b, verbose);
    if (over && verbose) cout << endl << "interaction is over" << endl;

#if 0
    for_all_daughters(hdyn2, b, bb) {
	cout << bb->get_index() << " " << bb->get_pos()
	     << " " << bb->get_vel() << endl << flush;
    }
#endif

    // Clean up and exit.

    rmtree(b);
    return over;
}

hdyn *get_tree(hdyn *bin)
{
    return (hdyn*)get_tree2(bin);	// cast OK?
}

#else

//#include <string>
//#include <unistd.h>

template <class Q>
bool get_quantity(string s, const char *q, Q& value)
{
    size_t i = s.find(q);
    if (i != string::npos) {
	i = s.find("=", i);
	if (i != string::npos) {
	    s = s.substr(i+1);
	    istringstream ss(s);
	    ss >> value;
	    return true;
	}
    }
    return false;
}

void initialize(hdyn *b)		// b = root node
{
    int n = 0;
    hdyn *bp = NULL, *bb = NULL;
    string s;

    while (getline(cin, s)) {
	if (s[0] != ';') {
	    istringstream ss(s);
	    int i = -1;
	    real m;
	    vec x, v;
	    ss >> i >> m >> x >> v;
	    if (i < 0) break;
	    bb = new hdyn();
	    n++;
	    if (bp) {
		bp->set_younger_sister(bb);
		bb->set_older_sister(bp);
	    } else
		b->set_oldest_daughter(bb);
	    bb->set_parent(b);
	    bb->set_index(i);
	    bb->set_mass(m);
	    bb->set_pos(x);
	    bb->set_vel(v);
	    bp = bb;
	} else {
	    real t;
	    if (get_quantity(s, "system_time", t)) b->set_system_time(t);
	    int seed;
	    if (get_quantity(s, "seed", seed)) {PRL(seed); b->set_seed(seed);}
	}
    }
    PRC(n); PRL(b->get_system_time());
}

int main(int argc, char *argv[])
{
    hdyn *b = new hdyn;
    initialize(b);
    check_structure(b, true);
    return 0;
}

#endif
