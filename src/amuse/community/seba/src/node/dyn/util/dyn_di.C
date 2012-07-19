
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// dyn_di.C: dyn-specific diagnostic functions.
//

#include "dyn.h"

void dbg_message(const char* s, dyn* b) {
#ifdef DEBUG    
    cerr << s << " ";
    b->pretty_print_node(cerr); cerr << "\n";
#else
    const char* dummy_char = s;   	// to keep the compiler happy
    dyn* dummy = b;       		// to keep the compiler happy
#endif    
}

void dbg_message(const char* s, dyn* bj, dyn *bi) {
#ifdef DEBUG    
    cerr << s << " ";
    bj->pretty_print_node(cerr); cerr << "->";
    bi->pretty_print_node(cerr); cerr << "\n";
#else
    const char* dummy_char = s;		// to keep the compiler happy
    dyn* dummy_i = bi;       		// to keep the compiler happy
    dyn* dummy_j = bj;       		// to keep the compiler happy
#endif    
}

local inline void accumulate_pot_on_work(real mass,
					 vec d_pos,
//					 vec d_vel,	// not used
					 real eps2,
					 real & p) {
    real r2inv = 1.0 / (d_pos*d_pos + eps2);
    real mrinv = mass * sqrt(r2inv);
    p -= mrinv;

//  silly trick to make the compiler shut up about "warning:  d_vel not used":
//    vec dummy = d_vel;

//    cout << "work p = " << p << " \n";    
}

local void tree_calculate_pot(dyn * bj,
			      dyn *mask,
			      dyn *bi,
			      vec offset_pos,
//			      vec offset_vel,
			      real eps2,
			      real & p)
{
    if (bj == mask) return;

    if (bj->get_oldest_daughter()) {
	for_all_daughters(dyn, bj, bb)
	    tree_calculate_pot(bb, mask, bi,
			       offset_pos + bb->get_pos(),
//			       offset_vel + bb->get_vel(),
			       eps2, p);			// recursion
    } else {
	if (bj != bi)
	    accumulate_pot_on_work(bj->get_mass(), offset_pos,
//				   offset_vel,
				   eps2, p);
    }
}


local void calculate_pot_on_leaf(dyn * bj,
				 dyn *mask,
				 dyn *bi,
				 real eps2,
				 real & p)
{
    // dbg_message("calculate_pot_on_leaf", bi);

    vec d_pos;
//    vec d_vel;
    for(dyn * b = bi; b != bj; b = b->get_parent()){
	d_pos -= b->get_pos();
//	d_vel -= b->get_vel();
    }
    tree_calculate_pot(bj, mask, bi,
		       d_pos,
//		       d_vel,
		       eps2, p);
}

local real pot_on_node(dyn * bj,		// root node for computation
		       dyn * mask,		// exclude all nodes under mask
		       dyn * bi,		// compute potential on node bi
		       real eps2)		// softening parameter
{
    // Compute the potential per unit mass on node bi due to all nodes
    // under node bj, excluding nodes under node mask.

    real p = 0;
    real mtot = 0;

    if (bi->get_oldest_daughter() == NULL) {

	calculate_pot_on_leaf(bj, mask, bi, eps2, p);

    } else {

	// Potential is sum of daughter potentials.

	real p_daughter;
	for_all_daughters(dyn, bi, bd) {
	    real m_daughter = bd->get_mass();
	    p_daughter = pot_on_node(bj, mask, bd, eps2);	// recursion
	    p += m_daughter * p_daughter;
	    mtot += m_daughter;
	}

	real minv = 1/mtot;
	p *= minv;
    }
    return p;
}

local real pot_on_low_level_node(dyn * bi, real eps2)
{
    // Compute the potential on bi due to its binary sister only.

    dyn * parent = bi->get_parent();

    if (parent->get_oldest_daughter()
	      ->get_younger_sister()->get_younger_sister())
	err_exit("pot_on_low_level_node: not a binary node!");

    dyn * sister = bi->get_elder_sister();
    if (sister == NULL) sister = bi->get_younger_sister();
    
    return pot_on_node(parent, bi, bi, eps2);
}

local inline real cm_pot_on_node(dyn * bj, dyn * bi, real eps2)
{
    // Compute the potential per unit mass on node bi due to all daughters
    // of node bj, using the center of mass approximation throughout.

    real p = 0;
    for_all_daughters(dyn, bj, bb)
	if (bb != bi)
	    accumulate_pot_on_work(bb->get_mass(),
				   bb->get_pos() - bi->get_pos(),
//				   bb->get_vel() - bi->get_vel(),
				   eps2, p);
    return p;
}

real pot_on_general_node(dyn * bj, dyn * bi, real eps2, bool cm = false)
{
    // Calculate potential (per unit mass) on top-level node bi due to
    // all nodes under node bj, or the potential on low-level node bi due
    // to its binary sister.  If cm = true, then bj is the root node, bi
    // is a top-level node, and the center of mass approximation is to be
    // used.

    if (bi->is_top_level_node()) {

	// Silly to carry the cm flag through the many levels of recursion
	// used in the general case.  Just compute the CM pot as a special
	// case.  Do this here, rather than in accumulate_energies, because
	// the hdyn version of this function will want to set the value
	// of pot on return.

	if (cm)
	    return cm_pot_on_node(bj, bi, eps2);
	else
	    return pot_on_node(bj, bi, bi, eps2);

    } else

	return pot_on_low_level_node(bi, eps2);
}

// Note that the kinetic energy is defined in the "root" frame
// -- i.e. it does not include any bulk motion associated with the
// root node itself.  (Usually what we want in practice, but don't
// forget the CM kinetic energy for diagnostics!)

local void accumulate_energies(dyn * root, dyn * b, real eps2,
			       real & epot, real & ekin, real & etot,
			       bool cm = false)
{
    if (b->get_oldest_daughter()) {

	// Start/continue the recursion.  In the cm = true case,
	// expand daughters only if root really is root.

	if (!cm || b->is_root())
	    for_all_daughters(dyn, b, bb)
		accumulate_energies(root, bb, eps2, epot, ekin, etot, cm);

	etot = epot + ekin;
    }

    // Do the work.

    if (!b->is_root()) {
	real m2 = 0.5 * b->get_mass();
	epot += m2 * pot_on_general_node(root, b, eps2, cm);
	ekin += m2 * square(b->get_vel());
    }
}

void calculate_energies(dyn * root, real eps2,
			real & epot, real & ekin, real & etot,
			bool cm)	// default = false
{
    epot = ekin = etot = 0;
    accumulate_energies(root, root, eps2, epot, ekin, etot, cm);
}

static real initial_etot = 0;

void print_recalculated_energies(dyn * b,
				 int mode,	// default = 0
				 real eps2,	// default = 0
				 real e_corr)	// default = 0
{
    real epot = 0;
    real ekin = 0;
    real etot = 0;

    if (!b->get_ignore_internal())
	accumulate_energies(b, b, eps2, epot, ekin, etot);

    cerr << "Energies: " << epot << " " << ekin << " " << etot
	 << " " << initial_etot;

    if (mode) {

	// real de = (etot - e_corr - initial_etot) / epot;

	// The above normalization fails if many stars have escaped
	// the system with HUGE kicks...

	real de = (etot - e_corr - initial_etot) / initial_etot;
	cerr << " de = " << de << endl;

    } else {
	initial_etot = etot;
	cerr << "\n";
    }	
}

void initialize_print_energies(dyn * b,
			       real eps2)	// default = 0
{
    real epot = 0;
    real ekin = 0;
    real etot = 0;

    if (!b->get_ignore_internal())
	accumulate_energies(b, b, eps2, epot, ekin, etot);

    initial_etot = etot;
}
