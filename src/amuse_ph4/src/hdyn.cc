
#include "hdyn.h"

// Initialization of hdyn static variables.

real hdyn::system_time = -1;
real hdyn::eta = -1;
int  hdyn::n_iter = 2;		// default of 2 seems to work well
int  hdyn::seed = 0;
real hdyn::dt_crit = 1.e-3;
real hdyn::r2_crit = 1.e-2;
bool hdyn::allow_full_unperturbed = true;
real hdyn::gamma2_unpert = 1.e-12;
real hdyn::gamma_inv3 = pow(fmax(gamma2_unpert, 1.e-120), -1./6);
int hdyn::cm_index = -1;

// Tree management:

hdyn* hdyn::next_node(hdyn* base)
{
    // Return the next node in the tree traversal.
    // Order: go down, right, or up.

    if (oldest_daughter) return oldest_daughter;	// Down.
    else if (this == base) return NULL;
    else if (younger_sister) return younger_sister;	// Right.
    else {

	// Up.

	if (parent == NULL) return NULL;

	hdyn* tmp = parent;
	while (tmp->get_younger_sister() == NULL) {
	    if (tmp == base) return NULL;
	    else tmp = tmp->get_parent();
	}

	// Now tmp is the lowest-level ancestor with a younger sister.

	if (tmp == base) return NULL;
	else return tmp->get_younger_sister();
    }

    return NULL;
}

// add_node: Insert node n into the tree before (i.e. as the older
//	     sister of) node m.  Set pointers only.  Don't adjust
//	     masses or dynamics.

void add_node(hdyn *n, hdyn *m)
{
    // const char *func = "add_node";
    if (!n || !m) return;
    
   hdyn *p = m->get_parent();
    if (!p) return;
    n->set_parent(p);
    if (p->get_oldest_daughter() == m) p->set_oldest_daughter(n);

   hdyn *older_sister = m->get_older_sister();
    n->set_older_sister(older_sister);
    m->set_older_sister(n);

    if (older_sister) older_sister->set_younger_sister(n);
    n->set_younger_sister(m);
}

// detach_node: Detach node n from the tree.  Do not check whether the
//		parent has more than two remaining daughters.  Set
//		pointers only.  Don't adjust dynamics, and don't
//		correct the parent mass.

void detach_node(hdyn *n)
{
    const char *func = "detach_node";

    if (!n) {
	cerr << func << ": n is NULL" << endl << flush;
	return;
    }

    hdyn *parent = n->get_parent();
    
    // Check if n is head without parent or sisters.

    if (!parent) {
	cerr << func << ": n has no parent" << endl << flush;
	return;
    }

    n->set_parent(NULL);

    hdyn *older_sister = n->get_older_sister();
    hdyn *younger_sister = n->get_younger_sister();

    if (parent->get_oldest_daughter() == n)
	parent->set_oldest_daughter(younger_sister);

    if (older_sister) older_sister->set_younger_sister(younger_sister);
    if (younger_sister)	younger_sister->set_older_sister(older_sister);

    n->set_older_sister(NULL);
    n->set_younger_sister(NULL);
}

void create_binary_node(hdyn *cm, hdyn *bi, hdyn *bj)
{
    // Create a binary consisting of nodes bi and bj.  CM node cm
    // replaces bi in the tree.  Set pointers only.  Don't adjust
    // dynamics or masses.

    add_node(cm, bi);
    detach_node(bi);
    detach_node(bj);
    cm->set_oldest_daughter(bi);
    bi->set_parent(cm);
    bj->set_parent(cm);
    bi->set_older_sister(NULL);
    bi->set_younger_sister(bj);
    bj->set_older_sister(bi);
    bj->set_younger_sister(NULL);
}

void rmtree(hdyn* b,
	    bool delete_b)	// default = true
{
    hdyn* d = b->get_oldest_daughter();
    while (d) {
	hdyn* tmp = d->get_younger_sister();
	rmtree(d);
	d = tmp;
    }

    if (delete_b) delete b;	// optionally leave node itself untouched
}

void pp2(hdyn *b,
	 int level)		// default = 0
{
    for (int i = 0; i<level*2; i++) {cout << " ";}
    cout << b->get_index() << " " << b->get_mass()
	 << " " << b->get_pos() << endl;
    for_all_daughters(hdyn, b, daughter)
	pp2(daughter, level + 1);	
}


// Prediction:

void predict_loworder_all(hdyn *b, real t) {
    if (t != b->get_t_pred()) {
	for_all_daughters(hdyn, b, bi) bi->predict_loworder(t);
	b->set_t_pred(t);
    }
}

void hdyn::copy_data_from(hdyn *from_b)
{
    if (from_b == NULL) return;

    // Copy all dynamic data from from_k to this instance.
    // Is there an easier way to do this...?

    index = from_b->index;
    mass = from_b->mass;
    radius = from_b->radius;
    pot = from_b->pot;
    t_pred = from_b->t_pred;
    fully_unperturbed = from_b->fully_unperturbed;
    pos = from_b->pos;
    vel = from_b->vel;
    pred_pos = from_b->pred_pos;
    pred_vel = from_b-> pred_vel;
    acc = from_b->acc;
    jerk = from_b->jerk;
    old_acc = from_b->old_acc;
    old_jerk = from_b->old_jerk;

    if (from_b->kep != NULL) {
	kepler *k = new kepler();
	k->copy_data_from(from_b->kep);
	kep = k;
    }
}
