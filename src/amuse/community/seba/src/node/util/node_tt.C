
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// node_tt.C : basic tree handling functions
//

#include "node.h"

void node::null_pointers()
{
    // Clear all pointers (don't touch what they point to)
    // -- for use in cleaning up temporary nodes...  Careful!!

    name = NULL;
    parent = oldest_daughter = elder_sister = younger_sister = NULL;
    log_story = dyn_story = NULL;
    hbase = NULL;
    sbase = NULL;
}

int node::n_daughters() const
{
    if(oldest_daughter == NULL){
	return 0;
    }else{
	int  n = 0;
	for ( node * daughter = get_oldest_daughter();
	      daughter != NULL;
	      daughter = daughter->get_younger_sister() )
	    n++;
	return n;
    }
}

int node::n_leaves() const
{
    if(oldest_daughter == NULL){
	return 1;
    }else{
	int  n = 0;
	for ( node * daughter = get_oldest_daughter();
	      daughter != NULL;
	      daughter = daughter->get_younger_sister() )
	    n += daughter->n_leaves();
	return n;
    }
}

bool node::is_grandparent() const
{
    for_all_daughters(node, const_cast<node*>(this), n)
	if (n->is_parent()) return true;

    return false;
}

// The following four functions are now inlined and defined in node.h
// (and use the new root member data) -- Steve 9/19/98

// bool node::is_top_level_node()
// {
//     if (parent == NULL) return FALSE;
//     if (parent->get_parent()) return FALSE;
//     return TRUE;
// }

// bool node::is_top_level_leaf()
// {
//     if (parent == NULL) return FALSE;
//     if (parent->get_parent()) return FALSE;
//     if (oldest_daughter) return FALSE;
//     return TRUE;
// }

// bool node::is_low_level_node()
// {
//     if (parent == NULL) return FALSE;
//     if (parent->get_parent() == NULL) return FALSE;
//     return TRUE;
// }

// bool node::is_low_level_leaf()
// {
//     if (parent == NULL) return FALSE;
//     if (parent->get_parent() == NULL) return FALSE;
//     if (oldest_daughter) return FALSE;
//     return TRUE;
// }

//node* node::get_root()
//{
//    if (root) return root;
//
//    // If root not already set, set it now for future use.
//
//    if (parent == NULL)
//	root = this;
//    else
//	root = get_top_level_node()->get_parent();
//
//    return root;
//}

node* node::get_binary_sister()
{
    bool err = false;
    node * sister;
    if (elder_sister) {
	if (elder_sister->get_elder_sister()) {
	    err = true;
	} else {
	    if (younger_sister) {
		err = true;
	    } else {
		sister = elder_sister;
	    }
	}
    } else {
	if (younger_sister) {
	    if (younger_sister->get_younger_sister()) {
	        err = true;
	    } else {
		sister = younger_sister;
	    }
	} else {
	    // err_exit("get_binary_sister: no sisters");
	    warning("get_binary_sister: no sisters!");
	    PRL(format_label());
	    PRL(get_parent()->format_label());
	    if (get_parent() != get_root()) pp2(get_parent(), cerr);
	    return NULL;
	}
	
    }
    if (err) {
	warning("get_binary_sister: too many sisters!");
	PRL(format_label());
	PRL(get_parent()->format_label());
	if (get_parent() != get_root()) pp2(get_parent(), cerr);
	return NULL;
    }
    return sister;
}

real total_mass(node * b)	// should be a member function?
{
    real total_mass_nodes = 0;
    for_all_daughters(node, b, bi)
	total_mass_nodes += bi->get_mass();
    return total_mass_nodes;
}

node * mk_flat_tree(int n, npfp new_node_type, hbpfp a_hbpfp, sbpfp a_sbpfp,
		    bool use_stories)
{
    node * root, * by, * bo;

    root = (* new_node_type)(a_hbpfp, a_sbpfp, use_stories);
    root->set_root(root);
    bo = (* new_node_type)(a_hbpfp, a_sbpfp, use_stories);
    root->set_oldest_daughter(bo);
    bo->set_parent(root);
    bo->set_label(1);

    for (int i = 1; i < n; i++) {

        by = (* new_node_type)(a_hbpfp, a_sbpfp, use_stories);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
	by->set_parent(root);
	by->set_label(i+1);
        bo = by;
    }

    return root;
}

char * string_index_of_node(node * n)
{
    char * tmp;
    tmp = n->get_name();
    if (tmp != NULL) {
	return tmp;
    } else {
	static char integer_id[20];
	sprintf(integer_id, "%d", n->get_index());
	n->set_label(integer_id);
	return n->get_name();
    }
}

char * construct_binary_label(node * ni, node * nj)
{
    static char new_name[256];
    sprintf(new_name,"(%s,%s)",string_index_of_node(ni),
	    string_index_of_node(nj));
    return new_name;
}

void label_binary_node(node* n)
{
    if (n->is_root()) return;

    for_all_daughters(node, n, nn) label_binary_node(nn);

    if (!n->is_leaf()) {
        node* od = n->get_oldest_daughter();
        node* yd = od->get_younger_sister();
        n->set_label(construct_binary_label(od, yd));
    }
}

// Default label for merger of "a" and "b" is simply "a+b".
// This may lead to unwieldy names in case of multiple mergers.
//
// Alternative is to construct a name that retains the identity
// of the massive component, and merely indicates how many other
// nodes have been destroyed -- really only works well in case
// of dominant motion... (Steve 12/98)
//
// Scheme:	1st merger (b)	-->  a+b (a more massive)
//		2nd merger (c)	-->  a<+2> (a+b more massive)
//				     c<+2> (c more massive)
//
// Make the alternative the default for now (12/98) -- return later
// to allow options to be passed as arguments.

local char * construct_merger_label(node * ni, node * nj)
{
    static char new_name[1024];
    sprintf(new_name, "%s+%s",
	    string_index_of_node(ni),
	    string_index_of_node(nj));
    return new_name;
}

local char* string_header(char* string)
{
    // Return everything preceding the "<", or everything preceding
    // the "+", or else the entire string.

    static char header[1024];
    strcpy(header, string);

    char *c = strchr(header, '<');	// seek "a<+n>"
    if (!c) {
	c = strchr(header, '+');		// seek "a+b"
	if (!c) return header;
    }

    *c = '\0';
    return header;
}

local int string_count(char* string)
{
    // Return the integer between the "<" and the ">", or else 1, if
    // the string is of the form "a+b", or else 0.

    static char header[1024];
    strcpy(header, string);

    char* c1 = strchr(header, '<');
    if (!c1) {
	if (strchr(header, '+')) {

	    // return 1;

	    // Would be OK to return 1 here if naming were consistent.
	    // Better to count "+" signs, just in case.

	    c1 = header;
	    int count = 0;
	    while (*c1 > '\0') {
		if (*c1 == '+') count++;
		c1++;
	    }

	    return count;

	} else
	    return 0;
    }

    char* c2 = strchr(header, '>');
    if (!c2) return 0;			// (actually an error)
    *c2 = '\0';

    return atoi(c1+1);
}

local char * alt_construct_merger_label(node * ni, node * nj)
{
    // Names are of the form:
    //
    //		ni:	"a" or "a+x" or "a<+n>"
    //		nj:	"b" or "b+y" or "b<+m>".
    //
    // Decode the names and combine them to make the new name.
    // Not very efficiently written, but... (Steve, 12/98)

    static char new_name[1024];
    int count_i = string_count(string_index_of_node(ni));
    int count_j = string_count(string_index_of_node(nj));

    if (count_i+count_j == 0)
	return(construct_merger_label(ni, nj));

    sprintf(new_name, "%s<+%d>",
	    string_header(string_index_of_node(ni)),
	    count_i + count_j + 1);
    return new_name;
}

void label_merger_node(node* n)
{
    node* od = n->get_oldest_daughter();
    node* yd = od->get_younger_sister();

    // Added by Steve, 12/98: new label lists more massive
    // component first.

    if (od->get_mass() < yd->get_mass()) {
	node* tmp = od;
	od = yd;
	yd = tmp;
    }

    n->set_label(alt_construct_merger_label(od, yd));

//               ^^^^	(note)

}

int depth_of_node(node * ni)
{
    int depth = 0;
    while((ni = ni->get_parent()) != NULL)depth ++;
    return depth;
}

// is_descendent_of: return TRUE if child is among the offspring of parent.
//                   mode = 1 ==> include parent in "offspring" list
//                   mode = 0 ==> exclude parent from consideration
//
int is_descendent_of(node *child, node *parent, int mode)
{
    if (mode == 0 && child == parent) return 0;

    while (child != NULL) {
	if(child == parent) {
	    return 1;
	} else {
	    child = child->get_parent();
	}
    }
    return 0;
}

node * common_ancestor(node * ni, node * nj)
{
    int i;
    real difference = depth_of_node(ni) - depth_of_node(nj);
    if(difference > 0){
	for(i = 0; i<difference; i++) ni = ni->get_parent();
    }else if (difference < 0){
	for(i = 0; i<-difference; i++) nj = nj->get_parent();
    }
    while(ni != nj){
	ni = ni->get_parent();
	nj = nj->get_parent();
    }
    return ni;
}

node * node_with_index(int i, node * top)	// recursive; default top = NULL
{
    if (!top) return NULL;
    node * n = top;

    if (n->get_index() == i) return n;

    if (n->get_oldest_daughter() != NULL) {
	for_all_daughters(node, n, nn) {
	    node * tmp = node_with_index(i, nn) ;
	    if (tmp != NULL) return tmp;
	}
    }
    return NULL;
}

node * node_with_name(char* s, node * top)	// recursive; default top = NULL
{
    if (!top) return NULL;
    node * n = top;

    // cerr << "node_with_name: "; PRC(s); PRL(n->format_label());

    if (streq(n->format_label(), s)) return n;

    if (n->get_oldest_daughter() != NULL) {
	for_all_daughters(node, n, nn) {
	    node * tmp = node_with_name(s, nn) ;
	    if (tmp != NULL) return tmp;
	}
    }
    return NULL;
}

// detach_node_from_general_tree
// Detach node n from the tree.  Do not check whether the parent has
// more than two remaining daughters.  Do not correct the parent mass.

void detach_node_from_general_tree(node *n)
{
    if (n == NULL) {
	cerr << "detach_node_from_general_tree: n is NULL" << endl;
	return;
    }

    node *parent = n->get_parent();
    
    // Check if n is head without parent or sisters.

    if (parent == NULL) {
	cerr << "detach_node_from_general_tree: n has no parent" << endl;
	return;
    }

    n->set_parent(NULL);

    node *elder_sister = n->get_elder_sister();
    node *younger_sister = n->get_younger_sister();

    if (parent->get_oldest_daughter() == n)
	parent->set_oldest_daughter(younger_sister);

    if (elder_sister) elder_sister->set_younger_sister(younger_sister);
    if (younger_sister)	younger_sister->set_elder_sister(elder_sister);
}

// remove_node_with_one_daughter
// remove the node from the tree
// and replace it by its daughter
// the node should be deleted afterwords
// In principle, virtual destructor should
// work, but I'm not sure what is actually
// done if delete is called in this
// function...
void remove_node_with_one_daughter(node *n)
{
    if (n->get_oldest_daughter()->get_younger_sister()) {
	cerr << "remove_node_with_one_daughter: #daughters of n is "
	     << "larger than 1\n";
	exit(1);
    }
    node *daughter = n->get_oldest_daughter();
    if (daughter->get_elder_sister() ||
	daughter->get_younger_sister()) {
	cerr << "remove_node_with_one_daughter: daughter of n has sisters\n";
	exit(1);
    }

    node *parent = n->get_parent();
    
    // check if n is head without parent or sisters

    if (parent == NULL) {
	cerr << "Warning: remove_node_with_one_daughter, n has no parent\n";
	return;
    }

    node *elder_sister = n->get_elder_sister();
    node *younger_sister = n->get_younger_sister();

    if (parent->get_oldest_daughter() == n) {
	parent->set_oldest_daughter(daughter);
    }

    if (elder_sister) {
	elder_sister->set_younger_sister(daughter);
    }
    if (younger_sister) {
	younger_sister->set_elder_sister(daughter);
    }

    daughter->set_parent(parent);
    daughter->set_younger_sister(younger_sister);
    daughter->set_elder_sister(elder_sister);

    // cerr << "remove_node_with_one_daughter: ";
    // cerr << "garbage collection not yet implemented\n";
}

// detach_node_from_binary_tree
// delete node n and repair the tree

void detach_node_from_binary_tree(node *n)
{
    if (n == NULL) {
	cerr << "Warning: detach_node_from_binary_tree, n is NULL\n";
	return;
    }

    node *parent = n->get_parent();

    detach_node_from_general_tree(n);
    remove_node_with_one_daughter(parent);
}
    
// extend_tree
// Extend a tree by inserting a new node to a location presently
// occupied by old_n. old_n becomes the only child of the newly
// inserted node.  This function returns the address of the
// newly created node
//

void extend_tree(node *old_n, node *new_n)
{
    if (old_n == NULL) {
	cerr << "Warning:extend_tree, old_n is NULL\n";
	return;
    }

    // check if old_n is root

    if (old_n->get_parent() == NULL) {
	cerr << "extend_tree, old_n is root, cannot insert\n";
	exit(1);
    }

    // set the pointer of ex-parent of old_n if she was the oldest daughter

    node *parent = old_n->get_parent();
    if (parent->get_oldest_daughter() == old_n) {
	// old_n is the oldest daughter of her parent
	parent->set_oldest_daughter(new_n);
    }
    
    // set the pointers of ex-sisters of old_n
    node *elder = old_n->get_elder_sister();
    if (elder) {
	elder->set_younger_sister(new_n);
    }
    node *younger = old_n->get_younger_sister();
    if (younger) {
	younger->set_elder_sister(new_n);
    }
    
    new_n->set_parent(old_n->get_parent());
    new_n->set_elder_sister(old_n->get_elder_sister());
    new_n->set_younger_sister(old_n->get_younger_sister());
    new_n->set_oldest_daughter(old_n);

    old_n->set_elder_sister(NULL);
    old_n->set_younger_sister(NULL);
    old_n->set_parent(new_n);
}

// add_node
// insert n into the tree as the oldest_daughter of the
// parent. The ex-oldest node of the parent becomes the
// younger sister of n.

void add_node(node *n, node *parent)
{
    if (n == NULL) {
	cerr << "Warning:add_node, n is NULL\n";
	return;
    }

    if (parent == NULL) {
	cerr << "Warning:add_node, parent is NULL\n";
	return;
    }

// The following part tests if n is completely isolated.

#if 0
    // Check if n is head without parent or sisters

    if (n->get_parent()) {
	cerr << "add_node, n has parent\n";
	exit(1);
    }
    if (n->get_elder_sister()) {
	cerr << "add_node, n has an elder sister\n";
	exit(1);
    }
    if (n->get_younger_sister()) {
	cerr << "add_node, n has an younger sister\n";
	exit(1);
    }
#endif

    // Set the pointers of ex-oldest-daughter of parent.

    node *ex = parent->get_oldest_daughter();
    if (ex) {
	ex->set_elder_sister(n);
    }

    parent->set_oldest_daughter(n);

    n->set_elder_sister(NULL);
    n->set_younger_sister(ex);
    n->set_parent(parent);
}

// add_node_before: insert the given node n into the tree before a
// specified node m.

void add_node_before(node *n, node *m)
{
    if (!n || !m) return;
    
    node *p = m->get_parent();
    if (!p) return;
    n->set_parent(p);
    if (p->get_oldest_daughter() == m) p->set_oldest_daughter(n);

    node *elder_sister = m->get_elder_sister();
    n->set_elder_sister(elder_sister);
    m->set_elder_sister(n);

    if (elder_sister) elder_sister->set_younger_sister(n);
    n->set_younger_sister(m);
}

void insert_node_into_binary_tree(node *n, node *old_n, node *new_n)
{
    extend_tree(old_n, new_n);
    add_node(n, new_n);
}

// next_node: return the next node to look at when traversing the tree
//            below node b

node* node::orig_next_node(node *b)
{

    if (oldest_daughter)
	return oldest_daughter;
    else if (younger_sister)
	return younger_sister;
    else {

	if (this == b) return NULL;		// In case b is a leaf...
	if (parent == b) return NULL;
	if (parent == NULL) return NULL;	// Just in case b is atomic...

	node* tmp = parent;
	while (tmp->get_younger_sister() == NULL)
	    if (tmp != b)
		tmp = tmp->get_parent();
	    else
		return NULL;
	
	return tmp->get_younger_sister();
    }
    return NULL;	// To keep some compilers happy... 
}

// Note from Steve 7/9/98.  This can fail if base is not the root node.
// An attempt to traverse only the tree below a top-level leaf will include
// all nodes to the right of base.

// This function is used in many places -- changes made here may have
// unexpected consequences elsewhere!  If so, fix is to rename and use
// the new version in in dstar_kira.C and elsewhere

#define DEBUG 0

node* node::next_node(node* base)
{
    if (DEBUG) {
	cerr << "next_node:  node "; print_label(cerr); cerr << endl;
	cerr << "            base "; base->print_label(cerr); cerr << endl;
    }

    if (oldest_daughter ) {

	if (DEBUG) {
	    cerr << "return od ";
	    oldest_daughter->print_label(cerr);
	    cerr << endl;
	}
	return oldest_daughter;


    // NEW!

    } else if (this == base) {		// In case base is a leaf...

	if (DEBUG) cerr << "return NULL 0\n";
	return NULL;


    } else if (younger_sister) {

	if (DEBUG) {
	    cerr << "return ys ";
	    younger_sister->print_label(cerr);
	    cerr << endl;
	}
	return younger_sister;

    } else {

	// Can't go down or right.  See if we can go up.

	if (this == base) {		// In case base is a leaf...
	    if (DEBUG) cerr << "return NULL 1\n";
	    return NULL;
	}

	if (parent == NULL) {		// In case b is atomic...
	    if (DEBUG) cerr << "return NULL 2\n";
	    return NULL;
	}

	node* tmp = parent;
	while (tmp->get_younger_sister() == NULL) {
	    if (tmp == base) {
		if (DEBUG) cerr << "return NULL 3\n";
		return NULL;
	    } else
		tmp = tmp->get_parent();
	}

	// Now tmp is the lowest-level ancestor with a younger sister.

	if (tmp == base) {
	    if (DEBUG) cerr << "return NULL 4\n";
	    return NULL;
	} else {
	    if (DEBUG) {
		cerr << "return tys ";
		tmp->get_younger_sister()->print_label(cerr);
		cerr << endl;
	    }
	    return tmp->get_younger_sister();
	}

    }

    if (DEBUG) cerr << "return NULL 5\n";
    return NULL;	// To keep some compilers happy... 
}

void  err_exit(char * line, node* n)
{
    cerr << "error: ";
    print_message(line);
    if (n) rmtree(n);
    exit(1);
}
