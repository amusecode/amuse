
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~


/// @file node.h  Defines base class for N-body systems.

//  version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
//  version 2:
//
//  This file includes:
//  1) definition of class node


#ifndef  STARLAB_NODE_H
#  define  STARLAB_NODE_H

#include  "starlab_vector.h"
#include  "story.h"
#include  "hydrobase.h"
#include  "starbase.h"

#define __VALID_NODE__		123456789
#define __INVALID_NODE__	-1

/// \a node: Base class for the nodes in a tree of dynamical particles.

/// Nodes have a name and a mass, together with pointers defining the
/// data tree describing the N-body system, starbase and hydrobase
/// pointers allowing extension to stellar and hydrodynamical
/// applications, and log and dyn stories for additional run-time
/// variables.

class node
    {
    protected:             // not private, since other classes, such as dyn
	                   // and hdyn should have access to these data.

        // Global variable!

    	static node* root;	  ///< Address of the root node.


	long int node_flag;	  ///< "Magic number" identifying valid node.
				  // Experimental attempt to check for
				  // and avoid deleted nodes (Steve, 8/98)
				  //  - set only in node constructor
				  //  - reset only in node destructor

	int  index;               ///< Nodes can be numbered
	char * name;              ///< Or nodes can receive individual names.

	real  mass;               ///< Node mass (= sum of daughter masses).

	node * parent;            // oya       | As you can see, Japanese is
	node * oldest_daughter;   // choujo    | a more natural language for
	node * elder_sister;      // oneesan   | building a tree structure;
	node * younger_sister;    // imouto    | the names are so much simpler!

	/// Class underlying all hydrodynamic classes.
	                  
        hydrobase * hbase;        // This is the only hydro class provided 
	                          // here, in order to allow separate
	                          // recompilation of the inherited classes
				  // based on hydro without necessitating
	                          // compilation of the dynamics part.

	/// Class underlying all stellar classes.
	                  
        starbase * sbase;         // This is the only star class provided 
	                          // here, in order to allow separate
	                          // recompilation of the inherited classes
				  // based on star without necessitating
	                          // compilation of the dynamics part.

        // The log story is a generalized scratchpad for notes and log entries.

	story * log_story;

	// The dyn story is a placeholder for unrecognized dynamical data.
                                  
	story * dyn_story;        // allows information to be preserved
				  // and passed down a pipe

    public:

	inline void clear_node() {
	    if (name) delete [] name;
	    parent = oldest_daughter = elder_sister = younger_sister = NULL;
	}

	inline void node_init() {
	    node_flag = __VALID_NODE__;
	    index = -1;     // < 0 to show that no number has yet been assigned
	    mass = 1;       // to allow star and hydro to rescale this value
	    name = NULL;
	    clear_node();
	}

	inline void node_set_stories(bool use_stories) {

	    // Potentially very dangerous to set any of the following NULL!

	    if (use_stories) {
		log_story = mk_story_chapter(LOG_ID);
		dyn_story = mk_story_chapter(DYNAMICS_ID);
	    } else {
		log_story = NULL;
		dyn_story = NULL;
	    }
	}

	inline void node_set_hbase(hbpfp the_hbpfp) {
	    if (the_hbpfp)
		hbase = (*the_hbpfp)();
	    else
		hbase = NULL;
	}

	inline void node_set_sbase(sbpfp the_sbpfp) {
	    if (the_sbpfp) {
		sbase = (*the_sbpfp)();
		sbase->set_node(this);
	    } else
		sbase = NULL;
	}

	/// Default node constructor includes minimal star and hydro support.

	node(hbpfp the_hbpfp = new_hydrobase,
	     sbpfp the_sbpfp = new_starbase,
	     bool use_stories = true) {
	    node_init();
	    node_set_stories(use_stories);
	    node_set_hbase(the_hbpfp);
	    node_set_sbase(the_sbpfp);
	}

	/// Delete the node's stories.

	inline void rmstory() {
	    if (log_story) {
		delete log_story;
		log_story = NULL;
	    }
	    if (dyn_story) {
		delete dyn_story;
		dyn_story = NULL;
	    }
	}

	/// Delete the star base.

	inline void rmstarbase() {
	    if (sbase) {
		delete sbase;
		sbase = NULL;
	    }
	}

	/// Delete the hydro base.

	inline void rmhydrobase() {
	    if (hbase) {
		delete hbase;
		hbase = NULL;
	    }
	}

	virtual ~node() {
	    node_flag = __INVALID_NODE__;
	    if (name) delete [] name;
	    rmstory();
	    rmhydrobase();
	    rmstarbase();
	    if (this == root) root = NULL;
	}

	/// Test the node's "magic number" to check validity.

	inline bool is_valid() const
	    {return (node_flag == __VALID_NODE__);}

	/// Flag the node as invalid.

	inline void set_invalid()
	    {node_flag = __INVALID_NODE__;}

	/// Set the node ID label.

	void  set_label(int number)             // to set the index to a number
	    {index = number;}

	/// Set the node ID label.

	void  set_label(const char * a_string)  // to set the name to a string
	    {
	    if(name != NULL)
	        delete [] name;
	    name = new char[strlen(a_string)+1];
	    strcpy(name, a_string);
	    }

	/// Set the node ID label (alternate name).

	void  set_index(int number)             // to set the index to a number
	    {index = number;}

	/// Set the node ID label (alternate name).

	void  set_name(const char * a_string)   // to set the name to a string
	{
	    if (name)
	        delete [] name;
	    if (!a_string)
		name = NULL;
	    else {
		name = new char[strlen(a_string)+1];
		strcpy(name, a_string);
	    }
	}

	/// Delete the node ID label.

	void clear_name()		    {set_name(NULL);}

	/// Delete the node ID label.

	void clear_label()		    {clear_name();}

	/// Set the node mass.

        void  set_mass(const real new_mass) {mass = new_mass;}

	/// Set the parent pointer.

	void  set_parent(node * b)          {parent = b;}

	/// Set the oldest daughter pointer.

	void  set_oldest_daughter(node * b) {oldest_daughter = b;}

	/// Set the elder sister pointer.

	void  set_elder_sister(node * b)    {elder_sister = b;}

	/// Set the younger sister pointer.

	void  set_younger_sister(node * b)  {younger_sister = b;}

#if 0
	void  set_log_story(story * s)      {log_story = s;}	
#else
	// note that some functions (e.g. node::scan_log_story)
	// simply set the story, and don't bother deleting
	// the old one, which is normally made with new_dyn()

	void  set_log_story(story * s) {
	    if (log_story != NULL) delete log_story;
	    log_story = s;
	}	
#endif
	void  set_dyn_story(story * s)      {dyn_story = s;}

	/// Add a comment to the log story.

        void  log_comment(char *);

	/// Add the argument list to the log story.

        void  log_history(int, char **);

	/// Increment the mass by the specified amount.

	void  inc_mass(const real d_mass)          {mass += d_mass;}

	/// Scale the mass by the specified factor.

	void  scale_mass(const real scale_factor)  {mass *= scale_factor;}

	int  get_index()		     const {return index;}
	char * get_name()		     const {return name;}

	inline real  get_mass()		     const {return mass;}

	inline node * get_parent()	     const {return parent;}
	inline node * get_oldest_daughter()  const {return oldest_daughter;}
	inline node * get_younger_sister()   const {return younger_sister;}
	inline node * get_elder_sister()     const {return elder_sister;}

	/// Return the top-level node of this node.

	inline node* get_top_level_node() const {

	    if (parent == NULL) return NULL;	// root node

	    node* n = const_cast<node*>(this);
	    node* g = parent->get_parent();

	    while (g) {
		n = n->get_parent();
		g = g->get_parent();
	    }
	    return n;
	}

	/// Set or find the root node pointer.

	inline void set_root(node * b = NULL) {
	    if (b)
		root = b;
	    else {
		if (parent == NULL)
		    root = this;
		else
		    root = get_top_level_node()->get_parent();
	    }
	}

	// As of 4/05, get_root() *never* sets the root node pointer.
	// We assume that set_root() has always been called first.
	// All functions that create or read data from a file should
	// call set_root() once the system has been created.  If for
	// some reason root is not set, this function reverts to the
	// original (slow) option of searching for the root node.

	/// Return (possibly find and return) the root node pointer.

	inline node* get_root() const
	{
	    if (root) return root;

	    if (parent == NULL)
		return const_cast<node*>(this);
	    else
		return get_top_level_node()->get_parent();
	}

	/// Find the binary sister of this node.

	node * get_binary_sister();

	void set_hydrobase(hydrobase * hb) {hbase = hb;}
	void set_starbase(starbase * sb)   {sbase = sb;}

	hydrobase * get_hydrobase()	const {return hbase;}
	starbase  * get_starbase()	const {return sbase;}

	story * get_log_story()		const {return log_story;}
	story * get_dyn_story()		const {return dyn_story;}
	story * get_hydro_story()	const {return hbase->get_hydro_story();}
	story * get_star_story()	const {return sbase->get_star_story();}

	/// Clear story and base pointers (dangerous!).

	virtual void null_pointers();

	/// Print the static members of this class.

	virtual void print_static(ostream &s = cerr);

	/// Read the log story from a stream.

	istream& scan_log_story(istream&, char *);

	/// Read the hydro story from a stream.

	istream& scan_hydro_story(istream&);

	/// Read the star story from a stream.

	virtual istream& scan_star_story(istream&, int level = 0);

	/// Read the dyn story from a stream.

	virtual istream& scan_dyn_story(istream&);

	/// Repair tree structure.

	virtual bool check_and_correct_node(bool verbose = false);

	/// Write the log story to a stream.

	ostream& print_log_story(ostream &s = cout);

	/// Write the hydro story to a stream.

	ostream& print_hydro_story(ostream &s = cout);

	/// Write the star story to a stream.

	ostream& print_star_story(ostream &s = cout,
				  int short_output = 0);

	/// Write the dyn story to a stream.

	virtual ostream& print_dyn_story(ostream &s = cout,
					 bool print_xreal = true,
					 int short_output = 0);

	/// Is this a solo node?

	inline bool is_isolated() const
	    {return (parent == NULL && oldest_daughter==NULL);}

	/// Is this the root node?

        inline bool is_root() const
	    {return (parent == NULL);}

	/// Is this a leaf node?

	inline bool is_leaf() const
	    {return (parent != NULL && oldest_daughter==NULL);}

	/// Is this a top-level node?

	inline bool is_top_level_node() const {
	    return (parent == get_root());
	}

	/// Is this a top-level leaf?

	inline bool is_top_level_leaf() const {
	    if (parent != get_root()) return false;
	    return (oldest_daughter == NULL);
	}

	/// Is this a low-level node?

	inline bool is_low_level_node() const {
	    return (parent != get_root());
	}

	/// Is this a low-level leaf?

	inline bool is_low_level_leaf() const {
	    if (parent == get_root()) return false;
	    return (oldest_daughter == NULL);
	}

	/// Is this a parent node?

	inline bool is_parent() const
	    {return oldest_daughter != NULL;}

	/// Is this a grandparent node?

        bool is_grandparent() const;

	/// Getthe next node in a tree-traversal loop.

        node* next_node(node*);
        node* orig_next_node(node*);

	/// Is this the name of this node?

        bool name_is(const char*) const;

	/// Return the label (index or name) of this node.

	char* format_label() const;

	/// Print the label (index or name) of this node.

	void print_label(ostream&) const;

	/// Pretty-print the label (index or name) of this node.

	void pretty_print_node(ostream& s = cerr) const;

	/// Recursively pretty-print the label of this node and thise below it.

	void pretty_print_tree(ostream& s = cerr, int level = 0);

	/// Compute the number of leaves under this node.

	int  n_leaves() const;

	/// Compute the number of daughters of this node.

	int  n_daughters() const;

	// SeBa counter access functions.
        //seba_counters* get_seba_counters() {
	//     return starbase->get_seba_counters();}
        //void set_seba_counters(seba_counters *sb) {
	//     starbase->set_seba_counters(sb);}

};

typedef  node *(*npfp)(hbpfp, sbpfp, bool);

typedef node * nodeptr;  // to enable dynamic array declarations such as
                         //    node** node_ptr_array = new nodeptr[n];
                         // (note that the following expression is illegal:
                         //    node** node_ptr_array = new (node *)[n];)

/// Create a new node, with standard starbase and hydrobase pointers.

inline  node * new_node(hbpfp the_hbpfp,
			sbpfp the_sbpfp ,
			bool use_stories)
    {return  new node(the_hbpfp, the_sbpfp, use_stories);}

/// Make a flat tree containing the specified number of nodes.

node * mk_flat_tree(int, npfp, hbpfp, sbpfp, bool use_stories = true);

/// Synonym for mk_flat_tree.

inline node * mknode(int n, hbpfp the_hbpfp = new_hydrobase,
	                    sbpfp the_sbpfp = new_starbase,
		     	    bool use_stories = true)
    {return  mk_flat_tree(n, new_node, the_hbpfp, the_sbpfp, use_stories);}

// Modified first argument (default) 5/03 (Steve):

/// Recursively read a node (or tree) from a stream.

node * get_node(istream &s = cin,
		npfp the_npfp = new_node,
		hbpfp the_hbpfp = new_hydrobase,
		sbpfp the_sbpfp = new_starbase,
		bool use_stories = true);

// Modified arguments 5/03 (Steve):
//
//  void put_node(ostream &s, node &b,
//  	      bool print_xreal = true,
//  	      int short_output = 0);
//  void put_single_node(ostream &s, node &b,
//  		     bool print_xreal = true,
//  		     int short_output = 0);

/// Recursively write a node to a stream.

void put_node(node *b,
	      ostream &s = cout,
	      bool print_xreal = true,
	      int short_output = 0);

/// Write a single node to a stream.

void put_single_node(node *b,
		     ostream &s = cout,
		     bool print_xreal = true,
		     int short_output = 0);

bool forget_node(istream &s = cin);

/// True if index i corresponds to node b or a descendent.

bool node_contains(node * b, int i);

/// True if string s corresponds to node b or a descendent.

bool node_contains(node * b, const char* s);

/// True if index i corresponds to b->get_top_level_node() or a descendent.

bool clump_contains(node * b, int i);

/// True if string s corresponds to b->get_top_level_node() or a descendent.

bool clump_contains(node * b, const char *s);

/// Recursively pretty-print a node.

void pp(const node *, ostream & s = cerr);

/// Recursively pretty-print a node, version 2.

void pp2(const node *, ostream & s = cerr, int level = 0);

#define  for_all_daughters(dyntype, mother_node, daughter_name)               \
         for (dyntype* daughter_name = mother_node->get_oldest_daughter();    \
	      daughter_name != NULL;                                          \
	      daughter_name = daughter_name->get_younger_sister())

// Note: for_all_nodes and for_all_leaves INCLUDE the base node.

#define  for_all_nodes(dyntype, base, node_name)                              \
         for (dyntype* node_name = base;                                      \
	      node_name != NULL;                                              \
	      node_name = (dyntype*) node_name->next_node(base))

#define  for_all_leaves(dyntype, base, node_name)                             \
         for (dyntype* node_name = base;                                      \
	      node_name != NULL;                                              \
	      node_name = (dyntype*) node_name->next_node(base))              \
	      if (node_name->get_oldest_daughter() == NULL)

node * mknode_mass(int n, real m = 1.0);

// Declaration of functions defined in node_tt.C

/// Compute total mass of the daughters of node n (should equal n->get_mass().

real total_mass(node *n);

/// Recursively delete node b and its descendents.

void rmtree(node *b, bool delete_b = true);

void detach_node_from_general_tree(node *n);

/// Remove node from the tree and replace it by its daughter.

void remove_node_with_one_daughter(node *n);
void detach_node_from_binary_tree(node *n);

/// Insert new node n at location of old, which becomes the only child of n.

void extend_tree(node *old, node *n);

/// Insert n into the tree as the oldest_daughter of the parent.

void add_node(node *n, node *parent);

/// Insert node n into the tree before specified node m.

void add_node_before(node *n, node *m);

/// Insert node n into a tree as the sister of node m.

void insert_node_into_binary_tree(node *n, node *m, node *new_node);

/// Is a among the offspring of b?

int is_descendent_of(node *a, node *b, int mode);

/// Return the common ancestor of a and b.

node * common_ancestor(node *a, node *b);

/// Return a pointer to the node with the specified index.

node * node_with_index(int i, node * top = NULL);

/// Return a pointer to the node with the specified name.

node * node_with_name(char* s, node * top = NULL);

/// Return the depth of node n in the tree.

int depth_of_node(node *n);

/// Construct a label for a binary from the component names.

char * construct_binary_label(node * ni, node * nj);

/// Label a binary node based on the component names (a,b).

void label_binary_node(node*);

/// Label a merger node based on the component names (a+b).

void label_merger_node(node*);

void print_normal_form(node*, ostream&);
char* get_normal_form(node*);

// From node/util:

void renumber(node* b, int istart, bool mass_order, bool name_nodes = false,
	      bool single_number = false);
void construct_node_name(node* b);

void makemass(node* b, mass_function mf,
	      real m_lower, real m_upper,
	      real exponent, real total_mass, bool renumber_stars);

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/node.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~
