
// Stripped-down and self-contained version of the Starlab hdyn
// class for use with the smallN integrator and kepler programs.
// No stories, bases, individual times or time steps...

#ifndef  HDYN_H
#  define  HDYN_H

#include "stdinc.h"
#include "vec.h"
#include "kepler.h"

class hdyn
{
  protected:

    // Static variables are initialized in hdyn.cc.

    static real system_time;

    static real eta;			// time step parameter
    static int  n_iter;			// symmetrization parameter
					// (0 ==> explicit, > 0 ==> implicit)
    static int  seed;			// random seed

    // Unperturbed motion:

    static bool allow_full_unperturbed;
    static real dt_crit;
    static real r2_crit;
    static real gamma2_unpert;		// threshold for unperturbed motion
    static real gamma_inv3;
    static int cm_index;

    // Particle properties:

    int  index;
    real mass, radius;
    real pot;			// pot for top-level node
				// tidal potential for low-level node
    real t_pred;		// pred time for root and top-level node
				// limit of kepler motion for low-level node
    bool fully_unperturbed;
    vec  pos, vel, pred_pos, pred_vel;
    vec  acc, jerk, old_acc, old_jerk;

    // Tree structure:

    hdyn *parent, *oldest_daughter, *older_sister, *younger_sister;
    kepler *kep;

  public:

    hdyn() {
	index = -1;
	mass = radius = t_pred = -1;
	pot = 0;
	fully_unperturbed = false;
	pos = vel = pred_pos = pred_vel = 0;
	acc = jerk = old_acc = old_jerk = 0;
	parent = oldest_daughter = older_sister = younger_sister = NULL;
	kep = NULL;
    }

    ~hdyn() {if (kep) delete kep;}

    // Setters:

    void set_system_time(const real t)  {system_time = t;}
    void set_eta(const real e)	        {eta = e;}
    void set_n_iter(const int n)        {n_iter = n;}
    void set_seed(const int s)	        {seed = s;}
    void set_allow_full_unperturbed(bool unp)
				  	{allow_full_unperturbed = unp;}
    void set_r2_crit(real r2)	  	{r2_crit = r2;}
    void set_dt_crit(real dt)	  	{dt_crit = dt;}
    void set_gamma(const real g)        {gamma2_unpert = g*g;
					 gamma_inv3 = pow(fmax(gamma2_unpert,
							  1.e-120), -1./6);}
    void set_cm_index(int i)		{cm_index = i;}

    void set_index(const int n)         {index = n;}
    void set_mass(const real m)         {mass = m;}
    void set_radius(const real r)       {radius = r;}
    void set_pot(const real p)          {pot = p;}
    void set_t_pred(const real t)       {t_pred = t;}
    void set_pos(const vec& p)          {pos = p;}
    void set_vel(const vec& v)          {vel = v;}
    void set_pred_pos(const vec& p)     {pred_pos = p;}
    void set_pred_vel(const vec& v)     {pred_vel = v;}
    void set_acc(const vec& a)          {acc = a;}
    void set_jerk(const vec& j)         {jerk = j;}
    void set_old_acc(const vec& a)      {old_acc = a;}
    void set_old_jerk(const vec& j)     {old_jerk = j;}
    void set_fully_unperturbed(const bool b)     {fully_unperturbed = b;}

    void set_parent(hdyn *b)            {parent = b;}
    void set_oldest_daughter(hdyn *b)   {oldest_daughter = b;}
    void set_older_sister(hdyn *b)      {older_sister = b;}
    void set_younger_sister(hdyn *b)    {younger_sister = b;}
    void set_kepler(kepler *k)	        {kep = k;}

    void inc_pos(vec dx)	        {pos += dx;}
    void inc_vel(vec dv)	        {vel += dv;}
    void inc_acc(vec da)	        {acc += da;}
    void inc_jerk(vec dj)	        {jerk += dj;}

    // Getters:

    inline real get_system_time()	const {return system_time;}
    inline real get_eta()		const {return eta;}
    inline int  get_n_iter()		const {return n_iter;}
    inline int  get_seed()		const {return seed;}
    inline real get_allow_full_unperturbed()
				  	const {return allow_full_unperturbed;}
    inline real get_r2_crit()	  	const {return r2_crit;}
    inline real get_dt_crit()	  	const {return dt_crit;}
    inline real get_gamma2_unpert()  	const {return gamma2_unpert;}
    inline real get_gamma_inv3()  	const {return gamma_inv3;}
    inline int  get_cm_index()		const {return cm_index;}

    inline int  get_index()	 	const {return index;}
    inline real get_mass()		const {return mass;}
    inline real get_radius()		const {return radius;}
    inline real get_pot()		const {return pot;}
    inline real get_t_pred()		const {return t_pred;}
    inline vec  get_pos()		const {return pos;}
    inline vec  get_vel()		const {return vel;}
    inline vec  get_pred_pos()		const {return pred_pos;}
    inline vec  get_pred_vel()		const {return pred_vel;}
    inline vec  get_acc()		const {return acc;}
    inline vec  get_jerk()		const {return jerk;}
    inline vec  get_old_acc()		const {return old_acc;}
    inline vec  get_old_jerk()		const {return old_jerk;}
    inline bool get_fully_unperturbed() const {return fully_unperturbed;}

    inline hdyn *get_parent()		const {return parent;}
    inline hdyn *get_oldest_daughter()	const {return oldest_daughter;}
    inline hdyn *get_younger_sister()	const {return younger_sister;}
    inline hdyn *get_older_sister()	const {return older_sister;}
    inline kepler *get_kepler()		const {return kep;}

    // Node accessors:

    inline hdyn* get_top_level_node() const {

	if (parent == NULL) return NULL;	// root node

	hdyn* n = (hdyn *)this;
	hdyn* g = parent->get_parent();
	while (g) {
	    n = n->get_parent();
	    g = g->get_parent();
	}
	return n;
    }

    inline hdyn *get_root() const {
      hdyn *n = (hdyn*)this;
	hdyn *p = get_parent();
	while (p) {
	    n = p;
	    p = p->get_parent();
	}
	return n;
    }

    inline hdyn * get_binary_sister() const {
	if (older_sister) return older_sister;
	if (younger_sister) return younger_sister;
	return NULL;
    }

    // Node classification:

    inline bool is_root() const {return !parent;}
    inline bool is_leaf() const {return parent && !oldest_daughter;}
    inline bool is_top_level_node() const {return parent && parent->is_root();}
    inline bool is_parent() const {return oldest_daughter;}

    // Get the next node in a tree-traversal loop.

    hdyn* next_node(hdyn*);			// in hdyn.cc

    // Prediction:

    void predict_loworder(real t) {
	if (t != t_pred) {
	    real dt = t - system_time;
	    real dt3 = dt/3;
	    real dt2 = dt/2;
	    pred_pos = ((old_jerk * dt3 + old_acc) * dt2 + vel) * dt + pos;
	    pred_vel = (old_jerk * dt2 + old_acc) * dt + vel;
	    t_pred = t;
	}
    }

    inline vec get_pred_pos(real t)   {predict_loworder(t);
				       return pred_pos;}
    inline vec get_pred_vel(real t)   {predict_loworder(t);
				       return pred_vel;}

    void  init_pred()	    {pred_pos = pos; pred_vel = vel;}
    void  store_old_force() {old_acc = acc; old_jerk = jerk;}

    // Correction:

    inline void correct_acc_and_jerk(const real new_dt,
				     const real prev_new_dt);
    inline void correct_pos_and_vel(const real new_dt);
};

typedef hdyn * hdynptr;  // to enable dynamic array declarations such as
                         //    hdyn** hdyn_ptr_array = new hdynptr[n];

#define  for_all_daughters(dyntype, mother_node, daughter_name)              \
         for (dyntype* daughter_name = mother_node->get_oldest_daughter();   \
	      daughter_name != NULL;                                         \
	      daughter_name = daughter_name->get_younger_sister())

// Note: for_all_nodes and for_all_leaves INCLUDE the base node.

#define  for_all_nodes(dyntype, base, node_name)                             \
         for (dyntype* node_name = base;                                     \
	      node_name != NULL;                                             \
	      node_name = (dyntype*) node_name->next_node(base))

#define  for_all_leaves(dyntype, base, node_name)                            \
         for (dyntype* node_name = base;                                     \
	      node_name != NULL;                                             \
	      node_name = (dyntype*) node_name->next_node(base))             \
	      if (node_name->get_oldest_daughter() == NULL)

// In hdyn.cc:

void predict_loworder_all(hdyn *b, real t);
void add_node(hdyn *n, hdyn *m);
void detach_node(hdyn *n);
void create_binary_node(hdyn *cm, hdyn *bi, hdyn *bj);
void rmtree(hdyn *b, bool delete_b = true);

// In smallN.cc:

real total_energy(hdyn *b);
void advance_components_to_time(hdyn *bi, real t);
void update_components_from_pred(hdyn *bi);
int smallN_evolve(hdyn *b,
		  real t_end = _INFINITY_,
		  real dt_log = _INFINITY_,
		  real dt_check = _INFINITY_,
		  real break_r2 = _INFINITY_,
		  int  verbose = 0);

// In smallN_unpert.cc:

real get_tidal_potential(hdyn *b, hdyn *bi, hdyn *bj, hdyn *cm,
			 bool absolute = true);
bool is_unperturbed_and_approaching(hdyn *b, real dt, hdyn *bi, hdyn *bj);
bool create_binary(hdyn *bi, hdyn *bj, bool verbose = false);
bool extend_or_end_binary(hdyn*& bi, bool verbose = false);

// In analyze.cc:

bool  check_structure(hdyn *bin, int verbose = 1);
hdyn* get_tree(hdyn *bin);

#endif
