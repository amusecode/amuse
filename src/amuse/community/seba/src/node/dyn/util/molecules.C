
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Print a hierarchical decomposition of a (small!) system.
////
//// Usage: molecules [OPTIONS] < input > output
////
//// Options:
//// None.
////
//// Written by Piet Hut.
//// 
//// Report bugs to starlab@sns.ias.edu.

//.............................................................................
//    version 1:  July 1989   Piet Hut               email: piet@iassns.bitnet
//                            Institute for Advanced Study, Princeton, NJ, USA
//    version 2:  Dec 1992   Piet Hut  --  adopted to the new C++-based Starlab
//.............................................................................
//     NOTE: this program is an adoptation of a Newton0 diagnostics piece.
//.............................................................................

#include "dyn.h"

#ifdef TOOLBOX

#define  BUFF_LENGTH    128

/*-----------------------------------------------------------------------------
 *  macros for the hierarchical decomposition of a few-body system:
 *      the following macro definitions govern the bookkeeping;
 *      they should migrate to their own header file in due time.
 *             NOTE: only one-digit star numbering implemented, i.e.
 *                   only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
#define  NBODY_MAX    10                     /* maximum number of particles  */
#define  NNODE_MAX    (2 * NBODY_MAX - 1)    /* maximum number of nodes      */

/*-----------------------------------------------------------------------------
 *  the following macros are used to access
 *                int  admin_table[NNODE_MAX][N_ADMIN];
 *  the entrees for the second variable are:
 *                                           number of daughter nodes   
 *                                           covered or visible         
 *                                           stablve, unstable, barely unbound
 *                                              or strongly unbound
 *                                           index of first daughter  
 *                                           index of second daughter 
 *                                           ...
 *                                           index of last daughter   
 *-----------------------------------------------------------------------------
 */
#define  N_ADMIN      (3 + NBODY_MAX)

#define  Ndaughter(ptr, offset)          ((ptr)[(offset)][0])

#define  Visibility(ptr, offset)         ((ptr)[(offset)][1])
#define  COVERED      1
#define  VISIBLE      2

#define  Stability(ptr, offset)          ((ptr)[(offset)][2])
#define  STABLE            1                          /* bound, and stable   */
#define  UNSTABLE          2                          /* bound, and unstable */
#define  BARELY_UNBOUND    3                          /* barely unbound      */
#define  STRONGLY_UNBOUND  4                          /* barely unbound      */
#define  UNKNOWN           5

#define  Daughter(ptr, offset, i)        ((ptr)[(offset)][3 + (i)])

/*-----------------------------------------------------------------------------
 *  the following macros are used to access
 *                real  struct_table[NNODE_MAX][N_STRUCT];
 *  the entrees of the second variable are:
 *                                          radius
 *                                          a: semimajor axis
 *                                          e: eccentricity
 *-----------------------------------------------------------------------------
 */
#define  N_STRUCT     3

#define  Radius(ptr, offset)    ((ptr)[(offset)][0])
#define  Apair(ptr, offset)     ((ptr)[(offset)][1])
#define  Epair(ptr, offset)     ((ptr)[(offset)][2])

local void  append_object(char hier_string[BUFF_LENGTH],
			  char tmp_string[BUFF_LENGTH]);
local int  count_visible_nodes(int admin_table[NNODE_MAX][N_ADMIN],
				   int nnode);
local void  dissolve_node(dyn * b, int * nnodeptr, int i,
			  int admin_table[NNODE_MAX][N_ADMIN],
			  real struct_table[NNODE_MAX][N_STRUCT],
			  real dist_table[NNODE_MAX][NNODE_MAX]);
local bool  dissolve_visible_strongly_unbound_nodes(dyn *b, int *nnodeptr,
					int admin_table[NNODE_MAX][N_ADMIN],
					real struct_table[NNODE_MAX][N_STRUCT],
				        real dist_table[NNODE_MAX][NNODE_MAX]);
local void  drop_tailnode(int * nnodeptr,
			      int admin_table[NNODE_MAX][N_ADMIN]);
local char *end_of_string(char * the_string);
local void  find_namenumbers(dyn *, char **);
local bool  find_new_node(dyn * b, int * nnodeptr,
			  int admin_table[NNODE_MAX][N_ADMIN],
			  real struct_table[NNODE_MAX][N_STRUCT],
			  real dist_table[NNODE_MAX][NNODE_MAX], int k);
local void get_binary_parameters(dyn *g1, dyn *g2, real * aptr, real * eptr);
local void  init_admin_table(int admin_table[NNODE_MAX][N_ADMIN]);
local void init_dist_table(dyn *b, real dist_table[NNODE_MAX][NNODE_MAX]);
local void  init_struct_table(real  struct_table[NNODE_MAX][N_STRUCT]);
local void  init_tuple(int tuple[NBODY_MAX],
		     int admin_table[NNODE_MAX][N_ADMIN], int n);
local void  install_com_dynamics(dyn * b, int nnode,
			       int tuple[NBODY_MAX], int k);
local void  install_node(int tuple[NBODY_MAX], dyn * b, int * nnodeptr,
		       int admin_table[NNODE_MAX][N_ADMIN],
		       real struct_table[NNODE_MAX][N_STRUCT],
		       real dist_table[NNODE_MAX][NNODE_MAX],
		       real radius, int k);
local bool is_approximately_isolated_tuple(int tuple[NBODY_MAX], int nnode,
					 int admin_table[NNODE_MAX][N_ADMIN],
					 real dist_table[NNODE_MAX][NNODE_MAX],
					 int k);
local bool  is_isolated_tuple(int tuple[NBODY_MAX], dyn * b, int nnode,
			      int admin_table[NNODE_MAX][N_ADMIN],
			      real struct_table[NNODE_MAX][N_STRUCT],
			      real * radiusptr, int k);
local bool  is_new_node(int tuple[NBODY_MAX], dyn * b, int * nnodeptr,
			int admin_table[NNODE_MAX][N_ADMIN],
			real struct_table[NNODE_MAX][N_STRUCT],
			real dist_table[NNODE_MAX][NNODE_MAX], int k);
local void  iterate_dissolve_visible_strongly_unbound_nodes(dyn * b,
							    int * nnodeptr,
					int admin_table[NNODE_MAX][N_ADMIN],
					real struct_table[NNODE_MAX][N_STRUCT],
					real dist_table[NNODE_MAX][NNODE_MAX]);
local void  iterate_find_new_node(dyn * b, int * nnodeptr,
				  int admin_table[NNODE_MAX][N_ADMIN],
				  real struct_table[NNODE_MAX][N_STRUCT],
				  real dist_table[NNODE_MAX][NNODE_MAX],int k);
local void  map_to_integers(char hier_string[BUFF_LENGTH],
			    int unordered[BUFF_LENGTH], int * nptr);
local int  next_element(int element, int tuple[NBODY_MAX], int n);
local bool  next_tuple(int tuple[NBODY_MAX], int k, int n,
		     int ordered_tuple[NBODY_MAX]);
local int  old_ranking(int old_array[BUFF_LENGTH],
		           int new_array[BUFF_LENGTH], int i, int n);
local void  print_binary_parameters(char node_report[BUFF_LENGTH],
				    real struct_table[NNODE_MAX][N_STRUCT],
				    int member);
local void  print_group_internal_structure(int nnode,
					int admin_table[NNODE_MAX][N_ADMIN],
					real struct_table[NNODE_MAX][N_STRUCT],
					char * name_table[NNODE_MAX]);
local void  print_node(char node_report[BUFF_LENGTH],
		       int admin_table[NNODE_MAX][N_ADMIN],
		       real struct_table[NNODE_MAX][N_STRUCT],
		       char * name_table[NNODE_MAX],
                       int member);
local void  print_radius(char node_report[BUFF_LENGTH],
		       real struct_table[NNODE_MAX][N_STRUCT], int member);
local void  select_object(char tmp_string[BUFF_LENGTH],
			  char hier_string[BUFF_LENGTH], int member);
local void  sort_int_array(int old_array[BUFF_LENGTH],
			   int new_array[BUFF_LENGTH], int n);
local void  swap_body_nodes(dyn *, dyn *);
local void  swap_nodes(dyn * b, int nnode, int i1, int i2,
		     int admin_table[NNODE_MAX][N_ADMIN],
		     real struct_table[NNODE_MAX][N_STRUCT],
		     real dist_table[NNODE_MAX][NNODE_MAX]);
local real  tuple_kin_energy(int tuple[NBODY_MAX], dyn * b, int k);
local real  tuple_pot_energy(int tuple[NBODY_MAX], dyn * b,
			   real dist_table[NNODE_MAX][NNODE_MAX], int k);
local real  tuple_radius(int tuple[NBODY_MAX], dyn * b,
			 real struct_table[NNODE_MAX][N_STRUCT],
			 vec & com_pos, int k);   
local void  unwrap_string(char substring[BUFF_LENGTH],
			  char * head_ptr, char * tail_ptr);
local void  update_admin_table(int tuple[NBODY_MAX], dyn * b, int nnode,
			     int admin_table[NNODE_MAX][N_ADMIN],
			     real struct_table[NNODE_MAX][N_STRUCT],
			     real dist_table[NNODE_MAX][NNODE_MAX], int k);
local void  update_dist_table(dyn * b, int nnode,
			    real dist_table[NNODE_MAX][NNODE_MAX]);
local void  update_struct_table(int nnode,
			      real struct_table[NNODE_MAX][N_STRUCT],
			      real radius, real a, real e, int k);
local void  wrap_string(char substring[BUFF_LENGTH],
			char head_char, char tail_char);

/*-----------------------------------------------------------------------------
 *  find_group_internal_structure  --  
 *-----------------------------------------------------------------------------
 */
void  find_group_internal_structure(dyn * b
//				    , bool s_flag
				    )
//bool  s_flag;                      /* if TRUE: simple, streamlined output  */
    {
//    int  pid;                           /* group id */
    int  nnode;
    int *nnodeptr;
    int  admin_table[NNODE_MAX][N_ADMIN];
    char *name_table[NNODE_MAX];
    real  struct_table[NNODE_MAX][N_STRUCT];
    real  dist_table[NNODE_MAX][NNODE_MAX];
    dyn *bi;

    for (nnode = 0, bi = b->get_oldest_daughter(); bi != NULL;
         bi = bi->get_younger_sister())
        nnode++;    

    if (nnode > NBODY_MAX)
	{
        cerr << "find_group_internal_structure: n = " << nnode
	     << " > NBODY_MAX = " << NBODY_MAX << "\n";
        exit(1);
	}

    find_namenumbers(b, name_table);

    init_dist_table(b, dist_table);
    init_admin_table(admin_table);
    init_struct_table(struct_table);

    nnodeptr = &nnode;

    iterate_find_new_node(b, nnodeptr, admin_table, struct_table, dist_table,
			  2);

    iterate_dissolve_visible_strongly_unbound_nodes(b, nnodeptr, admin_table,
						    struct_table, dist_table);

//    if (!s_flag && findiq(Pdata(Pdown(pn)), "group", &pid))
//        printf("group #%d\n", pid);

    print_group_internal_structure(nnode, admin_table, struct_table,
                                   name_table);
    for (int i = 0; i < NNODE_MAX; i++)
	if (name_table[i])
	    delete name_table[i];
    }

/*-----------------------------------------------------------------------------
 *  find_namenumbers  --  
 *                        NOTE: implemented temporarily through index numbers
 *-----------------------------------------------------------------------------
 */
local void  find_namenumbers(dyn * b, char * name_table[NNODE_MAX])
    {
    int  i;
    dyn * bi;

    for (i = 0, bi = b->get_oldest_daughter(); bi != NULL;
         i++, bi = bi->get_younger_sister())
	{
	if (bi->get_index() >= 0)
	    {
	    name_table[i] = new char[BUFF_LENGTH];
	    sprintf(name_table[i], "%d", bi->get_index());
	    }
	else
	    name_table[i] = NULL;
	}
    }

/*-----------------------------------------------------------------------------
 *  iterate_find_new_node  --  
 *-----------------------------------------------------------------------------
 */
local void  iterate_find_new_node(dyn * b, int * nnodeptr,
				  int admin_table[NNODE_MAX][N_ADMIN],
				  real struct_table[NNODE_MAX][N_STRUCT],
				  real dist_table[NNODE_MAX][NNODE_MAX], int k)
    {
    int  n_visible_nodes;

    n_visible_nodes = count_visible_nodes(admin_table, *nnodeptr);
    if (k > n_visible_nodes)
	return;

    if (find_new_node(b, nnodeptr, admin_table, struct_table, dist_table,
		      k) == TRUE)
	iterate_find_new_node(b, nnodeptr, admin_table, struct_table,
			      dist_table, 2);
    else
	iterate_find_new_node(b, nnodeptr, admin_table, struct_table,
			      dist_table, k+1);
    }

/*-----------------------------------------------------------------------------
 *  find_new_node  --  
 *-----------------------------------------------------------------------------
 */
local bool  find_new_node(dyn * b, int * nnodeptr,
			  int admin_table[NNODE_MAX][N_ADMIN],
			  real struct_table[NNODE_MAX][N_STRUCT],
			  real dist_table[NNODE_MAX][NNODE_MAX], int k)
    {
    int  i;
    int  n;                                       /* number of visible nodes */
    int  tuple[NBODY_MAX];
    int  ordered_tuple[NBODY_MAX];

    n = count_visible_nodes(admin_table, *nnodeptr);
    if (k > n)
	return(FALSE);

    init_tuple(tuple, admin_table, n);
    for (i = 0; i < n; i++)
        ordered_tuple[i] = tuple[i];

    do
        if (is_new_node(tuple, b, nnodeptr, admin_table, struct_table,
			dist_table, k) == TRUE)
	    return(TRUE);
    while
	(next_tuple(tuple, k, n, ordered_tuple) == TRUE);
	
    return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  is_new_node  --  checks wether a node is ilsolated
 *              note: old version:
 *                   if this particular tuple is not isolated, or does not form
 *                     a bound subsystem, then the value  FALSE  is returned;
 *                   if the tuple is isolated and bound,  TRUE  is returned
 *                     after the following actions have been carried out:
 *                        1) the members of this tuple are covered,
 *                        2) a new visible node is assigned to this tuple,
 *                        3) both tables are updated.
 *-----------------------------------------------------------------------------
 */
local bool  is_new_node(int tuple[NBODY_MAX], dyn * b, int * nnodeptr,
			int admin_table[NNODE_MAX][N_ADMIN],
			real struct_table[NNODE_MAX][N_STRUCT],
			real dist_table[NNODE_MAX][NNODE_MAX], int k)
    {
    real  radius;

    if (k < 2)
	{
        cerr << "is_new_node: k = " << k << " < 2\n";
        exit(1);
	}

    if (is_approximately_isolated_tuple(tuple, *nnodeptr, admin_table,
					dist_table, k) == FALSE)
	return(FALSE);

    if (is_isolated_tuple(tuple, b, *nnodeptr, admin_table, struct_table,
			  &radius, k) == FALSE)
	return(FALSE);

    install_node(tuple, b, nnodeptr, admin_table, struct_table, dist_table,
		 radius, k);

    return(TRUE);
    }

/*-----------------------------------------------------------------------------
 *  is_approximately_isolated_tuple  --  a tuple is defined to be approximately
 *                                       isolated if the minimum distance
 *                                       between a member of the tuple and a
 *                                       non-member exceeds the maximum
 *                                       distance between a pair of members.
 *-----------------------------------------------------------------------------
 */
local bool  is_approximately_isolated_tuple(int tuple[NBODY_MAX], int nnode,
					 int admin_table[NNODE_MAX][N_ADMIN],
					 real dist_table[NNODE_MAX][NNODE_MAX],
					 int k)
    {
    int  i, j, jt;
    real  max_internal_pair_distance;
    int  visible_non_tuple[NBODY_MAX];
    int  n_visible_non_tuple;

    if (k < 2)
	{
        cerr << "is_approximately_isolated_tuple: k = " << k << " < 2\n";
        exit(1);
	}

    max_internal_pair_distance = 0.0;

    for (i = 0; i < k - 1; i++)
	for (j = i + 1; j < k; j++)
	    if (dist_table[tuple[i]][tuple[j]] > max_internal_pair_distance)
		max_internal_pair_distance = dist_table[tuple[i]][tuple[j]];

    jt = 0;
    j = 0;
    for (i = 0; i < nnode; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    if (jt >= k)                            /* i.e.:  i > tuple[k-1] */
		visible_non_tuple[j++] = i;
	    else if (i < tuple[jt])
		visible_non_tuple[j++] = i;
	    else if (i == tuple[jt])
		jt++;
            else
		{
	        cerr << 
		      "is_approximately_isolated_tuple: non-visible tuple?\n";
	        exit(1);
		}

    n_visible_non_tuple = j;
    
    for (i = 0; i < n_visible_non_tuple; i++)
	for (j = 0; j < k; j++)
	    if (dist_table[visible_non_tuple[i]][tuple[j]]
		< max_internal_pair_distance)
		return(FALSE);

    return(TRUE);    
    }

/*-----------------------------------------------------------------------------
 *  is_isolated_tuple  -- two objects, with radius Ri and Rj, are considered
 *                        isolated, if each point inside a sphere belonging
 *                        to an object is closer to any point in its own sphere
 *                        than to any point in the other sphere. This implies
 *                        that the distance between the sphere surfaces should
 *                        be at least as large as the largest of the two 
 *                        diameters, i.e. 2.0*max(Ri, Rj), and the distance
 *                        between the two centers of mass of the objects should
 *                        be at least  
 *                                   Ri + Rj + 2.0*max(Ri, Rj)
 *-----------------------------------------------------------------------------
 */
local bool  is_isolated_tuple(int tuple[NBODY_MAX], dyn * b, int nnode,
			      int admin_table[NNODE_MAX][N_ADMIN],
			      real struct_table[NNODE_MAX][N_STRUCT],
			      real * radiusptr, int k)
    {
    int  i, j, jt;
    real  radius;
    real  radius_i;
    vec  com_pos;
    real  distance_to_com;
    int  visible_non_tuple[NBODY_MAX];
    int  n_visible_non_tuple;
    int  pi_offset;
    dyn * bi;

    radius = tuple_radius(tuple, b, struct_table, com_pos, k);
    *radiusptr = radius;

    jt = 0;
    j = 0;

    for (i = 0; i < nnode; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    if (jt >= k)                            /* i.e.:  i > tuple[k-1] */
                visible_non_tuple[j++] = i;
	    else if (i < tuple[jt])
		visible_non_tuple[j++] = i;
	    else if (i == tuple[jt])
		jt++;
            else
		{
	        cerr << "is_isolated_tuple: non-visible tuple?\n";
	        exit(1);
		}

    n_visible_non_tuple = j;
    for (i = 0; i < n_visible_non_tuple; i++)
	{
	bi = b->get_oldest_daughter();
	pi_offset = visible_non_tuple[i];
	while (pi_offset--)
	    bi = bi->get_younger_sister();
	distance_to_com = abs(bi->get_pos() - com_pos);
	radius_i = Radius(struct_table, visible_non_tuple[i]);
	if (distance_to_com < radius + radius_i + 2.0 * Starlab::max(radius, radius_i))
	    return(FALSE);
	}

    return(TRUE);
    }

/*-----------------------------------------------------------------------------
 *  how_bound_tuple  --  
 *-----------------------------------------------------------------------------
 */
local void  how_bound_tuple(int tuple[NBODY_MAX], dyn * b,
			    real dist_table[NNODE_MAX][NNODE_MAX],
			    int k, bool * isbound, bool * isbarelyunbound)
    {
    real  kinetic_energy;
    real  abs_potential_energy;

    kinetic_energy = tuple_kin_energy(tuple, b, k);
    abs_potential_energy = -tuple_pot_energy(tuple, b, dist_table, k);

    *isbound = *isbarelyunbound = FALSE;

    if (kinetic_energy < abs_potential_energy)
	*isbound = TRUE;
    else if (kinetic_energy < 1.5 * abs_potential_energy)
	*isbarelyunbound = TRUE;
    }

/*-----------------------------------------------------------------------------
 *  tuple_pot_energy  --  return the total potential energy of all the nodes in
 *                        the tuple, due to their mutual gravitational
 *                        attraction.
 *-----------------------------------------------------------------------------
 */
local real  tuple_pot_energy(int tuple[NBODY_MAX], dyn * b,
			     real dist_table[NNODE_MAX][NNODE_MAX], int k)
    {
    int  i, j;
    real  pot_energy;
    int  pi_offset;
    int  pj_offset;
    dyn * bi;
    dyn * bj;

    pot_energy = 0.0;

    for (i = 0; i < k - 1; i++)
        {
	bi = b->get_oldest_daughter();
	pi_offset = tuple[i];
	while (pi_offset--)
	    bi = bi->get_younger_sister();

	for (j = i + 1; j < k; j++)
	    {
	    bj = b->get_oldest_daughter();
	    pj_offset = tuple[j];
	    while (pj_offset--)
	        bj = bj->get_younger_sister();

	    pot_energy -= ( bi->get_mass() * bj->get_mass() )
		          / dist_table[ tuple[i] ][ tuple[j] ];
	    }
        }

    return(pot_energy);
    }

/*-----------------------------------------------------------------------------
 *  tuple_kin_energy  --  return the total kinetic energy of all the nodes in
 *                        the tuple, in the center of mass frame of the tuple.
 *-----------------------------------------------------------------------------
 */
local real  tuple_kin_energy(int tuple[NBODY_MAX], dyn * b, int k)
    {
    int  i;
    real  velocity_squared;
    real  kin_energy;
    real  total_mass;
    vec  mass_times_vel;
    vec  rel_vel;                  /* velocity with respect to c.o.m. */
    vec  com_vel;
    int  pi_offset;
    dyn * bi;

    total_mass = 0;

    com_vel = 0;
    for (i = 0; i < k; i++)
	{
	bi = b->get_oldest_daughter();
	pi_offset = tuple[i];
	while (pi_offset--)
	    bi = bi->get_younger_sister();

	total_mass += bi->get_mass();
	mass_times_vel = bi->get_mass() * bi->get_vel();
	com_vel += mass_times_vel;
	}
    com_vel /= total_mass;                       /* the real c.o.m. velocity */

    kin_energy = 0;

    for (i = 0; i < k; i++)
	{
	bi = b->get_oldest_daughter();
	pi_offset = tuple[i];
	while (pi_offset--)
	    bi = bi->get_younger_sister();

	rel_vel = bi->get_vel() - com_vel;
	velocity_squared = rel_vel * rel_vel;
	kin_energy += bi->get_mass() * velocity_squared;     /* 2 x too high */
	}

    kin_energy *= 0.5;                 /* corrects for the factor 2 left out */

    return(kin_energy);
    }

#define  LARGE_INTEGER   100000
/*-----------------------------------------------------------------------------
 *  representative_body  --  returns the lowest of the numbers of the atoms
 *                           which make up the node  i .
 *-----------------------------------------------------------------------------
 */
local int  representative_body(int admin_table[NNODE_MAX][N_ADMIN], int i)
    {
    int  j;
    int  number_of_daughters;
    int  daughter_j_representative;
    int  lowest_number;

    number_of_daughters = Ndaughter(admin_table, i);

    if (number_of_daughters != 0)
	{
	lowest_number = LARGE_INTEGER;
	for (j = 0; j < number_of_daughters; j++)
	    {
	    daughter_j_representative = representative_body(admin_table,
						    Daughter(admin_table,i,j));
	    if (lowest_number > daughter_j_representative)
		lowest_number = daughter_j_representative;
	    }
	return(lowest_number);
	}
    else
	return(i);
    }

/*-----------------------------------------------------------------------------
 *  init_dist_table  --  computes all node-node distances between their
 *                       centers of mass, and installs these values in the
 *                       distance table.
 *                       note: a moderate amount of inefficiency is allowed,
 *                             for the sake of simplicity.
 *                                 For example, the table is  n x n ,  rather 
 *                                 than n x (n-1) / 2  [symmetry is not used];
 *                                 and the distances, rather than the cheaper
 *                                 squared distances, are computed.
 *-----------------------------------------------------------------------------
 */
local void  init_dist_table(dyn * b, real dist_table[NNODE_MAX][NNODE_MAX])
    {
    int  i, j;
    dyn *bi, *bj;
    real  distance;

    for (i = 0, bi = b->get_oldest_daughter(); bi != NULL;
         i++, bi = bi->get_younger_sister())
	dist_table[i][i] = 0.0;

    for (i = 0, bi = b->get_oldest_daughter(); bi != NULL;
                 i++, bi = bi->get_younger_sister())
	for (j = i+1, bj = bi->get_younger_sister(); bj != NULL;
                     j++, bj = bj->get_younger_sister())
	    {
	    distance = abs(bi->get_pos() - bj->get_pos());
	    dist_table[i][j] = dist_table[j][i] = distance;
	    }
    }

/*-----------------------------------------------------------------------------
 *  init_admin_table  --  
 *-----------------------------------------------------------------------------
 */
local void  init_admin_table(int admin_table[NNODE_MAX][N_ADMIN])
    {
    int  i;

    for (i = 0; i < NNODE_MAX; i++)      /* initially the following holds:   */
	{
	Ndaughter(admin_table, i) = 0;         /* no daughter nodes present  */
	Visibility(admin_table, i) = VISIBLE;  /* and all nodes are visible  */
	Stability(admin_table, i) = STABLE;    /* and stable point particles */
	}
    }

/*-----------------------------------------------------------------------------
 *  init_struct_table  --  
 *-----------------------------------------------------------------------------
 */
local void  init_struct_table(real  struct_table[NNODE_MAX][N_STRUCT])
    {
    int  i;

    for (i = 0; i < NNODE_MAX; i++)
	Radius(struct_table, i) = 0.0;
    }

/*-----------------------------------------------------------------------------
 *  count_visible_nodes  --  
 *-----------------------------------------------------------------------------
 */
local int  count_visible_nodes(int admin_table[NNODE_MAX][N_ADMIN], int nnode)
    {
    int  i;
    int  n_visible_nodes;

    n_visible_nodes = 0;

    for (i = 0; i < nnode; i++)
	if (Visibility(admin_table, i) == VISIBLE)
	    n_visible_nodes++;

    return( n_visible_nodes );
    }

/*-----------------------------------------------------------------------------
 *  init_tuple  -  find the indices of the first  n  visible nodes, and list
 *                 these indices at the beginning of the  tuple[]  array.
 *-----------------------------------------------------------------------------
 */
local void  init_tuple(int tuple[NBODY_MAX],
		       int admin_table[NNODE_MAX][N_ADMIN], int n)
    {
    int  i, j;

    if (n > NBODY_MAX)
	{
        cerr << "init_tuple: n = " << n << " > NBODY_MAX = "
	     << NBODY_MAX << "\n";
        exit(1);
	}

    i = j = 0;

    while(i < NNODE_MAX && j < n)
	{
	if (Visibility(admin_table, i) == VISIBLE)
	    tuple[j++] = i;
	i++;
	}
    if (j != n)
	{
        cerr << "init_tuple: j = " << j << " != n = " << n << "\n";
        exit(1);
	}
    }

/*-----------------------------------------------------------------------------
 *  next_tuple  -  search for the next k-tuple from among the ordered list of
 *                 indices of visible nodes (provided in  ordered_tuple[] );
 *                 if there is a next k-tuple,
 *                     then replace the first k entrees in tuple[] with the
 *                     new k-tuple, and return TRUE;
 *                 else
 *                     return FALSE.
 *           note: the contents of tuple[] beyond the first k entrees, i.e.
 *                     tuple[k], tuple[k+1], ..., tuple[NBODY_MAX -1]
 *                 is undetermined, and should not be accessed.
 *-----------------------------------------------------------------------------
 */
local bool  next_tuple(int tuple[NBODY_MAX], int k, int n,
		       int ordered_tuple[NBODY_MAX])
    {
    if (k > n)
	{
        cerr << "next_tuple: k = " << k << " > n = " << n << "\n";
        exit(1);
	}

    if (k < 1)
	return(FALSE);

    if (tuple[k-1] < ordered_tuple[n-1])
	{
	tuple[k-1] = next_element(tuple[k-1], ordered_tuple, n);
	return(TRUE);
	}
    else if (next_tuple(tuple, k-1, n-1, ordered_tuple) == TRUE)
	{
        tuple[k-1] = next_element(tuple[k-2], ordered_tuple, n);
	return(TRUE);
	}
    else
	return(FALSE);
    }

/*-----------------------------------------------------------------------------
 *  next_element  -  return the next element from an array tuple[]
 *-----------------------------------------------------------------------------
 */
local int  next_element(int element, int tuple[NBODY_MAX], int n)
    {
    int  i;

    for (i = 0; i < n - 1; i++)
	if (element == tuple[i])
	    return(tuple[i+1]);

    int shhh = 1;           // to make the SG compiler shut up

    if (element == tuple[n-1])
	{
        cerr << 
	      "next_element: element provided is the last element of tuple\n";
        exit(1);
	}
    else if (shhh)           // to make the SG compiler shut up
	{
        cerr << "next_element: element provided is not an element of tuple\n";
        exit(1);
	}
    else                    // to make the SG compiler shut up
        return 0;           // to make the g++ compiler happy

    return 0;               // to make HP g++ happy
    }

/*-----------------------------------------------------------------------------
 *  tuple_radius  --  defined as the largest value of
 *                    (distance from the center of mass of tuple to a member
 *                     + internal radius of that member);
 *                    as a side effect, the center-of-mass positon is returned
 *                    in the argument com_pos[].
 *           new note:
 *                    for a k_tuple with  k = 2 but unbound, same treatment as
 *                    for k > 2.
 *           old note:
 *                    for a k_tuple with  k > 2  the instantaneous radius
 *                    is computed; 
 *                    for a pair (k = 2), the maximum radius is computed in the
 *                    approximation of an unperturbed ellipse,
 *                    i.e. the maximum value of
 *                    (distance from the center of mass of the binary to the
 *                    apocenter of a member + internal radius of that member).
 *-----------------------------------------------------------------------------
 */
local real  tuple_radius(int tuple[NBODY_MAX], dyn * b,
			 real struct_table[NNODE_MAX][N_STRUCT],
			 vec & com_pos, int k)
    {
    int  i;
    real  member_radius;
    real  max_member_radius;
    real  total_mass;
    real  lever_arm_factor;          /* in case of a binary, i.e. k = 2      */
    vec  mass_times_pos;
    real  a, e;                      /* binary parameters for the case k = 2 */
    int  member_offset;
    int  g0_offset;
    int  g1_offset;
    dyn * member;
    dyn * g0;
    dyn * g1;

    total_mass = 0.0;
    com_pos = 0;
    for (i = 0; i < k; i++)
	{
	member = b->get_oldest_daughter();
	member_offset = tuple[i];
	while (member_offset--)
	    member = member->get_younger_sister();

	total_mass += member->get_mass();
	mass_times_pos = member->get_mass() * member->get_pos();
	com_pos += mass_times_pos;
	}
    com_pos /= total_mass;                       /* the real c.o.m. position */

    if (k == 2)
        {
	g0 = g1 = b->get_oldest_daughter();
	g0_offset = tuple[0];
	g1_offset = tuple[1];
	while (g0_offset--)
	    g0 = g0->get_younger_sister();
	while (g1_offset--)
	    g1 = g1->get_younger_sister();

	get_binary_parameters(g0, g1, &a, &e);
        }

    max_member_radius = 0.0;

    for (i = 0; i < k; i++)
	{
	if (k == 2 && a > 0.0)
	    {
	    if (i == 0)
		lever_arm_factor = g1->get_mass() / total_mass;
	    else
		lever_arm_factor = g0->get_mass() / total_mass;
	    member_radius = a * (1.0 + e) * lever_arm_factor;
	    }
	else
	    {
	    member = b->get_oldest_daughter();
	    member_offset = tuple[i];
	    while (member_offset--)
	        member = member->get_younger_sister();

	    member_radius = abs(member->get_pos() - com_pos);
	    }
	member_radius += Radius(struct_table, tuple[i]);
	if (member_radius > max_member_radius)
	    max_member_radius = member_radius;
	}

    return(max_member_radius);
    }

/*-----------------------------------------------------------------------------
 *  install_node  --  does all the bookkeeping needed for the introduction of
 *                    a new node.
 *                    note: incrementing the node counter should be postponed
 *                          to the last line, since all functions called
 *                          presume that the new node is  nodes[*nnodeptr] .
 *                     note:
 *                          update_struct_table()  has to be invoked before
 *                          invoking  update_admin_table() , so that the
 *                          physical parameters of the new node are availalbe
 *                          to determine stability of the new node.
 *-----------------------------------------------------------------------------
 */
local void  install_node(int tuple[NBODY_MAX], dyn * b, int * nnodeptr,
			 int admin_table[NNODE_MAX][N_ADMIN],
			 real struct_table[NNODE_MAX][N_STRUCT],
			 real dist_table[NNODE_MAX][NNODE_MAX],
			 real radius, int k)
    {
    real  a;                                   /* semimajor axis of a binary */
    real  e;                                   /* eccentricity of a binary   */
    int  g0_offset;
    int  g1_offset;
    dyn * g0;
    dyn * g1;

    if (*nnodeptr >= NNODE_MAX)
	{
        cerr << "install_node: *nnodeptr = " << *nnodeptr
	     << " >= NNODE_MAX = " << NNODE_MAX << "\n";
        exit(1);
	}

    if (k == 2)
        { 
	g0 = g1 = b->get_oldest_daughter();
	g0_offset = tuple[0];
	g1_offset = tuple[1];
	while (g0_offset--)
	    g0 = g0->get_younger_sister();
	while (g1_offset--)
	    g1 = g1->get_younger_sister();

	get_binary_parameters(g0, g1, &a, &e);
        }

    install_com_dynamics(b, *nnodeptr, tuple, k);

    update_dist_table(b, *nnodeptr, dist_table);
    update_struct_table(*nnodeptr, struct_table, radius, a, e, k);
    update_admin_table(tuple, b, *nnodeptr, admin_table, struct_table,
		       dist_table, k);
    *nnodeptr += 1;
    }

/*-----------------------------------------------------------------------------
 *  update_dist_table  --  
 *-----------------------------------------------------------------------------
 */
local void  update_dist_table(dyn * b, int nnode,
			      real dist_table[NNODE_MAX][NNODE_MAX])
    {
    int  i;
    real  distance;
    dyn * bi;
    dyn * glast;

    dist_table[nnode][nnode] = 0.0;

    glast = b->get_oldest_daughter();
    for (i = 0; i < nnode; i++)
        glast = glast->get_younger_sister();

    for (i = 0, bi = b->get_oldest_daughter(); i < nnode;
         i++, bi = bi->get_younger_sister())
	{
	distance = abs(bi->get_pos() - glast->get_pos());
	dist_table[i][nnode] = dist_table[nnode][i] = distance;
	}
    }

#define  STABLE_SEPARATION_FACTOR     2.0              /* somewhat arbitrary */

/*-----------------------------------------------------------------------------
 *  update_admin_table  --  determines the values for the visibility and
 *                          stability of a new node, and enters those values
 *                          in the  admin_table[] .
 *                          Check the node for being barely or strongly 
 *                          unbound; if not, then bound, so check for
 *                          stable or unstable.
 *                     note:
 *                          the use of how_bound_tuple() may well be overkill,
 *                          but as long as it works, no problem for now.
 *    note: oldversion:
 *                     note:
 *                          a k-tuple with  k > 2  is considered to be always
 *                          unstable; a binary (k = 2) is stable if both
 *                          members are internally stable and in addition
 *                          the pericenter exceeds the sum of the radii of the
 *                          members by a factor  STABLE_SEPARATION_FACTOR .
 *-----------------------------------------------------------------------------
 */
local void  update_admin_table(int tuple[NBODY_MAX], dyn * b, int nnode,
			       int admin_table[NNODE_MAX][N_ADMIN],
			       real struct_table[NNODE_MAX][N_STRUCT],
			       real dist_table[NNODE_MAX][NNODE_MAX], int k)
    {
    int  i;
    bool  isbound;
    bool  isbarelyunbound;
    real  pericenter_distance;
    real  sum_of_radii;

    Ndaughter(admin_table, nnode) = k;
    Visibility(admin_table, nnode) = VISIBLE;

    for (i = 0; i < k; i++)
	{
        Daughter(admin_table, nnode, i) = tuple[i];
	Visibility(admin_table, tuple[i]) = COVERED;
	}

    how_bound_tuple(tuple, b, dist_table, k, &isbound, &isbarelyunbound);

    if (isbarelyunbound == TRUE)
	Stability(admin_table, nnode) = BARELY_UNBOUND;
    else if (isbound == FALSE)
	Stability(admin_table, nnode) = STRONGLY_UNBOUND;
    else if (k > 2)
	Stability(admin_table, nnode) = UNSTABLE;
    else
	{
	Stability(admin_table, nnode) = STABLE;
	for (i = 0; i < k; i++)
	    if (Ndaughter(admin_table, Daughter(admin_table, nnode, i)) > 0)
	        if (Stability(admin_table, Daughter(admin_table, nnode, i))
		    != STABLE)
		    Stability(admin_table, nnode) = UNSTABLE;

	pericenter_distance =
	    Apair(struct_table, nnode) * (1.0 - Epair(struct_table, nnode));
	sum_of_radii = Radius(struct_table, Daughter(admin_table, nnode, 0))
	             + Radius(struct_table, Daughter(admin_table, nnode, 1));

	if (pericenter_distance < STABLE_SEPARATION_FACTOR * sum_of_radii)
	    Stability(admin_table, nnode) = UNSTABLE;
	}
    }

/*-----------------------------------------------------------------------------
 *  update_struct_table  --  
 *-----------------------------------------------------------------------------
 */
local void  update_struct_table(int nnode,
				real struct_table[NNODE_MAX][N_STRUCT],
				real radius, real a, real e, int k)
    {
    Radius(struct_table, nnode) = radius;

    if (k == 2)
	{
	Apair(struct_table, nnode) = a;
	Epair(struct_table, nnode) = e;
	}
    }

/*-----------------------------------------------------------------------------
 *  get_binary_parameters  --  
 *-----------------------------------------------------------------------------
 */
local void get_binary_parameters(dyn *g1, dyn *g2, real * aptr, real * eptr)
//realptr  aptr;             /* pointer to a, the semimajor axis of a binary */
//realptr  eptr;             /* pointer to e, the eccentricity of a binary   */
    {
    real  delta_r, delta_v, m_sum;
    vec  r_rel, v_rel, r_out_v;
    real  r_out_v_squared;

    delta_r = abs(g1->get_pos() - g2->get_pos());
    delta_v = abs(g1->get_vel() - g2->get_vel());
    m_sum = g1->get_mass() + g2->get_mass();

    *aptr = 1.0 / (2.0/delta_r - delta_v*delta_v / m_sum);

    r_rel = g1->get_pos() - g2->get_pos();
    v_rel = g1->get_vel() - g2->get_vel();
    r_out_v = r_rel ^ v_rel;
    r_out_v_squared = r_out_v * r_out_v;

    *eptr = sqrt(Starlab::max(0.0, 1.0 - (r_out_v_squared / (m_sum * *aptr))));
    }

/*-----------------------------------------------------------------------------
 *  iterate_dissolve_visible_strongly_unbound_nodes  --  
 *-----------------------------------------------------------------------------
 */
local void  iterate_dissolve_visible_strongly_unbound_nodes(dyn * b,
							    int * nnodeptr,
					int admin_table[NNODE_MAX][N_ADMIN],
					real struct_table[NNODE_MAX][N_STRUCT],
					real dist_table[NNODE_MAX][NNODE_MAX])
    {
    while (dissolve_visible_strongly_unbound_nodes(b, nnodeptr,admin_table,
						   struct_table, dist_table))
	;
    }	  

/*-----------------------------------------------------------------------------
 *  dissolve_visible_strongly_unbound_nodes  --  
 *-----------------------------------------------------------------------------
 */
local bool  dissolve_visible_strongly_unbound_nodes(dyn * b, int * nnodeptr,
					int admin_table[NNODE_MAX][N_ADMIN],
					real struct_table[NNODE_MAX][N_STRUCT],
					real dist_table[NNODE_MAX][NNODE_MAX])
    {
    int  i;
    bool  new_activity;

    new_activity = FALSE;
    for (i = 0; i < *nnodeptr; i++)
	if (Visibility(admin_table, i) == VISIBLE &&
	    Stability(admin_table, i) == STRONGLY_UNBOUND)
	    {
	    dissolve_node(b, nnodeptr, i, admin_table, struct_table,
			  dist_table);
	    new_activity = TRUE;
	    break;                        /* only one dissolution at a time, */
	    }                             /* to retain integrity of nnodeptr */

    return(new_activity);
    }

/*-----------------------------------------------------------------------------
 *  dissolve_node  --  
 *-----------------------------------------------------------------------------
 */
local void  dissolve_node(dyn * b, int * nnodeptr, int i,
			  int admin_table[NNODE_MAX][N_ADMIN],
			  real struct_table[NNODE_MAX][N_STRUCT],
			  real dist_table[NNODE_MAX][NNODE_MAX])
    {
    if (i != *nnodeptr-1)
	swap_nodes(b, *nnodeptr, i, *nnodeptr-1, admin_table, struct_table,
		   dist_table);

    drop_tailnode(nnodeptr, admin_table);
    }

/*-----------------------------------------------------------------------------
 *  swap_nodes  --  
 *-----------------------------------------------------------------------------
 */
local void  swap_nodes(dyn * b, int nnode, int i1, int i2,
		       int admin_table[NNODE_MAX][N_ADMIN],
		       real struct_table[NNODE_MAX][N_STRUCT],
		       real dist_table[NNODE_MAX][NNODE_MAX])
    {
    int  tmp_int;
    real  tmp_real;
    int  j;
    int  i;
    dyn * bi1;
    dyn * bi2;

    if (i1 >= nnode)
	{
        cerr << "swap_nodes: i1 = " << i1 << " > nnode = " << nnode << "\n";
        exit(1);
	}
    if (i2 >= nnode)
	{
        cerr << "swap_nodes: i2 = " << i2 << " > nnode = " << nnode << "\n";
        exit(1);
	}

    bi1 = bi2 = b->get_oldest_daughter();
    while (i1--)
        bi1 = bi1->get_younger_sister();
    while (i2--)
        bi2 = bi2->get_younger_sister();

    swap_body_nodes(bi1, bi2);

    tmp_int = Ndaughter(admin_table, i1);
    Ndaughter(admin_table, i1) = Ndaughter(admin_table, i2);
    Ndaughter(admin_table, i2) = tmp_int;

    tmp_int = Visibility(admin_table, i1);
    Visibility(admin_table, i1) = Visibility(admin_table, i2);
    Visibility(admin_table, i2) = tmp_int;

    tmp_int = Stability(admin_table, i1);
    Stability(admin_table, i1) = Stability(admin_table, i2);
    Stability(admin_table, i2) = tmp_int;
    
    j = Starlab::max(Ndaughter(admin_table, i1), Ndaughter(admin_table, i2));
    while (j-- > 0)           /* this swaps also unused numbers, if the two  */
	{                     /* Ndaughter values are unequal, but who cares */
	tmp_int = Daughter(admin_table, i1, j);
	Daughter(admin_table, i1, j) = Daughter(admin_table, i2, j);
	Daughter(admin_table, i2, j) = tmp_int;
	}
/*
 * Now switch also the pointers to i1 and i2 from the outside,
 * including the cases where i = i1 and i = i2 !
 */
    for (i = 0; i < nnode; i++)
	if (Ndaughter(admin_table, i) > 0)
	    for (j = 0; j < Ndaughter(admin_table, i); j++)
		{
		if (Daughter(admin_table, i, j) == i1)
		    Daughter(admin_table, i, j) = i2;
		else if (Daughter(admin_table, i, j) == i2)
		    Daughter(admin_table, i, j) = i1;
		}
//
// NOTE: 921219: the body of the above loop used to be as follows
//       (obviously different from what originally was intended).
//       The question is now whether the above repair will introduce
//       new bugs or not! 
//
//		{
//		if (Daughter(admin_table, i, j) == i1)
//		    Daughter(admin_table, i, j) == i2;
//		else if (Daughter(admin_table, i, j) == i2)
//		    Daughter(admin_table, i, j) == i1;
//		}
//
    tmp_real = Radius(struct_table, i1);
    Radius(struct_table, i1) = Radius(struct_table, i2);
    Radius(struct_table, i2) = tmp_real;
    
    tmp_real = Apair(struct_table, i1);
    Apair(struct_table, i1) = Apair(struct_table, i2);
    Apair(struct_table, i2) = tmp_real;
    
    tmp_real = Epair(struct_table, i1);
    Epair(struct_table, i1) = Epair(struct_table, i2);
    Epair(struct_table, i2) = tmp_real;
    
    for (i = 0; i < nnode; i++)
	{
	tmp_real = dist_table[i][i1];
	dist_table[i][i1] = dist_table[i][i2];
	dist_table[i][i2] = tmp_real;
	}

    for (j = 0; j < nnode; j++)
	{
	tmp_real = dist_table[i1][j];
	dist_table[i1][j] = dist_table[i2][j];
	dist_table[i2][j] = tmp_real;
	}
    }

/*-----------------------------------------------------------------------------
 *  swap_body_nodes  --  swap two sister nodes
 *-----------------------------------------------------------------------------
 */
local void  swap_body_nodes(dyn * bi, dyn * bj)
    {
    dyn *bir, *bil, *bjr, *bjl, *bb;

    if ((bb = bi->get_parent()) == NULL)
	{
	cerr <<  "swap_body_nodes: node #" << bi << " has no parent\n";
	exit(1);
	}
    
    if (bi->get_parent() != bj->get_parent())
	{
	cerr << "nodes #" << bi << " and #" << bj
	     << " have different parents\n";
	exit(1);
	}

    if (bi == bj)
        return;

    bir = bi->get_younger_sister();
    bil = bi->get_elder_sister();
    bjr = bj->get_younger_sister();
    bjl = bj->get_elder_sister();

    if (bil == bj)
	{
	bj->set_elder_sister(bi);
	bi->set_younger_sister(bj);

	bi->set_elder_sister(bjl);
	if (bjl)
	    bjl->set_younger_sister(bi);
	else
	    bb->set_oldest_daughter(bi);

	bj->set_younger_sister(bir);
	if (bir)
	    bir->set_elder_sister(bj);
	}
    else if (bjl == bi)
	{
	bi->set_elder_sister(bj);
	bj->set_younger_sister(bi);

	bj->set_elder_sister(bil);
	if (bil)
	    bil->set_younger_sister(bj);
	else
	    bb->set_oldest_daughter(bj);

	bi->set_younger_sister(bjr);
	if (bjr)
	    bjr->set_elder_sister(bi);
	}
    else
        {
	bi->set_elder_sister(bjl);
	bi->set_younger_sister(bjr);
	bj->set_elder_sister(bil);
	bj->set_younger_sister(bir);

	if (bil)
	    bil->set_younger_sister(bj);
	else
	    bb->set_oldest_daughter(bj);

	if (bjl)
	    bjl->set_younger_sister(bi);
	else
	    bb->set_oldest_daughter(bi);
 
	if (bir)
	    bir->set_elder_sister(bj);

	if (bjr)
	    bjr->set_elder_sister(bi);
        }
    }

/*-----------------------------------------------------------------------------
 *  drop_tailnode  --  
 *-----------------------------------------------------------------------------
 */
local void  drop_tailnode(int * nnodeptr, int admin_table[NNODE_MAX][N_ADMIN])
    {
    int  j;

    if (Visibility(admin_table, *nnodeptr-1) != VISIBLE)
	{
        cerr << "drop_tailnode: can't drop an invisible tail\n";
        exit(1);
	}

    for (j = 0; j < Ndaughter(admin_table, *nnodeptr-1); j++)
	Visibility(admin_table, Daughter(admin_table, *nnodeptr-1, j))
	    = VISIBLE;

    (*nnodeptr)--;
    }

/*-----------------------------------------------------------------------------
 *  tail_insert  --  
 *-----------------------------------------------------------------------------
 */
local void tail_insert(dyn * b, dyn * new_node)
    {
    dyn *bi, *b_prev;

    b_prev = b->get_oldest_daughter();

    if (b_prev == NULL)
	{
	cerr << "tail_insert: b_prev == NULL\n";
	exit(1);
	}

    while (bi = b_prev->get_younger_sister())
	b_prev = bi;
	
    b_prev->set_younger_sister(new_node);
    }

/*-----------------------------------------------------------------------------
 *  install_com_dynamics  --  compute and install the values of the  mass,
 *                            position and velocity for a new node, with the
 *                            index  nnode , and having  k  daughters,
 *                            whose indices are contained in the first  k
 *                            entrees of  tuple[] .
 *-----------------------------------------------------------------------------
 */
local void  install_com_dynamics(dyn * b, int nnode,
				 int tuple[NBODY_MAX], int k)
    {
    int  i;
    real  total_mass;
    vec  com_pos;
    vec  com_vel;
    vec  mass_times_pos;
    vec  mass_times_vel;
    int  member_offset;
    dyn * member;
    dyn * new_node;

// just to keep the compiler happy:
    nnode++; nnode--;  

    total_mass = 0.0;
    com_pos = com_vel = 0;
    for (i = 0; i < k; i++)
	{
	member = b->get_oldest_daughter();
	member_offset = tuple[i];
	while (member_offset--)
	    member = member->get_younger_sister();

	total_mass += member->get_mass();
	mass_times_pos = member->get_mass() * member->get_pos();
	com_pos += mass_times_pos;
	mass_times_vel = member->get_mass() * member->get_vel();
	com_vel += mass_times_vel;
	}
    com_pos /= total_mass;                       /* the real c.o.m. position */
    com_vel /= total_mass;                       /* the real c.o.m. position */

    new_node = new dyn;
    new_node->set_mass(total_mass);
    new_node->set_pos(com_pos);
    new_node->set_vel(com_vel);

    tail_insert(b, new_node);
    }

/*-----------------------------------------------------------------------------
 *  print_group_internal_structure  --  
 *-----------------------------------------------------------------------------
 */
local void  print_group_internal_structure(int nnode,
					int admin_table[NNODE_MAX][N_ADMIN],
					real struct_table[NNODE_MAX][N_STRUCT],
					char * name_table[NNODE_MAX])
    {
    int  i;
    char  node_report[BUFF_LENGTH];

    for (i = 0; i < nnode; i++)
	if (Ndaughter(admin_table, i) > 0)
	    {
	    sprintf(node_report, "  ", i);

	    print_node(node_report, admin_table, struct_table, name_table, i);

	    if ((Ndaughter(admin_table, i) == 2))
		print_binary_parameters(node_report, struct_table, i);
	    else
		print_radius(node_report, struct_table, i);

	    printf(node_report);
	    }
    }

/*-----------------------------------------------------------------------------
 *  print_node  --  
 *-----------------------------------------------------------------------------
 */
local void  print_node(char node_report[BUFF_LENGTH],
		       int admin_table[NNODE_MAX][N_ADMIN],
		       real struct_table[NNODE_MAX][N_STRUCT],
		       char * name_table[NNODE_MAX],
                       int member)
    {
    int  i;

    if (Ndaughter(admin_table, member) == 0)
	{
	if (name_table[member])
            {
	    sprintf(end_of_string(node_report), "%s", name_table[member]);
//            *(end_of_string(node_report) - 1) = NULL;        /* erases '\n\' */
	    }
	}
    else
	{
	if (Stability(admin_table, member) == STABLE)
	    sprintf(end_of_string(node_report), "[");
	else if (Stability(admin_table, member) == UNSTABLE)
	    sprintf(end_of_string(node_report), "(");
	else if (Stability(admin_table, member) == BARELY_UNBOUND)
	    sprintf(end_of_string(node_report), "<");
	else if (Stability(admin_table, member) == STRONGLY_UNBOUND)
	    sprintf(end_of_string(node_report), "{");
	else if (Stability(admin_table, member) == UNKNOWN)
	    cerr << "print_node: stability of member " << member
		<< " is UNKNOWN\n";
	else
	    cerr << "print_node: invalid value "
		 << Stability(admin_table, member)
	         << " for stability of member " << member << "\n";
	
	for (i = 0; i < Ndaughter(admin_table, member); i++)
	    {
	    if (i > 0)
		sprintf(end_of_string(node_report), ", ");
	    print_node(node_report, admin_table, struct_table, name_table,
	    	       Daughter(admin_table, member, i));
	    }

	if (Stability(admin_table, member) == STABLE)
	    sprintf(end_of_string(node_report), "]");
	else if (Stability(admin_table, member) == UNSTABLE)
	    sprintf(end_of_string(node_report), ")");
	else if (Stability(admin_table, member) == BARELY_UNBOUND)
	    sprintf(end_of_string(node_report), ">");
	else if (Stability(admin_table, member) == STRONGLY_UNBOUND)
	    sprintf(end_of_string(node_report), "}");
	else if (Stability(admin_table, member) == UNKNOWN)
	    cerr << "print_node: stability of member " << member
		<< " is UNKNOWN\n";
	else
	    cerr << "print_node: invalid value "
		 << Stability(admin_table, member)
	         << " for stability of member " << member << "\n";
	}
    }

/*-----------------------------------------------------------------------------
 *  print_binary_parameters  --  
 *-----------------------------------------------------------------------------
 */
local void  print_binary_parameters(char node_report[BUFF_LENGTH],
				    real struct_table[NNODE_MAX][N_STRUCT],
				    int member)
    {
    int  offset_in_node_report;
    real  semimajor_axis;

    offset_in_node_report = end_of_string(node_report) - node_report;
/*
 * if the radius information starts in the 49th column, the output will be
 * nicely lined up with the distances.
 */
    while (offset_in_node_report < 49)
	node_report[offset_in_node_report++] = ' ';

    semimajor_axis = Apair(struct_table, member);

#if 0

    if (semimajor_axis < 10.0 && semimajor_axis > 0.0)
	sprintf(node_report + offset_in_node_report,"a = %.16f  ;  e = %.16f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 100.0 && semimajor_axis > -10.0)
	sprintf(node_report + offset_in_node_report,"a = %.15f  ;  e = %.16f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 1000.0 && semimajor_axis > -100.0)
	sprintf(node_report + offset_in_node_report,"a = %.14f  ;  e = %.16f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 10000.0 && semimajor_axis > -1000.0)
	sprintf(node_report + offset_in_node_report,"a = %.13f  ;  e = %.16f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 100000.0 && semimajor_axis > -10000.0)
	sprintf(node_report + offset_in_node_report,"a = %.12f  ;  e = %.16f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 1000000.0 && semimajor_axis > -100000.0)
	sprintf(node_report + offset_in_node_report,"a = %.11f  ;  e = %.16f\n",
		semimajor_axis, Epair(struct_table, member));
    else
	sprintf(node_report + offset_in_node_report,"a = %.10f  ;  e = %.16f\n",
		semimajor_axis, Epair(struct_table, member));

#else


    if (semimajor_axis < 10.0 && semimajor_axis > 0.0)
	sprintf(node_report + offset_in_node_report,"a = %.6f  ;  e = %8f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 100.0 && semimajor_axis > -10.0)
	sprintf(node_report + offset_in_node_report,"a = %.5f  ;  e = %8f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 1000.0 && semimajor_axis > -100.0)
	sprintf(node_report + offset_in_node_report,"a = %.4f  ;  e = %8f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 10000.0 && semimajor_axis > -1000.0)
	sprintf(node_report + offset_in_node_report,"a = %.3f  ;  e = %8f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 100000.0 && semimajor_axis > -10000.0)
	sprintf(node_report + offset_in_node_report,"a = %.2f  ;  e = %8f\n",
		semimajor_axis, Epair(struct_table, member));
    else if (semimajor_axis < 1000000.0 && semimajor_axis > -100000.0)
	sprintf(node_report + offset_in_node_report,"a = %.1f  ;  e = %8f\n",
		semimajor_axis, Epair(struct_table, member));
    else
	sprintf(node_report + offset_in_node_report,"a = %.0f  ;  e = %8f\n",
		semimajor_axis, Epair(struct_table, member));

#endif
    }

/*-----------------------------------------------------------------------------
 *  print_radius  --  
 * NOTE: CLEANUP THE SILLY OUTPUT SWITCHES WITH A NEW FIXED_WIDTH %f PROCEDURE
 *-----------------------------------------------------------------------------
 */
local void  print_radius(char node_report[BUFF_LENGTH],
			 real struct_table[NNODE_MAX][N_STRUCT], int member)
    {
    int  offset_in_node_report;
    real  r_node;

    offset_in_node_report = end_of_string(node_report) - node_report;
/*
 * if the radius information starts in the 49th column, the output will be
 * nicely lined up with the distances.
 */
    while (offset_in_node_report < 49)
	node_report[offset_in_node_report++] = ' ';

    r_node = Radius(struct_table, member);
    if (r_node < 10.0)
	sprintf(node_report + offset_in_node_report, "R = %.6f\n", r_node);
    else if (r_node < 100.0)
	sprintf(node_report + offset_in_node_report, "R = %.5f\n", r_node);
    else if (r_node < 1000.0)
	sprintf(node_report + offset_in_node_report, "R = %.4f\n", r_node);
    else if (r_node < 10000.0)
	sprintf(node_report + offset_in_node_report, "R = %.3f\n", r_node);
    else if (r_node < 100000.0)
	sprintf(node_report + offset_in_node_report, "R = %.2f\n", r_node);
    else if (r_node < 1000000.0)
	sprintf(node_report + offset_in_node_report, "R = %.1f\n", r_node);
    else
	sprintf(node_report + offset_in_node_report, "R = %.0f\n", r_node);
    }

/*-----------------------------------------------------------------------------
 *  how_much_offspring  --  computes the number of real particles below the
 *                          the node  member .
 *-----------------------------------------------------------------------------
 */
local int  how_much_offspring(int admin_table[NNODE_MAX][N_ADMIN], int member)
    {
    int  i, k, n_total;

    k = Ndaughter(admin_table, member);
    n_total = 0;
    
    if (k > 0)
	for (i = 0; i < k; i++)
	    n_total += how_much_offspring(admin_table,
					  Daughter(admin_table, member, i));
    else
	n_total = 1;

    return(n_total);
    }

#define  BUFF_SAFETY_MARGIN  28

/*-----------------------------------------------------------------------------
 *  end_of_string  --  
 *-----------------------------------------------------------------------------
 */
local char *end_of_string(char * the_string)
    {
    int  char_position;

    char_position = 0;
    while(the_string[char_position] != '\0')
	char_position++;

    if(char_position > BUFF_LENGTH - 1 - BUFF_SAFETY_MARGIN)
	{
        cerr << "end_of_string: too little room left in buffer\n";
        exit(1);
	}

    return(the_string + char_position);
    }

/*-----------------------------------------------------------------------------
 *  bring_to_normal_form  --  
 *                   example: 3,(7,[5,[4,8]],0),6,(2,1);  ==>
 *                            (0,[[4,8],5],7),(1,2),3,6;
 *                      NOTE: only one-digit star numbering implemented, i.e.
 *                            only valid for an N-body sub-system with N < 11
 *                  OH, WELL: a much cleaner way, of course, would be to first
 *                            map the string to an integer array of symbols,
 *                            where, e.g.,  ";" --> -1 , "(" --> -2 , 
 *                            ")" --> -3 , "[" --> -4 , "]" --> -5 , etc,
 *                            and the delimiting , is simply skipped since the
 *                            array elements are automatically separated.
 *-----------------------------------------------------------------------------
 */
local void  bring_to_normal_form(char old_hier_string[BUFF_LENGTH])
    {
    int  i, j;
    int  number_of_objects;
    int  ordered[BUFF_LENGTH];
    int  unordered[BUFF_LENGTH];
    char  new_hier_string[BUFF_LENGTH];
    char  substring[BUFF_LENGTH];
    char  head_char, tail_char;

    if (old_hier_string[1] == ';')
	return;                     /* a singleton is already in normal form */

    new_hier_string[0] = '\0';         /* initialize with a null_string, for */
    substring[0] = '\0';              /* end_of_string()  to work properly  */

    map_to_integers(old_hier_string, unordered, &number_of_objects);

    sort_int_array(unordered, ordered, number_of_objects);

    for (i = 0; i < number_of_objects; i++)
	{
	j = old_ranking(unordered, ordered, i, number_of_objects);
        select_object(substring, old_hier_string, j);
	unwrap_string(substring, &head_char, &tail_char);
	bring_to_normal_form(substring);
	wrap_string(substring, head_char, tail_char);
	append_object(new_hier_string, substring);
	}
/*
 * copy the new string back unto the old string:
 */
    i = 0;    
    while (new_hier_string[i] != '\0' && i < BUFF_LENGTH)
	{
	old_hier_string[i] = new_hier_string[i];
	i++;
	}
    if (i >= BUFF_LENGTH)
	{
        cerr << "bring_to_normal_form: buffer length exceeded\n";
        exit(1);
	}
    if (old_hier_string[i] != '\0')
	{
        cerr << 
	      "bring_to_normal_form: new and old hier. string diff. length\n";
        exit(1);
	}
    }

/*-----------------------------------------------------------------------------
 *  map_to_integers  --  
 *                 NOTE: only one-digit star numbering implemented, i.e.
 *                       only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  map_to_integers(char hier_string[BUFF_LENGTH],
			    int unordered[BUFF_LENGTH], int * nptr)
//int *nptr;                           /* will count total number of objects */
    {
    int  i;
    int  lowest_number;
    int  level;
    char  c;

    *nptr = 0;
    level = 0;
    lowest_number = 10;        /* higher than highest allowed number */
    i = 0;
    while (i < BUFF_LENGTH)
	{
	c = hier_string[i++];

	if (c >= '0' && c <= '9')
	    {
	    if (c - '0' < lowest_number)
		lowest_number = c - '0';
	    }
	else if (c == '(' || c == '[' || c == '{' || c == '<')
	    level++;
	else if (c == ')' || c == ']' || c == '}' || c == '>')
	    level--;
	else if (level == 0 && (c == ',' || c == ';'))
	    {
	    unordered[*nptr] = lowest_number;
            (*nptr)++;
	    lowest_number = 10;
	    }

	if (c == ';')
	    break;
	}
    if (i >= BUFF_LENGTH)
	{
        cerr << "map_to_integers: buffer length exceeded\n";
        exit(1);
	}
    }

/*-----------------------------------------------------------------------------
 *  sort_int_array  --  
 *                NOTE: only one-digit star numbering implemented, i.e.
 *                      only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  sort_int_array(int old_array[BUFF_LENGTH],
			   int new_array[BUFF_LENGTH], int n)
    {
    int  i, j;
    int  dummy;

    for (i = 0; i < n; i++)
	new_array[i] = old_array[i];

    for (i = 0; i < n - 1; i++)         /* lazy programmer's sorting method  */
	for (j = i+1; j < n; j++)
	    if (new_array[j] < new_array[i])
		{
		dummy = new_array[i];
		new_array[i] = new_array[j];
		new_array[j] = dummy;
		}
    }

/*-----------------------------------------------------------------------------
 *  old_ranking  --  
 *               NOTE: only one-digit star numbering implemented, i.e.
 *                     only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local int  old_ranking(int old_array[BUFF_LENGTH], int new_array[BUFF_LENGTH],
			int i, int n)
    {
    int  j;

    j = 0;
    while (j < n)
	if (old_array[j++] == new_array[i])
	    break;
    if (j > n)
	{
        cerr << "old_ranking: used buffer length exceeded\n";
        exit(1);
	}

    return(--j);
    }

/*-----------------------------------------------------------------------------
 *  select_object  --  
 *               NOTE: only one-digit star numbering implemented, i.e.
 *                     only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  select_object(char tmp_string[BUFF_LENGTH],
			  char hier_string[BUFF_LENGTH], int member)
    {
    int  i, j;
    int  level;
    char  c;

    level = 0;
    i = 0;
    while (member > 0)
	{
	c = hier_string[i];
	if (c == '(' || c == '[' || c == '{' || c == '<')
	    level++;
	else if (c == ')' || c == ']' || c == '}' || c == '>')
	    level--;
	else if (c == ',' && level == 0)
	    member--;
	else if (c == ';')
	    cerr << "select_object: reached end of string\n";
        if (++i >= BUFF_LENGTH)
	    cerr << "select_object: at position #1: buffer length exceeded\n";
	}

    if (level != 0)
	{
        cerr << "select_object: number of (,[,{,< and of ),],},> unbalanced\n";
        exit(1);
	}
    j = 0;
    while (c = hier_string[i])
	{
	if (c == '(' || c == '[' || c == '{' || c == '<')
	    level++;
	else if (c == ')' || c == ']' || c == '}' || c == '>')
	    level--;
	else if ((c == ',' && level == 0) || c == ';')
	    {
	    sprintf(tmp_string + j, ";");      /* end with ';' and  NULL  !! */
	    if (level != 0)
		{
	        cerr << "select_object: number of (,[, etc. unbalanced\n";
	        exit(1);
		}
	    break;
	    }
	tmp_string[j++] = c;

        if (++i >= BUFF_LENGTH)
	    cerr << "select_object: at position #2: buffer length exceeded\n";
	}
    }

/*-----------------------------------------------------------------------------
 *  append_object  --  
 *               NOTE: only one-digit star numbering implemented, i.e.
 *                     only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  append_object(char hier_string[BUFF_LENGTH],
			  char tmp_string[BUFF_LENGTH])
    {
    char *last_endptr;
    char *new_startingptr;
    char *end_of_string(char *);

    if (hier_string[0] == '\0')
	new_startingptr = hier_string;
    else
	{
	last_endptr = end_of_string(hier_string) - 1;
	if (*last_endptr == ';')
	    *last_endptr = ',';                      /* insert new delimiter */
	else	    
	    cerr << "append_object: hier_string not properly terminated\n";
	new_startingptr = last_endptr + 1;
	}

    sprintf(new_startingptr, tmp_string);
    }

/*-----------------------------------------------------------------------------
 *  unwrap_string  --  takes off the leading and trailing bracket,
 *                     if present, and only if the string 'substring' contains
 *                     only one object (the first two 'return;' statements
 *                     take care of these two exceptions).
 *                     If no enclosing brackets are present, the head and tail
 *                     characters '*head_ptr' and '*tail_ptr' are assigned the
 *                     value '*' to indicate the absence of brackets.
 *               NOTE: only one-digit star numbering implemented, i.e.
 *                     only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  unwrap_string(char substring[BUFF_LENGTH],
			  char * head_ptr, char * tail_ptr)
    {
    int  i;
    int  level;
    char  c;

    c = substring[0];

    if (c == '(' || c == '[' || c == '{' || c == '<')
	*head_ptr = c;
    else
	{
	*head_ptr = *tail_ptr = '*';
	return;                         /* no enclosing brackets of any kind */
	}

    level = 1;
    i = 1;
    for (;;)
	{
	c = substring[i];
	if (c == '(' || c == '[' || c == '{' || c == '<')
	    level++;
	else if (c == ')' || c == ']' || c == '}' || c == '>')
	    level--;

	if (level == 0)
	    {
	    if (substring[i+1] != ';')    /* more than one object,           */
		{                         /* therefore no enclosing brackets */
		*head_ptr = *tail_ptr = '*';
		return;
		}
	    else
		break;
	    }

        if (++i >= BUFF_LENGTH)
	    cerr << "unwrap_string: buffer length exceeded\n";
	}

    i = 1;
    while ((c = substring[i]) != ';')
	{
	substring[i-1] = c;
	i++;
	}
    *tail_ptr = substring[i-2];
    substring[i-2] = ';';
    substring[i-1] = '\0';                    /* proper string ending */
    }

/*-----------------------------------------------------------------------------
 *  wrap_string  --  puts the leading and trailing bracket in place again,
 *                   that is, if there are any brackets.
 *       convention: the character '*' for the value of 'head_char' and 
 *                   'tail_char' indicates the absence of brackets.
 *             NOTE: only one-digit star numbering implemented, i.e.
 *                   only valid for an N-body sub-system with N < 11
 *-----------------------------------------------------------------------------
 */
local void  wrap_string(char substring[BUFF_LENGTH],
			char head_char, char tail_char)
    {
    int  i;
    char  last_c;
    char  c;

    if (head_char == '*')
	return;

    last_c = head_char;
    i = 0;
    while((c = substring[i]) != ';')
	{
	substring[i] = last_c;
	last_c = c;
	i++;
	}
    substring[i] = last_c;
    substring[++i] = tail_char;
    substring[++i] = ';';
    substring[++i] = '\0';
    }

main(int argc, char ** argv)
{
    check_help();
    pgetopt(argc, argv, "", "$Revision: 1.8 $", _SRC_);

    if (argc != 1) {
	cerr << "molecules: no arguments allowed\n";
	exit(1);
    }

    dyn *b =  get_dyn();

    printf("\n");

    while (b) {			// return quietly if no new body is found

        dyn *bi;

        for (bi = b->get_oldest_daughter(); bi != NULL;
             bi = bi->get_younger_sister())
            if (bi->get_oldest_daughter() != NULL)
		{
		cerr << "molecules: flat tree required, sorry!\n";
		exit(1);
		}

        if (find_qmatch(b->get_dyn_story(), "t"))
	    printf("Time = %15g   \n", getrq(b->get_dyn_story(), "t"));

	find_group_internal_structure(b);

	rmtree(b);

        b = get_dyn();
    }
}

#endif

// endof: molecules.C
