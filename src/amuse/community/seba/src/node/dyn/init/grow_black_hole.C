
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Replace the specified inner mass in stars with a single black hole
//// of the same mass.
////                
//// Usage:  grow_black_hole [OPTIONS]
////
//// Options:
////        -M    select black hole mass (if <1 mass is read as a fraction)
////
//// Written by Simon Portegies Zwart
////
//// Report bugs to starlab@sns.ias.edu.

//                 Simon Portegies Zwart, MIT June 2000

#include "dyn.h"

#define  SEED_STRING_LENGTH  256
char  tmp_string[SEED_STRING_LENGTH];

#ifdef TOOLBOX

typedef  struct
{
    real  radius;
    real  mass;
} rm_pair, *rm_pair_ptr;

//-----------------------------------------------------------------------------
//  compare_radii  --  compare the radii of two particles
//-----------------------------------------------------------------------------

local int compare_radii(const void * pi, const void * pj)  // increasing radius
{
    if (((rm_pair_ptr) pi)->radius > ((rm_pair_ptr) pj)->radius)
        return +1;
    else if (((rm_pair_ptr)pi)->radius < ((rm_pair_ptr)pj)->radius)
        return -1;
    else
        return 0;
}

//-----------------------------------------------------------------------------
//  radius_containting_mass  --  Get the massradii for all particles.
//-----------------------------------------------------------------------------

real radius_containting_mass(dyn * b, real mass) {

    char lagr_string[64] = "central black hole";

    int n = b->n_daughters();

    // Use the density center if known and up to date (preferred).
    // Otherwise, use modified center of mass, if known and up to date.
    // Otherwise, use the geometric center.

    vec lagr_pos = 0, lagr_vel = 0;
    bool try_com = true;
    bool modify_com = false;

#if 0
    if (find_qmatch(b->get_dyn_story(), "density_center_pos")) {

	if (getrq(b->get_dyn_story(), "density_center_time")
		!= b->get_system_time())
	    warning("lagrad: neglecting out-of-date density center");
	else {
	    lagr_pos = getvq(b->get_dyn_story(), "density_center_pos");
	    strcpy(lagr_string, "density center");
	    try_com = false;

	    // Assume that density_center_vel exists if density_center_pos
	    // is OK.

	    lagr_vel = getvq(b->get_dyn_story(), "density_center_vel");
	}
    }

    if (try_com && find_qmatch(b->get_dyn_story(), "com_pos")) {

	if (getrq(b->get_dyn_story(), "com_time")
		!= b->get_system_time()) {
	    warning("lagrad: neglecting out-of-date center of mass");
	} else {
	    lagr_pos = getvq(b->get_dyn_story(), "com_pos");
	    lagr_vel = getvq(b->get_dyn_story(), "com_vel");
	    strcpy(lagr_string, "center of mass");
	    modify_com = true;
	}
    }
#endif

    rm_pair_ptr rm_table = new rm_pair[n];

    if (rm_table == NULL) {
	cerr << "radius_containing_mass: "
	     << "not enough memory left for rm_table\n";
	return -1;
    }

    // Set up an array of (radius, mass) pairs.  Also find the total
    // mass of all nodes under consideration.

    real total_mass = 0;
    int i = 0;

    for_all_daughters(dyn, b, bi) {
      total_mass += bi->get_mass();
      rm_table[i].radius = abs(bi->get_pos() - lagr_pos);
      rm_table[i].mass = bi->get_mass();
      i++;
    }

    // Sort the array by radius.
    qsort((void *)rm_table, (size_t)i, sizeof(rm_pair), compare_radii);
    
    real rlagr = 0;
    real cumulative_mass = 0.0;
    
    // cerr << "Determining Lagrangian radii 2" << endl << flush;
    
    i=0;
    while (cumulative_mass < mass)
      cumulative_mass += rm_table[i++].mass;

    rlagr = rm_table[i-1].radius;

    delete [] rm_table;

    PRL(rlagr);
    return rlagr;
}

void merge_bh_with_coll(dyn *bh, dyn *bcoll) {

    detach_node_from_general_tree(bcoll);
    bcoll->set_younger_sister(NULL);
    bcoll->set_elder_sister(NULL);

    real mass = bh->get_mass() + bcoll->get_mass();

    real f1 = bh->get_mass() / mass;
    real f2 = bcoll->get_mass() / mass;
    vec pos = f1 * bh->get_pos() + f2 * bcoll->get_pos();
    vec vel = f1 * bh->get_vel() + f2 * bcoll->get_vel();
    vec acc = f1 * bh->get_acc() + f2 * bcoll->get_acc();

    // Set up the physical coordinates of the daughters.

    bh->set_mass(mass);
    bh->set_pos(bh->get_pos() - pos);
    bh->set_vel(bh->get_vel() - vel);
    bh->set_acc(bh->get_acc() - acc);

    bcoll->set_pos(bcoll->get_pos() - pos);
    bcoll->set_vel(bcoll->get_vel() - vel);
    bcoll->set_acc(bcoll->get_acc() - acc);
  }


local void grow_black_hole(dyn* b, real m_bh) {

    PRL(m_bh);
    if(m_bh<=1) {
	cerr << "fractional bh mass" << endl;
	m_bh *= b->get_mass();
	PRL(m_bh);
    }

    real bh_radius = radius_containting_mass(b, m_bh);
    PRL(bh_radius);


    real sum = 0;

    int ibh = 0;
    dyn *bh = NULL;
    int n = b->n_daughters();
    do {
      for_all_daughters(dyn, b, bi) {
	if(abs(bi->get_pos())<=bh_radius) {
	  ibh++;
	  if (bh==NULL) 
	    bh = bi;
	  else if(bh!=bi) {
//	    cerr << "merge bh with " << bi->format_label() << endl;
	    merge_bh_with_coll(bh, bi);
//	    PRL(bh->get_mass());
//	    put_dyn(bh, cerr);
	    break;
	  }
	}
      }
    }
    while(bh->get_mass()<m_bh);

    putiq(bh->get_log_story(), "black_hole", 1);
  }

int main(int argc, char ** argv) {

    real m_bh = 0;

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "M:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.8 $", _SRC_)) != -1)
	switch(c) {
	    case 'M': m_bh = atof(poptarg);
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}


    dyn* b;
    b = get_dyn();
    b->log_history(argc, argv);

    grow_black_hole(b, m_bh);

    put_dyn(b);
    rmtree(b);
    return 0;
}
#endif
