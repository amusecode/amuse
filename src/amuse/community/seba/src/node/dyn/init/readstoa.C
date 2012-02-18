
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Convert NEMO "stoa" format:
////
////           (N, ndim, time,
////            mass[i], i = 1,...,N,
////            pos[i],  i = 1,...,N,
////            vel[i],  i = 1,...,N)
////
//// data into a Starlab snapshot.
////
//// Usage:  readstoa [OPTIONS]
////
//// Options:
////           -i    number the particles sequentially [don't number]
////           -w    write stoa insted of reading it
////
//// Written by Jun Makino.
////
//// Report bugs to starlab@sns.ias.edu.


//   Jun Makino, Aug 1996

#include "dyn.h"

#ifdef TOOLBOX

local dyn* read_stoa(bool i_flag) {

    dyn *root, *by, *bo;

    // Create root node.

    root = new dyn();
    if (i_flag) root->set_label(0);

    int n; scanf("%d", &n); PRL(n);
    int ndim; scanf("%d", &ndim); PRL(ndim);
    real time; scanf("%lf", &time); PRL(time);
    root->set_system_time(time);
    
    // Create first daughter node.

    bo = new dyn();
    root->set_oldest_daughter(bo);
    bo->set_parent(root);
    if (i_flag) bo->set_label(1);

    // Create other daughter nodes.

    for (int i = 1; i < n; i++) {
        by = new dyn();
	if (i_flag) by->set_label(i+1);
	by->set_parent(root);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
	by->set_mass(bo->get_mass());
        bo = by;
    }

    real total_mass = 0;

    for_all_daughters(dyn, root, b) {
	real mass; cin >> mass; 
	b->set_mass(mass);
	total_mass += mass;
    }

    root->set_mass(total_mass);

    for_all_daughters(dyn, root, b) {
	vec pos; cin >> pos; 
	b->set_pos(pos);
    }

    for_all_daughters(dyn, root, b){
	vec vel; cin >> vel; 
	b->set_vel(vel);
    }

    return root;
  }

local void write_stoa() {

  real ndim = 3; // x, y and z

  dyn *b;
  b = get_dyn();
  
  b->flatten_node();
  
  // for safetly check number of leaves.
  int n=0;
  for_all_daughters(dyn, b, bi) {
    n++;
  }
  
  //start writing NEMO stoa
  cout << n << endl;
  cout << ndim << endl;
  cout << b->get_system_time() << endl;
  for_all_daughters(dyn, b, bi) {
    cout << bi->get_mass() << endl;
  }
  
  for_all_daughters(dyn, b, bi) {
    cout << bi->get_pos() << endl;
  }
  
  for_all_daughters(dyn, b, bi){
    cout << bi->get_vel() << endl;
  }
  
}

int main(int argc, char ** argv)
{
    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "iw";

    bool i_flag = false;
    bool w_flag = false;

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.10 $", _SRC_)) != -1)
	switch(c) {

	    case 'i': i_flag = true;
		      break;
	    case 'w': w_flag = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}

    dyn *root;

    if(w_flag) {
        write_stoa();
    }
    else {
        root = read_stoa(i_flag);
        root->log_history(argc, argv);
        put_dyn(root);
    }
    return 0;
}

#endif

/* end of: readstoa.c */
