
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// addstoa.C: Applies nemo snapshot dynamics information to 
//		 an existing snapshot.  No scaling is performed.
//
//		 Steve McMillan, July 1996
//

#include "dyn.h"

#ifdef TOOLBOX


local void  addstoa(dyn * b, int m_flag)
{
    real time;
    int ndim;
    int n;
    int n_stoa;


    n = 0;
    for_all_daughters(dyn, b, bi) n++;

    PRL(n);
    cin >> n_stoa;
    cin >> ndim;
    cin >>  time;
    PRC(n_stoa); PRC(ndim); PRL(time);
    if (n != n_stoa){
	err_exit("n and n_stoa must be equal");
    }
    for_all_daughters(dyn, b, bi) {
	real mass;
	cin >> mass;
	if(m_flag)bi->set_mass(mass);
    }
    for_all_daughters(dyn, b, bi) {
	vec p;
	cin >> p;
	bi->set_pos(p);
    }
    for_all_daughters(dyn, b, bi) {
	vec v;
	cin >> v;
	bi->set_vel(v);
    }
}

main(int argc, char ** argv)
{
    int m_flag = false;
    extern char *poptarg;
    int c;
    const char *param_string = "m:";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.7 $", _SRC_)) != -1)
	switch(c) {

	    case 's': m_flag= true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
		      exit(1);
	}            
    
    dyn *b;
    b = get_dyn();

    b->log_history(argc, argv);

    addstoa(b, m_flag);

    put_dyn(b);
}

#endif
