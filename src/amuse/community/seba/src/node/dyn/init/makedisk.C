
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

//// Construct a near-keplerian disk, with N low-mass objects
//// orbiting a single massive object on almost circular paths.
////
//// Units:  masses are measured in millions of solar masses,
//// lengths are measured in parsecs.  Hence (with G = 1) the
//// time unitis 1.49e4 yr and the velocity unit is 65.7 km/s.
////
//// Usage:  makedisk [OPTIONS]
////
//// Options:
////              -c    add a comment to the output snapshot [false]
////              -C    output data in 'col' format [no]
////              -i    number the particles sequentially [don't number]
////              -l    specify particle radius [0]
////              -m    specify the total mass of the disk [no default]
////              -M    specify the mass of the central object [no default]
////              -n    specify number of disk particles [no default]
////              -o    echo value of random seed [don't echo]
////              -r    specify inner disk radius [no default]
////              -R    specify outer disk radius [no default]
////              -s    specify random seed [random from system clock]
////              -v    specify 3-D disk velocity dispersion [0]
////              -V    verbose mode [off]
////
//// Written by Steve McMillan.
////
//// Report bugs to starlab@sns.ias.edu.

#define Gcnst	6.67e-8		// cgs units
#define PCcm	3.086e18	// cm
#define MSUN	1.989e33	// g
#define YR	3.156e7		// s

#include "dyn.h"

#ifdef TOOLBOX

local void makedisk(dyn* b, int n,
		    real m_central, real m_disk,
		    real r_inner, real r_outer, real radius, real v_disp)
{
    real pmass = m_disk / n;
    v_disp *= sqrt(3.0);	// so a uniform distribution has proper v_disp;

    vec cmpos = 0, cmvel = 0;

    for_all_daughters(dyn, b, bi) {

	if (bi->get_elder_sister() == NULL) {

	    // Central object.

	    bi->set_mass(m_central);
	    bi->set_pos(0);
	    bi->set_vel(0);

	    // Particle 0 is a massive black hole.  We can include
	    // tidal destruction within the kira collision framework
	    // by giving the hole a radius such that it has the same
	    // "density" as the other particles in the system.  This
	    // works so long as the radius of the hole thus defined is
	    // much greater than the radius of any other particle in
	    // the system.

	    putrq(bi->get_dyn_story(), "R_eff",
		  radius * pow(m_central/pmass, 1.0/3));

	    // Note that R_eff is meaningful only to hdyn tools.

	} else {

	    // Disk particle; disk lies in the (x-y) plane, mean motion
	    // in the positive sense.

	    bi->set_mass(pmass);	

	    real r = sqrt(r_inner*r_inner
			     + randinter(0,1)
			  	  * (r_outer*r_outer - r_inner*r_inner));
	    real theta = randinter(0, 2*M_PI);
	    bi->set_pos(vec(r*cos(theta), r*sin(theta), 0));

	    real v_orb = sqrt(m_central/r);
	    bi->set_vel(vec(-v_orb*sin(theta) + v_disp*randinter(-1,1),
			        v_orb*cos(theta) + v_disp*randinter(-1,1),
						   v_disp*randinter(-1,1)));

	    // All particles have the same mass and radius, and hence
	    // the same density.

	    putrq(bi->get_dyn_story(), "R_eff", radius);

	    cmpos += pmass * bi->get_pos();
	    cmvel += pmass * bi->get_vel();
	}
    }

    b->set_mass(m_central + m_disk);

    // Force the system CM to be at rest at the origin.

    cmpos /= (m_central + m_disk);
    cmvel /= (m_central + m_disk);

    for_all_daughters(dyn, b, bi) {
	bi->inc_pos(-cmpos);
	bi->inc_vel(-cmvel);
    }
}

#define  SEED_STRING_LENGTH  60

main(int argc, char ** argv)
{
    bool c_flag = false;
    bool C_flag = false;
    bool n_flag = false;
    bool s_flag = false;
    bool i_flag = false;
    bool l_flag = false;
    bool o_flag = false;
    bool verbose = false;

    int  i;
    int  n = 0;
    int  input_seed, actual_seed;

    real m_central = 0, m_disk = 0;
    real r_inner = 0, r_outer = 0, radius = 0;
    real v_disp = 0;

    char  *comment;
    char  seedlog[SEED_STRING_LENGTH];

    check_help();

    extern char *poptarg;
    int c;
    const char *param_string = "c:Cil:m:M:n:or:R:s:v:V";

    while ((c = pgetopt(argc, argv, param_string,
		    "$Revision: 1.11 $", _SRC_)) != -1)
	switch(c)
	    {
	    case 'c': c_flag = true;
		      comment = poptarg;
		      break;
	    case 'C': C_flag = true;
		      break;
	    case 'i': i_flag = true;
		      break;
	    case 'l': radius = atof(poptarg);
		      break;
	    case 'm': m_disk = atof(poptarg);
		      break;
	    case 'M': m_central = atof(poptarg);
		      break;
	    case 'n': n = atoi(poptarg);
		      break;
	    case 'o': o_flag = true;
                      break;
	    case 'r': r_inner = atof(poptarg);
		      break;
	    case 'R': r_outer = atof(poptarg);
		      break;
	    case 's': s_flag = true;
		      input_seed = atoi(poptarg);
		      break;
	    case 'v': v_disp = atof(poptarg);
		      break;
	    case 'V': verbose = true;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	              get_help();
                      exit(1);
	    }            
    
    if (n <= 0) {
        cerr << "makedisk: must specify number of particles with -n#\n";
	exit(1);
    }
    
    if (m_disk <= 0) {
        cerr << "makedisk: must specify disk mass with -m\n";
	exit(1);
    }
    
    if (m_central <= 0) {
        cerr << "makedisk: must specify central mass with -M\n";
	exit(1);
    }
    
    if (r_inner <= 0) {
        cerr << "makedisk: must specify inner disk radius with -r\n";
	exit(1);
    }
    
    if (r_outer <= 0) {
        cerr << "makedisk: must specify outer disk radius with -R\n";
	exit(1);
    }
    
    if (r_outer <= r_inner) {
        cerr << "makedisk: must specify r_outer > r_inner\n";
	exit(1);
    }

    if (verbose) {

	real t_unit = sqrt(pow(PCcm,3)/(Gcnst*1.e6*MSUN));		// unit: s
	cerr << "[M] = 1.e6 solar masses" << endl;
	cerr << "[L] = 1 pc" << endl;
	cerr << "[T] = " << t_unit/YR << " yr" << endl;
	real v_unit = 1.e-5*PCcm/t_unit;
	cerr << "[v] = " << v_unit << " km/s" << endl;

	cerr << endl;

	PRL(n);
	PRC(m_central), PRL(m_disk);
	PRC(r_inner), PRC(r_outer), PRL(radius);
	cerr << "v_disp = " << v_disp
	     << " = " << v_disp*v_unit << " km/s" << endl;

	cerr << endl;

	real p_inner = 2*M_PI*sqrt(pow(r_inner,3)/m_central);
	cerr << "orbital period at disk inner edge = " << p_inner
	     << " = " << p_inner*t_unit/YR << " yr" << endl;

	real p_outer = 2*M_PI*sqrt(pow(r_outer,3)/m_central);
	cerr << "orbital period at disk outer edge = " << p_outer
	     << " = " << p_outer*t_unit/YR << " yr" << endl;

	cerr << endl;

	real v_inner = sqrt(m_central/r_inner);
	cerr << "orbital speed at disk inner edge = " << v_inner
	     << " = " << v_inner*v_unit << " km/s" << endl;

	real v_outer = sqrt(m_central/r_outer);
	cerr << "orbital speed at disk outer edge = " << v_outer
	     << " = " << v_outer*v_unit << " km/s" << endl;

	if (radius > 0) {
	    real v_esc = sqrt(2*(m_disk/n)/radius);
	    cerr << "cloud escape speed = " << v_esc
		 << " = " << v_esc*v_unit << " km/s" << endl;
	}
    }
    
    dyn *b, *by, *bo;
    b = new dyn();				// root node
    b->set_root(b);

    if (C_flag) b->set_col_output(true);

    bo = new dyn();				// central mass
    if (i_flag) 
        bo->set_label(0);
    b->set_oldest_daughter(bo);
    bo->set_parent(b);

    for (i = 1; i <= n; i++) {			// disk particles
        by = new dyn();
	if (i_flag) by->set_label(i);
	by->set_parent(b);
        bo->set_younger_sister(by);
        by->set_elder_sister(bo);
        bo = by;
    }

    if (c_flag == true) b->log_comment(comment);

    b->log_history(argc, argv);

    if (s_flag == false) input_seed = 0;
    actual_seed = srandinter(input_seed);

    if (o_flag) cerr << "makedisk: random seed = " << actual_seed << endl;
    sprintf(seedlog, "       random number generator seed = %d",actual_seed);
    b->log_comment(seedlog);

    makedisk(b, n,
	   m_central, m_disk,
	   r_inner, r_outer, radius, v_disp);

    put_dyn(b);
}

#endif
