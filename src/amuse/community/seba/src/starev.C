//// starev: evolve a single star.
////         creates a single star and evolves it in time.
////
//// Options:    -c    comment to put in the starbase log structure.
////             -M    initial mass of the star [in solar units].
////             -n    number of output timesteps (timesteps are taken
////                   with constant time intervals) 
////             -R    Dynamical size scaling for the star
////                   [in units of the virial radius].
////             -S    Random seed.
////             -s    Initial stellar type [default is main_sequence].
////             -T or -t end time of the stellar evolution [in Million year].
////
//// Latest version (SPZ:1.0) February 1993.

//++ the single stars are part of the dynamical tree.
//++ Since the assymetry in supernovae are taken care of by the 
//++   dynamical model no kicks are applied.

//   version 1.0:  Februari 1993   Simon F. Portegies Zwart
//                                 spz@grape.c.u-tokyo.ac.jp 
//
//#include "dyn.h"
#include "node.h"
#include "single_star.h"
#include "main_sequence.h"
//#include "sstar_to_dyn.h"
#define EPSILON 1.e-10

#ifdef TOOLBOX
local bool read_single_params(ifstream& in, real &mass, real &time, real &metal) {
    mass = 0;
    metal = 0;
    
    while (mass > 100 || metal < 0.0001 || metal > 0.03){
        if(in.eof())
            return false;
        
        // reading from input file
        in >> mass>>time>>metal;
    }

    PRC(mass);PRC(time);PRL(metal);
    return true;
}


/*-------------------------------------------------------------------------*/
local void evolve_star_until_next_time(node* bi, const real out_time, const int n_steps) {
    //FILE *f1;
    //f1=fopen("binev.data", "a");    
    //fprintf(f1, "something");
    //fclose(f1);
    ofstream starev("data/starev.data", ios::app|ios::out);
    bi->get_starbase()->dump(starev, false);  
    real current_time = ((star*)bi->get_starbase())->get_current_time();
    real time_step    =  bi->get_starbase()->get_evolve_timestep();

    while (out_time>current_time+time_step ) {
        bi->get_starbase()->evolve_element(current_time+time_step);
        bi->get_starbase()->dump(starev, false);                
        current_time = ((star*)bi->get_starbase())->get_current_time();
        time_step    =  bi->get_starbase()->get_evolve_timestep();
        
        star_state ss(dynamic_cast(star*, bi->get_starbase()));
    }
    
      
    bi->get_starbase()->evolve_element(out_time);
    bi->get_starbase()->dump(starev, false);
    bi->get_starbase()->dump(cerr, false);
    print_star(bi->get_starbase(), cerr);
    starev.close();
}

/*-----------------------------------------------------------------------------
 *  main  --
 *	usage:
 *		addstar -t # [options]  ,
 *
 *		where # is the initial age of the cluster.
 *	options:
 *            	The following options are allowed:
 *	cluster age:
 *		-t #	Where # stands for the initial age of the
 *			in Myear.
 *
 *		At present the running time of the integrator correspnds
 *		to the stellar age an a one by 10e6year basis.
 *		This however should be scaled to the cluster parameters.
 *-----------------------------------------------------------------------------
 */
int main(int argc, char ** argv)
    {
        stellar_type type = Main_Sequence;
    char * star_type_string;
    int  c;

    bool  t_flag = FALSE;
    bool  S_flag = FALSE;
    bool  c_flag = FALSE;
    bool  M_flag = FALSE;
    bool  n_flag = FALSE;
    bool  R_flag = FALSE;
    bool  I_flag = FALSE;
    real  m_tot;
    real  r_hm = 100;
    real  t_hc = 1;
    real  t_start = 0;           // default value;
    real  t_end;
    
    int n_steps = 1;
    int n_steps_per_phase = 10;
    int n_init = 0;
    int n =1;
    char* input_filename;
    real z;

    real mass=1;
    real endtime=100;
    real metal= cnsts.parameters(Zsun);
    
    char  *comment;
    int input_seed=0, actual_seed;
    extern char *poptarg;
    const char *param_string = "d:D:M:R:T:t:S:s:N:I:n:c:z:";

    //check_help();
    
    if (argc <= 1)
        {
            cerr <<"usage: starev -M # -R # -T # -t # -S # -s # -N # -I # -n # -z #[-c \"..\"]\n";
        exit(1);
        }
 
    while ((c = pgetopt(argc, argv, param_string)) != -1)
	switch(c)
	    {
            case 'd': n_steps = atoi(poptarg);//delta steps
            break;
            case 'D': n_steps_per_phase= atoi(poptarg);//delta steps
            break;
            case 'M': M_flag = TRUE;
	              mass = atof(poptarg);
                      break;
            case 'R': r_hm = atof(poptarg);
                      break;
            case 'T':                      
            case 't': endtime = atof(poptarg);
                      break;
            case 'S': S_flag = TRUE;
                      input_seed = atoi(poptarg);
                      break;
            case 's': 
                strcpy(star_type_string, poptarg);
		      type = extract_stellar_type_string(star_type_string);
                      break;
            case 'N': n_init = atoi(poptarg);
                      break;
            case 'I': I_flag = true;
                      input_filename = poptarg;
                      break;
            case 'n': n = atoi(poptarg);
                      break;
            case 'z': metal = atof(poptarg);
                      break;
            case 'c': c_flag = TRUE;
		      comment = poptarg;
		      break;
            case '?': params_to_usage(cerr, argv[0], param_string);
	      //get_help();
	    	      exit(1);
	    }
    ifstream infile(input_filename, ios::in);
    if(I_flag) {
        if (!infile) {
            cerr << "error: couldn't read file " << input_filename <<endl;
            exit(-1);
        }
        else cerr << "Reading input from file "<< input_filename <<endl;
    }

    cerr.precision(HIGH_PRECISION);

    if(!S_flag) actual_seed = 0;
    actual_seed = srandinter(input_seed);

    // make flat tree 
    node *root;
    root= mknode(1);
    root->log_history(argc, argv);
    
    for (int i=0;i<n;i++){
    if(I_flag) {
        if (read_single_params(infile, m_tot, t_end, z)){
            n=i+2;
        }
        else 
            break;
    }
    /*else if{ // random distribution of n stars
        
    }*/
    else { 
        if (mass<=100.0 && metal >= 0.0001 && metal <= 0.03 ){
            m_tot = mass;
            t_end = endtime;
            z = metal;
            n = 1;//if no input file, then evolve only 1 star
        }
        else{
            cerr<<"Parameters are not within valid range"<<endl;    
            cerr<<"M <= 100 "<<endl;
            cerr<<"0.0001 <= z <= 0.03"<<endl;
            return 0;
        }
     }
    root= mknode(1);
    //root->set_mass(1);
    root->get_starbase()->set_stellar_evolution_scaling(m_tot, r_hm, t_hc);

    //node *the_star = root->get_oldest_daughter();

    addstar(root, t_start, type, z, n_init+i, false);
    
    // Starev does not include hdyn.h nor the hdyn library
    // get_use_hdyn is therefore not defined.
    // The result is that kick velocities will be scaled
    // spuriously......
    root->get_starbase()->set_use_hdyn(false);
    cerr.precision(STD_PRECISION);

    //    put_node(root);

    real delta_t = t_end/((real)n_steps);
    real out_time; 
    for_all_daughters(node, root, bi) {
       out_time = 0;
       do {
            out_time = Starlab::min(out_time+delta_t, t_end);
            evolve_star_until_next_time(bi, out_time, n_steps_per_phase);
       }
       while(out_time < t_end);
    }

    for_all_daughters(node, root, bi) {
       cerr << "Time = " << bi->get_starbase()->get_current_time()
	    << " [Myear],  mass = " << bi->get_starbase()->get_total_mass() 
	    << " [Msun],  radius = " << bi->get_starbase()
                                          ->get_effective_radius() 
            << "   " << type_string(bi->get_starbase()->get_element_type())
	    << endl;
    }
   }//for loop
//    put_node(root);
    delete root;
    return 0;
}

#endif
