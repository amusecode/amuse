//
// nbody.C
//
// Version 1999/2/9 Jun Makino
//
// Change in the calling format for apply_vf, from function name to
// add explicit address operator. This was necessary to pass g++ 2.8.1
// 
// Version 1.1 Jun Makino : take parameters from command line argumensts
//                          Coding for 1.1 started on 1998/12/31
// Version 1.0 Jun Makino : First publisized version, parameters
//                          read in from the standard input

#include "BHTC.h"
using namespace std;

//#ifdef TOOLBOX

//------------------------------------------
//
// list of options
// -i        name of snapshot input file       (no default)
// -o        name of snapshot output file      (default: no output)
// -d        timestep                          (default: 0.015625)
// -D        time interval for snapshot output (default: 1)
// -T        time to stop integration          (default: 10)
// -e        softening parameter epsilon2      (default: 0.025)
// -t        opening angle theta               (default: 0.75)
// -n        ncrit for Barnes' vectrization    (default: 1024)
//           (ONLY used with GRAPE/HARP inplementation)
// -w        window size for PGPLOT snapshot plot (default: 10)
// -c        flag for collision run
// -x        relative position vec for collision run (no default)
// -v        relative velocity vec for collision run (no default)
// -s        scale factor for position scaling (default: 1)
// -S        scale factor for velocity scaling (default: 1)
//---------------------------------------------------------------------
int main(int argc, char ** argv)
{
    static real_system pb;
    ifstream fsnapin;
    int snapin_flag = 0;
    char fname[255];
    ofstream fsnapout;
    int snapout = 0;
    char foname[255];
    real  dt = 0.015625;
    real  dtsnapout = 1;
    real tend = 10;

    extern char *poptarg;
    int c;
    char* param_string = "i:o:d:D:T:e:t:n:w:cx:v:s:S:h";
    foname[0] = '?';
    foname[1] = '\0';
    real eps = 0.025;
    real theta = 0.75;
    int ncrit = 1024;
    int collision_flag = 0;
    int relx_flag = 0;
    int relv_flag = 0;
    vec relpos, relv;
    real pos_scale = 1;
    real vel_scale = 1;
    pb.plot_xmax = 10;
    while ((c = pgetopt(argc, argv, param_string)) != -1){
        switch(c) {
            case 'i': strcpy(fname,poptarg);
		      snapin_flag = 1;
		      break;
            case 'o': strcpy(foname,poptarg);
		      snapout = 1;
		      break;
            case 'd': dt = atof(poptarg);
		      break;
            case 'D': dtsnapout = atof(poptarg);
		      break;
            case 'T': tend = atof(poptarg);
		      break;
            case 'e': eps = atof(poptarg);
		      break;
            case 't': theta = atof(poptarg);
		      break;
            case 'n': ncrit = atoi(poptarg);
		      break;
            case 'w': pb.plot_xmax = atof(poptarg);
		      break;
            case 'c': collision_flag = 1;
		      break;
            case 'x': relx_flag = 1;
		      relpos = vec(atof(poptarg),
				      atof(poptarg+1),
				      atof(poptarg+2));
		      pskipopt();pskipopt();
		      break;
	    case 'v': relv_flag = 1;
		      relv = vec(atof(poptarg),
				      atof(poptarg+1),
				      atof(poptarg+2));
		      pskipopt();pskipopt();
		      break;
	    case 's': pos_scale = atof(poptarg);
		      break;
	    case 'S': vel_scale = atof(poptarg);
		      break;
            case 'h':		      
		      cerr << "list of options\n";
		      cerr << "-i        name of snapshot input file       (no default)\n";
		      cerr << "-o        name of snapshot output file      (default: no output)\n";
		      cerr << "-d        timestep                          (default: 0.015625)\n";
		      cerr << "-D        time interval for snapshot output (default: 1)\n";
		      cerr << "-T        time to stop integration          (default: 10)\n";
		      cerr << "-e        softening parameter epsilon2      (default: 0.025)\n";
		      cerr << "-t        opening angle theta               (default: 0.75)\n";
		      cerr << "-n        ncrit for Barnes' vectrization    (default: 1024)\n";
		      cerr << "          (ONLY used with GRAPE/HARP inplementation)\n";
		      cerr << "-w        window size for PGPLOT snapshot plot (default: 10)\n";
		      cerr << "-c        flag for collision run\n";
		      cerr << "-x        relative position vec for collision run (no default)\n";
		      cerr << "-v        relative velocity vec for collision run (no default)\n";
		      cerr << "-s        scale factor for position scaling (default: 1)\n";
		      cerr << "-S        scale factor for velocity scaling (default: 1)\n";
		      cerr << "-h        print this help\n";
		      exit(1);
		      break;
		  }
    }
    if (snapin_flag == 0){
	cerr << "Snapshot input file required (-i)\n";
	exit(1);
    }

    PR(dt); PR(dtsnapout); PRL(tend);
    cout << "dt= " << dt
         << " dtsnapout= " << dtsnapout
         << " tend= " <<tend<<endl;
    cerr << "outfile=<" <<foname<<">\n";
    cout << "snapin=" << fname << " snapout=" << foname <<endl;
    fsnapin.open(fname,ios::in);
    pb.atos(fsnapin);
    cout << "n= " << pb.n << endl;
    pb.use_self_gravity = 1;
    fsnapin.close();
    if (snapout) fsnapout.open(foname,ios::out); 
    PR(eps); PR(theta); PR(ncrit);
    cout << "eps= " << eps << " theta=" << theta << " ncrit=" <<ncrit <<endl;
    pb.eps2_for_gravity = eps*eps;
    pb.theta_for_tree = theta;
    pb.ncrit_for_tree = ncrit;
    if(collision_flag == 1){
	if (relx_flag == 0){
	    cerr << "relative position required (x option)\n";
	    exit(1);
	}
	if (relv_flag == 0){
	    cerr << "relative velocity required (v option)\n";
	    exit(1);
	}
	pb.make_collision(relpos, relv);
	cout << "Collision relp= " << relpos << " relv=" << relv <<endl;
    }else{
	pb.apply_vf(&real_particle::scale_pos,pos_scale);
	pb.apply_vf(&real_particle::scale_vel,vel_scale);
	cout << "posscale= " << pos_scale << " velscale=" << vel_scale <<endl;
    }

    cout << endl;
#ifdef GRAPHICS    
    initgraph();
#endif
    int nstep = (int)(tend/dt+0.1);
    dt = tend/nstep;
    int outsnapstep = (int) (dtsnapout/dt+0.1);
    pb.calculate_gravity();

    real E0 = pb.energy();
    PRL(E0);
    cout << "E0 = " <<E0<<endl;
    pb.plot(0.0);
    for(int i=0;i<nstep;i++){
	pb.time = (i+1)*dt;
	pb.integrate(dt);
	real KE = pb.kinetic_energy();
	real E = pb.energy();
	real Eerr = fabs((E0-E)/E0);
	real Q = KE/(KE-E);
	real T = pb.time;
	PR(T); PR(KE); PR(Q); PR(E); PRL(Eerr);
	cout << "T= " << T << " KE= " << KE << " Q= " << Q << " E= " << E << " Eerr= " << Eerr <<endl;
	pb.plot(pb.time);
	pb.calculate_cmterms();
	if(i % outsnapstep == outsnapstep - 1){
	    if (snapout)pb.stoa(fsnapout);
	}
	cout << "CPU sec = " <<cpusec() << endl <<endl;
    }
    if (snapout)fsnapout.close();
    
}


//#endif 

