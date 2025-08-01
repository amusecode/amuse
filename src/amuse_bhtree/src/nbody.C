
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
#include <iostream>
using namespace std;

void real_particle::read(istream& ss)
{
    ss 	>> index ;
    real dum;
    ss >>pos >>vel >> mass>> acc_gravity >> phi_gravity >> dum;
}

void real_system::read(istream& s)
{
    cerr << "enter read: ";PRL(nsize);
    s >> n >> eps2_for_gravity >> use_self_gravity
      >> theta_for_tree >> ncrit_for_tree;
    if (nsize < n){
	delete [] pb;
	pb = NULL;
    }

    if (pb == NULL){
	pb = new real_particle[n];
	nsize = n;
    }
    for(int i = 0; i<n; i++){
	(pb+i)->read(s);
    }
    cerr << "exit  read: ";PRL(nsize);
}

void real_system::atos(istream& s)
{
    cerr << "enter atos: ";PRL(nsize);
    int nddum;
    s>>n >> nddum >> time;
	
    if (nsize < n){
	delete [] pb;
	pb = NULL;
    }

    if (pb == NULL){
	pb = new real_particle[n];
	nsize = n;
    }
    for(int i = 0; i<n; i++){
	real mtmp;
	s>>mtmp;
	(pb+i)->set_mass(mtmp);
	(pb+i)->set_index(i);
	
    }
    for(int i = 0; i<n; i++){
	vec ptmp;
	s>>ptmp;
	(pb+i)->set_pos(ptmp);
    }
    for(int i = 0; i<n; i++){
	vec vtmp;
	s>>vtmp;
	(pb+i)->set_vel(vtmp);
    }
    cerr << "exit  atos: ";PRL(nsize);
}


void real_system::stoa(ostream& ss)
{
    ss.precision(12);
    ss <<n <<endl;
    ss << "3\n";
    ss <<  time << endl;
	
    for(int i = 0; i<n; i++){ss<<(pb+i)->get_mass() <<endl;}
    for(int i = 0; i<n; i++){ss<<(pb+i)->get_pos() <<endl;}
    for(int i = 0; i<n; i++){ss<<(pb+i)->get_vel() <<endl;}
    for(int i = 0; i<n; i++){ss<<(pb+i)->get_phi_gravity() <<endl;}
}



void real_particle::write(ostream& ss)
{
    ss << index << " " ;
    ss << pos << " "<< vel << " "<< mass << " "  <<acc_gravity << " " <<  phi_gravity <<endl;
}

void real_system::write(ostream& s)
{
    s<<n << " " << eps2_for_gravity <<" " << use_self_gravity
     << " " << theta_for_tree << " " <<ncrit_for_tree <<endl;;
    for(int i = 0; i<n; i++){
	(pb+i)->write(s);
    }
}



void real_particle::dump()
{
    PRL(index);
    PRC(pos);  PRC(vel);    PRL(acc_gravity);
    PRC(mass);  PRL(phi_gravity);
}

void real_system::dump()
{
    PRL(n);
    for(int i = 0; i<n; i++){
	(pb+i)->dump();
    }
}

const real piinv = 1.0/M_PI;
real real_particle::kinetic_energy()
{
    return (0.5*vel*vel)*mass;
}

real real_particle::energy()
{
    return kinetic_energy();
}


real real_system::kinetic_energy()
{
    real sum = 0;
    for(int i = 0; i<n; i++){
	sum+= (pb+i)->kinetic_energy();
    }
    return sum;
}




real real_system::energy()
{
    real sum = kinetic_energy();
    if(use_self_gravity){
	for(int i = 0; i<n; i++){
	    sum+= (pb+i)->get_phi_gravity()*(pb+i)->get_mass()
#ifdef REAL_GRAVITY		
		*0.5;
#else
	    ;
#endif	    
	}
    }
    return sum;
}


void real_system::apply_vf(real_VF_ptr f)
{
    for(int i = 0;i<n;i++){
	(pb[i].*f)();
    }
}
	
void real_system::apply_vf(real_RF_ptr f, real r)
{
    for(int i = 0;i<n;i++){
	(pb[i].*f)(r);
    }
}
	
void real_system::apply_vf(real_RRF_ptr f, real r, real r2)
{
    for(int i = 0;i<n;i++){
	(pb[i].*f)(r, r2);
    }
}

void real_system::integrate(real dt)
{
    apply_vf(&real_particle::predict,dt);
    calculate_gravity();
    apply_vf(&real_particle::correct,dt);
}    
    
const vec uniform_random_position_in_sphere(double power_index)
{
    vec x;
    do{
	//for(int i = 0; i<3;i++) x[i] = drand48()*2-1;
    }while (x*x >= 1);
    x *=  pow(x*x, 3.0/(power_index+3)-1);
    return x;
}

void real_system::calculate_cmterms()
{
    pos = 0.0;
    vel = 0.0;
    mass = 0.0;
    real_particle * p = pb;
    for(int i = 0; i<n;i++){
	pos += p->get_mass()*p->get_pos();
	vel += p->get_mass()*p->get_vel();
	mass += p->get_mass();
	p++;
    }

    vec CMpos = pos/mass, CMvel = vel/mass;
    // PRL(CMpos); PRL(CMvel);
}
    

void real_system::create_uniform_sphere(int nbody, real power_index, real r0)
{
    PRC(nbody); PRL(power_index);
    n = nbody;
    pb = new real_particle[n];
    nsize = n;
    real_particle * p = pb;
    for(int i = 0; i<n;i++){
	p->set_pos(uniform_random_position_in_sphere(power_index)*r0);
	p->set_vel(0);
	p->set_mass(1.0/n);
	p->set_index(i);
	p++;
    }
}



void copy_nbody_particles(real_system * source,
			real_system * destination)
{

    destination->pb = new real_particle[source->nsize];
    real_particle * pbs = source->pb;
    real_particle * pbd = destination->pb;
    for(int i = 0; i<source->n;i++){
	*pbd = *pbs; 
	pbd++;
	pbs++;
    }
}

int compare_particle_pointer(real_particle * * p1, real_particle * *p2)
{
    int i1 = (*p1)->get_index();
    int i2 = (*p2)->get_index();
    if (i1 > i2){
	return 1;
    }else if (i1 == i2){
	return 0;
    }else{
	return -1;
    }
}

void initgraph()
{
#ifdef GRAPHICS    
    if(cpgopen("") != 1) exit(EXIT_FAILURE);
    cpgask(0);
    cpgvsiz(1.0,6.0,1.0,6.0);
#endif    
}


void real_particle::plot(real parm)
{
#ifdef GRAPHICS    
    int iparm = (int)parm;
    cpgpt1(pos[0],pos[1],-1);
#endif
}
void real_system::plot(real t)
{
#ifdef GRAPHICS    
    real xmax = plot_xmax;
    cpgbbuf();
    cpgsci(1);
    cpgslw(2);
    cpgenv(-xmax, xmax, -xmax, xmax,  1, 0);
    char label[255];
    sprintf(label, "T = %g N=%d", t,n);
    cpglab("x","y", label);
    cpgsci(1);
    int psize = 4;
    if (n < 1000){
	psize *= (int)pow(1000.0/n,0.25);
    }
    cpgslw(psize);
    apply_vf(&real_particle::plot, 0.0);
    cpgebuf();
#endif
}

void real_system::make_collision(vec relpos, vec relv)
{
    real_particle *new_pb;
    new_pb = new real_particle[n*2];
    int i;
    for( i=0; i<n;i++){
	*(new_pb+i) = *(pb+i);
	*(new_pb+i+n) = *(pb+i);
    }
    delete [] pb;
    pb = new_pb;
    for( i=0; i<n;i++){
	(pb+i)->inc_pos(0.5*relpos);
	(pb+i)->inc_vel(0.5*relv);
	(pb+i)->set_index(i);
    }
    for( i=n; i<n*2;i++){
	(pb+i)->inc_pos(-0.5*relpos);
	(pb+i)->inc_vel(-0.5*relv);
	(pb+i)->set_index(i);
    }
    n *= 2;
    nsize *= 2;
    PR(n); PRL(nsize);
}


#if 0
#include<fstream>

//#ifndef TOOLBOX
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



#endif 
//#endif 

