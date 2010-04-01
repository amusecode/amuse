#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include "octgrav.h"

// shortcuts for printing parameters
#define PR(x)  cout << #x << " = " << x << " "
#define PRC(x) cout << #x << " = " << x << ",  "
#define PRL(x) cout << #x << " = " << x << "\n"

using namespace std;

struct momvar
{
  float x,y,z,abs,i;
};


vector<float4> pos;
vector<float4> acc;
vector<float3> vel,vel1;


// global variables with defaults
char  inp_file[255];  // input file name
float tend  = 10.;    // end time of simulation
float dtmax = 0.125;  // maximum time step (also used in case of fixed ts
float tnow;           // current time (read from input snapshot)
int   idt   = 0;      // fixed (=0) or variable (=1) ts
float eps=0.05,eps2=eps*eps; // softening
int   N;              // number of particles (read from input snapshot)
float dtlog = 0.125;
float dtout = 1.0;
float theta = 0.5;
float ekin,epot,etot,etot0,etot1,detot,detot1;
float tlog,tout;
int    isnapshot=0;
int    irestart=0,nrestart=10;
int    nsteps=0;
float3 pcom,vcom;
float  mtot;
momvar amom,amom0;


int init_parameters(int argc, char **argv)
{
  int c;
  int power;
  float dthelp;
  char* param_string = "i:o:d:D:l:T:e:E:t:n:w:s:S:r::R:h";
  while ((c = getopt(argc, argv, param_string)) != -1)
    {
      switch(c) 
	{
	case 'i':
	  strcpy(inp_file, optarg);
	  break;
	case 't': 
	  tend = atof(optarg);
	  break;
	case 'd':
	  dthelp = atof(optarg);
	  if (dthelp < 0) 
	    {
	      power = (int) dthelp;
	    }
	  else
	    {
	      power = (int)(log(dthelp)/log(2.));
	    }
	  dtmax = pow(2.,(float)power);
	  break;
	case 'l':
	  dthelp = atof(optarg);
	  if (dthelp < 0) 
	    {
	      power = (int) dthelp;
	    }
	  else
	    {
	      power = (int)(log(dthelp)/log(2.));
	    }
	  dtlog = pow(2.,(float)power);
	  break;
	case 'D':
	  dthelp = atof(optarg);
	  if (dthelp < 0) 
	    {
	      power = (int) dthelp;
	    }
	  else
	    {
	      power = (int)(log(dthelp)/log(2.));
	    }
	  dtout = pow(2.,(float)power);
	  break;
	case 'e':
	  eps  = atof(optarg);
	  eps2 = eps*eps;
	  break;
	case 'T':
	  theta = atof(optarg);
	  break;
	case 'R':
          nrestart = atoi(optarg);
	  break;
	case 'r':
	  irestart = 1;
	  break;
	case 'h':
	  cerr << "nbint:   simple n-body integrator" << endl;
	  cerr << endl;
	  cerr << "usage:   > ./nbint -i <input-file> [optional parameters]" 
	       << endl;
	  cerr << endl;
	  cerr << "options: -i <input-file>        give name of input file []*"
	       << endl;
	  cerr << "         -t                     t_end [10.]" << endl;
	  cerr << "         -d                     dtmax [0.125]" << endl;
	  cerr << "         -e                     softening [0.0]" << endl;
          cerr << "         -T                     opening angle theta [0.5]" << endl;
	  cerr << "         -l                     log output interval [0.125]" << endl;
	  cerr << "         -D                     snapshot interval [1.0]" << endl;
	  cerr << "         -R                     write restart file 'nbint.rst'" << endl;
	  cerr << "                                every n'th time step [10]" << endl;
	  cerr << "         -r                     indicate input file is restart file" << endl;
	  cerr << "         -h                     print this help" << endl;
          cerr << endl;
	  cerr << "         [] indicates the default * indicates mandatory arg"
	       << endl;
	  exit(1);
	}
    }
}

float set_dt()
{
  if (idt == 0) 
    {
      return dtmax;
    }
  else
    {
      cerr << "variable time step not yet implemented, using dtmax" << endl;
      return dtmax;
    }
}

int read_snapshot(char *fname)
{
  ifstream infile;
  int ndummy,i;

  infile.open(fname);

  if (!infile)
    {
      cerr << "Unable to open input file!" << endl;
      exit(1);
    }


  // reading header
  infile >> N;
  pos.resize(N);    // resize array to deal with actual number of particles
  vel.resize(N);
  vel1.resize(N);
  acc.resize(N);

  infile >> ndummy;
  infile >> tnow;
  
  // reading data
  for(i=0;i<N;i++)
    {
      infile >> pos[i].w;
    }
  for(i=0;i<N;i++)
    {
      infile >> pos[i].x >> pos[i].y >> pos[i].z;
    }
  for(i=0;i<N;i++)
    {
      infile >> vel[i].x >> vel[i].y >> vel[i].z;
    }
 
  cout << "Number of particles read: "; PRL(N);
  
  // debugging info
  /*  float mm;
  body  p;
  for(i=0;i<N;i++)
    {
      mm = m[i];
      //      x = pos[i].x; y = pos[i].y; z = pos[i].z;
      p = pos[i];
      cerr << "Particle: "; PR(i); PR(mm); PR(p.x); PR(p.y); PRL(p.z);
      }*/
  // end debug info

}

int write_snapshot(int iss)
{
  ofstream outfile;
  int i;
  char fname[255];

  // make file name
  sprintf(fname,"%10.10i.dat",iss);

  
  outfile.open(fname);

  // witing header
  outfile << N << endl;;
  outfile << "3" << endl;
  outfile << tnow << endl;
  
  // reading data
  for(i=0;i<N;i++)
    {
      outfile << pos[i].w << endl;
    }
  for(i=0;i<N;i++)
    {
      outfile << "  " <<  pos[i].x << "  " << pos[i].y << "  " << pos[i].z << endl;
    }
  for(i=0;i<N;i++)
    {
      outfile << "  " << vel[i].x << "  " << vel[i].y << "  " << vel[i].z << endl;
    }

  outfile.close();
  cout << "written snapshot at "; PR(tnow);
  cout << " to file "; PRL(fname);
  


}

int write_snapshot_fast(int iss)
{
  FILE *outfile;
  int i;
  char fname[255];

  // make file name
  sprintf(fname,"%10.10i.dat",iss);

  outfile = fopen(fname, "w");

  // witing header
  fprintf(outfile, "%d \n", N);
  fprintf(outfile, "3 \n");
  fprintf(outfile, "%g \n", tnow);
  
  
  // reading data
  for(i=0;i<N;i++)
    {
      fprintf(outfile, "%g \n", pos[i].w);
    }
  for(i=0;i<N;i++)
    {
      fprintf(outfile, "%g %g %g\n", pos[i].x, pos[i].y, pos[i].z);
    }
  for(i=0;i<N;i++)
    {
      fprintf(outfile, "%g %g %g\n", vel[i].x, vel[i].y, vel[i].z);
    }
  
  fclose(outfile);
  cout << "written snapshot at "; PR(tnow);
  cout << " to file "; PRL(fname);
  
}

int write_restart_file()
{
  FILE *outfile;
  int i;
  char fname[255];

  // use fixed file name
  sprintf(fname,"nbint.rst");

  outfile = fopen(fname, "w");

  // witing header
  fwrite(&N,sizeof(int),1,outfile);
  fwrite(&isnapshot,sizeof(int),1,outfile);
  fwrite(&nsteps,sizeof(int),1,outfile);
  fwrite(&nrestart,sizeof(int),1,outfile);
  fwrite(&idt,sizeof(int),1,outfile);
  fwrite(&tnow,sizeof(float),1,outfile);
  fwrite(&tend,sizeof(float),1,outfile);
  fwrite(&dtmax,sizeof(float),1,outfile);
  fwrite(&dtlog,sizeof(float),1,outfile);
  fwrite(&tlog,sizeof(float),1,outfile);
  fwrite(&dtout,sizeof(float),1,outfile);
  fwrite(&tout,sizeof(float),1,outfile);
  fwrite(&theta,sizeof(float),1,outfile);
  fwrite(&eps,sizeof(float),1,outfile);
  fwrite(&etot0,sizeof(float),1,outfile);
  

  // write particle data
  fwrite(&pos[0],sizeof(float4),N,outfile);  // mass and pos 
  fwrite(&vel[0],sizeof(float3),N,outfile);  // velocities
  fwrite(&acc[0],sizeof(float4),N,outfile);  // acc and pot

  fclose(outfile);
}

int read_restart_file(char *fname)
{
  FILE *infile;
  int i;

  // open restart file
  infile = fopen(fname, "r");

  if (!infile)
    {
      cerr << "Unable to open restart file!" << endl;
      exit(1);
    }
  
  // read header
  fread(&N,sizeof(int),1,infile);
  fread(&isnapshot,sizeof(int),1,infile);
  fread(&nsteps,sizeof(int),1,infile);
  fread(&nrestart,sizeof(int),1,infile);
  fread(&idt,sizeof(int),1,infile);
  fread(&tnow,sizeof(float),1,infile);
  fread(&tend,sizeof(float),1,infile);
  fread(&dtmax,sizeof(float),1,infile);
  fread(&dtlog,sizeof(float),1,infile);
  fread(&tlog,sizeof(float),1,infile);
  fread(&dtout,sizeof(float),1,infile);
  fread(&tout,sizeof(float),1,infile);
  fread(&theta,sizeof(float),1,infile);
  fread(&eps,sizeof(float),1,infile);
  fread(&etot0,sizeof(float),1,infile);

  // set some more vars
  eps2 = eps*eps;

  // resize array to deal with actual number of particles
  pos.resize(N);    
  vel.resize(N);
  vel1.resize(N);
  acc.resize(N);
  // read particle data
  fread(&pos[0],sizeof(float4),N,infile);  // mass and pos 
  fread(&vel[0],sizeof(float3),N,infile);  // velocities
  fread(&acc[0],sizeof(float4),N,infile);  // acc and pot

  fclose(infile);
  cout << "number of particles read from restart file ";
  PRL(N);

}


int gravity()
{
  int    i,j;
  float r2,r2inv,rinv,r3inv,dx,dy,dz;

  // reset  acc and pot
  for(i=0;i<N;i++)
    {
      acc[i].x = 0.0;
      acc[i].y = 0.0;
      acc[i].z = 0.0;
      acc[i].w = 0.0;
    }


  for(i=0;i<N-1;i++)
    {
      for(j=i+1;j<N;j++)
	{
	  dx = pos[j].x - pos[i].x;
	  dy = pos[j].y - pos[i].y;
	  dz = pos[j].z - pos[i].z;
	  r2 = eps2 + dx*dx + dy*dy + dz*dz;	  
	  r2inv = 1.0/r2;
	  rinv  = sqrt(r2inv);
	  r3inv = r2inv*rinv;
	  acc[i].x += pos[j].w*r3inv*dx; 
	  acc[i].y += pos[j].w*r3inv*dy; 
	  acc[i].z += pos[j].w*r3inv*dz; 
	  acc[j].x += -pos[i].w*r3inv*dx; 
	  acc[j].y += -pos[i].w*r3inv*dy; 
	  acc[j].z += -pos[i].w*r3inv*dz; 
	  acc[i].w += -pos[j].w*rinv;
	  acc[j].w += -pos[i].w*rinv;
	}
    }
}

int energy()
{
  int   i;
  float v2;
  ekin = 0.0;
  epot = 0.0;

  for(i=0;i<N;i++)
    {
      v2   = (vel[i].x*vel[i].x) + (vel[i].y*vel[i].y) + (vel[i].z*vel[i].z);
      ekin += pos[i].w*v2;
      epot += pos[i].w*acc[i].w;
//       fprintf(stderr, "i= %d  pot= %g [%g %g %g] \n",
// 	      i, acc[i].w, 
// 	      acc[i].x,
// 	      acc[i].y,
// 	      acc[i].z);
    }
  ekin *= 0.5;
  epot *= 0.5;
  etot  = ekin + epot; 
}


int center_of_mass()
{
  int i;

  mtot   = 0.;
  pcom.x = 0.;
  pcom.y = 0.;
  pcom.z = 0.;
  vcom.x = 0.;
  vcom.y = 0.;
  vcom.z = 0.;

  for(i=0;i<N;i++)
    {
      mtot   += pos[i].w;
      pcom.x += pos[i].w * pos[i].x;
      pcom.y += pos[i].w * pos[i].y;
      pcom.z += pos[i].w * pos[i].z;
      vcom.x += pos[i].w * vel[i].x;
      vcom.y += pos[i].w * vel[i].y;
      vcom.z += pos[i].w * vel[i].z;
    }

  pcom.x /= mtot;
  pcom.y /= mtot;
  pcom.z /= mtot;
  vcom.x /= mtot;
  vcom.y /= mtot;
  vcom.z /= mtot;


}

int momentum()
{
  int   i;
  float xy;

  amom.x = 0.;
  amom.y = 0.;
  amom.z = 0.;

  for(i=0;i<N;i++)
    {
      amom.x += pos[i].w * (pos[i].y*vel[i].z - pos[i].z*vel[i].y);
      amom.y += pos[i].w * (pos[i].z*vel[i].x - pos[i].x*vel[i].z);
      amom.z += pos[i].w * (pos[i].x*vel[i].y - pos[i].y*vel[i].x);
    }

  amom.abs = sqrt(amom.x*amom.x + amom.y*amom.y + amom.z*amom.z);
  xy       = sqrt(amom.x*amom.x + amom.y*amom.y);
  amom.i   = atan2(xy,amom.z+1.e-30);

}

int step_pos(float dt)
{
  int i;
  for(i=0;i<N;i++)
    {
      pos[i].x += vel1[i].x*dt;
      pos[i].y += vel1[i].y*dt;
      pos[i].z += vel1[i].z*dt;
    }
}

int step_vel(float dt)
{
  int i;
  for(i=0;i<N;i++)
    {
      // advance from staggered vel half a step to get unstaggered vel
      vel[i].x = vel1[i].x + acc[i].x*dt*0.5;
      vel[i].y = vel1[i].y + acc[i].y*dt*0.5;
      vel[i].z = vel1[i].z + acc[i].z*dt*0.5;

      // advance staggered vel by full step
      vel1[i].x += acc[i].x*dt;
      vel1[i].y += acc[i].y*dt;
      vel1[i].z += acc[i].z*dt;
    }
}


//======================================================================
//
// nbint: simple integrator to be linked with GPU-TreeCode octgrav
//
//======================================================================
int main(int argc, char *argv[])
{
  //cout << "running nbint" << endl;

  // some vars
  float destep=0.;
  
  // get input parameters from command line
  init_parameters(argc,argv);

  // write log to stdout
  PR(tend);PRL(dtmax);
  PR(irestart);PRL(inp_file);

  float dt,dt05,tnext;

  static octgrav system;

  

  if (irestart==0) 
    {
      // load initial conditions
      read_snapshot(inp_file);

      // initialize system
      system.set_softening(eps);
      system.set_opening_angle(theta);
      system.evaluate_gravity(pos, acc);
      // gravity();        // use this for force_on_host
      energy();
      etot0 = etot;
      etot1 = etot;
      detot = fabs(etot-etot0)/fabs(etot0);
      PR(tnow);PR(etot);PR(epot);PR(ekin);PRL(detot);
      center_of_mass();
      printf("pcom = %e %e %e  ",pcom.x,pcom.y,pcom.z);
      printf("vcom = %e %e %e\n",vcom.x,vcom.y,vcom.z);
      momentum();
      printf("amom = %e %e %e  amom.abs = %e  amom.i = %e\n",amom.x,amom.y,amom.z,amom.abs,amom.i); 

      // set output times
      tlog = (float)((int)(tnow/dtlog))*dtlog + dtlog;
      tout = (float)((int)(tnow/dtout))*dtout + dtout;
    }
  else
    {
      // all parameters are read from binary restart file
      read_restart_file(inp_file);
      system.set_softening(eps);
      system.set_opening_angle(theta);
    }


  // set time step
  dt   = set_dt();
  dt05 = 0.5*dt;

  // initialize staggered velocities
  int i;
  for(i=0;i<N;i++)
    {
      vel1[i].x = vel[i].x + dt05*acc[i].x;
      vel1[i].y = vel[i].y + dt05*acc[i].y;
      vel1[i].z = vel[i].z + dt05*acc[i].z;
    }


  // main loop begins here
  while (tnow<tend)
    {
      tnext = tnow + dt;

      // cout << "integrating to ";PRL(tnext);

      // leap-frog integration
      step_pos(dt);              // step postions by a step
      double tt1 = get_time();
      system.evaluate_gravity(pos, acc); // update force
      double tt2 = get_time();
      fprintf(stderr, "total_time= %lf sec\n", tt2 - tt1);
      // gravity();        // use this for force_on_host
      step_vel(dt);              // step staggerd vel by a full step
                                 // also update un-staggered vel 

      tnow  = tnext;
      nsteps++;

      // write restart file every n'th step
      if (nsteps%nrestart == 0)
	{
	  cerr << "writing restart file ...";
	  double t1  = get_time();
	  write_restart_file();
	  double t2 = get_time();
	  cerr << " ... done in " << t2 -  t1 << " sec " << endl;
	}

      // output 
      if (tnow>=tlog) 
	{
	  energy();
	  detot = fabs(etot-etot0)/fabs(etot0);
	  PR(tnow);PR(etot);PR(epot);PR(ekin);PRL(detot);
	  center_of_mass();
	  printf("pcom = %e %e %e  ",pcom.x,pcom.y,pcom.z);
	  printf("vcom = %e %e %e\n",vcom.x,vcom.y,vcom.z);
	  momentum();
	  printf("amom = %e %e %e  amom.abs = %e  amom.i = %e\n",amom.x,amom.y,amom.z,amom.abs,amom.i);  
	  cout.flush();
	  tlog += dtlog;
	}

      // write snapshot
       if (tnow>=tout) 
	{
	  isnapshot++;
	  write_snapshot_fast(isnapshot);
	  tout += dtout;
	}
    

    }
  // end of main loop

  destep = destep/(float)nsteps;
  cout << "end of run:   "; PR(dt);PR(nsteps);PRL(detot);

  return 0;

}
