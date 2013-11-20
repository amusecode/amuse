#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "forcetree.h"




#define KERN_LEN   10000


static struct NODE 
{ double center[3],len;                /* center and sidelength of treecubes */
  double mass,oc;                      /* mass and variable for opening criter*/
  double s[3];                         /* center of mass */
  double Q11,Q22,Q33,Q12,Q13,Q23,P;      /* quadrupol tensor */  
  struct NODE *next,*sibling,*father,*suns[8];
  int    partind;
  int    cost;
  int    count;          /* total # of particles in cell and below   */
} *nodes;




static struct NODE *last;
static struct NODE *trees[5];       /* we construct one tree for each particle type */
static int    ntype[5];             /* holds the number of particles of each type */
static int    treecost[5];          /* sums the cell-particle interactions  (diagnostic) */
static int    treecost_quadru[5];   /* sums the cell-particle interactions  (diagnostic) */

static int    numnodestotal;        /* total number of nodes */
static int    numnodestree[5];     /* number of nodes of each particle type */
static int    numnodes,MaxNodes;



static double xmin[3],xmax[3],len;




static struct COORDINATES
{
  double xyz[3];
  double mass;
  double accel[3];
  double pot;
  int    type;
} **Part;

static int    N;






struct xyz_data 
{ 
  double xyz[3]; 
}; 



static double *softTable;  /* points to softenings for different part. types */







static double  knlrad  [KERN_LEN+1],
               knlforce[KERN_LEN+1],
               knlpot  [KERN_LEN+1],
               knlW2   [KERN_LEN+1],
               knlW3   [KERN_LEN+1],
               knlW4   [KERN_LEN+1];










void add_particle_props_to_node(struct NODE *no,struct COORDINATES *pa)
{
  int i;

  for(i=0;i<3;i++)
    no->s[i] += pa->mass*pa->xyz[i];

  no->mass += pa->mass;

  no->Q11+=pa->mass*pa->xyz[0]*pa->xyz[0];
  no->Q22+=pa->mass*pa->xyz[1]*pa->xyz[1];
  no->Q33+=pa->mass*pa->xyz[2]*pa->xyz[2];
  no->Q12+=pa->mass*pa->xyz[0]*pa->xyz[1];
  no->Q13+=pa->mass*pa->xyz[0]*pa->xyz[2];
  no->Q23+=pa->mass*pa->xyz[1]*pa->xyz[2];

  no->P+= pa->mass*(pa->xyz[0]*pa->xyz[0] + pa->xyz[1]*pa->xyz[1] + pa->xyz[2]*pa->xyz[2]);
}







int force_treebuild(void **pospointer, int npart, double thetamax,int costresetflag)
/* packs the particles 0...Npart-1 into BH-trees, for each particle type 
   (different softening lengths) a different tree */
{
  int i,j,tr,n;
  struct NODE *nfree,*th;
  struct NODE *force_treebuild_single(struct NODE *startnode, int type, double thetamax);


  Part=(struct COORDINATES **)pospointer;
  N=npart;


  for(tr=0;tr<5;tr++)
    ntype[tr]=0;

  for(i=0;i<npart;i++)
    {
      ntype[ Part[i]->type ]++;

    }
  

  for(tr=0,nfree=nodes,numnodestotal=0; tr<5; tr++)
    {
      treecost[tr]=0;
      treecost_quadru[tr]=0;
      if(ntype[tr]>0)
	{
	  trees[tr]=nfree;
	  nfree = force_treebuild_single(trees[tr], tr, thetamax);
	}
      else
	{
	  trees[tr]=0;
	  numnodestree[tr]=0;
	}
    }

  for(tr=0;tr<5;tr++) 
    { 
      if(numnodestree[tr]) 
 	{ 
 	  for(i=0,th=trees[tr]; i<numnodestree[tr]; i++,th++) 
	    th->cost=0; 
 	} 
    } 

  return numnodestotal;
}





struct NODE *force_treebuild_single(struct NODE *startnode, int type, double thetamax)
/* packs the particles of the correct type in a BH-tree */
{
  int i,j,subp,subi,p,ni,subnode,fak,fp;
  double x;
  double dx,dy,dz;
  struct NODE *nfree,*th,*nn,*ff;
  void add_particle_props_to_node(struct NODE *no,struct COORDINATES *pa);
  void force_setupnonrecursive(struct NODE *no);


 

  for(fp=0;fp<N;fp++)          /* find first particle of this type */
    if(Part[fp]->type == type)
      break;


  for(j=0;j<3;j++)                        /* find enclosing rectangle */
    xmin[j]=xmax[j]=Part[fp]->xyz[j];

  for(i=fp+1;i<N;i++)
    {
      if(Part[i]->type == type)
	{
	  for(j=0;j<3;j++)
	    {
	      if(Part[i]->xyz[j]>xmax[j]) 
		xmax[j]=Part[i]->xyz[j];
	      if(Part[i]->xyz[j]<xmin[j]) 
		xmin[j]=Part[i]->xyz[j];
	    }
	}
    }


  for(j=1,len=xmax[0]-xmin[0];j<3;j++)  /* determine maxmimum externsion */
    if((xmax[j]-xmin[j])>len)
      len=xmax[j]-xmin[j];
  
  len*=1.0001;





  /* insert first particle of correct type in root node */
 
  nfree=startnode;
 
  for(j=0;j<3;j++)
    nfree->center[j]=(xmax[j]+xmin[j])/2;
  nfree->len=len;

  nfree->father=0;
  for(i=0;i<8;i++)
    nfree->suns[i]=0;
  nfree->partind=fp;

  nfree->mass=Part[fp]->mass;

  for(i=0;i<3;i++)
    nfree->s[i]= Part[fp]->mass*Part[fp]->xyz[i];
  
  nfree->Q11=Part[fp]->mass*Part[fp]->xyz[0]*Part[fp]->xyz[0];
  nfree->Q22=Part[fp]->mass*Part[fp]->xyz[1]*Part[fp]->xyz[1];
  nfree->Q33=Part[fp]->mass*Part[fp]->xyz[2]*Part[fp]->xyz[2];
  nfree->Q12=Part[fp]->mass*Part[fp]->xyz[0]*Part[fp]->xyz[1];
  nfree->Q13=Part[fp]->mass*Part[fp]->xyz[0]*Part[fp]->xyz[2];
  nfree->Q23=Part[fp]->mass*Part[fp]->xyz[1]*Part[fp]->xyz[2];
  nfree->P =  Part[fp]->mass*(Part[fp]->xyz[0]*Part[fp]->xyz[0]+
			     Part[fp]->xyz[1]*Part[fp]->xyz[1]+
			     Part[fp]->xyz[2]*Part[fp]->xyz[2]);
  nfree->sibling=0;
  nfree->count=1;

  numnodes=1;  numnodestotal++; nfree++;
  
  if(numnodestotal>=MaxNodes)
    {
      printf("maximum number %d of tree-nodes reached.\n",numnodestotal);
      exit(1);
    }





  /* insert all other particles */

  for(i=fp+1;i<N;i++)  
    {
      if(Part[i]->type != type)
	continue;

      th=startnode;
      
      while(1)
	{
	  th->count++;

	  add_particle_props_to_node(th,Part[i]);

	  if(th->partind>=0)
	    break;
	  
	  for(j=0,subnode=0,fak=1;j<3;j++,fak<<=1)
	    if(Part[i]->xyz[j]>th->center[j])
	      subnode+=fak;
	  
	  if(nn=th->suns[subnode])
	    th=nn;
	  else
	    break;
	}

      
      if(th->partind>=0)  /* cell is occcupied with one particle */
	{
	  while(1)
	    {
	      p=th->partind;

	      for(j=0,subp=0,fak=1;j<3;j++,fak<<=1)
		if(Part[p]->xyz[j]>th->center[j])
		  subp+=fak;

	      
	      nfree->father=th;
	      
	      for(j=0;j<8;j++)
		nfree->suns[j]=0;
	      nfree->sibling=0;
	      
	      nfree->len=th->len/2;
    
	      for(j=0;j<3;j++)
		nfree->center[j]=th->center[j];

	      for(j=0;j<3;j++)
		if(Part[p]->xyz[j]>nfree->center[j])
		  nfree->center[j]+=nfree->len/2;
		else
		  nfree->center[j]-=nfree->len/2;

	      nfree->partind=p;
	      nfree->count=1;

	      nfree->mass=Part[p]->mass;
	      for(j=0;j<3;j++)
		nfree->s[j]=Part[p]->mass*Part[p]->xyz[j];
	      
	      nfree->Q11=Part[p]->mass*Part[p]->xyz[0]*Part[p]->xyz[0];
	      nfree->Q22=Part[p]->mass*Part[p]->xyz[1]*Part[p]->xyz[1];
	      nfree->Q33=Part[p]->mass*Part[p]->xyz[2]*Part[p]->xyz[2];
	      nfree->Q12=Part[p]->mass*Part[p]->xyz[0]*Part[p]->xyz[1];
	      nfree->Q13=Part[p]->mass*Part[p]->xyz[0]*Part[p]->xyz[2];
	      nfree->Q23=Part[p]->mass*Part[p]->xyz[1]*Part[p]->xyz[2];
	      nfree->P  =Part[p]->mass*(Part[p]->xyz[0]*Part[p]->xyz[0]+
				       Part[p]->xyz[1]*Part[p]->xyz[1]+
				       Part[p]->xyz[2]*Part[p]->xyz[2]);
	  
	      th->partind=-1;
	      th->suns[subp]=nfree;
      
	      numnodes++; numnodestotal++; nfree++;

	      if(numnodestotal>=MaxNodes)
		{
		  printf("maximum number %d of tree-nodes reached.\n",numnodestotal);
		  exit(1);
		}

	      for(j=0,subi=0,fak=1;j<3;j++,fak<<=1)
		if(Part[i]->xyz[j]>th->center[j])
		  subi+=fak;
	      
	      if(subi==subp)   /* the new particle lies in the same sub-cube */
		{
		  th=nfree-1;
		  add_particle_props_to_node(th,Part[i]);		  
		  th->count++;
		}
	      else
		break;
	    }
	}
      

      
      for(j=0,subi=0,fak=1;j<3;j++,fak<<=1)
	if(Part[i]->xyz[j]>th->center[j])
	  subi+=fak;
      
      nfree->father=th;

      for(j=0;j<8;j++)
	nfree->suns[j]=0;
      nfree->sibling=0;

      nfree->len=th->len/2;
      for(j=0;j<3;j++)
	nfree->center[j]=th->center[j];

      for(j=0;j<3;j++)
	if(Part[i]->xyz[j]>nfree->center[j])
	  nfree->center[j]+=nfree->len/2;
	else
	  nfree->center[j]-=nfree->len/2;

      nfree->mass=Part[i]->mass;
      for(j=0;j<3;j++)
	nfree->s[j]=Part[i]->mass*Part[i]->xyz[j];

      nfree->Q11=Part[i]->mass*Part[i]->xyz[0]*Part[i]->xyz[0];
      nfree->Q22=Part[i]->mass*Part[i]->xyz[1]*Part[i]->xyz[1];
      nfree->Q33=Part[i]->mass*Part[i]->xyz[2]*Part[i]->xyz[2];
      nfree->Q12=Part[i]->mass*Part[i]->xyz[0]*Part[i]->xyz[1];
      nfree->Q13=Part[i]->mass*Part[i]->xyz[0]*Part[i]->xyz[2];
      nfree->Q23=Part[i]->mass*Part[i]->xyz[1]*Part[i]->xyz[2];
      nfree->P = Part[i]->mass*(Part[i]->xyz[0]*Part[i]->xyz[0]+
				Part[i]->xyz[1]*Part[i]->xyz[1]+
				Part[i]->xyz[2]*Part[i]->xyz[2]);

      nfree->partind=i;
      th->suns[subi]=nfree;
      
      nfree->count=1;

      numnodes++; numnodestotal++; nfree++;

      if(numnodestotal>=MaxNodes)
	{
	  printf("maximum number %d of tree-nodes reached.\n",numnodestotal);
	  exit(1);
	}
    }



 
  
  /* now finish-up center-of-mass and quadrupol computation */
  
  for(i=0,th=startnode; i<numnodes; i++,th++)
    {
      for(j=0;j<3;j++)
	th->s[j] /= th->mass;
      
      if(th->partind<0)   /* cell contains more than one particle */
	{
	  th->Q11 -= th->mass*th->s[0]*th->s[0];
	  th->Q22 -= th->mass*th->s[1]*th->s[1];
	  th->Q33 -= th->mass*th->s[2]*th->s[2];
	  th->Q12 -= th->mass*th->s[0]*th->s[1];
	  th->Q23 -= th->mass*th->s[1]*th->s[2];
	  th->Q13 -= th->mass*th->s[0]*th->s[2];

	  th->P   -= th->mass*(th->s[0]*th->s[0]+
			       th->s[1]*th->s[1]+
			       th->s[2]*th->s[2]); 
	  dx=th->s[0] - th->center[0];
	  dy=th->s[1] - th->center[1];
	  dz=th->s[2] - th->center[2];
	  
	  th->oc=sqrt(dx*dx+dy*dy+dz*dz);
	  th->oc += th->len/thetamax; 
	  th->oc *= th->oc;     /* used in cell-opening criterion */
	}
      
      for(j=7,nn=0;j>=0;j--)  	/* preparations for non-recursive walk */
	{
	  if(th->suns[j])
	    {
	      th->suns[j]->sibling=nn;
	      nn=th->suns[j];
	    }
	}
    }

  
  last=0;
  force_setupnonrecursive(startnode);  	/* set up non-recursive walk */
  last->next=0;

  
  
  for(i=0,th=startnode; i<numnodes; i++,th++)
    if(!(th->sibling))
      {
	ff=th;
	nn=ff->sibling;

	while(!nn)
	  {
	    ff=ff->father;
	    if(!ff)
	      break;
	    nn=ff->sibling;
	  }
	
	th->sibling=nn;
      }

	

/*   printf("Tree contruction finished. Number of nodes: %d\n",numnodes); */
  
  numnodestree[type]=numnodes;

  return nfree;
}



void force_setupnonrecursive(struct NODE *no)
{
  int i;
  struct NODE *nn;
  
  if(last)
    last->next=no;

  last=no;
  
  for(i=0;i<8;i++)
    if(nn=no->suns[i])
      force_setupnonrecursive(nn);
}
 







void force_treeevaluate(int target) 
{
  int    i,tr;
  double epsilon;
  void force_treeevaluate_single(int tree, int targetpart, double epsilon);


  Part[target]->pot=0;


  for(i=0;i<3;i++)
    Part[target]->accel[i]=0;
  
  for(tr=0;tr<5;tr++)
    {
      if(ntype[tr]>0)
	{
	  epsilon=0.5*(softTable[tr] + softTable[Part[target]->type]);

	  force_treeevaluate_single(tr,target,epsilon);
	}
    }
}




void force_treeevaluate_single(int tree, int targetpart, double epsilon)  /* non-recursive walk */
{
  struct NODE *no,*nn;
  int i,k,p,ii;
  double r2,r5,dx,dy,dz,r,fac,theta,u,h,h2_inv,h3_inv,h5_inv,h4_inv,h6_inv,ff;
  double wf,wp,w2,w3,w4,potq,r3,maxd;
  struct xyz_data *tppos, *tpaccel;   /* target particle for treewalk */
  double *tppot;
  double q11dx,q12dy,q13dz,q12dx,q22dy,q23dz,q13dx,q23dy,q33dz;
  double r_inv,r2_inv,r3_inv,r5_inv;
  double h_inv;

  h = 2.8*epsilon;
  h_inv=1/h;
  h2_inv=h_inv*h_inv;
  h3_inv=h2_inv*h_inv;
  h4_inv=h2_inv*h2_inv;
  h5_inv=h2_inv*h3_inv;
  h6_inv=h3_inv*h3_inv;



  tppos   = (struct xyz_data *)Part[targetpart]->xyz;
  tpaccel = (struct xyz_data *)Part[targetpart]->accel;
  tppot   = &Part[targetpart]->pot;


  no=trees[tree];

  while(no)
    {
      dx=no->s[0]- tppos->xyz[0];     /* observe the sign ! */
      dy=no->s[1]- tppos->xyz[1];     /* this vector is -y in my thesis notation */
      dz=no->s[2]- tppos->xyz[2];

      r2=dx*dx+dy*dy+dz*dz;


      if((p=no->partind)>=0)   /* single particle */
	{
	  r=sqrt(r2);  
	   
	  u=r*h_inv;


	  treecost[tree]+=1;	  	  

	  if(u>=1)
	    {
	      r_inv=1/r;

	      fac=no->mass*r_inv*r_inv*r_inv;
	  
	      tpaccel->xyz[0]+=dx*fac;
	      tpaccel->xyz[1]+=dy*fac;
	      tpaccel->xyz[2]+=dz*fac;
#ifdef POTENTIAL_SIMULTO_FORCE
	      *tppot-=no->mass*r_inv;
#endif
	    }
	  else
	    {
	      ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
	      wf=knlforce[ii]+(knlforce[ii+1]-knlforce[ii])*ff;
	      wp=knlpot[ii]+(knlpot[ii+1]-knlpot[ii])*ff;
	      
	      if(r>1.0e-20)
		{
		  fac=no->mass*h_inv*h_inv/r*wf;
		  
		  tpaccel->xyz[0]+=dx*fac;
		  tpaccel->xyz[1]+=dy*fac;
		  tpaccel->xyz[2]+=dz*fac;
		}
#ifdef POTENTIAL_SIMULTO_FORCE
	      *tppot+=no->mass*h_inv*wp;
#endif
	    }
	  no=no->sibling;
	}
      else
	{
	  if(r2 < no->oc)
	    {
	      no=no->next;  /* open cell */
	    }
	  else
	    {
	      /*	      treecost[tree]+=1; */
  treecost_quadru[tree]+=1;	  	  
	      
	      no->cost += 1;

	      r=sqrt(r2);  
	  
	      u=r*h_inv;
	  
	      if(u>=1)  /* ordinary quadrupol moment */
		{
		  r_inv=1/r;
		  r2_inv=r_inv*r_inv;
		  r3_inv=r2_inv*r_inv;
		  r5_inv=r2_inv*r3_inv;

                  q11dx=no->Q11*dx;
                  q12dy=no->Q12*dy;
		  q13dz=no->Q13*dz;
		  q12dx=no->Q12*dx;
		  q22dy=no->Q22*dy;
		  q23dz=no->Q23*dz;
		  q13dx=no->Q13*dx;
                  q23dy=no->Q23*dy;
		  q33dz=no->Q33*dz;

		  potq=0.5*(q11dx*dx+              /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;
#ifdef POTENTIAL_SIMULTO_FORCE 
		  *tppot += -no->mass*r_inv  /* monopole */
                           +r3_inv*( -3*potq*r2_inv + 0.5*no->P);  /* quadrupole */
#endif
		  
		  fac=no->mass*r3_inv  /* monopole force*/
		      +(15*potq*r2_inv -1.5*no->P)*r5_inv;  /* radial quadrupole part */

		  tpaccel->xyz[0]+=dx*fac;
		  tpaccel->xyz[1]+=dy*fac;
		  tpaccel->xyz[2]+=dz*fac;

		  /* add tensor contribution */

		  ff=-3*r5_inv;
		  tpaccel->xyz[0] += ff*(q11dx + q12dy + q13dz);
		  tpaccel->xyz[1] += ff*(q12dx + q22dy + q23dz);
		  tpaccel->xyz[2] += ff*(q13dx + q23dy + q33dz);

		}
	      else    /* softened quadrupol moment */
		{
		  ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
		  wf=knlforce[ii] + (knlforce[ii+1]-knlforce[ii])*ff;
		  wp=knlpot[ii]   + (knlpot[ii+1]-knlpot[ii])*ff;
		  w2=knlW2[ii]    + (knlW2[ii+1]-knlW2[ii])*ff;
		  w3=knlW3[ii]    + (knlW3[ii+1]-knlW3[ii])*ff;
		  w4=knlW4[ii]    + (knlW4[ii+1]-knlW4[ii])*ff;

		  r_inv=1/r;

                  q11dx=no->Q11*dx;
                  q12dy=no->Q12*dy;
		  q13dz=no->Q13*dz;
		  q12dx=no->Q12*dx;
		  q22dy=no->Q22*dy;
		  q23dz=no->Q23*dz;
		  q13dx=no->Q13*dx;
                  q23dy=no->Q23*dy;
		  q33dz=no->Q33*dz;

		  potq=0.5*(q11dx*dx+      /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;
#ifdef POTENTIAL_SIMULTO_FORCE
		  *tppot += no->mass*h_inv*wp +  /* monopole */   
		            potq*w2*h5_inv + 0.5*no->P*wf*h2_inv*r_inv ; /* quadru contribution */
#endif
		  /* note: observe definition (sign!) of dx,dy,dz */

		  fac=no->mass*h2_inv*r_inv*wf +  /* monopole force */ 
                     +potq*h6_inv * w3*r_inv   + 0.5*no->P * w4 *h4_inv*r_inv; /* radial contribution
									     of quadrupole */
		  tpaccel->xyz[0]+=dx*fac;
		  tpaccel->xyz[1]+=dy*fac;
		  tpaccel->xyz[2]+=dz*fac;
		  
		  /* add tensor contribution */
		  ff=w2*h5_inv;
		  tpaccel->xyz[0] += ff*(q11dx + q12dy + q13dz);
		  tpaccel->xyz[1] += ff*(q12dx + q22dy + q23dz);
		  tpaccel->xyz[2] += ff*(q13dx + q23dy + q33dz);
		}

	      no=no->sibling;
	    }
	}
    }

}



























void force_treeevaluate_potential(int target) 
{
  int    i,tr;
  double epsilon;
  void force_treeevaluate_potential_single(int tree, int targetpart, double epsilon);


  Part[target]->pot=0;

  for(tr=0;tr<5;tr++)
    {
      if(ntype[tr]>0)
	{
	  epsilon=0.5*(softTable[tr] + softTable[Part[target]->type]);

	  force_treeevaluate_potential_single(tr,target,epsilon);
	}
    }
}





void force_treeevaluate_potential_single(int tree, int targetpart, double epsilon)  /* non-recursive walk */
{
  struct NODE *no,*nn;
  int i,k,p,ii;
  double r2,r5,dx,dy,dz,r,fac,theta,u,h,h2_inv,h3_inv,h5_inv,h4_inv,h6_inv,ff;
  double wf,wp,w2,w3,w4,potq,r3,maxd;
  struct xyz_data *tppos;   /* target particle for treewalk */
  double *tppot;
  double q11dx,q12dy,q13dz,q12dx,q22dy,q23dz,q13dx,q23dy,q33dz;
  double r_inv,r2_inv,r3_inv,r5_inv;
  double h_inv;

  h = 2.8*epsilon;
  h_inv=1/h;
  h2_inv=h_inv*h_inv;
  h3_inv=h2_inv*h_inv;
  h4_inv=h2_inv*h2_inv;
  h5_inv=h2_inv*h3_inv;
  h6_inv=h3_inv*h3_inv;



  tppos   = (struct xyz_data *)Part[targetpart]->xyz;
  tppot   = &Part[targetpart]->pot;


  no=trees[tree];

  while(no)
    {
      dx=no->s[0]- tppos->xyz[0];     /* observe the sign ! */
      dy=no->s[1]- tppos->xyz[1];     /* this vector is -y in my thesis notation */
      dz=no->s[2]- tppos->xyz[2];

      r2=dx*dx+dy*dy+dz*dz;


      if((p=no->partind)>=0)   /* single particle */
	{
	  r=sqrt(r2);  
	   
	  u=r*h_inv;

	  if(u>=1)
	    {
	      *tppot-=no->mass/r;
	    }
	  else
	    {
	      ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
	      wp=knlpot[ii]+(knlpot[ii+1]-knlpot[ii])*ff;
	      
	      *tppot+=no->mass*h_inv*wp;
	    }
	  no=no->sibling;
	}
      else
	{
	  if(r2 < no->oc)
	    {
	      no=no->next;  /* open cell */
	    }
	  else
	    {
	      r=sqrt(r2);  
	  
	      u=r*h_inv;
	  
	      if(u>=1)  /* ordinary quadrupol moment */
		{
		  r_inv=1/r;
		  r2_inv=r_inv*r_inv;
		  r3_inv=r2_inv*r_inv;
		  r5_inv=r2_inv*r3_inv;

                  q11dx=no->Q11*dx;
                  q12dy=no->Q12*dy;
		  q13dz=no->Q13*dz;
		  q12dx=no->Q12*dx;
		  q22dy=no->Q22*dy;
		  q23dz=no->Q23*dz;
		  q13dx=no->Q13*dx;
                  q23dy=no->Q23*dy;
		  q33dz=no->Q33*dz;

		  potq=0.5*(q11dx*dx+              /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;
 
		  *tppot += -no->mass*r_inv  /* monopole */
                           +r3_inv*( -3*potq*r2_inv + 0.5*no->P);  /* quadrupole */

		}
	      else    /* softened quadrupol moment */
		{
		  ii = (int)(u*KERN_LEN); ff=(u-knlrad[ii])*KERN_LEN;
		  wf=knlforce[ii] + (knlforce[ii+1]-knlforce[ii])*ff;
		  wp=knlpot[ii]   + (knlpot[ii+1]-knlpot[ii])*ff;
		  w2=knlW2[ii]    + (knlW2[ii+1]-knlW2[ii])*ff;


		  r_inv=1/r;

                  q11dx=no->Q11*dx;
                  q12dy=no->Q12*dy;
		  q13dz=no->Q13*dz;
		  q12dx=no->Q12*dx;
		  q22dy=no->Q22*dy;
		  q23dz=no->Q23*dz;
		  q13dx=no->Q13*dx;
                  q23dy=no->Q23*dy;
		  q33dz=no->Q33*dz;

		  potq=0.5*(q11dx*dx+      /* 1/2 y^T Q y */
			    q22dy*dy+
			    q33dz*dz)+
		       q12dx*dy+
		       q13dx*dz+
		       q23dy*dz;

		  *tppot += no->mass*h_inv*wp +  /* monopole */   
		            potq*w2*h5_inv + 0.5*no->P*wf*h2_inv*r_inv ; /* quadru contribution */

		}

	      no=no->sibling;
	    }
	}
    }

}





















void force_setkernel(void) 
{
  int i;
  double u;

  for(i=0;i<=KERN_LEN;i++)
    {
      u=((double)i)/KERN_LEN;

      knlrad[i] = u;

      if(u<=0.5)
	{
	  knlforce[i]=32 * (u/3 -6.0/5*pow(u,3) + pow(u,4));
	  knlpot[i]=16.0/3*pow(u,2)-48.0/5*pow(u,4)+32.0/5*pow(u,5)-14.0/5;

	  knlW2[i]= -384.0/5 +96.0*u;
	  knlW3[i]= 96.0;
	  knlW4[i]= 96.0/5*u*(5*u-4);
	}
      else
	{
	  knlforce[i]=64*(u/3 -3.0/4*u*u + 3.0/5*pow(u,3)-pow(u,4)/6) - 1.0/15/pow(u,2);
	  knlpot[i]=1.0/15/u +32.0/3*pow(u,2)-16.0*pow(u,3)+48.0/5*pow(u,4)-32.0/15*pow(u,5)-16.0/5;

	  knlW2[i]= 384.0/5 + 1/(5.0*pow(u,5)) -48.0/u -32*u;
	  knlW3[i]= -32-1/pow(u,6)+48/pow(u,2);
	  knlW4[i]= -48 +1/(5*pow(u,4)) +384.0/5*u -32*pow(u,2);
	}
    }
}




void force_treeallocate(int maxnodes)  /* usually maxnodes=2*npart is sufficient */
{
  long bytes;

  MaxNodes=maxnodes;

  if(!(nodes=malloc(bytes=MaxNodes*sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%d Mbytes).\n",MaxNodes,(int) (bytes/1024.0*1024.0));
      exit(3);
    }

  printf("\nAllocated %g MByte for BH-tree.\n\n",bytes/(1024.0*1024.0));
    
   force_setkernel();
}



void force_treefree(void)
{
  free(nodes);
}


void force_setpointers(double *_softTable)
{
  softTable=_softTable;
}




int force_getcost(void)
{
  int tr,cost;

  for(tr=cost=0;tr<5;tr++) 
    cost += treecost[tr];
  
  return cost;
}

int force_getcost_quadru(void)
{
  int tr,cost;

  for(tr=cost=0;tr<5;tr++) 
    cost += treecost_quadru[tr];
  
  return cost;
}
