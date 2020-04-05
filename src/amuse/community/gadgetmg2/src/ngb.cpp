#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#ifndef NOMPI
#include <mpi.h>
#endif

//#include "allvars.hpp"
#include "proto.hpp"


/*! \file ngb.c
 *  \brief neighbour search by means of the tree
 *
 *  This file contains routines for neighbour finding.  We use the
 *  gravity-tree and a range-searching technique to find neighbours.
 */


/*! these macros maps a coordinate difference to the nearest periodic
 * image
 */

#define NGB_PERIODIC_X(x) (xtmp=(x),(xtmp>boxHalf_X)?(xtmp-boxSize_X):((xtmp<-boxHalf_X)?(xtmp+boxSize_X):xtmp))
#define NGB_PERIODIC_Y(x) (xtmp=(x),(xtmp>boxHalf_Y)?(xtmp-boxSize_Y):((xtmp<-boxHalf_Y)?(xtmp+boxSize_Y):xtmp))
#define NGB_PERIODIC_Z(x) (xtmp=(x),(xtmp>boxHalf_Z)?(xtmp-boxSize_Z):((xtmp<-boxHalf_Z)?(xtmp+boxSize_Z):xtmp))



/*! This routine finds all neighbours `j' that can interact with the
 *  particle `i' in the communication buffer.
 *
 *  Note that an interaction can take place if
 *  \f$ r_{ij} < h_i \f$  OR if  \f$ r_{ij} < h_j \f$.
 *
 *  In the range-search this is taken into account, i.e. it is guaranteed that
 *  all particles are found that fulfil this condition, including the (more
 *  difficult) second part of it. For this purpose, each node knows the
 *  maximum h occuring among the particles it represents.
 */
int gadgetmg2::ngb_treefind_pairs(double searchcenter[3], double hsml, int *startnode)
{
  int k, no, p, numngb;
  double hdiff;
  double searchmin[3], searchmax[3];
  struct NODE *thiis;


  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - hsml;
      searchmax[k] = searchcenter[k] + hsml;
    }

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;

	  hdiff = SphP[p].Hsml - hsml;
	  if(hdiff < 0)
	    hdiff = 0;


	  if(P[p].Pos[0] < (searchmin[0] - hdiff))
	    continue;
	  if(P[p].Pos[0] > (searchmax[0] + hdiff))
	    continue;
	  if(P[p].Pos[1] < (searchmin[1] - hdiff))
	    continue;
	  if(P[p].Pos[1] > (searchmax[1] + hdiff))
	    continue;
	  if(P[p].Pos[2] < (searchmin[2] - hdiff))
	    continue;
	  if(P[p].Pos[2] > (searchmax[2] + hdiff))
	    continue;
	  Ngblist[numngb++] = p;

	  if(numngb == MAX_NGB)
	    {
	      printf
		("ThisTask=%d: Need to do a second neighbour loop in hydro-force for (%g|%g|%g) hsml=%g no=%d\n",
		 ThisTask, searchcenter[0], searchcenter[1], searchcenter[2], hsml, no);
	      *startnode = no;
	      return numngb;
	    }
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  thiis = &Nodes[no];
	  hdiff = Extnodes[no].hmax - hsml;
	  if(hdiff < 0)
	    hdiff = 0;

	  no = thiis->u.d.sibling;	/* in case the node can be discarded */


	  if((thiis->center[0] + 0.5 * thiis->len) < (searchmin[0] - hdiff))
	    continue;
	  if((thiis->center[0] - 0.5 * thiis->len) > (searchmax[0] + hdiff))
	    continue;
	  if((thiis->center[1] + 0.5 * thiis->len) < (searchmin[1] - hdiff))
	    continue;
	  if((thiis->center[1] - 0.5 * thiis->len) > (searchmax[1] + hdiff))
	    continue;
	  if((thiis->center[2] + 0.5 * thiis->len) < (searchmin[2] - hdiff))
	    continue;
	  if((thiis->center[2] - 0.5 * thiis->len) > (searchmax[2] + hdiff))
	    continue;

	  no = thiis->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}



/*! This function returns neighbours with distance <= hsml and returns them in
 *  Ngblist. Actually, particles in a box of half side length hsml are
 *  returned, i.e. the reduction to a sphere still needs to be done in the
 *  calling routine.
 */
int gadgetmg2::ngb_treefind_variable(double searchcenter[3], double hsml, int *startnode)
{
  int k, numngb;
  int no, p;
  struct NODE *thiis;
  double searchmin[3], searchmax[3];

  for(k = 0; k < 3; k++)	/* cube-box window */
    {
      searchmin[k] = searchcenter[k] - hsml;
      searchmax[k] = searchcenter[k] + hsml;
    }

  numngb = 0;
  no = *startnode;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  p = no;
	  no = Nextnode[no];

	  if(P[p].Type > 0)
	    continue;


	  if(P[p].Pos[0] < searchmin[0])
	    continue;
	  if(P[p].Pos[0] > searchmax[0])
	    continue;
	  if(P[p].Pos[1] < searchmin[1])
	    continue;
	  if(P[p].Pos[1] > searchmax[1])
	    continue;
	  if(P[p].Pos[2] < searchmin[2])
	    continue;
	  if(P[p].Pos[2] > searchmax[2])
	    continue;
	  Ngblist[numngb++] = p;

	  if(numngb == MAX_NGB)
	    {
	      numngb = ngb_clear_buf(searchcenter, hsml, numngb);
	      if(numngb == MAX_NGB)
		{
		  printf("ThisTask=%d: Need to do a second neighbour loop for (%g|%g|%g) hsml=%g no=%d\n",
			 ThisTask, searchcenter[0], searchcenter[1], searchcenter[2], hsml, no);
		  *startnode = no;
		  return numngb;
		}
	    }
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  thiis = &Nodes[no];

	  no = thiis->u.d.sibling;	/* in case the node can be discarded */

	  if((thiis->center[0] + 0.5 * thiis->len) < (searchmin[0]))
	    continue;
	  if((thiis->center[0] - 0.5 * thiis->len) > (searchmax[0]))
	    continue;
	  if((thiis->center[1] + 0.5 * thiis->len) < (searchmin[1]))
	    continue;
	  if((thiis->center[1] - 0.5 * thiis->len) > (searchmax[1]))
	    continue;
	  if((thiis->center[2] + 0.5 * thiis->len) < (searchmin[2]))
	    continue;
	  if((thiis->center[2] - 0.5 * thiis->len) > (searchmax[2]))
	    continue;
	  no = thiis->u.d.nextnode;	/* ok, we need to open the node */
	}
    }

  *startnode = -1;
  return numngb;
}




/*! The buffer for the neighbour list has a finite length MAX_NGB. For a large
 *  search region, this buffer can get full, in which case this routine can be
 *  called to eliminate some of the superfluous particles in the "corners" of
 *  the search box - only the ones in the inscribed sphere need to be kept.
 */
int gadgetmg2::ngb_clear_buf(double searchcenter[3], double hsml, int numngb)
{
  int i, p;
  double dx, dy, dz, r2;


  for(i = 0; i < numngb; i++)
    {
      p = Ngblist[i];

      dx = P[p].Pos[0] - searchcenter[0];
      dy = P[p].Pos[1] - searchcenter[1];
      dz = P[p].Pos[2] - searchcenter[2];

      r2 = dx * dx + dy * dy + dz * dz;

      if(r2 > hsml * hsml)
	{
	  Ngblist[i] = Ngblist[numngb - 1];
	  i--;
	  numngb--;
	}
    }

  return numngb;
}



/*! Allocates memory for the neighbour list buffer.
 */
void gadgetmg2::ngb_treeallocate(int npart)
{
  double totbytes = 0;
  size_t bytes;
  if(!(Ngblist = (int*)malloc(bytes = npart * (long) sizeof(int))))
    {
      printf("Failed to allocate %g MB for ngblist array\n", bytes / (1024.0 * 1024.0));
      endrun(78);
    }
  totbytes += bytes;

  if(ThisTask == 0)
    printf("allocated %g Mbyte for ngb search.\n", totbytes / (1024.0 * 1024.0));
}


/*! free memory allocated for neighbour list buffer.
 */
void gadgetmg2::ngb_treefree(void)
{
  free(Ngblist);
}

/*! This function constructs the neighbour tree. To this end, we actually need
 *  to construct the gravitational tree, because we use it now for the
 *  neighbour search.
 */
void gadgetmg2::ngb_treebuild(void)
{
  if(ThisTask == 0)
    printf("Begin Ngb-tree construction.\n");

  force_treebuild(N_gas);

  if(ThisTask == 0)
    printf("Ngb-Tree contruction finished \n");
}

