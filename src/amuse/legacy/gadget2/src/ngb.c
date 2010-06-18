#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file ngb.c
 *  \brief neighbour search by means of the tree
 *
 *  This file contains routines for neighbour finding.  We use the
 *  gravity-tree and a range-searching technique to find neighbours.
 */

#ifdef PERIODIC
static double boxSize, boxHalf;

#ifdef LONG_X
static double boxSize_X, boxHalf_X;
#else
#define boxSize_X boxSize
#define boxHalf_X boxHalf
#endif
#ifdef LONG_Y
static double boxSize_Y, boxHalf_Y;
#else
#define boxSize_Y boxSize
#define boxHalf_Y boxHalf
#endif
#ifdef LONG_Z
static double boxSize_Z, boxHalf_Z;
#else
#define boxSize_Z boxSize
#define boxHalf_Z boxHalf
#endif
#endif


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
int ngb_treefind_pairs(FLOAT searchcenter[3], FLOAT hsml, int *startnode)
{
  int k, no, p, numngb;
  FLOAT hdiff;
  FLOAT searchmin[3], searchmax[3];
  struct NODE *this;

#ifdef PERIODIC
  double xtmp;
#endif

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

#ifdef PERIODIC
	  if(NGB_PERIODIC_X(P[p].Pos[0] - searchcenter[0]) < (-hsml - hdiff))
	    continue;
	  if(NGB_PERIODIC_X(P[p].Pos[0] - searchcenter[0]) > (hsml + hdiff))
	    continue;
	  if(NGB_PERIODIC_Y(P[p].Pos[1] - searchcenter[1]) < (-hsml - hdiff))
	    continue;
	  if(NGB_PERIODIC_Y(P[p].Pos[1] - searchcenter[1]) > (hsml + hdiff))
	    continue;
	  if(NGB_PERIODIC_Z(P[p].Pos[2] - searchcenter[2]) < (-hsml - hdiff))
	    continue;
	  if(NGB_PERIODIC_Z(P[p].Pos[2] - searchcenter[2]) > (hsml + hdiff))
	    continue;
#else
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
#endif
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

	  this = &Nodes[no];
	  hdiff = Extnodes[no].hmax - hsml;
	  if(hdiff < 0)
	    hdiff = 0;

	  no = this->u.d.sibling;	/* in case the node can be discarded */

#ifdef PERIODIC
	  if((NGB_PERIODIC_X(this->center[0] - searchcenter[0]) + 0.5 * this->len) < (-hsml - hdiff))
	    continue;
	  if((NGB_PERIODIC_X(this->center[0] - searchcenter[0]) - 0.5 * this->len) > (hsml + hdiff))
	    continue;
	  if((NGB_PERIODIC_Y(this->center[1] - searchcenter[1]) + 0.5 * this->len) < (-hsml - hdiff))
	    continue;
	  if((NGB_PERIODIC_Y(this->center[1] - searchcenter[1]) - 0.5 * this->len) > (hsml + hdiff))
	    continue;
	  if((NGB_PERIODIC_Z(this->center[2] - searchcenter[2]) + 0.5 * this->len) < (-hsml - hdiff))
	    continue;
	  if((NGB_PERIODIC_Z(this->center[2] - searchcenter[2]) - 0.5 * this->len) > (hsml + hdiff))
	    continue;
#else
	  if((this->center[0] + 0.5 * this->len) < (searchmin[0] - hdiff))
	    continue;
	  if((this->center[0] - 0.5 * this->len) > (searchmax[0] + hdiff))
	    continue;
	  if((this->center[1] + 0.5 * this->len) < (searchmin[1] - hdiff))
	    continue;
	  if((this->center[1] - 0.5 * this->len) > (searchmax[1] + hdiff))
	    continue;
	  if((this->center[2] + 0.5 * this->len) < (searchmin[2] - hdiff))
	    continue;
	  if((this->center[2] - 0.5 * this->len) > (searchmax[2] + hdiff))
	    continue;
#endif
	  no = this->u.d.nextnode;	/* ok, we need to open the node */
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
int ngb_treefind_variable(FLOAT searchcenter[3], FLOAT hsml, int *startnode)
{
  int k, numngb;
  int no, p;
  struct NODE *this;
  FLOAT searchmin[3], searchmax[3];

#ifdef PERIODIC
  double xtmp;
#endif

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

#ifdef PERIODIC
	  if(NGB_PERIODIC_X(P[p].Pos[0] - searchcenter[0]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_X(P[p].Pos[0] - searchcenter[0]) > hsml)
	    continue;
	  if(NGB_PERIODIC_Y(P[p].Pos[1] - searchcenter[1]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_Y(P[p].Pos[1] - searchcenter[1]) > hsml)
	    continue;
	  if(NGB_PERIODIC_Z(P[p].Pos[2] - searchcenter[2]) < -hsml)
	    continue;
	  if(NGB_PERIODIC_Z(P[p].Pos[2] - searchcenter[2]) > hsml)
	    continue;
#else
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
#endif
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

	  this = &Nodes[no];

	  no = this->u.d.sibling;	/* in case the node can be discarded */
#ifdef PERIODIC
	  if((NGB_PERIODIC_X(this->center[0] - searchcenter[0]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_X(this->center[0] - searchcenter[0]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_Y(this->center[1] - searchcenter[1]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_Y(this->center[1] - searchcenter[1]) - 0.5 * this->len) > hsml)
	    continue;
	  if((NGB_PERIODIC_Z(this->center[2] - searchcenter[2]) + 0.5 * this->len) < -hsml)
	    continue;
	  if((NGB_PERIODIC_Z(this->center[2] - searchcenter[2]) - 0.5 * this->len) > hsml)
	    continue;
#else
	  if((this->center[0] + 0.5 * this->len) < (searchmin[0]))
	    continue;
	  if((this->center[0] - 0.5 * this->len) > (searchmax[0]))
	    continue;
	  if((this->center[1] + 0.5 * this->len) < (searchmin[1]))
	    continue;
	  if((this->center[1] - 0.5 * this->len) > (searchmax[1]))
	    continue;
	  if((this->center[2] + 0.5 * this->len) < (searchmin[2]))
	    continue;
	  if((this->center[2] - 0.5 * this->len) > (searchmax[2]))
	    continue;
#endif
	  no = this->u.d.nextnode;	/* ok, we need to open the node */
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
int ngb_clear_buf(FLOAT searchcenter[3], FLOAT hsml, int numngb)
{
  int i, p;
  FLOAT dx, dy, dz, r2;

#ifdef PERIODIC
  double xtmp;
#endif

  for(i = 0; i < numngb; i++)
    {
      p = Ngblist[i];
#ifdef PERIODIC
      dx = NGB_PERIODIC_X(P[p].Pos[0] - searchcenter[0]);
      dy = NGB_PERIODIC_Y(P[p].Pos[1] - searchcenter[1]);
      dz = NGB_PERIODIC_Z(P[p].Pos[2] - searchcenter[2]);
#else
      dx = P[p].Pos[0] - searchcenter[0];
      dy = P[p].Pos[1] - searchcenter[1];
      dz = P[p].Pos[2] - searchcenter[2];
#endif
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
void ngb_treeallocate(int npart)
{
  double totbytes = 0;
  size_t bytes;

#ifdef PERIODIC
  boxSize = All.BoxSize;
  boxHalf = 0.5 * All.BoxSize;
#ifdef LONG_X
  boxHalf_X = boxHalf * LONG_X;
  boxSize_X = boxSize * LONG_X;
#endif
#ifdef LONG_Y
  boxHalf_Y = boxHalf * LONG_Y;
  boxSize_Y = boxSize * LONG_Y;
#endif
#ifdef LONG_Z
  boxHalf_Z = boxHalf * LONG_Z;
  boxSize_Z = boxSize * LONG_Z;
#endif
#endif

  if(!(Ngblist = malloc(bytes = npart * (long) sizeof(int))))
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
void ngb_treefree(void)
{
  free(Ngblist);
}

/*! This function constructs the neighbour tree. To this end, we actually need
 *  to construct the gravitational tree, because we use it now for the
 *  neighbour search.
 */
void ngb_treebuild(void)
{
  if(ThisTask == 0)
    printf("Begin Ngb-tree construction.\n");

  force_treebuild(N_gas);

  if(ThisTask == 0)
    printf("Ngb-Tree contruction finished \n");
}

