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


/*! \file forcetree.c
 *  \brief gravitational tree and code for Ewald correction
 *
 *  This file contains the computation of the gravitational force by means
 *  of a tree. The type of tree implemented is a geometrical oct-tree,
 *  starting from a cube encompassing all particles. This cube is
 *  automatically found in the domain decomposition, which also splits up
 *  the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. Tree nodes can
 *  be dynamically updated in drift/kick operations to avoid having to
 *  reconstruct the tree every timestep.
 */

/*! auxialiary variable used to set-up non-recursive walk */
static int last;



/*! length of lock-up table for short-range force kernel in TreePM algorithm */
#define NTAB 1000
/*! variables for short-range lookup table */
static double tabfac, shortrange_table[NTAB], shortrange_table_potential[NTAB];

/*! toggles after first tree-memory allocation, has only influence on log-files */
static int first_flag = 0;


/*! This function is a driver routine for constructing the gravitational
 *  oct-tree, which is done by calling a small number of other functions.
 */
int gadgetmg2::force_treebuild(int npart)
{
  Numnodestree = force_treebuild_single(npart);

  force_update_pseudoparticles();

  force_flag_localnodes();

  TimeOfLastTreeConstruction = All.Time;

  return Numnodestree;
}



/*! Constructs the gravitational oct-tree.
 *
 *  The index convention for accessing tree nodes is the following: the
 *  indices 0...NumPart-1 reference single particles, the indices
 *  All.MaxPart.... All.MaxPart+nodes-1 reference tree nodes. `Nodes_base'
 *  points to the first tree node, while `nodes' is shifted such that
 *  nodes[All.MaxPart] gives the first tree node. Finally, node indices
 *  with values 'All.MaxPart + MaxNodes' and larger indicate "pseudo
 *  particles", i.e. multipole moments of top-level nodes that lie on
 *  different CPUs. If such a node needs to be opened, the corresponding
 *  particle must be exported to that CPU. The 'Extnodes' structure
 *  parallels that of 'Nodes'. Its information is only needed for the SPH
 *  part of the computation. (The data is split onto these two structures
 *  as a tuning measure.  If it is merged into 'Nodes' a somewhat bigger
 *  size of the nodes also for gravity would result, which would reduce
 *  cache utilization slightly.
 */
int gadgetmg2::force_treebuild_single(int npart)
{
  int i, j, subnode = 0, parent, numnodes;
  int nfree, th, nn, no;
  struct NODE *nfreep;
  double lenhalf, epsilon;
  peanokey key;


  /* create an empty root node  */
  nfree = All.MaxPart;		/* index of first free node */
  nfreep = &Nodes[nfree];	/* select first node */

  nfreep->len = DomainLen;
  for(j = 0; j < 3; j++)
    nfreep->center[j] = DomainCenter[j];
  for(j = 0; j < 8; j++)
    nfreep->u.suns[j] = -1;


  numnodes = 1;
  nfreep++;
  nfree++;

  /* create a set of empty nodes corresponding to the top-level domain
   * grid. We need to generate these nodes first to make sure that we have a
   * complete top-level tree which allows the easy insertion of the
   * pseudo-particles at the right place
   */

  force_create_empty_nodes(All.MaxPart, 0, 1, 0, 0, 0, &numnodes, &nfree);


  /* if a high-resolution region in a global tree is used, we need to generate
   * an additional set empty nodes to make sure that we have a complete
   * top-level tree for the high-resolution inset
   */

  nfreep = &Nodes[nfree];
  parent = -1;			/* note: will not be used below before it is changed */


  /* now we insert all particles */
  for(i = 0; i < npart; i++)
    {

      /* the softening is only used to check whether particles are so close
       * that the tree needs not to be refined further
       */
      epsilon = All.ForceSoftening[P[i].Type];

      key = peano_hilbert_key((P[i].Pos[0] - DomainCorner[0]) * DomainFac,
			      (P[i].Pos[1] - DomainCorner[1]) * DomainFac,
			      (P[i].Pos[2] - DomainCorner[2]) * DomainFac, BITS_PER_DIMENSION);

      no = 0;
      while(TopNodes[no].Daughter >= 0)
	no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Leaf;
      th = DomainNodeIndex[no];

      while(1)
	{
	  if(th >= All.MaxPart)	/* we are dealing with an internal node */
	    {
	      subnode = 0;
	      if(P[i].Pos[0] > Nodes[th].center[0])
		subnode += 1;
	      if(P[i].Pos[1] > Nodes[th].center[1])
		subnode += 2;
	      if(P[i].Pos[2] > Nodes[th].center[2])
		subnode += 4;

	      nn = Nodes[th].u.suns[subnode];

	      if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		{
		  parent = th;
		  th = nn;
		}
	      else
		{
		  /* here we have found an empty slot where we can attach
		   * the new particle as a leaf.
		   */
		  Nodes[th].u.suns[subnode] = i;
		  break;	/* done for this particle */
		}
	    }
	  else
	    {
	      /* We try to insert into a leaf with a single particle.  Need
	       * to generate a new internal node at this point.
	       */
	      Nodes[parent].u.suns[subnode] = nfree;

	      nfreep->len = 0.5 * Nodes[parent].len;
	      lenhalf = 0.25 * Nodes[parent].len;

	      if(subnode & 1)
		nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
	      else
		nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

	      if(subnode & 2)
		nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
	      else
		nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

	      if(subnode & 4)
		nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
	      else
		nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

	      nfreep->u.suns[0] = -1;
	      nfreep->u.suns[1] = -1;
	      nfreep->u.suns[2] = -1;
	      nfreep->u.suns[3] = -1;
	      nfreep->u.suns[4] = -1;
	      nfreep->u.suns[5] = -1;
	      nfreep->u.suns[6] = -1;
	      nfreep->u.suns[7] = -1;


	      subnode = 0;
	      if(P[th].Pos[0] > nfreep->center[0])
		subnode += 1;
	      if(P[th].Pos[1] > nfreep->center[1])
		subnode += 2;
	      if(P[th].Pos[2] > nfreep->center[2])
		subnode += 4;
#ifndef NOTREERND
	      if(nfreep->len < 1.0e-3 * epsilon)
		{
		  /* seems like we're dealing with particles at identical (or extremely close)
		   * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
		   * of tree are still correct, but this will only happen well below gravitational softening
		   * length-scale anyway.
		   */
		  subnode = (int) (8.0 * get_random_number((0xffff & P[i].ID) + P[i].GravCost));
		  P[i].GravCost += 1;
		  if(subnode >= 8)
		    subnode = 7;
		}
#endif
	      nfreep->u.suns[subnode] = th;

	      th = nfree;	/* resume trying to insert the new particle at
				 * the newly created internal node
				 */

	      numnodes++;
	      nfree++;
	      nfreep++;

	      if((numnodes) >= MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, MaxNodes);
		  printf("for particle %d\n", i);
		  dump_particles();
		  endrun(1);
		}
	    }
	}
    }


  /* insert the pseudo particles that represent the mass distribution of other domains */
  force_insert_pseudo_particles();


  /* now compute the multipole moments recursively */
  last = -1;

  force_update_node_recursive(All.MaxPart, -1, -1);

  if(last >= All.MaxPart)
    {
      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
	Nextnode[last - MaxNodes] = -1;
      else
	Nodes[last].u.d.nextnode = -1;
    }
  else
    Nextnode[last] = -1;

  return numnodes;
}



/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can
 *  easily associate the pseudo-particles of other CPUs with tree-nodes at
 *  a given level in the tree, even when the particle population is so
 *  sparse that some of these nodes are actually empty.
*/
void gadgetmg2::force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount,
			      int *nextfree)
{
  int i, j, k, n, sub, count;

  if(TopNodes[topnode].Daughter >= 0)
    {
      for(i = 0; i < 2; i++)
	for(j = 0; j < 2; j++)
	  for(k = 0; k < 2; k++)
	    {
	      sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

	      count = i + 2 * j + 4 * k;

	      Nodes[no].u.suns[count] = *nextfree;


	      Nodes[*nextfree].len = 0.5 * Nodes[no].len;
	      Nodes[*nextfree].center[0] = Nodes[no].center[0] + (2 * i - 1) * 0.25 * Nodes[no].len;
	      Nodes[*nextfree].center[1] = Nodes[no].center[1] + (2 * j - 1) * 0.25 * Nodes[no].len;
	      Nodes[*nextfree].center[2] = Nodes[no].center[2] + (2 * k - 1) * 0.25 * Nodes[no].len;

	      for(n = 0; n < 8; n++)
		Nodes[*nextfree].u.suns[n] = -1;

	      if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
		DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = *nextfree;

	      *nextfree = *nextfree + 1;
	      *nodecount = *nodecount + 1;

	      if((*nodecount) >= MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, MaxNodes);
		  printf("in create empty nodes\n");
		  dump_particles();
		  endrun(11);
		}

	      force_create_empty_nodes(*nextfree - 1, TopNodes[topnode].Daughter + sub,
				       bits + 1, 2 * x + i, 2 * y + j, 2 * z + k, nodecount, nextfree);
	    }
    }
}



/*! this function inserts pseudo-particles which will represent the mass
 *  distribution of the other CPUs. Initially, the mass of the
 *  pseudo-particles is set to zero, and their coordinate is set to the
 *  center of the domain-cell they correspond to. These quantities will be
 *  updated later on.
 */
void gadgetmg2::force_insert_pseudo_particles(void)
{
  int i, index, subnode, nn, th;

  for(i = 0; i < NTopleaves; i++)
    {
      index = DomainNodeIndex[i];

      DomainMoment[i].mass = 0;
      DomainMoment[i].s[0] = Nodes[index].center[0];
      DomainMoment[i].s[1] = Nodes[index].center[1];
      DomainMoment[i].s[2] = Nodes[index].center[2];
    }

  for(i = 0; i < NTopleaves; i++)
    {
      if(i < DomainMyStart || i > DomainMyLast)
	{
	  th = All.MaxPart;	/* select index of first node in tree */

	  while(1)
	    {
	      if(th >= All.MaxPart)	/* we are dealing with an internal node */
		{
		  if(th >= All.MaxPart + MaxNodes)
		    endrun(888);	/* this can't be */

		  subnode = 0;
		  if(DomainMoment[i].s[0] > Nodes[th].center[0])
		    subnode += 1;
		  if(DomainMoment[i].s[1] > Nodes[th].center[1])
		    subnode += 2;
		  if(DomainMoment[i].s[2] > Nodes[th].center[2])
		    subnode += 4;

		  nn = Nodes[th].u.suns[subnode];

		  if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		    {
		      th = nn;
		    }
		  else
		    {
		      /* here we have found an empty slot where we can
		       * attach the pseudo particle as a leaf
		       */
		      Nodes[th].u.suns[subnode] = All.MaxPart + MaxNodes + i;

		      break;	/* done for this pseudo particle */
		    }
		}
	      else
		{
		  endrun(889);	/* this can't be */
		}
	    }
	}
    }
}


/*! this routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *
 *  Note that the bitflags-variable for each node is used to store in the
 *  lowest bits some special information: Bit 0 flags whether the node
 *  belongs to the top-level tree corresponding to the domain
 *  decomposition, while Bit 1 signals whether the top-level node is
 *  dependent on local mass.
 *
 *  If UNEQUALSOFTENINGS is set, bits 2-4 give the particle type with
 *  the maximum softening among the particles in the node, and bit 5
 *  flags whether the node contains any particles with lower softening
 *  than that.
 */
void gadgetmg2::force_update_node_recursive(int no, int sib, int father)
{
  int j, jj, p, pp, nextsib, suns[8];
  double hmax;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int maxsofttype, diffsoftflag;
#else
  double maxsoft;
#endif
#endif
  struct particle_data *pa;
  double s[3], vs[3], mass;

  if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
    {
      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will
					   overwrite one element (union!) */
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      mass = 0;
      s[0] = 0;
      s[1] = 0;
      s[2] = 0;
      vs[0] = 0;
      vs[1] = 0;
      vs[2] = 0;
      hmax = 0;
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      maxsofttype = 7;
      diffsoftflag = 0;
#else
      maxsoft = 0;
#endif
#endif

      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      /* check if we have a sibling on the same level */
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      force_update_node_recursive(p, nextsib, no);


	      if(p >= All.MaxPart)	/* an internal node or pseudo particle */
		{
		  if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
		    {
		      /* nothing to be done here because the mass of the
		       * pseudo-particle is still zero. This will be changed
		       * later.
		       */
		    }
		  else
		    {
		      mass += Nodes[p].u.d.mass;
		      s[0] += Nodes[p].u.d.mass * Nodes[p].u.d.s[0];
		      s[1] += Nodes[p].u.d.mass * Nodes[p].u.d.s[1];
		      s[2] += Nodes[p].u.d.mass * Nodes[p].u.d.s[2];
		      vs[0] += Nodes[p].u.d.mass * Extnodes[p].vs[0];
		      vs[1] += Nodes[p].u.d.mass * Extnodes[p].vs[1];
		      vs[2] += Nodes[p].u.d.mass * Extnodes[p].vs[2];

		      if(Extnodes[p].hmax > hmax)
			hmax = Extnodes[p].hmax;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
		      diffsoftflag |= (Nodes[p].u.d.bitflags >> 5) & 1;

		      if(maxsofttype == 7)
			{
			  maxsofttype = (Nodes[p].u.d.bitflags >> 2) & 7;
			}
		      else
			{
			  if(((Nodes[p].u.d.bitflags >> 2) & 7) != 7)
			    {
			      if(All.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] >
				 All.ForceSoftening[maxsofttype])
				{
				  maxsofttype = ((Nodes[p].u.d.bitflags >> 2) & 7);
				  diffsoftflag = 1;
				}
			      else
				{
				  if(All.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] <
				     All.ForceSoftening[maxsofttype])
				    diffsoftflag = 1;
				}
			    }
			}
#else
		      if(Nodes[p].maxsoft > maxsoft)
			maxsoft = Nodes[p].maxsoft;
#endif
#endif
		    }
		}
	      else		/* a particle */
		{
		  pa = &P[p];

		  mass += pa->Mass;
		  s[0] += pa->Mass * pa->Pos[0];
		  s[1] += pa->Mass * pa->Pos[1];
		  s[2] += pa->Mass * pa->Pos[2];
		  vs[0] += pa->Mass * pa->Vel[0];
		  vs[1] += pa->Mass * pa->Vel[1];
		  vs[2] += pa->Mass * pa->Vel[2];

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
		  if(maxsofttype == 7)
		    {
		      maxsofttype = pa->Type;
		    }
		  else
		    {
		      if(All.ForceSoftening[pa->Type] > All.ForceSoftening[maxsofttype])
			{
			  maxsofttype = pa->Type;
			  diffsoftflag = 1;
			}
		      else
			{
			  if(All.ForceSoftening[pa->Type] < All.ForceSoftening[maxsofttype])
			    diffsoftflag = 1;
			}
		    }
#else
		  if(pa->Type == 0)
		    {
		      if(SphP[p].Hsml > maxsoft)
			maxsoft = SphP[p].Hsml;
		    }
		  else
		    {
		      if(All.ForceSoftening[pa->Type] > maxsoft)
			maxsoft = All.ForceSoftening[pa->Type];
		    }
#endif
#endif
		  if(pa->Type == 0)
		    if(SphP[p].Hsml > hmax)
		      hmax = SphP[p].Hsml;
		}
	    }
	}


      if(mass)
	{
	  s[0] /= mass;
	  s[1] /= mass;
	  s[2] /= mass;
	  vs[0] /= mass;
	  vs[1] /= mass;
	  vs[2] /= mass;
	}
      else
	{
	  s[0] = Nodes[no].center[0];
	  s[1] = Nodes[no].center[1];
	  s[2] = Nodes[no].center[2];
	}

      Nodes[no].u.d.s[0] = s[0];
      Nodes[no].u.d.s[1] = s[1];
      Nodes[no].u.d.s[2] = s[2];
      Nodes[no].u.d.mass = mass;


#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      Nodes[no].u.d.bitflags = 4 * maxsofttype + 32 * diffsoftflag;
#else
      Nodes[no].u.d.bitflags = 0;
      Nodes[no].maxsoft = maxsoft;
#endif
#else
      Nodes[no].u.d.bitflags = 0;
#endif


      Extnodes[no].vs[0] = vs[0];
      Extnodes[no].vs[1] = vs[1];
      Extnodes[no].vs[2] = vs[2];
      Extnodes[no].hmax = hmax;

      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father = father;
    }
  else				/* single particle or pseudo particle */
    {
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      if(no < All.MaxPart)	/* only set it for single particles */
	Father[no] = father;
    }

}



/*! This function updates the multipole moments of the pseudo-particles
 *  that represent the mass distribution on different CPUs. For that
 *  purpose, it first exchanges the necessary data, and then updates the
 *  top-level tree accordingly. The detailed implementation of these two
 *  tasks is done in separate functions.
 */
void gadgetmg2::force_update_pseudoparticles(void)
{
  force_exchange_pseudodata();

  force_treeupdate_pseudos();
}



/*! This function communicates the values of the multipole moments of the
 *  top-level tree-nodes of the domain grid.  This data can then be used to
 *  update the pseudo-particles on each CPU accordingly.
 */
void gadgetmg2::force_exchange_pseudodata(void)
{
  int i, no;
#ifndef NOMPI
  MPI_Status status;
#endif
  int level, sendTask, recvTask;

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      /* read out the multipole moments from the local base cells */
      DomainMoment[i].s[0] = Nodes[no].u.d.s[0];
      DomainMoment[i].s[1] = Nodes[no].u.d.s[1];
      DomainMoment[i].s[2] = Nodes[no].u.d.s[2];
      DomainMoment[i].vs[0] = Extnodes[no].vs[0];
      DomainMoment[i].vs[1] = Extnodes[no].vs[1];
      DomainMoment[i].vs[2] = Extnodes[no].vs[2];
      DomainMoment[i].mass = Nodes[no].u.d.mass;
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      DomainMoment[i].bitflags = Nodes[no].u.d.bitflags;
#else
      DomainMoment[i].maxsoft = Nodes[no].maxsoft;
#endif
#endif
    }

  /* share the pseudo-particle data accross CPUs */

#ifndef NOMPI
  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&DomainMoment[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(struct DomainNODE),
		     MPI_BYTE, recvTask, TAG_DMOM,
		     &DomainMoment[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(struct DomainNODE),
		     MPI_BYTE, recvTask, TAG_DMOM, GADGET_WORLD, &status);
    }
#endif
}

/*! This function updates the top-level tree after the multipole moments of
 *  the pseudo-particles have been updated.
 */
void gadgetmg2::force_treeupdate_pseudos(void)
{
  int i, k, no;
  double sold[3], vsold[3], snew[3], vsnew[3], massold, massnew, mm;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int maxsofttype, diffsoftflag;
#else
  double maxsoft;
#endif
#endif

  for(i = 0; i < NTopleaves; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];

	for(k = 0; k < 3; k++)
	  {
	    sold[k] = Nodes[no].u.d.s[k];
	    vsold[k] = Extnodes[no].vs[k];
	  }
	massold = Nodes[no].u.d.mass;

	for(k = 0; k < 3; k++)
	  {
	    snew[k] = DomainMoment[i].s[k];
	    vsnew[k] = DomainMoment[i].vs[k];
	  }
	massnew = DomainMoment[i].mass;


#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	maxsofttype = (DomainMoment[i].bitflags >> 2) & 7;
	diffsoftflag = (DomainMoment[i].bitflags >> 5) & 1;
#else
	maxsoft = DomainMoment[i].maxsoft;
#endif
#endif
	do
	  {
	    mm = Nodes[no].u.d.mass + massnew - massold;
	    for(k = 0; k < 3; k++)
	      {
		if(mm > 0)
		  {
		    Nodes[no].u.d.s[k] =
		      (Nodes[no].u.d.mass * Nodes[no].u.d.s[k] + massnew * snew[k] - massold * sold[k]) / mm;
		    Extnodes[no].vs[k] =
		      (Nodes[no].u.d.mass * Extnodes[no].vs[k] + massnew * vsnew[k] -
		       massold * vsold[k]) / mm;
		  }
	      }
	    Nodes[no].u.d.mass = mm;


#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	    diffsoftflag |= (Nodes[no].u.d.bitflags >> 5) & 1;

	    if(maxsofttype == 7)
	      maxsofttype = (Nodes[no].u.d.bitflags >> 2) & 7;
	    else
	      {
		if(((Nodes[no].u.d.bitflags >> 2) & 7) != 7)
		  {
		    if(All.ForceSoftening[((Nodes[no].u.d.bitflags >> 2) & 7)] >
		       All.ForceSoftening[maxsofttype])
		      {
			maxsofttype = ((Nodes[no].u.d.bitflags >> 2) & 7);
			diffsoftflag = 1;
		      }
		    else
		      {
			if(All.ForceSoftening[((Nodes[no].u.d.bitflags >> 2) & 7)] <
			   All.ForceSoftening[maxsofttype])
			  diffsoftflag = 1;
		      }
		  }
	      }

	    Nodes[no].u.d.bitflags = (Nodes[no].u.d.bitflags & 3) + 4 * maxsofttype + 32 * diffsoftflag;
#else
	    if(Nodes[no].maxsoft < maxsoft)
	      Nodes[no].maxsoft = maxsoft;
	    maxsoft = Nodes[no].maxsoft;
#endif
#endif
	    no = Nodes[no].u.d.father;

	  }
	while(no >= 0);
      }
}



/*! This function flags nodes in the top-level tree that are dependent on
 *  local particle data.
 */
void gadgetmg2::force_flag_localnodes(void)
{
  int no, i;

  /* mark all top-level nodes */

  for(i = 0; i < NTopleaves; i++)
    {
      no = DomainNodeIndex[i];

      while(no >= 0)
	{
	  if((Nodes[no].u.d.bitflags & 1))
	    break;

	  Nodes[no].u.d.bitflags |= 1;

	  no = Nodes[no].u.d.father;
	}
    }

  /* mark top-level nodes that contain local particles */

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      /*
         if(DomainMoment[i].mass > 0)
       */
      {
	no = DomainNodeIndex[i];

	while(no >= 0)
	  {
	    if((Nodes[no].u.d.bitflags & 2))
	      break;

	    Nodes[no].u.d.bitflags |= 2;

	    no = Nodes[no].u.d.father;
	  }
      }
    }
}



/*! This function updates the side-length of tree nodes in case the tree is
 *  not reconstructed, but only drifted.  The grouping of particles to tree
 *  nodes is not changed in this case, but some tree nodes may need to be
 *  enlarged because particles moved out of their original bounds.
 */
void gadgetmg2::force_update_len(void)
{
  int i, no;
#ifndef NOMPI
  MPI_Status status;
#endif
  int level, sendTask, recvTask;

  force_update_node_len_local();

  /* first update the side-lengths of all local nodes */
  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      DomainTreeNodeLen[i] = Nodes[no].len;
    }

#ifndef NOMPI
  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&DomainTreeNodeLen[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(double),
		     MPI_BYTE, recvTask, TAG_NODELEN,
		     &DomainTreeNodeLen[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(double),
		     MPI_BYTE, recvTask, TAG_NODELEN, GADGET_WORLD, &status);
    }
#endif

  /* Finally, we update the top-level tree. */
  force_update_node_len_toptree();
}


/*! This function recursively enlarges nodes such that they always contain
 *  all their daughter nodes and daughter particles.
 */
void gadgetmg2::force_update_node_len_local(void)
{
  int i, p, k, no;
  double dist, distmax;

  for(i = 0; i < NumPart; i++)
    {
      no = Father[i];

      for(k = 0, distmax = 0; k < 3; k++)
	{
	  dist = P[i].Pos[k] - Nodes[no].center[k];
	  if(dist < 0)
	    dist = -dist;
	  if(dist > distmax)
	    distmax = dist;
	}

      if(distmax + distmax > Nodes[no].len)
	{
	  Nodes[no].len = distmax + distmax;
	  p = Nodes[no].u.d.father;

	  while(p >= 0)
	    {
	      distmax = Nodes[p].center[0] - Nodes[no].center[0];
	      if(distmax < 0)
		distmax = -distmax;
	      distmax = distmax + distmax + Nodes[no].len;

	      if(0.999999 * distmax > Nodes[p].len)
		{
		  Nodes[p].len = distmax;
		  no = p;
		  p = Nodes[p].u.d.father;
		}
	      else
		break;
	    }
	}
    }
}


/*! This function recursively enlarges nodes of the top-level tree such
 *  that they always contain all their daughter nodes.
 */
void gadgetmg2::force_update_node_len_toptree(void)
{
  int i, no, p;
  double distmax;

  for(i = 0; i < NTopleaves; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];

	if(Nodes[no].len < DomainTreeNodeLen[i])
	  Nodes[no].len = DomainTreeNodeLen[i];

	p = Nodes[no].u.d.father;

	while(p >= 0)
	  {
	    distmax = Nodes[p].center[0] - Nodes[no].center[0];
	    if(distmax < 0)
	      distmax = -distmax;
	    distmax = distmax + distmax + Nodes[no].len;

	    if(0.999999 * distmax > Nodes[p].len)
	      {
		Nodes[p].len = distmax;
		no = p;
		p = Nodes[p].u.d.father;
	      }
	    else
	      break;
	  }
      }
}




/*! This function updates the hmax-values in tree nodes that hold SPH
 *  particles. These values are needed to find all neighbors in the
 *  hydro-force computation.  Since the Hsml-values are potentially changed
 *  in the SPH-denity computation, force_update_hmax() should be carried
 *  out just before the hydrodynamical SPH forces are computed, i.e. after
 *  density().
 */
void gadgetmg2::force_update_hmax(void)
{
  int i, no;
#ifndef NOMPI
  MPI_Status status;
#endif
  int level, sendTask, recvTask;

  force_update_node_hmax_local();

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      DomainHmax[i] = Extnodes[no].hmax;
    }

  /* share the hmax-data of the pseudo-particles accross CPUs */

#ifndef NOMPI
  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;
      if(recvTask < NTask)
	MPI_Sendrecv(&DomainHmax[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(double),
		     MPI_BYTE, recvTask, TAG_HMAX,
		     &DomainHmax[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(double),
		     MPI_BYTE, recvTask, TAG_HMAX, GADGET_WORLD, &status);
    }
#endif


  force_update_node_hmax_toptree();
}

/*! This routine updates the hmax-values of local tree nodes.
 */
void gadgetmg2::force_update_node_hmax_local(void)
{
  int i, p, no;

  for(i = 0; i < N_gas; i++)
    {

      no = Father[i];

      if(SphP[i].Hsml > Extnodes[no].hmax)
	{

	  Extnodes[no].hmax = SphP[i].Hsml;
	  p = Nodes[no].u.d.father;

	  while(p >= 0)
	    {
	      if(Extnodes[no].hmax > Extnodes[p].hmax)
		{
		  Extnodes[p].hmax = Extnodes[no].hmax;
		  no = p;
		  p = Nodes[p].u.d.father;
		}
	      else
		break;
	    }
	}

    }
}




/*! This function recursively sets the hmax-values of the top-level tree.
 */
void gadgetmg2::force_update_node_hmax_toptree(void)
{

  int i, no, p;


  for(i = 0; i < NTopleaves; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];

	if(Extnodes[no].hmax < DomainHmax[i])
	  Extnodes[no].hmax = DomainHmax[i];

	p = Nodes[no].u.d.father;

	while(p >= 0)
	  {
	    if(Extnodes[no].hmax > Extnodes[p].hmax)
	      {
		Extnodes[p].hmax = Extnodes[no].hmax;
		no = p;
		p = Nodes[p].u.d.father;
	      }
	    else
	      break;
	  }
      }
}





/*! This routine computes the gravitational force for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 */
int gadgetmg2::force_treeevaluate(int target, int mode, double *ewaldcountsum)
{
  struct NODE *nop = 0;
  int no, ninteractions, ptype;
  double r2, dx, dy, dz, mass, r, fac, u, h, h_inv, h3_inv;
  double acc_x, acc_y, acc_z, pos_x, pos_y, pos_z, aold;
#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  ninteractions = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].Type;
#else
      ptype = P[0].Type;
#endif
      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }



#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
  h_inv = 1.0 / h;
  h3_inv = h_inv * h_inv * h_inv;
#endif
  no = All.MaxPart;		/* root node */

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign */

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;

	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }
	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;

	  mass = nop->u.d.mass;
	}
      r2 = dx * dx + dy * dy + dz * dz;

      if(no < All.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an  internal node. Need to check opening criterion */
	{
	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node
						 * which does not contain
						 * local particles we can
						 * continue to do a short-cut */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }


	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      /* check in addition whether we lie inside the cell */

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
          maxsofttype = (nop->u.d.bitflags >> 2) & 7;
          if(maxsofttype == 7) /* may only occur for zero mass top-level nodes */
            {
              if(mass > 0)
                endrun(986);
              no = nop->u.d.nextnode;
              continue;
            }
          else
            {
              if(h < All.ForceSoftening[maxsofttype])
                {
                  h = All.ForceSoftening[maxsofttype];
                  if(r2 < h * h)
                    {
                      if(((nop->u.d.bitflags >> 5) & 1))	/* bit-5 signals that there are particles of different softening in the node */
                        {
                          no = nop->u.d.nextnode;
                          continue;
                        }
                    }
                }
            }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif

	  no = nop->u.d.sibling;	/* ok, node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      r = sqrt(r2);

      if(r >= h)
	fac = mass / (r2 * r);
      else
	{
#ifdef UNEQUALSOFTENINGS
	  h_inv = 1.0 / h;
	  h3_inv = h_inv * h_inv * h_inv;
#endif
	  u = r * h_inv;
	  if(u < 0.5)
	    fac = mass * h3_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      mass * h3_inv * (21.333333333333 - 48.0 * u +
			       38.4 * u * u - 10.666666666667 * u * u * u - 0.066666666667 / (u * u * u));
	}

      acc_x += dx * fac;
      acc_y += dy * fac;
      acc_z += dz * fac;

      ninteractions++;
    }


  /* store result at the proper place */
  if(mode == 0)
    {
      P[target].GravAccel[0] = acc_x;
      P[target].GravAccel[1] = acc_y;
      P[target].GravAccel[2] = acc_z;
      P[target].GravCost = ninteractions;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
      GravDataResult[target].w.Ninteractions = ninteractions;
    }

  return ninteractions;
}



/*! This routine computes the gravitational potential by walking the
 *  tree. The same opening criteria is used as for the gravitational force
 *  walk.
 */
void gadgetmg2::force_treeevaluate_potential(int target, int mode)
{
  struct NODE *nop = 0;
  int no, ptype;
  double r2, dx, dy, dz, mass, r, u, h, h_inv, wp;
  double pot, pos_x, pos_y, pos_z, aold;
#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif

  pot = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].Type;
#else
      ptype = P[0].Type;
#endif
      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }


#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
  h_inv = 1.0 / h;
#endif
  no = All.MaxPart;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign */

	  dx = P[no].Pos[0] - pos_x;
	  dy = P[no].Pos[1] - pos_y;
	  dz = P[no].Pos[2] - pos_z;
	  mass = P[no].Mass;
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];
	  dx = nop->u.d.s[0] - pos_x;
	  dy = nop->u.d.s[1] - pos_y;
	  dz = nop->u.d.s[2] - pos_z;
	  mass = nop->u.d.mass;
	}
      r2 = dx * dx + dy * dy + dz * dz;

      if(no < All.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an internal node. Need to check opening criterion */
	{
	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node
						 * which does not contain
						 * local particles we can make
						 * a short-cut
						 */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2 * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(mass * nop->len * nop->len > r2 * r2 * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
          maxsofttype = (nop->u.d.bitflags >> 2) & 7;
          if(maxsofttype == 7) /* may only occur for zero mass top-level nodes */
            {
              if(mass > 0)
                endrun(988);
              no = nop->u.d.nextnode;
              continue;
            }
          else
            {
              if(h < All.ForceSoftening[maxsofttype])
                {
                  h = All.ForceSoftening[maxsofttype];
                  if(r2 < h * h)
                    {
                      if(((nop->u.d.bitflags >> 5) & 1))	/* bit-5 signals that there are particles of different softening in the node */
                        {
                          no = nop->u.d.nextnode;
                          continue;
                        }
                    }
                }
            }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2 < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif

	  no = nop->u.d.sibling;	/* node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      r = sqrt(r2);

      if(r >= h)
	pot -= mass / r;
      else
	{
#ifdef UNEQUALSOFTENINGS
	  h_inv = 1.0 / h;
#endif
	  u = r * h_inv;

	  if(u < 0.5)
	    wp = -2.8 + u * u * (5.333333333333 + u * u * (6.4 * u - 9.6));
	  else
	    wp =
	      -3.2 + 0.066666666667 / u + u * u * (10.666666666667 +
						   u * (-16.0 + u * (9.6 - 2.133333333333 * u)));

	  pot += mass * h_inv * wp;
	}
    }

  /* store result at the proper place */

  if(mode == 0)
    P[target].Potential = pot;
  else
    GravDataResult[target].u.Potential = pot;
}



/*! This function allocates the memory used for storage of the tree and of
 *  auxiliary arrays needed for tree-walk and link-lists.  Usually,
 *  maxnodes approximately equal to 0.7*maxpart is sufficient to store the
 *  tree for up to maxpart particles.
 */
void gadgetmg2::force_treeallocate(int maxnodes, int maxpart)
{
  int i;
  size_t bytes;
  double allbytes = 0;
  double u;

  MaxNodes = maxnodes;

  if(!(Nodes_base = (NODE*)malloc(bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
      endrun(3);
    }
  allbytes += bytes;

  if(!(Extnodes_base = (extNODE*)malloc(bytes = (MaxNodes + 1) * sizeof(struct extNODE))))
    {
      printf("failed to allocate memory for %d tree-extnodes (%g MB).\n", MaxNodes,
	     bytes / (1024.0 * 1024.0));
      endrun(3);
    }
  allbytes += bytes;

  Nodes = Nodes_base - All.MaxPart;
  Extnodes = Extnodes_base - All.MaxPart;

  if(!(Nextnode = (int*)malloc(bytes = (maxpart + MAXTOPNODES) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n", maxpart + MAXTOPNODES,
	     bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;

  if(!(Father = (int*)malloc(bytes = (maxpart) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;

  if(first_flag == 0)
    {
      first_flag = 1;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for BH-tree. %ld\n\n", allbytes / (1024.0 * 1024.0),
	       sizeof(struct NODE) + sizeof(struct extNODE));

      tabfac = NTAB / 3.0;

      for(i = 0; i < NTAB; i++)
	{
	  u = 3.0 / NTAB * (i + 0.5);
	  shortrange_table[i] = erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
	  shortrange_table_potential[i] = erfc(u);
	}
    }
}


/*! This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function force_treeallocate().
 */
void gadgetmg2::force_treefree(void)
{
  free(Father);
  free(Nextnode);
  free(Extnodes_base);
  free(Nodes_base);
}




/*! This function does the force computation with direct summation for the
 *  specified particle in the communication buffer. This can be useful for
 *  debugging purposes, in particular for explicit checks of the force
 *  accuracy.
 */
#ifdef FORCETEST
int gadgetmg2::force_treeevaluate_direct(int target, int mode)
{
  double epsilon;
  double h, h_inv, dx, dy, dz, r, r2, u, r_inv, fac;
  int i, ptype;
  double pos_x, pos_y, pos_z;
  double acc_x, acc_y, acc_z;

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
#ifdef UNEQUALSOFTENINGS
      ptype = GravDataGet[target].Type;
#else
      ptype = P[0].Type;
#endif
    }

  for(i = 0; i < NumPart; i++)
    {
      epsilon = dmax(All.ForceSoftening[P[i].Type], All.ForceSoftening[ptype]);

      h = epsilon;
      h_inv = 1 / h;

      dx = P[i].Pos[0] - pos_x;
      dy = P[i].Pos[1] - pos_y;
      dz = P[i].Pos[2] - pos_z;
      r2 = dx * dx + dy * dy + dz * dz;

      r = sqrt(r2);

      u = r * h_inv;

      if(u >= 1)
	{
	  r_inv = 1 / r;

	  fac = P[i].Mass * r_inv * r_inv * r_inv;
	}
      else
	{
	  if(u < 0.5)
	    fac = P[i].Mass * h_inv * h_inv * h_inv * (10.666666666667 + u * u * (32.0 * u - 38.4));
	  else
	    fac =
	      P[i].Mass * h_inv * h_inv * h_inv * (21.333333333333 -
						   48.0 * u + 38.4 * u * u -
						   10.666666666667 * u * u *
						   u - 0.066666666667 / (u * u * u));
	}

      acc_x += dx * fac;
      acc_y += dy * fac;
      acc_z += dz * fac;
    }


  if(mode == 0)
    {
      P[target].GravAccelDirect[0] = acc_x;
      P[target].GravAccelDirect[1] = acc_y;
      P[target].GravAccelDirect[2] = acc_z;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
    }


  return NumPart;
}
#endif


/*! This function dumps some of the basic particle data to a file. In case
 *  the tree construction fails, it is called just before the run
 *  terminates with an error message. Examination of the generated file may
 *  then give clues to what caused the problem.
 */
void gadgetmg2::dump_particles(void)
{
  FILE *fd;
  char buffer[200];
  int i;

  sprintf(buffer, "particles%d.dat", ThisTask);
  fd = fopen(buffer, "wb");
  my_fwrite(&NumPart, 1, sizeof(int), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Pos[0], 3, sizeof(double), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Vel[0], 3, sizeof(double), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].ID, 1, sizeof(int), fd);

  fclose(fd);
}

