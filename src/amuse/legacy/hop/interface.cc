#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <map>
#include <vector>
#include "src/kd.h"
#include "src/kd.c"
#include "src/smooth.h"
#include "src/smooth.c"
#include "src/hop.cc"
#include "worker_code.h"
#include <iostream>

#define INFORM(string) printf(string); fflush(stdout)

class AmuseParticle{

public:
    int index;
    double mass,x,y,z, density;
    int group, neighbor;
    
    AmuseParticle(int index, double mass, double x, double y, double z):index(index), mass(mass),x(x), y(y), z(z), density(-1.0), group(-1), neighbor(-1) {
    }
    AmuseParticle(const AmuseParticle & original):
    index(original.index), mass(original.mass), x(original.x), y(original.y), z(original.z), density(original.density), group(original.group), neighbor(original.neighbor) {
    }
};

typedef std::map<int, AmuseParticle *> ParticlesMap;
typedef std::map<int, AmuseParticle *>::iterator ParticlesMapIterator;

typedef std::vector<AmuseParticle *> ParticlesList;
typedef std::vector<AmuseParticle *>::iterator ParticlesListIterator;

int highest_index = 0;
ParticlesMap particlesMap;

/* hop structs */
KD gkd;
SMX gsmx;

/* hop parameters */
int nBucket,nSmooth;
float fPeriod[3];
int nDens,nHop,nMerge;
float fDensThresh;
int nMethod;

/* control variables */
int bParamInit = 1;
int bHopInit = 1;
int bDensity = 0;
int bHopDone = 0;
int bNeighborsFound = 0;

int ReadPositions(KD &kd){
  std::size_t n = particlesMap.size();
  kd->nActive = n;
  
  if (n == 0) return -1;
  
  kd->p = (PARTICLE *)malloc(kd->nActive*sizeof(PARTICLE));

  ParticlesMapIterator i;
  int c = 0;
  INFORM("Reading Positions...\n");
  for (i = particlesMap.begin(); i != particlesMap.end(); i++) {
    AmuseParticle * p = (*i).second;
    kd->p[c].fMass = p->mass; 
    kd->p[c].r[0] = p->x; 
    kd->p[c].r[1] = p->y;
    kd->p[c].r[2] = p->z;
    c++;
  }
  
  kd->fMass = 1.0/kd->nActive; /* particles have equal mass */
  return 0;
}

int ReadDensities(SMX &smx){
  ParticlesMapIterator i;
  int c = 0;
  double dens;
    
  INFORM("Reading Densities...\n");
  for (i = particlesMap.begin(); i != particlesMap.end(); i++) {
    AmuseParticle * p = (*i).second;
    dens = p->density;
    if (dens < 0) return -1;
    smx->kd->p[c].fDensity = p->density;
    c++;
  }
  return 0;
}  

int InitParam(){ /* used to set the hop parameters to default values */
  INFORM("\nInitialising Parameters...\n");
  nBucket = 16;
  nSmooth = 64;
  nDens = 64;
  nHop = -1;
  fDensThresh = -1.0;
  
  nMethod = 0;

  for (int j=0;j<3;++j) fPeriod[j] = HUGE;
  nMerge = 4;
  
  bParamInit = 0;
  return 0;
}

int InitHop(KD &kd, SMX &smx, int bDens){
  if (bParamInit) InitParam(); /* make sure the parameters have values */
  INFORM("\nInitialising Hop...\n");
  if (nHop<0) nHop=nDens;
  if (bDens==0) nSmooth = nHop+1;
  else nSmooth = nDens+1;
	/* When smSmooth() is asked for nSmooth particles, it seems to
	generally return nSmooth-1 particles, including primary itself.
	Hence, when we want nDens or nSmooth particles (including the
	primary) we should ask for one more.  By convention, I've chosen
	to have nMerge reflect the number *not* including the primary,
	so in this case we need to ask for two more! */

  kdInit(&kd,nBucket);
  if (ReadPositions(kd) == -1) return -1;
  if (nMerge > kd->nActive) return -1;
  if (nHop > kd->nActive) return -1;
  if (nDens > kd->nActive) return -1;
  if (nSmooth > kd->nActive) return -1;
  if (nBucket > kd->nActive) return -1;
  
  PrepareKD(kd);
  
  smInit(&smx,kd,nSmooth,fPeriod);
  smx->nHop = nHop;
  smx->nDens = nDens;
  smx->nMerge = nMerge;
  smx->nGroups = 0;
  smx->fDensThresh = fDensThresh;
  
  if (bDens == 0) {
    if(ReadDensities(smx) == -1) return -1;
  }
  
  INFORM("Building Tree...\n");
  kdBuildTree(kd);

  bHopDone = 0;
  return 0;
}


/* parameter setters */
int set_nBucket(int value){
  if (bParamInit) InitParam();
  nBucket = value;
  bHopInit = 1;
  return 0;
}
int set_nDens(int value){
  if (bParamInit) InitParam();
  nDens = value;
  bHopInit = 1;
  return 0;
}
int set_nHop(int value){
  if (bParamInit) InitParam();
  if (value < nMerge +1) return -1;
  nHop = value;
  bHopInit = 1;
  return 0;
}
int set_fDensThresh(double value){
  if (bParamInit) InitParam();
  fDensThresh = value;
  bHopInit = 1;
  return 0;
}
int set_fPeriod(double x, double y, double z){
  if (bParamInit) InitParam();
  fPeriod[0] = x;
  fPeriod[1] = y;
  fPeriod[2] = z;
  bHopInit = 1;
  return 0;
}
int set_nMerge(int value){
  if (bParamInit) InitParam();
  if (value > nHop -1) return -1;
  nMerge = value;
  bHopInit = 1;
  return 0;
}
int set_density_method(int value){
  if (bParamInit) InitParam();
  switch (value) {
    case 0:
    case 1: 
    case 2:
      nMethod = value;
      return 0;
    default:
      return -1;
   }
}

/* densities */
int calculate_densities(){
  KD kd;
  SMX smx;
  ParticlesList particles_as_list;
  
  if(InitHop(kd, smx, 1) == -1) return -1;

  INFORM("Finding Densities...\n");
  switch (nMethod) {
    case 1: smSmooth(smx,smDensitySym); break;
    case 2: smSmooth(smx,smDensityTH); break;
    default: smSmooth(smx,smDensity);
  }
  INFORM("Storing Results...");
  ParticlesMapIterator i;
  int c = 0;
  for (i = particlesMap.begin(); i != particlesMap.end(); i++) {
      AmuseParticle * p = (*i).second;
      particles_as_list.push_back(p);
  }
  
  for (c = 0; c < particles_as_list.size(); c++) {
      AmuseParticle * p = particles_as_list[kd->p[c].iOrder];
      p->mass=kd->p[c].fMass;       
      p->x=kd->p[c].r[0]; 
      p->y=kd->p[c].r[1];
      p->z=kd->p[c].r[2];
      p->density = kd->p[c].fDensity;
  }

  gkd = kd;
  gsmx = smx;

  INFORM("Done!\n");
  
  return 0;
}

/* actual hop */
int do_hop(){
  KD kd;
  SMX smx;
  if(InitHop(kd, smx, 0) == -1) return -1;

  INFORM("Finding Densest Neighbors...\n");
  if (nHop>=nSmooth) {
    nSmooth = nHop+1;
    ReSizeSMX(smx,nSmooth);
  }
  smSmooth(smx,smHop);

  int c = 0;
  ParticlesMapIterator i;
  for (i = particlesMap.begin(); i != particlesMap.end(); i++) {
    AmuseParticle * p = (*i).second;
    p->neighbor = -1 - smx->kd->p[c].iHop;
    c++;
  }

	INFORM("Grouping...\n");
	FindGroups(smx);
  SortGroups(smx);

  INFORM("Merging Groups...\n");
  MergeGroupsHash(smx);

	kdOrder(kd);
	
	INFORM("Storing Results...");
  c = 0;
  for (i = particlesMap.begin(); i != particlesMap.end(); i++) {
      AmuseParticle * p = (*i).second;
      p->group = smx->kd->p[c].iHop;
      c++;
  }
  
  gkd = kd;
  gsmx = smx;
  
  bHopDone = 1;
  
  INFORM("Done!\n");
  
  return 0;
}

int new_particle(int * index_of_the_particle, double mass, double x, double y, double z) {
  *index_of_the_particle = highest_index;
  AmuseParticle * p = new AmuseParticle(highest_index, mass, x, y, z);
  particlesMap[highest_index] = p;
  highest_index++;
  bHopInit = 1;
  return 0;
}

int delete_particle(int index_of_the_particle) {
  if(index_of_the_particle > highest_index) return -1;
  particlesMap.erase(index_of_the_particle);
  bHopInit = 1;
  return 0;
}

int get_number_of_particles(int * value) {
  *value = (int) particlesMap.size();
  return 0;
}

int set_density(int index_of_the_particle, double density) {
  if(index_of_the_particle > highest_index) return -1;
  AmuseParticle * p = particlesMap[index_of_the_particle];
  p->density = density;
  bHopInit = 1;
  return 0;
}

int get_density(int index_of_the_particle, double * density) {
  if(index_of_the_particle > highest_index) return -1;
  AmuseParticle * p = particlesMap[index_of_the_particle];
  if (p->density < 0.0) return -1;
  *density = p->density;
  return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z) {
  if(index_of_the_particle > highest_index) return -1;
  AmuseParticle * p = particlesMap[index_of_the_particle];
  p->x = x;
  p->y = y;
  p->z = z;
  bHopInit = 1;
  return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, double * z) {
  if(index_of_the_particle > highest_index) return -1; 
  AmuseParticle * p = particlesMap[index_of_the_particle];
  *x = p->x;
  *y = p->y;
  *z = p->z;
  return 0;
}

int get_mass(int index_of_the_particle, double * mass) {
  if(index_of_the_particle > highest_index) return -1; 
  AmuseParticle * p = particlesMap[index_of_the_particle];
  *mass = p->mass;
  return 0;
}

int get_densest_neighbor(int index_of_the_particle, int * index_of_densest_neighbor) {
  if(index_of_the_particle > highest_index) return -1;
  AmuseParticle * p = particlesMap[index_of_the_particle];
  if (p->neighbor == -1 ) return -1;
  *index_of_densest_neighbor = p->neighbor;
  return 0;
}

int get_group_id(int index_of_the_particle, int * group_id) {
  if (!bHopDone) return -1;
  if(index_of_the_particle > highest_index) return -1;
  AmuseParticle * p = particlesMap[index_of_the_particle];
  if (p->group == -1){ *group_id = p->group; return -1; }
  *group_id = p->group;
  return 0;
}

int get_densest_particle_in_group(int group_id, int * index_of_the_particle) {
  if (!bHopDone) return -1;
  if(group_id > gsmx->nGroups) return -1;
  *index_of_the_particle = gsmx->densestingroup[group_id];
  return 0;
}

int get_number_of_particles_in_group(int group_id, int * number_of_particles) {
  if (!bHopDone) return -1;
  if(group_id >= gsmx->nGroups) return -1;
  *number_of_particles = gsmx->nmembers[group_id];
  return 0;
}

int get_average_boundary_density_of_groups(int first_group_id, int second_group_id, double * boundary_density) {
  if (!bHopDone) return -1;
  if(first_group_id > gsmx->nGroups) return -1;
  if(second_group_id > gsmx->nGroups) return -1;

  Boundary *hp;
  hp=gsmx->hash;
  for (int j=0;j<gsmx->nHashLength; j++) {
    if (hp->nGroup1 == first_group_id && hp->nGroup2 == second_group_id) {
      *boundary_density = hp->fDensity;
      return 0;
		}
		else if (hp->nGroup1 == second_group_id && hp->nGroup2 == first_group_id) {
		  *boundary_density = hp->fDensity;
		  return 0;
		}
		hp++;
  }
//  printf("No boundary exists between groups %4d and %4d.\n", first_group_id, second_group_id);
  return -2;
}

int get_number_of_groups(int * number_of_groups) {
  if (!bHopDone) return -1;
  *number_of_groups = gsmx->nGroups;
  return 0;
}

int get_number_of_particles_outside_groups(int * number_of_particles) {
  if (!bHopDone) return -1;
  *number_of_particles = gsmx->nmembers[gsmx->nGroups];
  return 0;
}

int show_parameters(){
  if (bParamInit) InitParam();
  printf("Printing Parameters...\n");
  printf("nBucket:\t%d\nfPeriod x:\t%f\nfPeriod y:\t%f\nfPeriod z:\t%f\nnDens:\t%d\nnHop:\t%d\nnMerge:\t%d\nfDensThresh:\t%f\nnMethod:\t%d\n",
  nBucket, fPeriod[0],fPeriod[1],fPeriod[2], nDens, nHop, nMerge, fDensThresh, nMethod);
  return 0;
}








