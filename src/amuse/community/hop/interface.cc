#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <map>
#include <vector>
#include <algorithm>
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


bool debug = true;

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

double outer_densthresh = -1.0;
double saddle_densthresh = -1.0;
double peak_densthresh = -1.0;
int relative_saddle_density_threshold = 0;
double saddle_density_threshold_factor = 0.80;

int weakly_connected; // used when merging to store which of the two is less strong connected to another one
int nGroups_before_regroup;
int nGroups_after_regroup = -1;

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


/* parameters */
int set_nBucket(int value){
  if (bParamInit) InitParam();
  
  nBucket = value;
  bHopInit = 1;
  return 0;
}
int get_nBucket(int * value){
  if (bParamInit) InitParam();
  
  *value = nBucket;
  return 0;
}

int set_nDens(int value){
  if (bParamInit) InitParam();
  
  nDens = value;
  bHopInit = 1;
  return 0;
}
int get_nDens(int * value){
  if (bParamInit) InitParam();
  
  *value = nDens;
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
int get_nHop(int * value){

  *value = nHop;
  
  return 0;
}

int set_fDensThresh(double value){
  if (bParamInit) InitParam();
  fDensThresh = value;
  bHopInit = 1;
  return 0;
}
int get_fDensThresh(double * value){
  *value = fDensThresh;
  return 0;
}

int set_saddle_densthresh(double value){
  saddle_densthresh = value;
  return 0;
}
int get_saddle_densthresh(double *value){
  *value = saddle_densthresh;
  return 0;
}

int set_peak_densthresh(double value){
  peak_densthresh = value;
  return 0;
}
int get_peak_densthresh(double *value){
  *value = peak_densthresh;
  return 0;
}

int set_saddle_density_threshold_factor(double value){
  saddle_density_threshold_factor = value;
  return 0;
}
int get_saddle_density_threshold_factor(double *value){
  *value = saddle_density_threshold_factor;
  return 0;
}
int set_relative_saddle_density_threshold(int value){
    relative_saddle_density_threshold = value;
    return 0;
}
int get_relative_saddle_density_threshold(int *value){
    *value = relative_saddle_density_threshold;
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

int get_fPeriod(double * x, double * y, double * z){
  *x = fPeriod[0];
  *y = fPeriod[1];
  *z = fPeriod[2];
  return 0;
}

int set_nMerge(int value){
  if (bParamInit) InitParam();
  if (value > nHop -1) return -1;
  nMerge = value;
  bHopInit = 1;
  return 0;
}
int get_nMerge(int *value){
  *value = nMerge;
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

int get_density_method(int * value){
   *value = nMethod;
   return 0;
}

int initialize_code()
{
    return InitParam();
}

int cleanup_code()
{
    return 0;
}

int commit_parameters() {
    if (outer_densthresh < 0.0) {
        outer_densthresh = fDensThresh;
    }
    if (saddle_densthresh < 0.0) {
        saddle_densthresh = 1.75 * outer_densthresh;
    }
    if (peak_densthresh < 0.0) {
        peak_densthresh = saddle_densthresh > 2.0*outer_densthresh ? saddle_densthresh : 2.0*outer_densthresh;
    }
    return 0;
}

int recommit_parameters() {
    return commit_parameters();
}

/* densities */
int calculate_densities() {
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

bool compare_boundaries(const Boundary &a, const Boundary &b) {
    return a.fDensity < b.fDensity;
}
bool is_weakly_connected(const Boundary &boundary) {
    return (boundary.nGroup1 == weakly_connected) || (boundary.nGroup2 == weakly_connected);
}

void merge_groups(int g1, int g2, 
        std::vector<std::vector<Boundary>*> &merge_priorities, 
        std::vector<double> &group_peak_density, 
        std::vector<int> &id_from_index, int index_from_id[]) {
    std::vector<Boundary>::iterator it1, it2;
    int target_group1, target_group2;
    bool have_common_neighbour;
    int index1 = index_from_id[g1];
    int index2 = index_from_id[g2];
    
    if (debug) {
        for (int j=0; j < merge_priorities.size(); j++) {
            printf("\nGROUP %d, boundaries:\n", id_from_index[j]);
            for (it1 = merge_priorities[j]->begin(); it1 != merge_priorities[j]->end(); it1++) {
                printf("%d %d %f\n", it1->nGroup1, it1->nGroup2, it1->fDensity);
            }
        }
    }
    for (it2 = merge_priorities[index2]->begin(); it2 < merge_priorities[index2]->end(); it2++) {
        have_common_neighbour = false;
        target_group2 = (it2->nGroup1 == g2) ? it2->nGroup2 : it2->nGroup1;
        if (target_group2 == g1) { // This is the connection between the two merging groups: no longer needed
            continue;
        }
        
        // Check if group 1 already is connected to target_group2, i.e. they have a neighbour in common:
        it1 = merge_priorities[index1]->begin();
        while (!have_common_neighbour && (it1 < merge_priorities[index1]->end())) {
            target_group1 = (it1->nGroup1 == g1) ? it1->nGroup2 : it1->nGroup1;
            if (target_group1 == target_group2) { // merging groups share a neighbour
                have_common_neighbour = true;
            } else {
                it1++;
            }
        }
        if (have_common_neighbour) { // Yes: replace with best connection
            if (it1->fDensity < it2->fDensity) {
                it1->fDensity = it2->fDensity;
                // make sure the ordering is maintained:
                while ((it1+1 < merge_priorities[index1]->end()) && (it1->fDensity > (it1+1)->fDensity)) {
                    std::iter_swap(it1, it1+1);
                    it1++;
                }
                weakly_connected = g1;
            } else {
                weakly_connected = g2;
            }
            // also remove the weaker connection from the merge priority vector of target_group:
            merge_priorities[index_from_id[target_group1]]->erase(
                std::remove_if(
                    merge_priorities[index_from_id[target_group1]]->begin(),
                    merge_priorities[index_from_id[target_group1]]->end(),
                    is_weakly_connected
                ), merge_priorities[index_from_id[target_group1]]->end()
            );
        } else { // No: insert this new connection
            merge_priorities[index1]->insert(
                std::lower_bound(
                    merge_priorities[index1]->begin(),
                    merge_priorities[index1]->end(),
                    *it2,
                    compare_boundaries
                ), 
                *it2
            );
        }
    }
    // All connections are now stored on merge_priorities[index1], so we can remove merge_priorities[index2]
    merge_priorities.erase(merge_priorities.begin() + index2);
    // make sure the other vectors still match
    if (group_peak_density[index2] > group_peak_density[index1]) {
        group_peak_density[index1] = group_peak_density[index2];
    }
    group_peak_density.erase(group_peak_density.begin() + index2);
    id_from_index.erase(id_from_index.begin() + index2);
    index_from_id[g2] = -1;
    for (int i = g2 + 1; i < nGroups_before_regroup; i++) {
        index_from_id[i]--;
    }
    
    // Now update all boundaries in merge_priorities from g2 to g1
    std::vector<std::vector<Boundary>*>::iterator it_group;
    for (it_group = merge_priorities.begin(); it_group < merge_priorities.end(); it_group++) {
        for (it1 = (*it_group)->begin(); it1 < (*it_group)->end(); it1++) {
            if (it1->nGroup1 == g2) {
                it1->nGroup1 = g1;
            } else if (it1->nGroup2 == g2) {
                it1->nGroup2 = g1;
            }
        }
    }
}

int regroup() {
    int nGroups = gsmx->nGroups;
    
    std::vector<std::vector<Boundary>*> merge_priorities(nGroups);
    std::vector<Boundary>::iterator it;
    std::vector<double> group_peak_density(nGroups);
    
    // to go from index in the vectors (which decrease in size due to merging) to group id
    std::vector<int> id_from_index(nGroups);
    int *index_from_id = new int[nGroups]; // and vice versa
    int *new_group_id = new int[nGroups]; // to go from old to new group id 
    Boundary *boundary;
    double d_tmp, d_saddle_tmp;
    int j, i;
    
    nGroups_before_regroup = nGroups;
    
    // Create a vector of (pointers to) vectors containing, for each group, the boundaries 
    // it has with other groups, sorted in ascending order of saddle density.
    for (j=0; j < nGroups; j++) {
        merge_priorities[j] = new std::vector<Boundary>;
        
        group_peak_density[j] = gsmx->kd->p[gsmx->densestingroup[j]].fDensity;
        index_from_id[j] = j;
        id_from_index[j] = j;
        new_group_id[j] = j;
    }
    for (j=0, boundary=gsmx->hash; j < gsmx->nHashLength; j++, boundary++){
        if (boundary->nGroup1 >= 0) {
            merge_priorities[boundary->nGroup1]->insert(
                std::lower_bound(
                    merge_priorities[boundary->nGroup1]->begin(),
                    merge_priorities[boundary->nGroup1]->end(),
                    *boundary,
                    compare_boundaries
                ), 
                *boundary
            );
            merge_priorities[boundary->nGroup2]->insert(
                std::lower_bound(
                    merge_priorities[boundary->nGroup2]->begin(),
                    merge_priorities[boundary->nGroup2]->end(),
                    *boundary,
                    compare_boundaries
                ), 
                *boundary
            );
        }
    }
    if (debug) {
        for (j=0; j < nGroups; j++) {
            printf("\nGROUP %d, boundaries:\n", j);
            for (it = merge_priorities[j]->begin(); it != merge_priorities[j]->end(); it++) {
                printf("%d %d %f\n", it->nGroup1, it->nGroup2, it->fDensity);
            }
        }
    }
    
    // Merge proper groups
    int this_group, target_group, index;
    j=0;
    while (j < nGroups) {
    //~for (j=0; j < nGroups; j++) { // For each group...
        this_group = id_from_index[j];
        if (debug) printf("\nMerge proper groups, with group1: %d\n", this_group);
        if (group_peak_density[j] > peak_densthresh) { // is this group a proper group?
            it = merge_priorities[j]->begin();
            while (it != merge_priorities[j]->end()) {
                target_group = (it->nGroup1 == this_group) ? it->nGroup2 : it->nGroup1;
                index = index_from_id[target_group];
                if (group_peak_density[index] > peak_densthresh) { // is the target group also a proper group?
                    d_tmp = it->fDensity;
                    it = merge_priorities[j]->erase(it);
                    if (relative_saddle_density_threshold) {
                        d_saddle_tmp = saddle_density_threshold_factor * std::min(group_peak_density[j], group_peak_density[index]);
                    } else {
                        d_saddle_tmp = saddle_densthresh;
                    }
                    if (d_tmp > d_saddle_tmp) { // and boundary density high enough? ==> merge groups:
                        if (debug) printf("merging proper groups %d and %d (original size %d)\n", this_group, target_group, merge_priorities[j]->size());
                        merge_groups(
                            this_group, target_group, 
                            merge_priorities, group_peak_density, id_from_index, index_from_id);
                        if (debug) printf("merged proper groups %d and %d (new size %d)\n", this_group, target_group, merge_priorities[j]->size());
                        nGroups--;
                        for (i = 0; i < nGroups_before_regroup; i++) {
                            if (new_group_id[i] == target_group) {
                                new_group_id[i] = this_group;
                            }
                        }
                        it = merge_priorities[j]->begin(); // merge_priorities[j] has changed -> start over
                    }
                } else { // the target group is not a proper group. Deal with fringe groups later:
                    it++;
                }
            }
        }
        j++;
    }
    
    // Merge fringe groups
    bool changes = true;
    Boundary next_boundary;
    while (changes) {
        changes = false;
        for (j=0; j < nGroups; j++) {
            if (merge_priorities[j]->size() > 0) {
                this_group = id_from_index[j];
                next_boundary = merge_priorities[j]->back(); // top merge priority for this group
                target_group = (next_boundary.nGroup1 == this_group) ? next_boundary.nGroup2 : next_boundary.nGroup1;
                next_boundary = merge_priorities[index_from_id[target_group]]->back(); // top merge priority for this group's top target
                int targets_target_group = (next_boundary.nGroup1 == target_group) ? next_boundary.nGroup2 : next_boundary.nGroup1;
                
                // Merging not allowed between proper groups separated by density below saddle density
                if ((group_peak_density[index_from_id[target_group]] > peak_densthresh) && 
                    (group_peak_density[j] > peak_densthresh)) {
                    merge_priorities[j]->pop_back();
                } else if (targets_target_group == this_group) { // Groups exchanged vows...
                    merge_priorities[j]->pop_back();
                        if (debug) printf("merging fringe groups %d and %d\n", this_group, target_group);
                    merge_groups( // and shall be one flesh
                        this_group, target_group, 
                        merge_priorities, group_peak_density, id_from_index, index_from_id);
                    for (i = 0; i < nGroups_before_regroup; i++) {
                        if (new_group_id[i] == target_group) {
                            new_group_id[i] = this_group;
                        }
                    }
                    nGroups--;
                    changes = true;
                }
            }
        }
    }
    
    // Reindex the new group ids to ids starting at 0, increment 1
    int *reindexed_new_group_id = new int[nGroups_before_regroup+1]; // +1 for particles not in any group
    int *reindexer = new int[nGroups_before_regroup]; // to go from new group id to reindexed new group id
    for (i = 0; i < nGroups_before_regroup; i++) {
        reindexer[i] = -1; // default, means will be dropped (for fringe groups that do not belong to any proper group)
    }
    int tmp_id, id_counter = 0;
    for (i = 0; i < nGroups_before_regroup; i++) {
        tmp_id = new_group_id[i];
        if ((reindexer[tmp_id] == -1) && (group_peak_density[index_from_id[tmp_id]] > peak_densthresh)) {
            reindexer[tmp_id] = id_counter++; // Post increment, next new group will get this id incremented by one
        }
    }
    reindexed_new_group_id[0] = -1; // for particles not in any group
    for (i = 0; i < nGroups_before_regroup; i++) {
        reindexed_new_group_id[i+1] = reindexer[new_group_id[i]];
    }
    
    // Apply results
    ParticlesMapIterator it_p;
    AmuseParticle *p;
    for (it_p = particlesMap.begin(); it_p != particlesMap.end(); it_p++) {
        p = (*it_p).second;
        p->group = reindexed_new_group_id[p->group+1];
    }
    nGroups_after_regroup = id_counter;
    
    delete[] index_from_id;
    delete[] new_group_id;
    delete[] reindexed_new_group_id;
    delete[] reindexer;
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
  
  regroup();
  
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
  if (p->neighbor == -1 ) return 0;
  *index_of_densest_neighbor = p->neighbor;
  return 0;
}

int get_group_id(int index_of_the_particle, int * group_id) {
  if (!bHopDone) return -1;
  if(index_of_the_particle > highest_index) return -2;
  AmuseParticle * p = particlesMap[index_of_the_particle];
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
  *number_of_groups = nGroups_after_regroup;
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








