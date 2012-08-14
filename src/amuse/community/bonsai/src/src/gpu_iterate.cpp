#include "octree.h"

#include <iostream>
#include <algorithm>
#include <iomanip>

#ifdef _AMUSE_STOPPING_CONDITIONS_
// AMUSE STOPPING CONDITIONS SUPPORT
#include <stopcond.h>
#include <map>
#endif

using namespace std;


static double de_max = 0;
static double dde_max = 0;  
extern bool maxlevels_exceeded;

void octree::makeLET()
{
   //LET code test
  double tTest = get_time();

  my_dev::dev_stream memCpyStream;
  
  localTree.bodies_Ppos.d2h(false, memCpyStream.s());
  localTree.bodies_Pvel.d2h(false, memCpyStream.s());
  localTree.multipole.d2h(false, memCpyStream.s());   
  localTree.boxSizeInfo.d2h(false, memCpyStream.s());
  localTree.boxCenterInfo.d2h(false, memCpyStream.s());
  
  //Exchange domain boundaries, while memory copies take place

  rMinLocalTreeGroups.w = this->maxLocalEps;
  sendCurrentRadiusInfo(rMinLocalTreeGroups,rMaxLocalTreeGroups);  

  memCpyStream.sync();  //Sync otherwise we are not certain the required data is on the host
    
  //Exchange particles and start LET kernels
  vector<real4> LETParticles;
  essential_tree_exchange(LETParticles, localTree, remoteTree);
  fprintf(stderr, "LET Exchange took (%d): %g \n", mpiGetRank(), get_time() - tTest);
  
  letRunning = false;
  execStream->sync();  //Sync LET execution

  
  #if 0
    //Write the LET tree structure to file

    char fileName[256];
    sprintf(fileName, "LETTreeStructure-%d.txt", mpiGetRank());
    ofstream nodeFile;
    //nodeFile.open(nodeFileName.c_str());
    nodeFile.open(fileName);

    //The texture offsets used:
    int particleTexOffset = (remoteTree.remoteTreeStruct.z >> 16);
    int nodeTexOffset     = (remoteTree.remoteTreeStruct.z & 0xFFFF);
 
    //Number of particles and number of nodes in the remote tree
    int remoteP = remoteTree.remoteTreeStruct.x;
    int remoteN = remoteTree.remoteTreeStruct.y;
    
    double massTEST = 0;
    
    for(int i=0; i < remoteN; i++)
    {
      int centIdx = (2*(particleTexOffset+remoteP) + remoteN + nodeTexOffset) + i;
      int sizeIdx = (2*(particleTexOffset+remoteP)) + i;
      int multIdx = (2*(particleTexOffset+remoteP) + 2*(remoteN+nodeTexOffset)) + i*3 + 0 ;
      
      nodeFile << i << "\t" << remoteTree.fullRemoteTree[centIdx].x-remoteTree.fullRemoteTree[sizeIdx].x << "\t" << remoteTree.fullRemoteTree[centIdx].y-remoteTree.fullRemoteTree[sizeIdx].y;
      
      if(remoteTree.fullRemoteTree[sizeIdx].x == 0 && remoteTree.fullRemoteTree[sizeIdx].y == 0)
      {
        nodeFile << "\t" << remoteTree.fullRemoteTree[centIdx].x+0.01 << "\t" << remoteTree.fullRemoteTree[centIdx].y+0.01 << "\t";
      }
      else
      {
        nodeFile << "\t" << remoteTree.fullRemoteTree[centIdx].x+remoteTree.fullRemoteTree[sizeIdx].x << "\t" << remoteTree.fullRemoteTree[centIdx].y+remoteTree.fullRemoteTree[sizeIdx].y << "\t";
      }
     

      nodeFile <<              remoteTree.fullRemoteTree[multIdx].x << "\t" 
               << remoteTree.fullRemoteTree[multIdx].y << "\t" << remoteTree.fullRemoteTree[multIdx].w << "\n";
      
      massTEST +=  remoteTree.fullRemoteTree[multIdx].w;
    }

    nodeFile.close();
    

    sprintf(fileName, "LETTreeStructureParticles-%d.txt", mpiGetRank());
    ofstream partFile;
    partFile.open(fileName);

                                                                                                                    
    for(int i=0; i < remoteP; i++)                                                                                     
    {                                                                                                                 
      float4  pos =  remoteTree.fullRemoteTree[i];                                                                              
      partFile << i << "\t" << pos.x << "\t" << pos.y << "\t" << pos.z << endl;       
      
      massTEST += pos.w;
    }                                                                                                                 
    partFile.close(); 
    
    cerr << "Mass test (rank= " << mpiGetRank() << " ) : " << massTEST << std::endl;
  
  exit(0);

  
  #endif

}


void octree::iterate() {
  
  #ifdef _AMUSE_STOPPING_CONDITIONS_
    //Reset stopping condition settings
    this->stopping_condition_found = 0;      
    reset_stopping_conditions();
  #endif    
      
  int Nact_since_last_tree_rebuild = 0;
  real4 r_min, r_max;
  
  if(execStream == NULL)
    execStream = new my_dev::dev_stream(0);
  
  letRunning = false;
  
  double totalGravTime = 0;
  double lastGravTime  = 0;
  double totalBuildTime = 0;
  double lastBuildTime  = 0;  
  double totalDomTime = 0;
  double lastDomTime  = 0;  
  double totalWaitTime = 0;
  double lastWaitTime  = 0;    
    
  double t1;

  //Initial prediction/acceleration to setup the system
  //Will be at time 0
  //predict localtree
  predict(this->localTree);
  this->getBoundaries(localTree, r_min, r_max);
  //Build the tree using the predicted positions  
  //Compute the (new) node properties
  compute_properties(this->localTree);
  
  t1 = get_time();
   
  //Approximate gravity
//   devContext.startTiming();
  approximate_gravity(this->localTree);
//   devContext.stopTiming("Approximation", 4);

  if(nProcs > 1)  makeLET();

  execStream->sync();  

  
  lastGravTime   = get_time() - t1;
  totalGravTime += lastGravTime;
  
  correct(this->localTree);
  compute_energies(this->localTree);

  //Print time 0 snapshot
  if(snapshotIter > 0 )
  {
      int time = t_current;
      if((time >= nextSnapTime))
      {
        nextSnapTime += snapshotIter;
        string fileName; fileName.resize(256);
        sprintf(&fileName[0], "%s_%06d", snapshotFile.c_str(), time + snapShotAdd);

        localTree.bodies_pos.d2h();
        localTree.bodies_vel.d2h();
        localTree.bodies_ids.d2h();

        write_dumbp_snapshot_parallel(&localTree.bodies_pos[0], &localTree.bodies_vel[0],
                                      &localTree.bodies_ids[0], localTree.n, fileName.c_str()) ;
      }
  }

  double t0 = get_time();
  for(int i=0; i < 10000000; i++) //Large number, limit
  {
    printf("At the start of iterate:\n");


    //predict localtree
    devContext.startTiming();
    predict(this->localTree);
    devContext.stopTiming("Predict", 9);
    
    bool needDomainUpdate = true;
    
   //Redistribute the particles
    if(1)
    {      
      if(nProcs > 1)
      { 
//       if(iter % 5 == 0)
//        if(0)
        if(1)
        {     
          //If we do a redistribution we _always_ have to do 
          //an update of the particle domain, otherwise the boxes 
          //do not match and we get errors of particles outside
          //domains
          t1 = get_time();
          
          devContext.startTiming();
          gpu_updateDomainDistribution(lastGravTime);          
          devContext.stopTiming("DomainUpdate", 6);
          
          devContext.startTiming();
          gpuRedistributeParticles();
          devContext.stopTiming("Exchange", 6);
          
          needDomainUpdate = false;
          
          lastDomTime   = get_time() - t1;
          totalDomTime += lastDomTime;          
        }
        else
        {
          //Only send new box sizes, incase we do not exchange particles
          //but continue with the current tree_structure
          gpu_updateDomainOnly();
          
          needDomainUpdate = false;
        } //if (iter % X )
//         else
//         {
//           //Only exchange, do not update domain decomposition
//         }
      } //if nProcs > 1
    }//if (0)        
    
    
    //Build the tree using the predicted positions
#ifdef ADAPTIVE_TIMESTEP
    bool rebuild_tree = Nact_since_last_tree_rebuild > 4*this->localTree.n;   
#else
    bool rebuild_tree = true;
#endif
   // bool rebuild_tree = true;
    
    if(rebuild_tree)
    {
      t1 = get_time();
      //Rebuild the tree
      this->sort_bodies(this->localTree, needDomainUpdate);

      devContext.startTiming();
      this->build(this->localTree);
      if (maxlevels_exceeded) return;
      devContext.stopTiming("Tree-construction", 2);

      devContext.startTiming();
      this->allocateTreePropMemory(this->localTree);
      devContext.stopTiming("Memory", 11);      

      devContext.startTiming();
      this->compute_properties(this->localTree);
      devContext.stopTiming("Compute-properties", 3);

      devContext.startTiming();
      setActiveGrpsFunc(this->localTree);
      devContext.stopTiming("setActiveGrpsFunc", 10);      
      Nact_since_last_tree_rebuild = 0;
      
      lastBuildTime   = get_time() - t1;
      totalBuildTime += lastBuildTime;  
    }
    else
    {
      devContext.startTiming();
      this->compute_properties(this->localTree);
      devContext.stopTiming("Compute-properties", 3);
    }//end rebuild tree

    //Approximate gravity
    t1 = get_time();
//     devContext.startTiming();
    approximate_gravity(this->localTree);
//     devContext.stopTiming("Approximation", 4);
    
    
    if(nProcs > 1)  makeLET();
    
    execStream->sync();
    
    lastGravTime   = get_time() - t1;
//     totalGravTime += lastGravTime;
    totalGravTime += lastGravTime - thisPartLETExTime;
//     lastGravTime -= thisPartLETExTime;
    
    fprintf(stderr,"APPTIME [%d]: Iter: %d\t%g \n", procId, iter, lastGravTime);
    
    //Corrector
    devContext.startTiming();
    correct(this->localTree);
    devContext.stopTiming("Correct", 8);
    
    t1 = get_time();
    devContext.startTiming();
    mpiSync();
    devContext.stopTiming("Unbalance", 12);
    lastWaitTime  += get_time() - t1;
    totalWaitTime += lastWaitTime;
    
    Nact_since_last_tree_rebuild += this->localTree.n_active_particles;

    //Compute energies
    devContext.startTiming();
    double de = compute_energies(this->localTree); de=de;
    devContext.stopTiming("Energy", 7);
    
    
    #ifdef _AMUSE_STOPPING_CONDITIONS_
      //Check if a stopping condition was triggered
      if(this->stopping_condition_found)
      {
        return;
      }
    #endif
      
    
        

    if(snapshotIter > 0)
    {
      int time = t_current;
      if((time >= nextSnapTime))
      {
        nextSnapTime += snapshotIter;
        string fileName; fileName.resize(256);
        sprintf(&fileName[0], "%s_%06d", snapshotFile.c_str(), time + snapShotAdd);

        localTree.bodies_pos.d2h();
        localTree.bodies_vel.d2h();
        localTree.bodies_ids.d2h();

        write_dumbp_snapshot_parallel(&localTree.bodies_pos[0], &localTree.bodies_vel[0],
                                      &localTree.bodies_ids[0], localTree.n, fileName.c_str()) ;
      }
    }

  //  checkMergingDistance(this->localTree, iter, de);    

    if(t_current >= tEnd)
    {
      compute_energies(this->localTree);
      cout << " Finished: "  << t_current << "  > "  << tEnd << " loop alone took: " << get_time() -t0 <<  endl;
      fprintf(stderr,"TIME [%02d] TOTAL: %g\t GRAV: %g\tBUILD: %g\tCOMM: %g\t WAIT: %g\n", 
              procId, get_time() -t0, totalGravTime, totalBuildTime, totalDomTime, lastWaitTime);      
     
      my_dev::base_mem::printMemUsage();

      if(execStream != NULL)
      {
        delete execStream;
        execStream = NULL;
      }
      
      return;
    }
    
    
    if((iter % 50) == 0)
    {
//       if(removeDistance > 0) checkRemovalDistance(this->localTree);
    }

    iter++;
    
  } //end for i
  
  fprintf(stderr,"TOTAL GRAV TIME [%d ] \t %g \n", procId, totalGravTime);
  
  if(execStream != NULL)
  {
    delete execStream;
    execStream = NULL;
  }
  
} //end iterate


void octree::predict(tree_structure &tree)
{
  //Functions that predicts the particles to the next timestep

//   tend is time per particle
//   tnext is reduce result

  //First we get the minimum time, which is the next integration time
  int blockSize = NBLOCK_REDUCE ;
  getTNext.set_arg<int>(0,    &tree.n);
  getTNext.set_arg<cl_mem>(1, tree.bodies_time.p());
  getTNext.set_arg<cl_mem>(2, tnext.p());
  getTNext.set_arg<double>(3,  NULL, 128); //Dynamic shared memory
  getTNext.setWork(-1, 128, blockSize);
  getTNext.execute();

  //Reduce the last parts on the host
  tnext.d2h();
  t_previous = t_current;
  t_current  = tnext[0];
  for (int i = 1; i < blockSize ; i++)
  {
      t_current = fmin(t_current, tnext[i]);
  }

/*
  bool timeFix = false;
  if(iter > 0){
  if(t_current == t_previous){
	  fprintf(stderr, "TIMEFIX: iter %d  t_current: %lg (%.15f)   t_prev: %lg (%.15f)  \n", iter, t_current, t_current, t_previous, t_previous);
	  t_current += 1./16384;
	  fprintf(stderr, "TIMEFIX2: iter %d  t_current: %lg (%.15f)   t_prev: %lg (%.15f) test: %lg  \n", iter, t_current, t_current, t_previous, t_previous, 1./16384);
	  t_current += 1.0/16384.0;
	  fprintf(stderr, "TIMEFIX3: iter %d  t_current: %lg (%.15f)   t_prev: %lg (%.15f) test: %lg  \n", iter, t_current, t_current, t_previous, t_previous, 1./16384.0);
	  timeFix = true;

	  fprintf(stderr, "TIMEFIX4: iter %d  t_current: %lg (%.15f)  (%.15f) (%.15f)  \n", iter, t_current, t_current+1./16384.0, t_current+1./8192., ((double)t_current)+1.0/16384.0);
  }}
*/

  tree.activeGrpList.zeroMem();      //Reset the active grps

  //Set valid list to zero
  predictParticles.set_arg<int>(0,    &tree.n);
  predictParticles.set_arg<double>(1,  &t_current);
  predictParticles.set_arg<double>(2,  &t_previous);
  predictParticles.set_arg<cl_mem>(3, tree.bodies_pos.p());
  predictParticles.set_arg<cl_mem>(4, tree.bodies_vel.p());
  predictParticles.set_arg<cl_mem>(5, tree.bodies_acc0.p());
  predictParticles.set_arg<cl_mem>(6, tree.bodies_time.p());
  predictParticles.set_arg<cl_mem>(7, tree.body2group_list.p());
  predictParticles.set_arg<cl_mem>(8, tree.activeGrpList.p());
  predictParticles.set_arg<cl_mem>(9, tree.bodies_Ppos.p());
  predictParticles.set_arg<cl_mem>(10, tree.bodies_Pvel.p());  

  predictParticles.setWork(tree.n, 128);
  predictParticles.execute();
  clFinish(devContext.get_command_queue());
  
  //Compact the valid list to get a list of valid groups
  gpuCompact(devContext, tree.activeGrpList, tree.active_group_list,
             tree.n_groups, &tree.n_active_groups);

  printf("t_previous: %lg t_current: %lg dt: %lg Active groups: %d \n",
         t_previous, t_current, t_current-t_previous, tree.n_active_groups);
  fprintf(stderr,"t_previous: %lg (%.15f) t_current: %lg (%.15f) dt: %lg (%.15f) Active groups: %d \n",
         t_previous, t_previous, t_current, t_current, t_current-t_previous, t_current-t_previous, tree.n_active_groups);



/*
  if(timeFix){
	  tree.bodies_time.d2h();
	  tree.body2group_list.d2h();
	  
	  for(int i=0; i < tree.n; i++)
	  {
		  fprintf(stderr,"%d\t t.x: %lg (%.15f) \tt.y: %lg (%.15f) \t %d  \n",
				  i, tree.bodies_time[i].x, tree.bodies_time[i].x,
				  tree.bodies_time[i].y, tree.bodies_time[i].y, tree.body2group_list[i]);
	  }

	 tree.active_group_list.d2h();
	 for(int i=0; i < tree.n_active_groups; i++)
	 {
		 fprintf(stderr, "%d\tgrp: %d \n", i, tree.active_group_list[i]);
	 }
	exit(0);
  }
*/



}
//End predict


void octree::setActiveGrpsFunc(tree_structure &tree)
{
  tree.activeGrpList.zeroMem();      //Reset the active grps

  //Set valid list to zero
  setActiveGrps.set_arg<int>(0,    &tree.n);
  setActiveGrps.set_arg<double>(1,  &t_current);
  setActiveGrps.set_arg<cl_mem>(2, tree.bodies_time.p());
  setActiveGrps.set_arg<cl_mem>(3, tree.body2group_list.p());
  setActiveGrps.set_arg<cl_mem>(4, tree.activeGrpList.p());

  setActiveGrps.setWork(tree.n, 128);
  setActiveGrps.execute();
  clFinish(devContext.get_command_queue());

  //Compact the valid list to get a list of valid groups
  gpuCompact(devContext, tree.activeGrpList, tree.active_group_list,
             tree.n_groups, &tree.n_active_groups);

  printf("t_previous: %lg t_current: %lg dt: %lg Active groups: %d \n",
         t_previous, t_current, t_current-t_previous, tree.n_active_groups);

}

void octree::approximate_gravity(tree_structure &tree)
{ 
  uint2 node_begend;
  int level_start = 2;
  node_begend.x   = tree.level_list[level_start].x;
  node_begend.y   = tree.level_list[level_start].y;

  printf("node begend: %d %d iter-> %d\n", node_begend.x, node_begend.y, iter);

  //Reset the active particles
  tree.activePartlist.zeroMem();

//   int grpOffset = 0;

  //Set the kernel parameters, many!
  approxGrav.set_arg<int>(0,    &tree.n_active_groups);
  approxGrav.set_arg<int>(1,    &tree.n);
  approxGrav.set_arg<float>(2,  &(this->eps2));
  approxGrav.set_arg<uint2>(3,  &node_begend);
  approxGrav.set_arg<cl_mem>(4, tree.active_group_list.p());
  approxGrav.set_arg<cl_mem>(5, tree.bodies_Ppos.p());
  approxGrav.set_arg<cl_mem>(6, tree.multipole.p());
  approxGrav.set_arg<cl_mem>(7, tree.bodies_acc1.p());
  approxGrav.set_arg<cl_mem>(8, tree.ngb.p());
  approxGrav.set_arg<cl_mem>(9, tree.activePartlist.p());
  approxGrav.set_arg<cl_mem>(10, tree.interactions.p());
  approxGrav.set_arg<cl_mem>(11, tree.boxSizeInfo.p());
  approxGrav.set_arg<cl_mem>(12, tree.groupSizeInfo.p());
  approxGrav.set_arg<cl_mem>(13, tree.boxCenterInfo.p());
  approxGrav.set_arg<cl_mem>(14, tree.groupCenterInfo.p());
  approxGrav.set_arg<cl_mem>(15, tree.bodies_Pvel.p());
  approxGrav.set_arg<cl_mem>(16,  tree.generalBuffer1.p()); //Instead of using Local memory
  
  approxGrav.set_arg<real4>(17, tree.boxSizeInfo, 4, "texNodeSize");
  approxGrav.set_arg<real4>(18, tree.boxCenterInfo, 4, "texNodeCenter");
  approxGrav.set_arg<real4>(19, tree.multipole, 4, "texMultipole");
  approxGrav.set_arg<real4>(20, tree.bodies_Ppos, 4, "texBody");
  
  approxGrav.setWork(-1, NTHREAD, nBlocksForTreeWalk);    
  approxGrav.execute(execStream->s());  //First half

  //Print interaction statistics
  #if 0
  
  tree.body2group_list.d2h();
  tree.interactions.d2h();
    long long directSum = 0;
    long long apprSum = 0;
    long long directSum2 = 0;
    long long apprSum2 = 0;
    
    
    int maxDir = -1;
    int maxAppr = -1;

    for(int i=0; i < tree.n; i++)
    {
      apprSum     += tree.interactions[i].x;
      directSum   += tree.interactions[i].y;
      
      maxAppr = max(maxAppr,tree.interactions[i].x);
      maxDir  = max(maxDir,tree.interactions[i].y);
      
      apprSum2     += tree.interactions[i].x*tree.interactions[i].x;
      directSum2   += tree.interactions[i].y*tree.interactions[i].y;      
    }
  
    //cerr << "Interaction at iter: " << iter << "\tdirect: " << directSum << "\tappr: " << apprSum << "\t";
    //cerr << "avg dir: " << directSum / tree.n << "\tavg appr: " << apprSum / tree.n << endl;

    cout << "Interaction at (rank= " << mpiGetRank() << " ) iter: " << iter << "\tdirect: " << directSum << "\tappr: " << apprSum << "\t";
    cout << "avg dir: " << directSum / tree.n << "\tavg appr: " << apprSum / tree.n << "\tMaxdir: " << maxDir << "\tmaxAppr: " << maxAppr <<  endl;
    cout << "sigma dir: " << sqrt((directSum2  - directSum)/ tree.n) << "\tsigma appr: " << std::sqrt((apprSum2 - apprSum) / tree.n)  <<  endl;    
 
    #if 0
      //Histogram of number of interactions
      const int bins = 256;
      const int jump = 15;
      int histoIDX[bins+1];
      for(int i=0; i < bins; i++)
        histoIDX[i] = 0;
      
      for(int i=0; i < tree.n; i++)
      {
          int idx = tree.interactions[i].x / jump;
          if(idx >= bins)
            idx = bins;
          histoIDX[idx]++;  
      }
      for(int i=0; i < bins; i++)
      {
        if(histoIDX[i] == 0)
          fprintf(stderr, "HISTO %d\t-\n", i*jump, histoIDX[i]);
        else
          fprintf(stderr, "HISTO %d\t%d\n", i*jump, histoIDX[i]);
      }     
    #endif
  #endif
  
  #if 0
  tree.bodies_acc1.d2h();
  
  for(int i=0; i < tree.n; i++)
  {
      fprintf(stderr, "%d\t Acc: %f %f %f %f \n", 
              i, tree.bodies_acc1[i].x,tree.bodies_acc1[i].y,
              tree.bodies_acc1[i].w,tree.bodies_acc1[i].z);
  }
  #endif

  if(mpiGetNProcs() == 1) //Only do it here if there is only one process
  {
    
    //Reduce the number of valid particles    
    getNActive.set_arg<int>(0,    &tree.n);
    getNActive.set_arg<cl_mem>(1, tree.activePartlist.p());
    getNActive.set_arg<cl_mem>(2, this->nactive.p());
    getNActive.set_arg<int>(3,    NULL, 128); //Dynamic shared memory , equal to number of threads
    getNActive.setWork(-1, 128,   NBLOCK_REDUCE);
    
    CU_SAFE_CALL(cuCtxSynchronize()); //Synchronize all streams, makes sure that the approx stream is finished
    getNActive.execute();
    
    //Reduce the last parts on the host
    this->nactive.d2h();
    tree.n_active_particles = this->nactive[0];
    for (int i = 1; i < NBLOCK_REDUCE ; i++)
        tree.n_active_particles += this->nactive[i];

    printf("Active particles: %d \n", tree.n_active_particles);
  }
}
//end approximate


void octree::approximate_gravity_let(tree_structure &tree, tree_structure &remoteTree, int bufferSize, bool doActiveParticles)
{

  //Start and end node of the remote tree structure
  uint2 node_begend;  
  node_begend.x =  0;
  node_begend.y =  remoteTree.remoteTreeStruct.w;
  
  //The texture offset used:
  int nodeTexOffset     = remoteTree.remoteTreeStruct.z ;
  
  //The start and end of the top nodes:
  node_begend.x = (remoteTree.remoteTreeStruct.w >> 16);
  node_begend.y = (remoteTree.remoteTreeStruct.w & 0xFFFF);  
 
  //Number of particles and number of nodes in the remote tree
  int remoteP = remoteTree.remoteTreeStruct.x;
  int remoteN = remoteTree.remoteTreeStruct.y;

  printf("LET node begend [%d]: %d %d iter-> %d\n", procId, node_begend.x, node_begend.y, iter);
  fflush(stderr);
  fflush(stdout);


  //Set the kernel parameters, many!
  approxGravLET.set_arg<int>(0,    &tree.n_active_groups);
  approxGravLET.set_arg<int>(1,    &tree.n);
  approxGravLET.set_arg<float>(2,  &(this->eps2));
  approxGravLET.set_arg<uint2>(3,  &node_begend);
  approxGravLET.set_arg<cl_mem>(4, tree.active_group_list.p());
  approxGravLET.set_arg<cl_mem>(5, remoteTree.fullRemoteTree.p());

  void *multiLoc = remoteTree.fullRemoteTree.a(2*(remoteP) + 2*(remoteN+nodeTexOffset));
  approxGravLET.set_arg<cl_mem>(6, &multiLoc);  

  approxGravLET.set_arg<cl_mem>(7, tree.bodies_acc1.p());
  approxGravLET.set_arg<cl_mem>(8, tree.ngb.p());
  approxGravLET.set_arg<cl_mem>(9, tree.activePartlist.p());
  approxGravLET.set_arg<cl_mem>(10, tree.interactions.p());
  
  void *boxSILoc = remoteTree.fullRemoteTree.a(2*(remoteP));
  approxGravLET.set_arg<cl_mem>(11, &boxSILoc);  

  approxGravLET.set_arg<cl_mem>(12, tree.groupSizeInfo.p());

  void *boxCILoc = remoteTree.fullRemoteTree.a(2*(remoteP) + remoteN + nodeTexOffset);
  approxGravLET.set_arg<cl_mem>(13, &boxCILoc);  

  approxGravLET.set_arg<cl_mem>(14, tree.groupCenterInfo.p());  
  
  void *bdyVelLoc = remoteTree.fullRemoteTree.a(1*(remoteP));
  approxGravLET.set_arg<cl_mem>(15, &bdyVelLoc);  //<- Remote bodies velocity
  
  approxGravLET.set_arg<cl_mem>(16, tree.bodies_Ppos.p()); //<- Predicted local body positions
  approxGravLET.set_arg<cl_mem>(17, tree.bodies_Pvel.p()); //<- Predicted local body velocity
  approxGravLET.set_arg<cl_mem>(18, tree.generalBuffer1.p()); //<- Predicted local body velocity
  
  approxGravLET.set_arg<real4>(19, remoteTree.fullRemoteTree, 4, "texNodeSize",
                               2*(remoteP), remoteN );
  approxGravLET.set_arg<real4>(20, remoteTree.fullRemoteTree, 4, "texNodeCenter",
                               2*(remoteP) + (remoteN + nodeTexOffset),
                               remoteN);
  approxGravLET.set_arg<real4>(21, remoteTree.fullRemoteTree, 4, "texMultipole",
                               2*(remoteP) + 2*(remoteN + nodeTexOffset), 
                               3*remoteN);
  approxGravLET.set_arg<real4>(22, remoteTree.fullRemoteTree, 4, "texBody", 0, remoteP);  
    
  approxGravLET.setWork(-1, NTHREAD, nBlocksForTreeWalk);
 
  printf("LET Approx config: "); approxGravLET.printWorkSize();
    
  if(letRunning)
  {
    //dont want to overwrite the data of previous LET tree
    execStream->sync();
  }
  
  remoteTree.fullRemoteTree.h2d(bufferSize); //Only copy required data
  tree.activePartlist.zeroMem();
//   devContext.startTiming();  
  approxGravLET.execute(execStream->s());
//   devContext.stopTiming("Approximation_let", 5);   
  
  letRunning = true;

 //Print interaction statistics
  #if 0
    tree.interactions.d2h();
//     tree.body2group_list.d2h();
    
    long long directSum = 0;
    long long apprSum = 0;
    
    int maxDir = -1;
    int maxAppr = -1;
    
    long long directSum2 = 0;
    long long apprSum2 = 0;
    
    
    for(int i=0; i < tree.n; i++)
    {
      apprSum     += tree.interactions[i].x;
      directSum   += tree.interactions[i].y;
      
      maxAppr = max(maxAppr,tree.interactions[i].x);
      maxDir  = max(maxDir, tree.interactions[i].y);
      
      apprSum2     += (tree.interactions[i].x*tree.interactions[i].x);
      directSum2   += (tree.interactions[i].y*tree.interactions[i].y);    
    }

    cout << "Interaction (LET) at (rank= " << mpiGetRank() << " ) iter: " << iter << "\tdirect: " << directSum << "\tappr: " << apprSum << "\t";
    cout << "avg dir: " << directSum / tree.n << "\tavg appr: " << apprSum / tree.n  << "\tMaxdir: " << maxDir << "\tmaxAppr: " << maxAppr <<  endl;
    cout << "sigma dir: " << sqrt((directSum2  - directSum)/ tree.n) << "\tsigma appr: " << std::sqrt((apprSum2 - apprSum) / tree.n)  <<  endl;
  #endif


  if(doActiveParticles)
  {
    //Reduce the number of valid particles
    getNActive.set_arg<int>(0,    &tree.n);
    getNActive.set_arg<cl_mem>(1, tree.activePartlist.p());
    getNActive.set_arg<cl_mem>(2, this->nactive.p());
    getNActive.set_arg<int>(3, NULL, 128); //Dynamic shared memory , equal to number of threads
    getNActive.setWork(-1, 128, NBLOCK_REDUCE);
    
    CU_SAFE_CALL(cuCtxSynchronize()); //Syncrhonize all streams, ensures that the approx streams are finished
    getNActive.execute();
    
    //Reduce the last parts on the host
    this->nactive.d2h();
    tree.n_active_particles = this->nactive[0];
    for (int i = 1; i < NBLOCK_REDUCE ; i++)
        tree.n_active_particles += this->nactive[i];

    printf("LET Active particles: %d (Process: %d ) \n",tree.n_active_particles, mpiGetRank());
  }
}
//end approximate



void octree::correct(tree_structure &tree)
{
  correctParticles.set_arg<int   >(0, &tree.n);
  correctParticles.set_arg<double>(1, &t_current);
  correctParticles.set_arg<cl_mem>(2, tree.bodies_time.p());
  correctParticles.set_arg<cl_mem>(3, tree.activePartlist.p());
  correctParticles.set_arg<cl_mem>(4, tree.bodies_vel.p());
  correctParticles.set_arg<cl_mem>(5, tree.bodies_acc0.p());
  correctParticles.set_arg<cl_mem>(6, tree.bodies_acc1.p());
  correctParticles.set_arg<cl_mem>(7, tree.bodies_pos.p());
  correctParticles.set_arg<cl_mem>(8, tree.bodies_Ppos.p());
  correctParticles.set_arg<cl_mem>(9, tree.bodies_Pvel.p());
#ifdef _AMUSE_STOPPING_CONDITIONS_
  //Stopping conditions requires extra buffers and tests
  
  //The buffers to hold the results
  my_dev::dev_mem<uint> pairDetectionBuffer(devContext);
  my_dev::dev_mem<uint> compactPairDetectionBuffer(devContext);  
  pairDetectionBuffer.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[0], 0,  
                         2*tree.n, getAllignmentOffset(0));  

  compactPairDetectionBuffer.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                         tree.generalBuffer1.get_flags(), 
                         tree.generalBuffer1.get_devMem(),
                         &tree.generalBuffer1[2*tree.n], 2*tree.n,  
                         2*tree.n, getAllignmentOffset(2*tree.n));   
  
  
  pairDetectionBuffer.zeroMem();

  correctParticles.set_arg<cl_mem>(11, tree.ngb.p());
  correctParticles.set_arg<cl_mem>(12, pairDetectionBuffer.p());
#endif   
     

  correctParticles.setWork(tree.n, 128);
  correctParticles.execute();
  clFinish(devContext.get_command_queue());
  
#ifdef _AMUSE_STOPPING_CONDITIONS_
  
  int is_collision_detection_enabled;
  int error = is_stopping_condition_enabled(COLLISION_DETECTION, &is_collision_detection_enabled);  
  if(is_collision_detection_enabled)
  {
    //Compact the detection list to check if there were collisions
    int collisionPairs = 0;
    gpuCompact(devContext, pairDetectionBuffer, compactPairDetectionBuffer,
              tree.n*2, &collisionPairs);
    
    if(collisionPairs > 0)
    {
      //We store keys with values to check on double occurances
      std::map<int, std::vector<int> > doubleCheck; 
      //Set the stopping info
      //Pairs are in : compactPairDetectionBuffer
      //Pair 0: compactPairDetectionBuffer[2*0+0],compactPairDetectionBuffer[2*0+1]
      //Pair 1: compactPairDetectionBuffer[2*1+0],compactPairDetectionBuffer[2*1+1]
      //etc.
      compactPairDetectionBuffer.d2h();
      tree.bodies_ids.d2h();
      for(int i =0; i < collisionPairs; i +=2)
      {
        int val1   = compactPairDetectionBuffer[i];
        int val2   = compactPairDetectionBuffer[i+1];
        int key    = std::min(val1, val2);
        int val    = std::max(val1, val2);
        bool isNew = true;
        
        if(doubleCheck.find( key ) != doubleCheck.end())
        {
          std::vector<int> temp = doubleCheck[key];
          for(uint j=0; j < temp.size(); j++)
          {
            if(temp[j] == val) isNew = false;
          }
        }
        
        if(isNew){
          //Add this pair
          std::vector<int> temp =  doubleCheck[key];
          temp.push_back(val);
          doubleCheck[key] = temp;
        
          int stopping_index  = next_index_for_stopping_condition();        //get the next storage location
          if(stopping_index  >= 0)
          {
            set_stopping_condition_info(stopping_index, COLLISION_DETECTION);
            set_stopping_condition_particle_index(stopping_index, 0, tree.bodies_ids[compactPairDetectionBuffer[i]]);        set_stopping_condition_particle_index(stopping_index, 1, tree.bodies_ids[compactPairDetectionBuffer[i+1]]);
            //fprintf(stderr, "STOP found: %d \t %d   %d \n", i, compactPairDetectionBuffer[i], compactPairDetectionBuffer[i+1]);
          }
        }
        else {
          //Double item
        }
        
      }
      //Set error
      this->stopping_condition_found = 1;
    }
  }
  
  
#endif  //Ifdef _AMUSE_STOPPING_CONDITIONS_
  
  
  

  computeDt.set_arg<int>(0,    &tree.n);
  computeDt.set_arg<double>(1,  &t_current);
  computeDt.set_arg<float>(2,  &(this->eta));
  computeDt.set_arg<int>(3,    &(this->dt_limit));
  computeDt.set_arg<float>(4,  &(this->eps2));
  computeDt.set_arg<cl_mem>(5, tree.bodies_time.p());
  computeDt.set_arg<cl_mem>(6, tree.bodies_vel.p());
  computeDt.set_arg<cl_mem>(7, tree.ngb.p());
  computeDt.set_arg<cl_mem>(8, tree.bodies_pos.p());
  computeDt.set_arg<cl_mem>(9, tree.bodies_acc0.p());
  computeDt.set_arg<cl_mem>(10, tree.activePartlist.p());
  computeDt.set_arg<float >(11, &timeStep);

  computeDt.setWork(tree.n, 128);
  computeDt.execute();
  clFinish(devContext.get_command_queue());
}

int octree::checkMergingDistance(tree_structure &tree, int iter, double dE)
{                                                                                                                                                       
  //Allocate memory to store the current positions                                                                                                          
  //of the blackholes  
  
  float starSoftening = 0.0;
  
  int nBH = NThird+1;   
  
  if(NThird == 1)                                                                                                                                              
   return 0;                                                                                                                                                
                                                                                                                                                            
  my_dev::dev_mem<real4>  bhDistance(devContext, nBH*2);                                                                                                    
  bhDistance.zeroMem();                                                                                                                                     
                                                                                                                                                            
                                                                                                                                                            
                                                                                                                                                            
  distanceCheck.set_arg<int>(0,    &tree.n);                                                                                                                
  distanceCheck.set_arg<cl_mem>(1, tree.bodies_pos.p());                                                                                                    
  distanceCheck.set_arg<cl_mem>(2, tree.bodies_ids.p());                                                                                                    
  distanceCheck.set_arg<cl_mem>(3, bhDistance.p());                                                                                                         
  distanceCheck.set_arg<cl_mem>(4, &nBH);        //Number of black holes                                                                                    
  distanceCheck.set_arg<cl_mem>(5, tree.bodies_vel.p());                                                                                                    
                                                                                                                                                            
  distanceCheck.setWork(tree.n, 128);                                                                                                                       
  distanceCheck.execute();                                                                                                                                  
                                                                                                                                                            
  bhDistance.d2h();                                                                                                                                         
                                                                                                                                                            
  //Calculate the distance between the black-holes                                                                                                          
  //Simple N*N loop                                                                                                                                         
//  for(int i=1; i < nBH; i++)                                                                                                                              
//  {                                                                                                                                                       
  //  for(int j=1; j < nBH; j++)                                                                                                                            
    {                                                                                                                                                       
      int i=1; int j=2;                                                                                                                                     
      if(i != j)                                                                                                                                            
      {                                                                                                                                                     
        real4 posi = bhDistance[i*2+0];                                                                                                                     
        real4 posj = bhDistance[j*2+0];                                                                                                                     
                                                                                                                                                            
        real4 veli = bhDistance[i*2+1];                                                                                                                     
        real4 velj = bhDistance[j*2+1];                                                                                                                     
                                                                                                                                                            
        real4 r;                                                                                                                                            
        r.x = (posi.x-posj.x); r.y = (posi.y-posj.y);r.z = (posi.z-posj.z);                                                                                 
                                                                                                                                                            
        float dist = (r.x*r.x) + (r.y*r.y) + (r.z*r.z);                                                                                                     
        dist = sqrt(dist);                                                                                                                                  
                                                                                                                                                            
        //cerr << "BH_Distance: " << dist << " bh mass: " << posi.w << " and: " << posj.w << "\t kill distance: " << starSoftening*2 << endl;               
        cerr << (iter-1) << "\t"  <<  posi.x << "\t" << posi.y << "\t" << posi.z << "\t" <<  posj.x << "\t" << posj.y << "\t" << posj.z << "\t" << dist << "\t" << starSoftening;                                                                                                                                       
        cerr << "\t" <<  veli.x << "\t" << veli.y << "\t" << veli.z << "\t" <<  velj.x << "\t" << velj.y << "\t" << velj.z << "\t" << dE <<  "\t" << de_max <<                
               "\t" << t_current << "\t" << "BHDIST" << endl;                                                                                                                   
                                                                                                                                                            
        //Some test here to see if they are close enough,
        //use softening of the star particles
       if(dist < (starSoftening))
          return 1; 
      }
    }                                                                                                                                                       
  return 0;                                                                                                                                                 
}    


void octree::checkRemovalDistance(tree_structure &tree)                                                                                                     
{                                                                                                                                                           
  //Download all particle properties to the host                                                                                                            
                                                                                                                                                            
  tree.bodies_pos.d2h();    //The particles positions                                                                                                       
  tree.bodies_key.d2h();    //The particles keys                                                                                                            
  tree.bodies_vel.d2h();    //Velocities                                                                                                                    
  tree.bodies_acc0.d2h();    //Acceleration                                                                                                                 
  tree.bodies_acc1.d2h();    //Acceleration                                                                                                                 
  tree.bodies_time.d2h();  //The timestep details (.x=tb, .y=te                                                                                             
  tree.bodies_ids.d2h();                                                                                                                                    
                                                                                                                                                            
  bool modified = false;                                                                                                                                    
                                                                                                                                                            
  tree.multipole.d2h();                                                                                                                                     
  real4 com = tree.multipole[0];                                                                                                                            
                                                                                                                                                            
  int storeIdx = 0;               
  
  int NTotalT = 0, NFirstT = 0, NSecondT = 0, NThirdT = 0;
                                                                                                                                                            
  for(int i=0; i < tree.n ; i++)                                                                                                                            
  {                                                                                                                                                         
    real4 posi = tree.bodies_pos[i];                                                                                                                        
                                                                                                                                                            
    real4 r;                                                                                                                                                
    r.x = (posi.x-com.x); r.y = (posi.y-com.y);r.z = (posi.z-com.z);                                                                                        
    float dist = (r.x*r.x) + (r.y*r.y) + (r.z*r.z);                                                                                                         
    dist = sqrt(dist);                                                                                                                                      
                                                                                                                                                            
    tree.bodies_pos[storeIdx] = tree.bodies_pos[i];                                                                                                         
    tree.bodies_key[storeIdx] = tree.bodies_key[i];                                                                                                         
    tree.bodies_vel[storeIdx] = tree.bodies_vel[i];                                                                                                         
    tree.bodies_acc0[storeIdx] = tree.bodies_acc0[i];                                                                                                       
    tree.bodies_acc1[storeIdx] = tree.bodies_acc1[i];                                                                                                       
    tree.bodies_time[storeIdx] = tree.bodies_time[i];                                                                                                       
    tree.bodies_ids[storeIdx] = tree.bodies_ids[i];                

    if(dist > removeDistance)                                                                                                                               
    {                                                                                                                                                       
        //Remove this particle                                                                                                                              
        cerr << "Removing particle: " << i << " distance is: " << dist;                                                                                     
        cerr << "\tPOSM: " << posi.x << " " << posi.y << " " << posi.z << " " << posi.w;                                                                    
        cerr << "\tCOM: " << com.x << " " << com.y << " " << com.z << " " << com.w << endl;                                                                 
                                                                                                                                                            
        //Add this particles potential energy to the sum                                                                                                    
//         removedPot += hostbodies[i].w*0.5*hostacc0[i].w;                                                                                                 
        modified =  true;                                                                                                                                   
    }                                                                                                                                                       
    else                                                                                                                                                    
    {                                                                                                                                                       
      storeIdx++; //Increase the store position           

      NTotalT++;
      NFirstT = 0, NSecondT = 0, NThirdT = 0;    
      
      //Specific for Jeroens files
      if(tree.bodies_ids[i] >= 0 && tree.bodies_ids[i] < 100000000) NThirdT++;
      if(tree.bodies_ids[i] >= 100000000 && tree.bodies_ids[i] < 200000000) NSecondT++;
      if(tree.bodies_ids[i] >= 200000000 && tree.bodies_ids[i] < 300000000) NFirstT++;      
    }                                                                                                                                                       
  } //end for loop           


  NTotal  = NTotalT;
  NFirst  = NFirstT;
  NSecond = NSecondT;
  NThird  = NThirdT;

                                                                                                                                                            
  if(modified)                                                                                                                                              
  {                                                                                                                                                         
    tree.setN(storeIdx);                                                                                                                                    
                                                                                                                                                            
    //Now copy them back!!! Duhhhhh                                                                                                                         
    tree.bodies_pos.h2d();    //The particles positions                                                                                                     
    tree.bodies_key.h2d();    //The particles keys                                                                                                          
    tree.bodies_vel.h2d();    //Velocities                                                                                                                  
    tree.bodies_acc0.h2d();    //Acceleration                                                                                                               
    tree.bodies_acc1.h2d();    //Acceleration                                                                                                               
    tree.bodies_time.h2d();  //The timestep details (.x=tb, .y=te                                                                                           
    tree.bodies_ids.h2d();                                                                                                                                  
                                                                                                                                                            
    //Compute the energy!                                                                                                                                   
    store_energy_flag = true;                                                                                                                               
    compute_energies(tree);                                                                                                                                
  }//end if modified                                                                                                                                        
  else                                                                                                                                                      
  {                                                                                                                                                         
        cerr << "Nothing removed! :-) \n";                                                                                                                  
  }                   
  
  //TODO sync the number of particles with process 0 for correct header file
  
  
}                                                                                                                                                           
     





#if 1
 //Double precision
double octree::compute_energies(tree_structure &tree)
{
  Ekin = 0.0; Epot = 0.0;

  #if 0
  double hEkin = 0.0;
  double hEpot = 0.0;

  tree.bodies_pos.d2h();
  tree.bodies_vel.d2h();
  tree.bodies_acc0.d2h();
  for (int i = 0; i < tree.n; i++) {
    float4 vel = tree.bodies_vel[i];
    hEkin += tree.bodies_pos[i].w*0.5*(vel.x*vel.x +
                               vel.y*vel.y +
                               vel.z*vel.z);
    hEpot += tree.bodies_pos[i].w*0.5*tree.bodies_acc0[i].w;
    
    fprintf(stderr, "%d\t Vel: %f %f %f Mass: %f Pot: %f \n", 
            i,vel.x, vel.y, vel.z,tree.bodies_pos[i].w, tree.bodies_acc0[i].w);

  }
  MPI_Barrier(MPI_COMM_WORLD);
  double hEtot = hEpot + hEkin;
  printf("Energy (on host): Etot = %.10lg Ekin = %.10lg Epot = %.10lg \n", hEtot, hEkin, hEpot);
  #endif

  //float2 energy : x is kinetic energy, y is potential energy
  int blockSize = NBLOCK_REDUCE ;
  my_dev::dev_mem<double2>  energy(devContext);
  energy.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                      tree.generalBuffer1.get_flags(), 
                      tree.generalBuffer1.get_devMem(),
                      &tree.generalBuffer1[0], 0,  
                      blockSize, getAllignmentOffset(0));    
    
  computeEnergy.set_arg<int>(0,    &tree.n);
  computeEnergy.set_arg<cl_mem>(1, tree.bodies_pos.p());
  computeEnergy.set_arg<cl_mem>(2, tree.bodies_vel.p());
  computeEnergy.set_arg<cl_mem>(3, tree.bodies_acc0.p());
  computeEnergy.set_arg<cl_mem>(4, energy.p());
  computeEnergy.set_arg<double>(5, NULL, 128*2); //Dynamic shared memory, equal to number of threads times 2

  computeEnergy.setWork(-1, 128, blockSize);
  computeEnergy.execute();

  //Reduce the last parts on the host
  energy.d2h();
  Ekin = energy[0].x;
  Epot = energy[0].y;
  for (int i = 1; i < blockSize ; i++)
  {
      Ekin += energy[i].x;
      Epot += energy[i].y;
  }
  
  //Sum the values / energies of the system using MPI
  AllSum(Epot); AllSum(Ekin);
  
  Etot = Epot + Ekin;

  if (store_energy_flag) {
    Ekin0 = Ekin;
    Epot0 = Epot;
    Etot0 = Etot;
    Ekin1 = Ekin;
    Epot1 = Epot;
    Etot1 = Etot;
    tinit = get_time();
    store_energy_flag = false;
  }

  
  double de = (Etot - Etot0)/Etot0;
  double dde = (Etot - Etot1)/Etot1;

  if(tree.n_active_particles == tree.n)
  {
    de_max  = std::max( de_max, std::abs( de));
    dde_max = std::max(dde_max, std::abs(dde));
  }  
  
  Ekin1 = Ekin;
  Epot1 = Epot;
  Etot1 = Etot;
  
  if(mpiGetRank() == 0)
  {
  printf("iter=%d : time= %lg  Etot= %.10lg  Ekin= %lg   Epot= %lg : de= %lg ( %lg ) d(de)= %lg ( %lg ) t_sim=  %lg sec\n",
		  iter, this->t_current, Etot, Ekin, Epot, de, de_max, dde, dde_max, get_time() - tinit);  
  fprintf(stderr,"iter=%d : time= %lg  Etot= %.10lg  Ekin= %lg   Epot= %lg : de= %lg ( %lg ) d(de)= %lg ( %lg ) t_sim=  %lg sec\n", 
		  iter, this->t_current, Etot, Ekin, Epot, de, de_max, dde, dde_max, get_time() - tinit);          
  }
  
  return de;
}
#endif

#if 0
//Single precision version
void octree::compute_energies(tree_structure &tree)
{
  Ekin = 0.0; Epot = 0.0;

 /* double hEkin = 0.0;
  double hEpot = 0.0;

  tree.bodies_pos.d2h();
  tree.bodies_vel.d2h();
  tree.bodies_acc0.d2h();
  for (int i = 0; i < tree.n; i++) {
    float4 vel = tree.bodies_vel[i];
    hEkin += tree.bodies_pos[i].w*0.5*(vel.x*vel.x +
                               vel.y*vel.y +
                               vel.z*vel.z);
    hEpot += tree.bodies_pos[i].w*0.5*tree.bodies_acc0[i].w;
  }
  double hEtot = hEpot + hEkin;
  printf("Energy (on host): Etot = %.10lg Ekin = %.10lg Epot = %.10lg \n", hEtot, hEkin, hEpot);
  */

  //float2 energy : x is kinetic energy, y is potential energy
  int blockSize = NBLOCK_REDUCE ;
//   my_dev::dev_mem<float2>  energy(devContext, blockSize);
  my_dev::dev_mem<float2>  energy(devContext);
  energy.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                      tree.generalBuffer1.get_flags(), 
                      tree.generalBuffer1.get_devMem(),
                      &tree.generalBuffer1[0], 0,  
                      blockSize, getAllignmentOffset(0));    
                      
  computeEnergy.set_arg<int>(0,    &tree.n);
  computeEnergy.set_arg<cl_mem>(1, tree.bodies_pos.p());
  computeEnergy.set_arg<cl_mem>(2, tree.bodies_vel.p());
  computeEnergy.set_arg<cl_mem>(3, tree.bodies_acc0.p());
  computeEnergy.set_arg<cl_mem>(4, energy.p());
  computeEnergy.set_arg<float>(5, NULL, 2*128); //Dynamic shared memory

  computeEnergy.setWork(-1, 128, blockSize);
  computeEnergy.execute();

  //Reduce the last parts on the host
  energy.d2h();
  Ekin = energy[0].x;
  Epot = energy[0].y;
  for (int i = 1; i < blockSize ; i++)
  {
      Ekin += energy[i].x;
      Epot += energy[i].y;
  }

  Etot = Epot + Ekin;

  if (store_energy_flag) {
    Ekin0 = Ekin;
    Epot0 = Epot;
    Etot0 = Etot;
    Ekin1 = Ekin;
    Epot1 = Epot;
    Etot1 = Etot;
    tinit = get_time();
    store_energy_flag = false;
  }
  double de = (Etot - Etot0)/Etot0;
  double dde = (Etot - Etot1)/Etot1;
  Ekin1 = Ekin;
  Epot1 = Epot;
  Etot1 = Etot;
  printf("iter=%d : time= %lg  Etot= %.10lg  Ekin= %lg   Epot= %lg : de= %lg d(de)= %lg t_sim= %lg sec\n",
          iter, this->t_current, Etot, Ekin, Epot, de, dde, get_time() - tinit);

}
#endif
