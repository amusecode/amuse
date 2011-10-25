#include "octree.h"

#include <iostream>
#include <algorithm>
#include <iomanip>
using namespace std;

static double de_max = 0;
static double dde_max = 0;  


my_dev::dev_stream *execStream;

void octree::makeLET()
{
   //LET code test
  double tTest = get_time();

//   int prevN            = remoteTree.n;
//   int prevNodes        = remoteTree.n_nodes;
 
  my_dev::dev_stream memCpyStream;
  
  localTree.bodies_Ppos.d2h(false, memCpyStream.s());
  localTree.bodies_Pvel.d2h(false, memCpyStream.s());
  localTree.multipole.d2h(false, memCpyStream.s());   
  localTree.boxSizeInfo.d2h(false, memCpyStream.s());
  localTree.boxCenterInfo.d2h(false, memCpyStream.s());
  
  
//   Hier gebleven, checken wat er precies met xlow en xhigh en dergelijke gebeurd
//   kijken of die wel hetzelfde zijn als we LET exchange doen met en zonder radius find

  //Exchange domain boundaries, while memory copies take place
  rMinLocalTree.w = this->maxLocalEps;
  getAllLETBoxes(rMinLocalTree,rMaxLocalTree);

  memCpyStream.sync();
  
//   fprintf(stderr, "LET (%d) 3: %g \n", mpiGetRank(), get_time() - tTest);

  vector<real4> LETParticles;
  essential_tree_exchange(LETParticles, localTree, remoteTree);

  fprintf(stderr, "LET Exchange took (%d): %g \n", mpiGetRank(), get_time() - tTest);
      

 /*  if(mpiGetRank() == 0)
   {
//      remote.remoteTreeStruct.x = Nparticles;
//     remote.remoteTreeStruct.y = Nnodes;
//     remote.remoteTreeStruct.z = 0;
//     remote.remoteTreeStruct.w = totalTopNodes;
    FILE *norm = fopen("norm.txt", "w");
    FILE *newf = fopen("new.txt", "w");     
     
    for(int i=0; i < remoteTree.remoteTreeStruct.x; i++)        //Particles
    {
       fprintf(norm,"Part %d: %f %f %f %f\n", i,
               remoteTree.bodies_Ppos[i].x, remoteTree.bodies_Ppos[i].y,
               remoteTree.bodies_Ppos[i].z, remoteTree.bodies_Ppos[i].w);
       fprintf(newf,"Part %d: %f %f %f %f\n", i,               
               remoteTree.fullRemoteTest[i].x, remoteTree.fullRemoteTest[i].y,
               remoteTree.fullRemoteTest[i].z, remoteTree.fullRemoteTest[i].w);   
    }
    
    for(int i=0; i < remoteTree.remoteTreeStruct.y; i++)        //Nodes
    {
       fprintf(norm,"BSI %d: %f %f %f %f\n", i,
               remoteTree.boxSizeInfo[i].x, remoteTree.boxSizeInfo[i].y,
               remoteTree.boxSizeInfo[i].z, remoteTree.boxSizeInfo[i].w);
       fprintf(newf,"BSI %d: %f %f %f %f\n", i,               
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + i].x, remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + i].y,
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + i].z, remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + i].w);   
    }    
     
    for(int i=0; i < remoteTree.remoteTreeStruct.y; i++)        //Nodes
    {
       fprintf(norm,"BCI %d: %f %f %f %f\n", i,
               remoteTree.boxCenterInfo[i].x, remoteTree.boxCenterInfo[i].y,
               remoteTree.boxCenterInfo[i].z, remoteTree.boxCenterInfo[i].w);
       fprintf(newf,"BCI %d: %f %f %f %f\n", i,               
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + remoteTree.remoteTreeStruct.y + i].x,
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + remoteTree.remoteTreeStruct.y + i].y,               
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + remoteTree.remoteTreeStruct.y + i].z,
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + remoteTree.remoteTreeStruct.y + i].w);   
    }
    
    for(int i=0; i < remoteTree.remoteTreeStruct.y; i++)        //Nodes
    {
       fprintf(norm,"MP %d: %f %f %f %f\n", i,
               remoteTree.multipole[i*3 + 0].x, remoteTree.multipole[i*3 + 0].y,
               remoteTree.multipole[i*3 + 1].z, remoteTree.multipole[i*3 + 2].w);
       fprintf(newf,"MP %d: %f %f %f %f\n", i,               
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + 2*remoteTree.remoteTreeStruct.y + i*3 + 0].x,
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + 2*remoteTree.remoteTreeStruct.y + i*3 + 0].y,               
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + 2*remoteTree.remoteTreeStruct.y + i*3 + 1].z,
               remoteTree.fullRemoteTest[remoteTree.remoteTreeStruct.x + 2*remoteTree.remoteTreeStruct.y + i*3 + 2].w);   
    }    
    
    fclose(norm);
    fclose(newf);
     
     
//      for(int i=0; i < remoteTree.n_nodes; i++)
//      for(int i=0; i < 40; i++)
//      {
//       union{float f; int i;} u; //__float_as_int
//       u.f           = remoteTree.boxSizeInfo[i].w;      
//       int childinfo = u.i;
// 
//        fprintf(stderr,"Part %d: %f\t%f\t%f\t%d \n", i,
//                remoteTree.boxSizeInfo[i].x, remoteTree.boxSizeInfo[i].y,
//                remoteTree.boxSizeInfo[i].z, childinfo);
// //                remoteTree.multipole[i*3].x, remoteTree.multipole[i*3].y,
// //                remoteTree.multipole[i*3].y, remoteTree.multipole[i*3].w);
//      }
   }
   
   fprintf(stderr,"Total nodes: %d \n", remoteTree.n_nodes);                      
   
   
   MPI_Barrier(MPI_COMM_WORLD);
  exit(0);*/
 
  remoteTree.fullRemoteTest.h2d(false, memCpyStream.s());//Async
  memCpyStream.sync();
  
//   devContext.startTiming();
  approximate_gravity_let(this->localTree, this->remoteTree);
//   devContext.stopTiming("Approximation_let", 5);


  #if 1
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
      
      nodeFile << i << "\t" << remoteTree.fullRemoteTest[centIdx].x-remoteTree.fullRemoteTest[sizeIdx].x << "\t" << remoteTree.fullRemoteTest[centIdx].y-remoteTree.fullRemoteTest[sizeIdx].y;
      
      if(remoteTree.fullRemoteTest[sizeIdx].x == 0 && remoteTree.fullRemoteTest[sizeIdx].y == 0)
      {
        nodeFile << "\t" << remoteTree.fullRemoteTest[centIdx].x+0.01 << "\t" << remoteTree.fullRemoteTest[centIdx].y+0.01 << "\t";
      }
      else
      {
        nodeFile << "\t" << remoteTree.fullRemoteTest[centIdx].x+remoteTree.fullRemoteTest[sizeIdx].x << "\t" << remoteTree.fullRemoteTest[centIdx].y+remoteTree.fullRemoteTest[sizeIdx].y << "\t";
      }
      
            
      
      
//       nodeFile << i << "\t" << remoteTree.fullRemoteTest[centIdx].x << "\t" << remoteTree.fullRemoteTest[centIdx].y;
//       nodeFile << "\t" <<      remoteTree.fullRemoteTest[sizeIdx].x << "\t" << remoteTree.fullRemoteTest[sizeIdx].y << "\t";
      nodeFile <<              remoteTree.fullRemoteTest[multIdx].x << "\t" 
               << remoteTree.fullRemoteTest[multIdx].y << "\t" << remoteTree.fullRemoteTest[multIdx].w << "\n";
      
      massTEST +=  remoteTree.fullRemoteTest[multIdx].w;
    }

    nodeFile.close();
    

    sprintf(fileName, "LETTreeStructureParticles-%d.txt", mpiGetRank());
    ofstream partFile;
    partFile.open(fileName);

                                                                                                                    
    for(int i=0; i < remoteP; i++)                                                                                     
    {                                                                                                                 
      float4  pos =  remoteTree.fullRemoteTest[i];                                                                              
      partFile << i << "\t" << pos.x << "\t" << pos.y << "\t" << pos.z << endl;       
      
      massTEST += pos.w;
    }                                                                                                                 
    partFile.close(); 
    
    cerr << "Mass test (rank= " << mpiGetRank() << " ) : " << massTEST << std::endl;
  
//     exit(0);

  
  #endif

}


void octree::iterate() {

  int Nact_since_last_tree_rebuild = 0;
  real4 r_min, r_max;
  
  execStream = new my_dev::dev_stream(0);

  //Initial prediction/acceleration to setup the system
  //Will be at time 0
  //predict localtree
  predict(this->localTree);
  this->getBoundaries(localTree, r_min, r_max);
  //Build the tree using the predicted positions  
  //Compute the (new) node properties
  compute_properties(this->localTree);

  //Approximate gravity
//   double tTest = get_time();

  devContext.startTiming();
  approximate_gravity(this->localTree);
  devContext.stopTiming("Approximation", 4);
//   fprintf(stderr, "Launch of approx took: %g \n", get_time() - tTest);  


  if(nProcs > 1)  makeLET();

//   fprintf(stderr, "LET + approx took: %g \n", get_time() - tTest);  
  //Corrector
  correct(this->localTree);
  
//   exit(0);

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

        //write_dumbp_snapshot(&localTree.bodies_pos[0], &localTree.bodies_vel[0], &localTree.bodies_ids[0], localTree.n, fileName.c_str()) ;
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
    //Get updated domain boundaries, for constructing the LET tre
    this->getBoundaries(localTree, r_min, r_max);
    devContext.stopTiming("Predict", 9);
    
    //Build the tree using the predicted positions
   // bool rebuild_tree = Nact_since_last_tree_rebuild > 4*this->localTree.n;

    if(1)
    //if(rebuild_tree)
    {
      //Rebuild the tree
      this->sort_bodies(this->localTree);

      devContext.startTiming();
      this->build(this->localTree);
      devContext.stopTiming("Tree-construction", 2);

      devContext.startTiming();
      this->allocateTreePropMemory(this->localTree);
      devContext.stopTiming("Memory", 11);      

      devContext.startTiming();
      this->compute_properties(this->localTree);
      devContext.stopTiming("Compute properties", 3);

      devContext.startTiming();
      setActiveGrpsFunc(this->localTree);
      devContext.stopTiming("setActiveGrpsFunc", 10);      
      Nact_since_last_tree_rebuild = 0;
    }
    else
    {
      devContext.startTiming();
      this->compute_properties(this->localTree);
      devContext.stopTiming("Compute properties", 3);
    }//end rebuild tree

    //Approximate gravity
    devContext.startTiming();
    approximate_gravity(this->localTree);
    devContext.stopTiming("Approximation", 4);
    
    if(nProcs > 1)  makeLET();
    
    //Corrector
    devContext.startTiming();
    correct(this->localTree);
    devContext.stopTiming("Correct", 8);

    Nact_since_last_tree_rebuild += this->localTree.n_active_particles;

    //Compute energies
    devContext.startTiming();
    double de = compute_energies(this->localTree); de=de;
    devContext.stopTiming("Energy", 7);
    
    //Redistribute the particles
    if(1)
    {      
      if(nProcs > 1)
      { 
//         if(iter % 5 == 0)
        if(1)
        {       
          devContext.startTiming();

          localTree.bodies_pos.d2h();
          localTree.bodies_vel.d2h();
          localTree.bodies_acc0.d2h();
          localTree.bodies_acc1.d2h();
          localTree.bodies_time.d2h();
          localTree.bodies_ids.d2h();
          
          devContext.stopTiming("Exchange_copy_d2h", 6);
          
          devContext.startTiming();
          //Update the distribution and then exchange particles        
          updateDistribution(&localTree.bodies_pos[0], localTree.n);        
          while(exchange_particles_with_overflow_check(localTree));
          devContext.stopTiming("Exchange", 6);
          
          devContext.startTiming();
          //Send the new and old particles to the device
          localTree.bodies_pos.h2d();
          localTree.bodies_vel.h2d();
          localTree.bodies_acc0.h2d();
          localTree.bodies_acc1.h2d();
          localTree.bodies_time.h2d();
          localTree.bodies_ids.h2d();
          
          localTree.body2group_list.zeroMem();
        
          devContext.stopTiming("Exchange_copy_h2d", 6);
        }
        else
        {
          //Get the new domain boundaries
//           real4 r_min, r_max;
//           this->getBoundaries(localTree, r_min, r_max);
        } //if (iter % X )
      } //if nProcs > 1
    }//if (0)    
    

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

        //write_dumbp_snapshot(&localTree.bodies_pos[0], &localTree.bodies_vel[0], &localTree.bodies_ids[0], localTree.n, fileName.c_str()) ;
        write_dumbp_snapshot_parallel(&localTree.bodies_pos[0], &localTree.bodies_vel[0],
                                      &localTree.bodies_ids[0], localTree.n, fileName.c_str()) ;
      }
    }

  //  checkMergingDistance(this->localTree, iter, de);    

    if(t_current >= tEnd)
    {
      compute_energies(this->localTree);
      cout << " Finished: "  << t_current << "  > "  << tEnd << " loop alone took: " << get_time() -t0 <<  endl;
      return;
    }
    
    
    if((iter % 50) == 0)
    {
//       if(removeDistance > 0) checkRemovalDistance(this->localTree);
    }
    
    
    

    iter++;
    
  } //end for i
} //end iterate


void octree::predict(tree_structure &tree)
{
  //Functions that predicts the particles to the next timestep

//   tend is time per particle
//   tnext is reduce result

  //First we get the minimum time, which is the next integration
  //time
  int blockSize = NBLOCK_REDUCE ;
  getTNext.set_arg<int>(0,    &tree.n);
  getTNext.set_arg<cl_mem>(1, tree.bodies_time.p());
  getTNext.set_arg<cl_mem>(2, tnext.p());
  getTNext.set_arg<float>(3, NULL, 128); //Dynamic shared memory

  getTNext.setWork(-1, 128, blockSize);
  getTNext.execute();

  //Reduce the last parts on the host
  tnext.d2h();
  t_previous = t_current;
  t_current = tnext[0];
  for (int i = 1; i < blockSize ; i++)
  {
      t_current = fmin(t_current, tnext[i]);
  }

  tree.activeGrpList.zeroMem();      //Reset the active grps


  //Set valid list to zero
  predictParticles.set_arg<int>(0,    &tree.n);
  predictParticles.set_arg<float>(1,  &t_current);
  predictParticles.set_arg<float>(2,  &t_previous);
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
  
  
//   if(mpiGetRank() == 1)
//   {
//     tree.body2group_list.d2h();
//     for(int i=0; i < tree.n; i++)
//     {
//        fprintf(stderr, "Body to grp: %d\t %d \n", i, tree.body2group_list[i]);
//     }
//   }
  
  

  //Compact the valid list to get a list of valid groups
  gpuCompact(devContext, tree.activeGrpList, tree.active_group_list,
             tree.n_groups, &tree.n_active_groups);

   printf("t_previous: %lg t_current: %lg dt: %lg Active groups: %d \n",
         t_previous, t_current, t_current-t_previous, tree.n_active_groups);
/*  printf("t_previous: %lg t_current: %lg dt: %lg Active groups: %d Active particles: %d \n",
         t_previous, t_current, t_current-t_previous, tree.n_active_groups, tree.n_active_particles);    */

}
//End predict


void octree::setActiveGrpsFunc(tree_structure &tree)
{
  tree.activeGrpList.zeroMem();      //Reset the active grps

  //Set valid list to zero
  setActiveGrps.set_arg<int>(0,    &tree.n);
  setActiveGrps.set_arg<float>(1,  &t_current);
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


#if 1

/*

void create_local_essential_tree_host(int node, 
                                      const real4 *bodies, 
                                      real4 *nodeCenterInfo,
                                      real4 *nodeSizeInfo, 
                                      real4 *multipole,                                       
                                      vector<real4> &essentialList, 
                                      int iter, 
                                      double &mass)
{

  //Read node data
  
    real4 nodeCenter = nodeCenterInfo[node];
    real4 nodeSize   = nodeSizeInfo[node];
    bool leaf        = nodeCenter.w <= 0;
    
    union{float f; int i;} u; //__float_as_int
    u.f           = nodeSize.w;      
    int childinfo = u.i;
    
    int child, nchild;
    if(!leaf)
    {
      //Node
      child    =    childinfo & 0x0FFFFFFF;                         //Index to the first child of the node
      nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              
    }
    else
    {
      //Leaf
      child   =   childinfo & BODYMASK;                                     //thre first body in the leaf
      nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag        
    }
    
   double massT = multipole[node*3 + 0].w;
   double idxT  = multipole[node*3 + 0].x;
   fprintf(stderr, "Node: %d leaf: %d  child: %d nchild: %d (%d, %f, %f) (%f %f) tm: %f\n",
                    node, leaf, child, nchild, childinfo, nodeCenter.w, nodeCenter.x,
                    idxT, massT, mass);
        

//    if(1) if(0)
    if(childinfo == 0) //All particles, eg always split
    {
      //Node is far enough away from the box we can use its cm (node and box are well separated)
//       essentialList.push_back(multipole[3*node + 0]); //push back the com, not the node position
      mass += multipole[3*node + 0].w;
    }
    else
    {
      if (leaf)
      {
        for (int i = child; i < child+nchild; i++)
        {
          const real4 p = bodies[i];     //Body position
          essentialList.push_back(p);

          mass += p.w;
        }
    }
    else
    {
      //Continue the tree-walk

      //Do for each child
      for(int ic = child; ic < child+nchild; ic++)
      {
        create_local_essential_tree_host(ic, bodies, nodeCenterInfo, nodeSizeInfo, multipole, essentialList, iter, mass);
      }
      fprintf(stderr, " Temp sum: %f \n", mass);
    } //end if is leaf
  } //end if sep
}//end create_local_essential_tree
*/

//Version that prints dot structure
void create_local_essential_tree_host(int node, 
                                      const real4 *bodies, 
                                      real4 *nodeCenterInfo,
                                      real4 *nodeSizeInfo, 
                                      real4 *multipole,                                       
                                      vector<real4> &essentialList, 
                                      int iter, 
                                      double &mass)
{

  //Read node data
  
    real4 nodeCenter = nodeCenterInfo[node];
    real4 nodeSize   = nodeSizeInfo[node];
    bool leaf        = nodeCenter.w <= 0;
    
//     int nchildren = 0;
    
    union{float f; int i;} u; //__float_as_int
    u.f           = nodeSize.w;      
    int childinfo = u.i;
    
    int child, nchild;
    if(!leaf)
    {
      //Node
      child    =    childinfo & 0x0FFFFFFF;                         //Index to the first child of the node
      nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              
    }
    else
    {
      //Leaf
      child   =   childinfo & BODYMASK;                                     //thre first body in the leaf
      nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag        
    }
    
//    double massT = multipole[node*3 + 0].w;
//    double idxT  = multipole[node*3 + 0].x;
   
//    if(1) if(0)
    if(childinfo == 0) //All particles, eg always split
    {
      //Node is far enough away from the box we can use its cm (node and box are well separated)
//       essentialList.push_back(multipole[3*node + 0]); //push back the com, not the node position
      mass += multipole[3*node + 0].w;
    }
    else
    {
      if (leaf)
      {
        for (int i = child; i < child+nchild; i++)
        {
          const real4 p = bodies[i];     //Body position
          essentialList.push_back(p);

          mass += p.w;
          
          fprintf(stderr, "%d -> %d ; \n", node, -1*i );  
          fprintf(stderr, "%d [shape=box]; \n", -1*i);  
        }
    }
    else
    {
      //Continue the tree-walk

      //Do for each child
      for(int ic = child; ic < child+nchild; ic++)
      {
        create_local_essential_tree_host(ic, bodies, nodeCenterInfo, nodeSizeInfo, multipole, essentialList, iter, mass);
        
        fprintf(stderr, "%d -> %d ; \n", node, ic );        
        
      }
//       fprintf(stderr, " Temp sum: %f \n", mass);
    } //end if is leaf
  } //end if sep
}//end create_local_essential_tree

//Version that prints dot structure andkeeps track of number of children
void create_local_essential_tree_host_childr(int node, 
                                      const real4 *bodies, 
                                      real4 *nodeCenterInfo,
                                      real4 *nodeSizeInfo, 
                                      real4 *multipole,                                       
                                      vector<real4> &essentialList, 
                                      int iter, 
                                      double &mass,
                                      int &nchildren)
{

  //Read node data
  
    real4 nodeCenter = nodeCenterInfo[node];
    real4 nodeSize   = nodeSizeInfo[node];
    bool leaf        = nodeCenter.w <= 0;
    
    int nchildren_loc = 0;
    
    union{float f; int i;} u; //__float_as_int
    u.f           = nodeSize.w;      
    int childinfo = u.i;
    
    int child, nchild;
    if(!leaf)
    {
      //Node
      child    =    childinfo & 0x0FFFFFFF;                         //Index to the first child of the node
      nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              
    }
    else
    {
      //Leaf
      child   =   childinfo & BODYMASK;                                     //thre first body in the leaf
      nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag        
    }
    
//    double massT = multipole[node*3 + 0].w;
//    double idxT  = multipole[node*3 + 0].x;
   
//    if(1) if(0)
    if(childinfo == 0) //All particles, eg always split
    {
      //Node is far enough away from the box we can use its cm (node and box are well separated)
//       essentialList.push_back(multipole[3*node + 0]); //push back the com, not the node position
      mass += multipole[3*node + 0].w;
    }
    else
    {
      if (leaf)
      {
        nchildren_loc += nchild;
        for (int i = child; i < child+nchild; i++)
        {
          const real4 p = bodies[i];     //Body position
          essentialList.push_back(p);

          mass += p.w;
          
          fprintf(stderr, "%d -> %d ; \n", node, -1*i );  
          fprintf(stderr, "%d [shape=box]; \n", -1*i);  
        }
    }
    else
    {
      //Continue the tree-walk
          
      //Do for each child
      for(int ic = child; ic < child+nchild; ic++)
      {
        create_local_essential_tree_host_childr(ic, bodies, nodeCenterInfo, nodeSizeInfo, multipole, essentialList, iter, mass, nchildren_loc);
        
        fprintf(stderr, "%d -> %d ; \n", node, ic );                
      }
    } //end if is leaf
  } //end if sep
  
  char buff[256];
  sprintf(buff, "%d ( %d )", node, nchildren_loc);
  fprintf(stderr, "%d [label=\"%s\"]; \n", node, buff);  
//   a [label="Foo"];fprintf(stderr, "%d [shape=box]; \n", -1*i);  

  
  
  nchildren += nchildren_loc;
}//end create_local_essential_tree


void create_local_essential_tree_count_modify(real4* bodies, real4* multipole, real4* nodeSizeInfo, real4* nodeCenterInfo,   
                                              int start, int end, 
                                              int &particles, int &nodes)
{
    //Walk the tree as is done on device, level by level
    vector<int> curLevel;
    vector<int> nextLevel;
    
    vector<int> newTreeSizeInfo;
    vector<int> newTreeCenterInfo;

    int particleCount   = 0;
    int nodeCount       = 0;

    start = 0;
    end = 1;
    
    //Add the initial nodes to the curLevel list
    for(int i=start; i < end; i++)
    {
      curLevel.push_back(i);
    }
    
    //Add the nodes before the start and end to the node list
    for(int i=0; i < start; i++)
    {
      nodeCount++;
    }
      
    //Start the tree-walk       
    while(curLevel.size() > 0)
    {
      for(unsigned int i=0; i < curLevel.size(); i++)
      {
        //Read node data
        int node         = curLevel[i];
        real4 nodeCenter = nodeCenterInfo[node];
        real4 nodeSize   = nodeSizeInfo[node];
        bool leaf        = nodeCenter.w <= 0;
        
        union{float f; int i;} u; //__float_as_int
        u.f           = nodeSize.w;      
        int childinfo = u.i;
        
        int child, nchild;
        if(!leaf)
        {
          //Node
          child    =    childinfo & 0x0FFFFFFF;                         //Index to the first child of the node
          nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              
        }
        else
        {
          //Leaf
          child   =   childinfo & BODYMASK;                                     //thre first body in the leaf
          nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag        
        }
        
        bool split = true;

        //if split & node add children to next lvl stack
        if(split && !leaf)
        {       
          
          //Check if we can compress the child nodes
          int childSum = 0;
          
          int newChildStart = -1;
//           int nNewChild     = -1;
                          
          for(int i=child; i < child+nchild; i++)
          {
            real4 childCenter = nodeCenterInfo[i];
            real4 childSize   = nodeSizeInfo[i];
            bool  leaf        = childCenter.w <= 0;
            
            if(leaf)
            {
              //Leaf
              int grandchildinfo = host_float_as_int(childSize.w);
              int grandchild   =   grandchildinfo & BODYMASK;                                     //thre first body in the leaf
              int grandnchild  = (((grandchildinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag                
              
              //Check if we can merge this child with the previous child
              if((childSum + grandnchild) <= NLEAF)
              {
                //We can!
                if(newChildStart < 0) newChildStart = grandchild;   //First time we set it otherwise no need
              }
              else
              {
                //Sum becomes too large, write the new node if there is any
                
              }
              childSum += grandnchild;              
            } 
            else
            {
              //Check if we merged any previous nodes
              //if so finish them otherwise
            }
            nextLevel.push_back(i);            
          }
          
          fprintf(stderr, "%d has %d childs and %d grandchilds\n", node, nchild, childSum);
          
        }
        
        //if split & leaf add particles to particle list
        if(split && leaf)
        { 
          for(int i=child; i < child+nchild; i++)
          {
            particleCount++;
          }
        }

        //Increase the nodeCount, since this node will be part of the tree-structure
        nodeCount++;        
      } //end for curLevel.size
      
      
//       newTreeCenterInfo.push_back(nodeCenter);
//       newTreeSizeInfo.push_back(nodeSize);
      
      //Put next level stack into current level and continue
      curLevel.clear();
      
//       cout << "Next level: " << nextLevel.size() << endl;
      curLevel.assign(nextLevel.begin(), nextLevel.end());
      nextLevel.clear();
    }//end while
    
    particles = particleCount;
    nodes = nodeCount;  
    
    fprintf(stderr, "Count found: %d particles and %d nodes \n", particles, nodes);   
}

void walkTreeOnHostModify(tree_structure &tree, int iter)
{
  tree.bodies_Ppos.d2h();
  tree.boxSizeInfo.d2h();
  tree.boxCenterInfo.d2h();  
  tree.multipole.d2h();



  int nparticles = 0;
  int nnodes = 0;
//   for(int i=start; i < end; i++)
  {
    
//     create_local_essential_tree_host(0, &tree.bodies_Ppos[0], &tree.boxCenterInfo[0],
  create_local_essential_tree_count_modify(&tree.bodies_Ppos[0], &tree.multipole[0],
                                           &tree.boxSizeInfo[0], &tree.boxCenterInfo[0],
                                           0, 1, nparticles, nnodes);
  }

  fprintf(stderr,"Modify done (own tree)!!! selected particles: %d \tNodes: %d\n", 
          nparticles, nnodes);
  

  
}

void walkTreeOnHost(tree_structure &tree, int iter)
{
  tree.bodies_Ppos.d2h();
  tree.boxSizeInfo.d2h();
  tree.boxCenterInfo.d2h();  
  tree.multipole.d2h();

  vector<real4> essentialList;

  double mass = 0;

//   printf("Starting walk! From: %d to %d\n");
//   start = tree.level_list[level_start].x;
//   end   = tree.level_list[level_start].y;  

int nchildren = 0;
//   for(int i=start; i < end; i++)
  {
//     create_local_essential_tree_host(0, &tree.bodies_Ppos[0], &tree.boxCenterInfo[0],
  create_local_essential_tree_host_childr(0, &tree.bodies_Ppos[0], &tree.boxCenterInfo[0],
                                   &tree.boxSizeInfo[0],
                                   &tree.multipole[0], essentialList, 
                                   iter, mass, nchildren);
  }

  fprintf(stderr,"Walk done (own tree)!!! selected particles: %ld  Mss sum: %f\tIter: %d  \tChildren: %d\n", 
          essentialList.size(), mass, iter, nchildren);
  
  
  
//   if(1)
//   {
//     for(int node=0; node < 40; node++){
//       real4 nodeCenter = tree.boxCenterInfo[node];
//       real4 nodeSize   = tree.boxSizeInfo[node];
//       bool leaf        = nodeCenter.w <= 0;
//       
//       union{float f; int i;} u; //__float_as_int
//       u.f           = nodeSize.w;      
//       int childinfo = u.i;
//       
//       int child, nchild;
//       if(!leaf)
//       {
//         //Node
//         child    =    childinfo & 0x0FFFFFFF;                         //Index to the first child of the node
//         nchild   = (((childinfo & 0xF0000000) >> 28)) ;         //The number of children this node has              
//       }
//       else
//       {
//         //Leaf
//         child   =   childinfo & BODYMASK;                                     //thre first body in the leaf
//         nchild  = (((childinfo & INVBMASK) >> LEAFBIT)+1);     //number of bodies in the leaf masked with the flag        
//       }
//       
//       fprintf(stderr, "Node: %d leaf: %d  child: %d nchild: %d (%d, %f, %f)\n",
//               node, leaf, child, nchild, childinfo, nodeCenter.w, nodeCenter.x);
//     }
//         
//   }
    
  
}

void walkRemoteTreeOnHost(tree_structure &tree, int iter)
{
  tree.bodies_Ppos.d2h();
  tree.boxSizeInfo.d2h();
  tree.boxCenterInfo.d2h();  
  tree.multipole.d2h();
  
// remoteTree.fullRemoteTest  

  vector<real4> essentialList;

  double mass = 0;
  
  int start =  tree.remoteTreeStruct.z;
  int end =  tree.remoteTreeStruct.w;  

  fprintf(stderr,"Starting walk LET! start: %d end: %d\n", start, end);
  
  double massTEST = 0;
  
//  TODO Doe een loop over de begin en end nodes
 for(int i=start; i < end; i++)
 {
   double massTwo = 0;
   create_local_essential_tree_host(i, &tree.bodies_Ppos[0], &tree.boxCenterInfo[0],
                                   &tree.boxSizeInfo[0],
                                   &tree.multipole[0], essentialList, 
                                   iter, massTwo);
   mass         += massTwo;                                   
   massTEST     += tree.multipole[i*3+0].w;
   fprintf(stderr, "MASS Compare %d\t: %f %f \n", i, massTwo, tree.multipole[i*3+0].w);
 }

  fprintf(stderr,"Walk done (remote tree)!!!11 selected particles: %ld  Mass sum: %f\tIter: %d MASSTEST: %f\n", 
          essentialList.size(), mass, iter, massTEST);
}

#endif


void octree::approximate_gravity(tree_structure &tree)
{
//  if(mpiGetRank() == 0)
//       walkTreeOnHost(tree, iter);
//   
//   MPI_Barrier(MPI_COMM_WORLD);
//    exit(0);  
//    

  //Memory needed: sizeof(int)*2048*NTHREAD*number of blocks
  
//   unsigned long sizeLMEM = sizeof(int)*2048*NTHREAD*tree.n_groups;  
//   printf("Bytes needed: %ld  kb: %ld  mbyte: %ld \n", 
//          sizeLMEM, sizeLMEM / 1024, sizeLMEM / (1024*1024));

 
  uint2 node_begend;
  int level_start = 2;
  node_begend.x = tree.level_list[level_start].x;
  node_begend.y = tree.level_list[level_start].y;

  printf("node begend: %d %d iter-> %d\n", node_begend.x, node_begend.y, iter);

  //Reset the active particles
  tree.activePartlist.zeroMem();

//   int grpOffset = 0;
  
//   my_dev::dev_mem<uint>  MEM_BUF(devContext, 120*2048*64); //blocks*length*threads
  
//
  //Set the kernel parameters, many!
  approxGrav.set_arg<int>(0,    &tree.n_active_groups);
  approxGrav.set_arg<float4>(1, &tree.corner);
  approxGrav.set_arg<float>(2,  &(this->inv_theta));
  approxGrav.set_arg<float>(3,  &(this->eps2));
  approxGrav.set_arg<uint2>(4,  &node_begend);
  approxGrav.set_arg<cl_mem>(5, tree.active_group_list.p());
//   approxGrav.set_arg<cl_mem>(6, tree.node_data.p());
//   approxGrav.set_arg<cl_mem>(7, tree.bodies_pos.p());
  approxGrav.set_arg<cl_mem>(6, tree.bodies_Ppos.p());
  approxGrav.set_arg<cl_mem>(7, tree.multipole.p());
  approxGrav.set_arg<cl_mem>(8, tree.bodies_acc1.p());
  approxGrav.set_arg<cl_mem>(9, tree.ngb.p());
  approxGrav.set_arg<cl_mem>(10, tree.activePartlist.p());
//   approxGrav.set_arg<cl_mem>(12, tree.group_data.p());
  approxGrav.set_arg<cl_mem>(11, tree.interactions.p());
  approxGrav.set_arg<cl_mem>(12, tree.boxSizeInfo.p());
  approxGrav.set_arg<cl_mem>(13, tree.groupSizeInfo.p());
  approxGrav.set_arg<cl_mem>(14, tree.boxCenterInfo.p());
  approxGrav.set_arg<cl_mem>(15, tree.groupCenterInfo.p());
  approxGrav.set_arg<cl_mem>(16, tree.bodies_Pvel.p());
  approxGrav.set_arg<cl_mem>(17, tree.group_list.p()); //Check why/if we need this!
//   approxGrav.set_arg<cl_mem>(20, tree.peanoOrder.p()); //Check why/if we need this!
//   approxGrav.set_arg<int>(21,  &grpOffset); //Check why/if we need this!
//   approxGrav.set_arg<cl_mem>(22,  MEM_BUF.p()); //Check why/if we need this!

  approxGrav.set_arg<real4>(18, tree.boxSizeInfo, 4, "texNodeSize");
  approxGrav.set_arg<real4>(19, tree.boxCenterInfo, 4, "texNodeCenter");
  approxGrav.set_arg<real4>(20, tree.multipole, 4, "texMultipole");
  approxGrav.set_arg<real4>(21, tree.bodies_Ppos, 4, "texBody");


/*  approxGrav.set_arg<real4>(23, tree.boxSizeInfo, 4, "texNodeSize");
  approxGrav.set_arg<real4>(24, tree.boxCenterInfo, 4, "texNodeCenter");
  approxGrav.set_arg<real4>(25, tree.multipole, 4, "texMultipole");
  approxGrav.set_arg<real4>(26, tree.bodies_Ppos, 4, "texBody");*/
  
//Trick to get timing of big moddels
//   tree.bodies_key.free_mem();
//   tree.bodies_ids.free_mem();
//
//   tree.bodies_vel.free_mem();
//   tree.bodies_acc0.free_mem();
//   tree.bodies_time.free_mem();
//
//
//   tree.node_key.free_mem();
//   tree.n_children.free_mem();
//   tree.node_bodies.free_mem();
//
//    //Allocate memory
//   tree.lowerBounds.free_mem();
//   tree.upperBounds.free_mem();
//   tree.nodeLowerBounds.free_mem();
//   tree.nodeUpperBounds.free_mem();

//    Build-in usage: free: 215203840 bytes ( 205 MB , total: 1535)
// Interaction at iter: 0  direct: 1951154601      appr: 6126429632        avg dir: 232    avg appr: 730
// Active particles: 8388608
// Approximation took:     1832.446777      millisecond
// iter=0 : time= 0  Etot= -0.2459316656  Ekin= 0.249861   Epot= -0.495792 : de= -0 d(de)= -0 t_sim= 0 sec
// At the start of iterate:

//   approxGrav.setWork(-1, 64, tree.n_active_groups);
 //  approxGrav.setWork(-1, NCRIT, tree.n_active_groups);
//   approxGrav.setWork(-1, NCRIT, tree.n_active_groups);

//     approxGrav.setWork(-1, NTHREAD, tree.n_active_groups / 2);
// approxGrav.setWork(-1, NTHREAD, 1);
//   approxGrav.setWork(-1, NTHREAD, tree.n_active_groups);


//   int nBlocks = tree.n_active_groups / 2;

  approxGrav.setWork(-1, NTHREAD, tree.n_active_groups);    
  approxGrav.execute(execStream->s());  //First half

/*
 for(int i=0; i <= 100; i++)
{
  devContext.startTiming();

  int activeTemp = tree.n_active_groups * (0.01*i);
  cerr << "ACTIVE: " << activeTemp << endl;
  
  if(activeTemp == 0) activeTemp = 1;

//   approxGrav.setWork(-1, NTHREAD, 120);
 // approxGrav.setWork(-1, NTHREAD, tree.n_active_groups);
  approxGrav.setWork(-1, NTHREAD, activeTemp);
    
  printf("Approx config: "); approxGrav.printWorkSize();
  my_dev::base_mem::printMemUsage();
  fflush(stderr);
//   fprintf(stderr, "Before thread launch! %g \n", get_time() - tTest);
//   approxGrav.execute();

//   devContext.startTiming();
  approxGrav.execute(execStream->s());  //First half
  
    devContext.stopTiming("Approximation", activeTemp);
}
*/
 
//   tree.activePartlist.d2h();
//   fprintf(stderr ,"TEST 1e HELFT: %d \n", tree.activePartlist[tree.n]);  
  
//   grpOffset     = nBlocks;
//   nBlocks       = tree.n_active_groups-nBlocks;
//   approxGrav.set_arg<int>(21,  &grpOffset); //Check why/if we need this!  
//   approxGrav.setWork(-1, NTHREAD, nBlocks);  
//   approxGrav.execute(execStream->s());  //Second half
//   
//   tree.activePartlist.d2h();
//   fprintf(stderr ,"TEST 2e HELFT: %d \n", tree.activePartlist[1024*1024]);
//   
//   
//   
//   grpOffset = 0;
//   approxGrav.set_arg<int>(21,  &grpOffset); //Check why/if we need this!  
//   approxGrav.setWork(-1, NTHREAD, tree.n_active_groups);
//   approxGrav.execute(execStream->s());  //Full
// //   approxGrav.execute(execStream->s());  
// //   devContext.stopTiming("Approximation_detail", 6);
// 
//   tree.activePartlist.d2h();
//   fprintf(stderr ,"TEST FULL: %d \n", tree.activePartlist[1024*1024]);

  
//   fprintf(stderr, "After thread launch! %g\n", get_time() - tTest);

 /*if(test == 0)
 {
  tree.bodies_acc1.d2h();
  tree.bodies_pos.d2h();

  tree.bodies_ids.d2h();
   for(int i=0; i < tree.n; i++)
      {
  	cerr << std::setprecision(7) << tree.bodies_ids[i] << "\tacc x: " << tree.bodies_acc1[i].x << "\tacc y: " << tree.bodies_acc1[i].y <<"\tacc z: " << tree.bodies_acc1[i].z <<"\t" << tree.bodies_acc1[i].w << "\tpos: " << tree.bodies_pos[i].x << " " <<  tree.bodies_pos[i].y << " " << tree.bodies_pos[i].z << " " << tree.bodies_ids[i] << endl;
      }
        //if(iter>0)
      for(int i=tree.n-10; i < tree.n; i++)
      {
      // cerr << i << "\tacc x: " << tree.bodies_acc1[i].x << "\tacc y: " << tree.bodies_acc1[i].y <<"\tacc z: " << tree.bodies_acc1[i].z << endl;
      }
 */
 
  //Print interaction statistics
  #if 0
  
  tree.body2group_list.d2h();
  tree.interactions.d2h();
    long long directSum = 0;
    long long apprSum = 0;
    
    int maxDir = -1;
    int maxAppr = -1;

    for(int i=0; i < tree.n; i++)
    {
      apprSum     += tree.interactions[i].x;
      directSum   += tree.interactions[i].y;
      
      maxAppr = max(maxAppr,tree.interactions[i].x);
      maxDir  = max(maxDir,tree.interactions[i].y);
      
      
//       if(tree.interactions[i].y > 5000)        cerr << i << "( " << tree.body2group_list[i] << " )\t dir: " << tree.interactions[i].y << " avg: " << tree.interactions[i].x << endl;
      
//       if(i == 0) cout << "id :\t " << i << "\tdirect: " << tree.interactions[i].y <<"\tapprox: " << tree.interactions[i].x  << endl;
      
    }
  
    //cerr << "Interaction at iter: " << iter << "\tdirect: " << directSum << "\tappr: " << apprSum << "\t";
    //cerr << "avg dir: " << directSum / tree.n << "\tavg appr: " << apprSum / tree.n << endl;

    cout << "Interaction at (rank= " << mpiGetRank() << " ) iter: " << iter << "\tdirect: " << directSum << "\tappr: " << apprSum << "\t";
    cout << "avg dir: " << directSum / tree.n << "\tavg appr: " << apprSum / tree.n << "\tMaxdir: " << maxDir << "\tmaxAppr: " << maxAppr <<  endl;
    
 //   cerr << "avg dir: " << directSum / tree.n << "\tavg appr: " << apprSum / tree.n << "\tMaxdir: " << maxDir << "\tmaxAppr: " << maxAppr <<  endl;
  
    #if 0
      //Histogram of number of interactions
      const int bins = 256;
      const int jump = 15;
      int histoIDX[bins+1];
      for(int i=0; i < bins; i++)
        histoIDX[i] = 0;
      
/*      for(int i=0; i < tree.n; i++)
      {
          int idx = tree.interactions[i].y / jump;
          if(idx >= bins)
            idx = bins;
          histoIDX[idx]++;  
      }*/

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
      
    
    
      exit(0);
    
    
    #endif
    
    
    
    
    
    
    #endif
//  }
//  exit(0);


  if(mpiGetNProcs() == 1) //Only do it here if there is only one process
  {
    
    //Reduce the number of valid particles    
    getNActive.set_arg<int>(0,    &tree.n);
    getNActive.set_arg<cl_mem>(1, tree.activePartlist.p());
    getNActive.set_arg<cl_mem>(2, this->nactive.p());
    getNActive.set_arg<int>(3, NULL, 128); //Dynamic shared memory , equal to number of threads
    getNActive.setWork(-1, 128, NBLOCK_REDUCE);
    
    CU_SAFE_CALL(cuCtxSynchronize()); //Syncrhonize all streams, makes sure that the approx stream is finished
    getNActive.execute();
    
    //Reduce the last parts on the host
    this->nactive.d2h();
    tree.n_active_particles = this->nactive[0];
    for (int i = 1; i < NBLOCK_REDUCE ; i++)
        tree.n_active_particles += this->nactive[i];

    printf("Active particles: %d \n",tree.n_active_particles);
  }
}
//end approximate


void octree::approximate_gravity_let(tree_structure &tree, tree_structure &remoteTree)
{

//    if(mpiGetRank() == 0)
//     walkRemoteTreeOnHost(remoteTree, iter);
  
//   MPI_Barrier(MPI_COMM_WORLD);
//    exit(0);  
  
//   my_dev::dev_mem<int2> interactions(devContext,tree.n);



  //Start and end node of the remote tree structure
  uint2 node_begend;  
  node_begend.x =  0;
  node_begend.y =  remoteTree.remoteTreeStruct.w;
  
  //The texture offsets used:
  int particleTexOffset = (remoteTree.remoteTreeStruct.z >> 16);
  int nodeTexOffset     = (remoteTree.remoteTreeStruct.z & 0xFFFF);
 
  //Number of particles and number of nodes in the remote tree
  int remoteP = remoteTree.remoteTreeStruct.x;
  int remoteN = remoteTree.remoteTreeStruct.y;

  printf("LET node begend: %d %d iter-> %d\n", node_begend.x, node_begend.y, iter);
  fflush(stderr);
  fflush(stdout);
  //Reset the active particles
//   tree.activePartlist.zeroMem();

//
  //Set the kernel parameters, many!
  approxGravLET.set_arg<int>(0,    &tree.n_active_groups);
  approxGravLET.set_arg<float4>(1, &tree.corner);
  approxGravLET.set_arg<float>(2,  &(this->inv_theta));
  approxGravLET.set_arg<float>(3,  &(this->eps2));
  approxGravLET.set_arg<uint2>(4,  &node_begend);
  approxGravLET.set_arg<cl_mem>(5, tree.active_group_list.p());
  approxGravLET.set_arg<cl_mem>(6, remoteTree.fullRemoteTest.p());

  void *multiLoc = remoteTree.fullRemoteTest.a(2*(particleTexOffset+remoteP) + 2*(remoteN+nodeTexOffset));
  approxGravLET.set_arg<cl_mem>(7, &multiLoc);  

  approxGravLET.set_arg<cl_mem>(8, tree.bodies_acc1.p());
  approxGravLET.set_arg<cl_mem>(9, tree.ngb.p());
  approxGravLET.set_arg<cl_mem>(10, tree.activePartlist.p());
  approxGravLET.set_arg<cl_mem>(11, tree.interactions.p());
  
  void *boxSILoc = remoteTree.fullRemoteTest.a(2*(particleTexOffset+remoteP));
  approxGravLET.set_arg<cl_mem>(12, &boxSILoc);  

  approxGravLET.set_arg<cl_mem>(13, tree.groupSizeInfo.p());

  void *boxCILoc = remoteTree.fullRemoteTest.a(2*(particleTexOffset+remoteP) + remoteN + nodeTexOffset);
  approxGravLET.set_arg<cl_mem>(14, &boxCILoc);  

  approxGravLET.set_arg<cl_mem>(15, tree.groupCenterInfo.p());  
  
  void *bdyVelLoc = remoteTree.fullRemoteTest.a(1*(particleTexOffset+remoteP));
  approxGravLET.set_arg<cl_mem>(16, &bdyVelLoc);  //<- Remote bodies velocity
  
  approxGravLET.set_arg<cl_mem>(17, tree.group_list.p()); //Check why/if we need this!  
  approxGravLET.set_arg<cl_mem>(18, tree.bodies_Ppos.p()); //<- Predicted local body positions
  approxGravLET.set_arg<cl_mem>(19, tree.bodies_Pvel.p()); //<- Predicted local body velocity
  
  approxGravLET.set_arg<real4>(20, remoteTree.fullRemoteTest, 4, "texNodeSize",
                               2*(remoteP + particleTexOffset), remoteN );
  approxGravLET.set_arg<real4>(21, remoteTree.fullRemoteTest, 4, "texNodeCenter",
                               2*(particleTexOffset + remoteP) + (remoteN + nodeTexOffset),
                               remoteN);
  approxGravLET.set_arg<real4>(22, remoteTree.fullRemoteTest, 4, "texMultipole",
                               2*(particleTexOffset + remoteP) + 2*(remoteN + nodeTexOffset), 
                               3*remoteN);
  approxGravLET.set_arg<real4>(23, remoteTree.fullRemoteTest, 4, "texBody", 0, remoteP);  
    
 //   approxGrav.setWork(-1, 64, tree.n_active_groups);
 //  approxGrav.setWork(-1, NCRIT, tree.n_active_groups);
//   approxGrav_let.setWork(-1, NCRIT, tree.n_active_groups);
  approxGravLET.setWork(-1, NTHREAD, tree.n_active_groups);
 
  printf("LET Approx config: "); approxGravLET.printWorkSize();

  my_dev::base_mem::printMemUsage();
 
  devContext.startTiming();
  approxGravLET.execute(execStream->s());
  devContext.stopTiming("Approximation_let", 5);   
  
  
  #if 0
    tree.bodies_acc1.d2h();
//     //for(int i=0; i < tree.n ; i++)
//     for(int i=0; i < 10 ; i++)
//     {
//         printf("TEST acc.w: %f \n", tree.bodies_acc1[i].w);      
//     }
//   
    double test = 0;
    for(int i=0; i < tree.n ; i++)    
    {
        test += tree.bodies_acc1[i].w;
    }
    
    printf("TEST2 %d \t %f \n", tree.n, test);
  
  
  
  exit(0);
  #endif
  
    #if 0
    tree.bodies_acc1.d2h();
    tree.bodies_ids.d2h();
    for(int i=0; i < tree.n ; i++)
//     for(int i=0; i < 10 ; i++)
    {
      if(tree.bodies_ids[i] < 10)
        printf("TEST (%d) %d %d acc: %f %f %f %f \n", mpiGetRank(), i, tree.bodies_ids[i],
               tree.bodies_acc1[i].x, tree.bodies_acc1[i].y, tree.bodies_acc1[i].z, tree.bodies_acc1[i].w);      
    }
//   


  
  
  
  exit(0);
  #endif

// if(mpiGetRank() == 0)
// for(int i=0; i < 15; i++){
// printf("MP1: %d %f %f %f %f \n" , i, remoteTree.multipole[i*3+0].x, remoteTree.multipole[i*3+0].y, remoteTree.multipole[i*3+0].z, remoteTree.multipole[i*3+0].w);
// printf("MP2: %f %f %f %f \n" , remoteTree.multipole[i*3+1].x, remoteTree.multipole[i*3+1].y, remoteTree.multipole[i*3+1].z, remoteTree.multipole[i*3+1].w); 
// printf("MP3: %f %f %f %f \n" , remoteTree.multipole[i*3+2].x, remoteTree.multipole[i*3+2].y, remoteTree.multipole[i*3+2].z, remoteTree.multipole[i*3+2].w); 
// }
//   
//  if(test == 0)
//  {
//   tree.bodies_acc1.d2h();
//   tree.bodies_pos.d2h();
// 
//   tree.bodies_ids.d2h();
//    for(int i=0; i < tree.n; i++)
//       {
// //      cerr << std::setprecision(7) << tree.bodies_ids[i] << "\tacc x: " << tree.bodies_acc1[i].x << "\tacc y: " << tree.bodies_acc1[i].y <<"\tacc z: " << tree.bodies_acc1[i].z <<"\t" << tree.bodies_acc1[i].w << "\tpos: " << tree.bodies_pos[i].x << " " <<  tree.bodies_pos[i].y << " " << tree.bodies_pos[i].z << " " << tree.bodies_ids[i] << endl;
//       }
//         //if(iter>0)
//       for(int i=tree.n-10; i < tree.n; i++)
//       {
//       // cerr << i << "\tacc x: " << tree.bodies_acc1[i].x << "\tacc y: " << tree.bodies_acc1[i].y <<"\tacc z: " << tree.bodies_acc1[i].z << endl;
//       }
// //exit(0);
 
  //Print interaction statistics
 //Print interaction statistics
  #if 1
    tree.interactions.d2h();
    long long directSum = 0;
    long long apprSum = 0;
    
    int maxDir = -1;
    int maxAppr = -1;


    for(int i=0; i < tree.n; i++)
    {
      apprSum     += tree.interactions[i].x;
      directSum   += tree.interactions[i].y;
      
      maxAppr = max(maxAppr,tree.interactions[i].x);
      maxDir  = max(maxDir, tree.interactions[i].y);
            
//       if(i == 0)    cout << "id :\t " << i << "\tdirect: " << interactions[i].y <<"\tapprox: " << interactions[i].x  << endl;
    }
  
    //cerr << "Interaction at iter: " << iter << "\tdirect: " << directSum << "\tappr: " << apprSum << "\t";
    //cerr << "avg dir: " << directSum / tree.n << "\tavg appr: " << apprSum / tree.n << endl;

    cout << "Interaction (LET) at (rank= " << mpiGetRank() << " ) iter: " << iter << "\tdirect: " << directSum << "\tappr: " << apprSum << "\t";
    cout << "avg dir: " << directSum / tree.n << "\tavg appr: " << apprSum / tree.n  << "\tMaxdir: " << maxDir << "\tmaxAppr: " << maxAppr <<  endl;
  #endif


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
//end approximate



void octree::correct(tree_structure &tree)
{
  correctParticles.set_arg<int   >(0, &tree.n);
  correctParticles.set_arg<float >(1, &t_current);
  correctParticles.set_arg<cl_mem>(2, tree.bodies_time.p());
  correctParticles.set_arg<cl_mem>(3, tree.activePartlist.p());
  correctParticles.set_arg<cl_mem>(4, tree.bodies_vel.p());
  correctParticles.set_arg<cl_mem>(5, tree.bodies_acc0.p());
  correctParticles.set_arg<cl_mem>(6, tree.bodies_acc1.p());
  correctParticles.set_arg<cl_mem>(7, tree.bodies_pos.p());
  correctParticles.set_arg<cl_mem>(8, tree.bodies_Ppos.p());
  correctParticles.set_arg<cl_mem>(9, tree.bodies_Pvel.p());
     

  correctParticles.setWork(tree.n, 128);
  correctParticles.execute();
  clFinish(devContext.get_command_queue());


  computeDt.set_arg<int>(0,    &tree.n);
  computeDt.set_arg<float>(1,  &t_current);
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
//   tree.bodies_ids.d2h();
//   tree.body2group_list.d2h();
  for (int i = 0; i < tree.n; i++) {
    float4 vel = tree.bodies_vel[i];
    hEkin += tree.bodies_pos[i].w*0.5*(vel.x*vel.x +
                               vel.y*vel.y +
                               vel.z*vel.z);
    hEpot += tree.bodies_pos[i].w*0.5*tree.bodies_acc0[i].w;

//      if(i > 350 && mpiGetRank() == 0){
// //     fprintf(stderr, "%d\t%f\t%f\n", i, tree.bodies_pos[i].w, tree.bodies_acc0[i].w);
// //    fprintf(stderr, "%d\t%f\t%f\n", i, hEkin, hEpot);
//   fprintf(stderr, "%d\t%f\t%f\t%f\t%f\t%f\t%f\n", i, tree.bodies_acc0[i].x, tree.bodies_acc0[i].y, tree.bodies_acc0[i].z, tree.bodies_acc0[i].w, hEkin, hEpot);
//    fprintf(stderr, "%d\t(id: %d  grp: %d )\t%f\t%f\t%f\t%f\n", i, tree.bodies_ids[i], tree.body2group_list[i], vel.x, vel.y, vel.z, tree.bodies_pos[i].w);
//    if(tree.bodies_ids[i] == 222) fprintf(stderr, "\nTEST TEST TEST TEST TEST TEST \n");
//    
//    if(tree.body2group_list[i] > 14) exit(0);
//    
//      }

  }
  MPI_Barrier(MPI_COMM_WORLD);
  double hEtot = hEpot + hEkin;
  printf("Energy (on host): Etot = %.10lg Ekin = %.10lg Epot = %.10lg \n", hEtot, hEkin, hEpot);
  #endif

  //float2 energy : x is kinetic energy, y is potential energy
  int blockSize = NBLOCK_REDUCE ;
//   my_dev::dev_mem<double2>  energy(devContext, blockSize);
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
  
//   mpiSync();
//   if(this->t_current > 0.01) exit(0);

  
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
