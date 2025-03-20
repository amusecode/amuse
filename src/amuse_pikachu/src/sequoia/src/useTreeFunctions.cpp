
#include "octree.h"

/*
void octree::approximate_gravity(tree_structure &tree,
                                 my_dev::dev_mem<float4>  &j_bodies_pos,  //Bodies that are part of the tree-structure
                                 my_dev::dev_mem<float4>  &i_bodies_pos,  //Bodies that are part of the groups
                                 int n_groupBodies,                     //Number of bodies that are part of the groups
                                 my_dev::dev_mem<real4>  &i_bodies_acc,
                                 my_dev::dev_mem<real>   &i_bodies_ds2,
                                 my_dev::dev_mem<int>    &i_bodies_ngb
){
*/
// by M.I.
void octree::approximate_gravity(tree_structure &tree,
                                 my_dev::dev_mem<float4>  &j_bodies_pos,  //Bodies that are part of the tree-structure
                                 my_dev::dev_mem<float4>  &i_bodies_pos,  //Bodies that are part of the groups
                                 int n_groupBodies,                     //Number of bodies that are part of the groups
                                 my_dev::dev_mem<real4>  &i_bodies_acc,
                                 my_dev::dev_mem<real>   &i_bodies_ds2,
                                 my_dev::dev_mem<int>    &i_bodies_ngb,
                                 my_dev::dev_mem<int>    &i_bodies_Nngb){ 
  uint2 node_begend;
  int level_start = 2;
  node_begend.x   = tree.level_list[level_start].x;
  node_begend.y   = tree.level_list[level_start].y;

  //Reset the active particles
  tree.activePartList.zeroMem();

  float eps2 = this->eps2;
  float rsearch_sq = this->rsearch_sq; // modifyied by M.I.

  //Set the kernel parameters, many!
  int argIdx = 0;
  approxGrav.set_arg<int>(argIdx++,    &tree.n_groups);
  approxGrav.set_arg<int>(argIdx++,    &n_groupBodies);
  approxGrav.set_arg<float>(argIdx++,  &eps2);
  approxGrav.set_arg<uint2>(argIdx++,  &node_begend);
  approxGrav.set_arg<cl_mem>(argIdx++, j_bodies_pos.p());
  approxGrav.set_arg<cl_mem>(argIdx++, i_bodies_acc.p());
  approxGrav.set_arg<cl_mem>(argIdx++, i_bodies_pos.p());
  approxGrav.set_arg<cl_mem>(argIdx++, i_bodies_ds2.p());
  approxGrav.set_arg<cl_mem>(argIdx++, i_bodies_ngb.p());
  approxGrav.set_arg<cl_mem>(argIdx++, tree.activePartList.p());
  approxGrav.set_arg<cl_mem>(argIdx++, tree.interactions.p());
  approxGrav.set_arg<cl_mem>(argIdx++, tree.group_list.p());
  approxGrav.set_arg<cl_mem>(argIdx++, tree.multipole.p());
  approxGrav.set_arg<cl_mem>(argIdx++, tree.boxSizeInfo.p());
  approxGrav.set_arg<cl_mem>(argIdx++, tree.boxCenterInfo.p());
  approxGrav.set_arg<cl_mem>(argIdx++, tree.generalBuffer1.p()); //Instead of using Local memory
  approxGrav.set_arg<float>(argIdx++,  &rsearch_sq); // modifyied by M.I.
  approxGrav.set_arg<cl_mem>(argIdx++,  i_bodies_Nngb.p()); // modifyied by M.I.
  approxGrav.set_arg<real4>(argIdx++, tree.boxSizeInfo, 4, "texNodeSize");
  approxGrav.set_arg<real4>(argIdx++, tree.boxCenterInfo, 4, "texNodeCenter");
  approxGrav.set_arg<real4>(argIdx++, tree.multipole, 4, "texMultipole");
  approxGrav.set_arg<real4>(argIdx++, j_bodies_pos, 4, "texBody");
  approxGrav.setWork(-1, NTHREAD, nBlocksForTreeWalk);
//   approxGrav.setWork(-1, NTHREAD, 1);
  approxGrav.execute(execStream->s());  //First half

  execStream->sync();
  
  //Print interaction statistics
  #if 0
  tree.n = n_groupBodies;
  
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

    cout << "Interaction at (rank= " << 0 << " ) iter: " << 0 << "\tdirect: " << directSum << "\tappr: " << apprSum << "\t";
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
    i_bodies_ngb.d2h();
    i_bodies_acc.d2h();
    i_bodies_ds2.d2h();
    
    for(int i=0; i < n_groupBodies; i++)
    {
        fprintf(stderr, "%d\t Acc: %f %f %f %f  \t Ds2: %f \t NGB: %d \n", 
                i, i_bodies_acc[i].x, i_bodies_acc[i].y,
                i_bodies_acc[i].z, i_bodies_acc[i].w,
                i_bodies_ds2[i],
                i_bodies_ngb[i]);
    }
  #endif
  
  
  
  
  
/*
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
    */

  my_dev::base_mem::printMemUsage();
}
//end approximate


