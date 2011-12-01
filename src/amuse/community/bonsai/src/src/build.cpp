#include "octree.h"
#include "peano.h"

/* Peano Hilbert Sorting! */

struct cmp_peanokey_index{
        bool operator() (const peano_hilbert::peano_struct &lhs, const peano_hilbert::peano_struct &rhs){
                return lhs.key < rhs.key;
        }
};

void sort_local_data_ph(real4* bodies_pos,
                        int local_n, 
                        real4 corner_info,
                        uint *newOrder
) {
  if (local_n == 0) return;
  
  static std::vector<int> order;
/*  static std::vector<int> order;
  order.resize(local_n);
  order.clear()*/;
  {
      static std::vector<peano_hilbert::peano_struct> keys;
      keys.resize(local_n);
      keys.clear();
      keys.resize(local_n);
//       tree.domain_fac   = size/(1 << MAXLEVELS);
//       float idomain_fac = 1.0/tree.domain_fac;
//       float domain_fac = tree.domain_fac;
//       tree.corner.w = domain_fac;        
      
      //Get the size back from the domain_fac
      const float size = corner_info.w* (1 << MAXLEVELS);
      
//       const vec3 domain_len(global_domain.hsize() * 2.0);
//       const float size = std::max(domain_len.x,
//                       std::max(domain_len.y, domain_len.z));

      const int domain_fac = 1.0f / size * (((peano_hilbert::peanokey)1) << (BITS_PER_DIMENSION));

//       const double xmin = global_domain.rmin.x;
//       const double ymin = global_domain.rmin.y;
//       const double zmin = global_domain.rmin.z;
      const double xmin = corner_info.x;
      const double ymin = corner_info.y;
      const double zmin = corner_info.z;

      for (int i = 0; i < local_n; i++) {
              const int x = (int)((bodies_pos[i].x - xmin) * domain_fac);
              const int y = (int)((bodies_pos[i].y - ymin) * domain_fac);
              const int z = (int)((bodies_pos[i].z - zmin) * domain_fac);
              keys[i].idx = i;
              keys[i].key = peano_hilbert::peano_hilbert_key(x, y, z, BITS_PER_DIMENSION);
      }

      cerr << "keys size: " << keys.size() << endl;
      
      std::sort(keys.begin(), keys.end(), cmp_peanokey_index());

      for (int i = 0; i < local_n; i++)
              newOrder[i] = keys[i].idx;
      
      
      
      //Now with the new key method

     order.resize(local_n);

      for (int i = 0; i < local_n; i++) {
              const int x = (int)((bodies_pos[i].x - xmin) * domain_fac);
              const int y = (int)((bodies_pos[i].y - ymin) * domain_fac);
              const int z = (int)((bodies_pos[i].z - zmin) * domain_fac);
              keys[i].idx = i;
              keys[i].key = peano_hilbert::peano_hilbert_key_new_thesis2(x, y, z, BITS_PER_DIMENSION);
              
              if(i < 10)
                fprintf(stderr, "%d\t%d\t%d \n", x, y, z);
      }
      

// for(int i=0; i < 10; i++)
// {
//     fprintf(stderr, "%d  \t (host)key: %016llX \t %lld\n", i, keys[i].key, keys[i].key);
// }


      
      std::sort(keys.begin(), keys.end(), cmp_peanokey_index());

// for(int i=0; i < 10; i++)
// {
//     fprintf(stderr, "%d  \t (host after sort)key: %016llX \t --> %d \n", i, keys[i].key, keys[i].idx);
// }

      
      for (int i = 0; i < local_n; i++)
              newOrder[i] = keys[i].idx;

	     
//      2       4       Approximation   461.501
// 7       4       Approximation   455.728
/*Interaction at iter: 0  direct: 886232903       appr: 2713473664        avg dir: 845    avg appr: 2587
avg dir: 845    avg appr: 2587  Maxdir: 26539   maxAppr: 9664*/
/*Interaction at iter: 0  direct: 746757852       appr: 2626137472        avg dir: 712    avg appr: 2504
avg dir: 712    avg appr: 2504  Maxdir: 1488    maxAppr: 3200
  */     
      
  }
  
//   for (int i = 0; i < local_n; i++) 
//   {
//     if(newOrder[i] != i)
//       fprintf(stderr, "Na reorder: %d \t-> %d \t %d \n", i, newOrder[i],  i);
//   }
  
  //                 reorder(order, ptcl);
}



//     Hier gebleven
//     
//     het moet anders, misschien gewoon in de basis een node 
//     alvast opbouwen
//     of op level X beginnen en dan recursief checken, maar dan kloppen de offsets niet 
//     meer. 
//     Of tijdens tree-build op een gegeven moment overgaan op een dense tree structure
//     
//     als nchilds < 512 -> maak 512/NLEAF nodes aan, dat is wel een idee kijken of dat makkelijk
//     kan tijdens het defineren van de nodes, want dan doen we al iets met het instellen van de
//     start eind particle. Misschien dat we dat daar al kunnen berekenen?

// uint setNewChildren(uint curValue, int removed)
// {
//   //Node
//   int child    =   curValue & 0x0FFFFFFF;         //TODO make this name/define?
//   int nchild   = ((curValue & 0xF0000000) >> 28); //TODO make this name/define?
//   
//   child -= removed;
//   uint newValue = (curValue & 0xF0000000) | child;
//   
//   return newValue;
// }

void create_dot(uint* n_children, uint2* node_bodies, uint* validList)
{
     //Walk the tree as is done on device, level by level
    vector<int> curLevel;
    vector<int> nextLevel;

    //Add the top node
    curLevel.push_back(0);
    //Start the tree-walk       
    while(curLevel.size() > 0)
    {      
      for(unsigned int i=0; i < curLevel.size(); i++)
      {
        int node = curLevel[i];
        //Read node data        
        bool leaf =  (validList[node] >> 31);
        
        if(leaf)
        {
          //Dont have to do anything          
          //Count the number of particles
          //Leaf
          uint2 bij = node_bodies[node];          
          int child     = bij.x & ILEVELMASK; 
          int nchild    = bij.y - (bij.x & ILEVELMASK);    
                            

          char buff[256];
          sprintf(buff, "%d ( L %d -> %d )", node, child, nchild);
          fprintf(stderr, "%d [label=\"%s\"]; \n", node, buff);            
          
          continue;
        }
        
         //Node
        int child    =   n_children[node] & 0x0FFFFFFF;         //TODO make this name/define?
        int nchild   = ((n_children[node] & 0xF0000000) >> 28); //TODO make this name/define?
    

        //Process the children of the non-leaf node  and check if we can compress the child nodes
        for(int j=child; j < child+nchild; j++)
        {
          fprintf(stderr, "%d -> %d ; \n", node, j);     
          
          nextLevel.push_back(j);
        }
        
        uint2 bij = node_bodies[node];          
        int allchild    = bij.y - (bij.x & ILEVELMASK);           
        
        char buff[256];
        sprintf(buff, "%d ( %d )", node, allchild);
        fprintf(stderr, "%d [label=\"%s\"]; \n", node, buff);  
                      
        
      } //end for curLevel.size
      
      //Put next level stack into current level and continue
      curLevel.clear();

      curLevel.assign(nextLevel.begin(), nextLevel.end());
      nextLevel.clear();
    }//end while
  
  
}

void create_local_essential_tree_count_modify(uint* n_children, uint2* node_bodies, uint* validList,
                                              int &particles, int &nodes, uint *removeOffsets, int &nRemovedNodes)
{
    //Walk the tree as is done on device, level by level
    vector<int> curLevel;
    vector<int> nextLevel;
    
//     vecotr<int> processed(nodes);
//     invalidList.assign(0,0, nodes.size());
   
    int particleCount   = 0;
    int nodeCount       = 0;    
    int removed         = 0;
    
    //Add the top node
    curLevel.push_back(0);
    nodeCount++;
    
    //Start the tree-walk       
    while(curLevel.size() > 0)
    {      
      for(unsigned int i=0; i < curLevel.size(); i++)
      {
        int node = curLevel[i];
        
        removeOffsets[node] = removed;
        //Read node data        
        bool leaf =  (validList[node] >> 31);
        
        if(leaf)
        {
          //Dont have to do anything          
          //Count the number of particles
          //Leaf
          uint2 bij = node_bodies[node];          
//           int child     = bij.x & ILEVELMASK; 
          int nchild    = bij.y - (bij.x & ILEVELMASK);    
          
//           fprintf(stderr,"Leaf: %d  %d %d \n", node, child, nchild);
        
          particleCount += nchild;
          
          continue;
        }
        
          //Node
        int child    =   n_children[node] & 0x0FFFFFFF;         //TODO make this name/define?
        int nchild   = ((n_children[node] & 0xF0000000) >> 28); //TODO make this name/define?
    
        int childSum      = 0;          
        int newChildStart = -1;
//         int nCount        = 0;
        int newChildIdx   = -1;

        //Process the children of the non-leaf node  and check if we can compress the child nodes
        for(int j=child; j < child+nchild; j++)
        {
          //Get info about the child
          bool childIsLeaf = (validList[j] >> 31);
          
          int childlevel = 0;

          if(childIsLeaf)
          {
            //Leaf
            uint2 grandchildinfo = node_bodies[j];          
            int   grandchild     = grandchildinfo.x & ILEVELMASK;                          //the first body in the leaf
            int   grandnchild    = grandchildinfo.y - (grandchildinfo.x & ILEVELMASK);     //number of bodies in the leaf 
            
            childlevel = (grandchildinfo.x &  LEVELMASK) >> BITLEVELS;
            
            if((childSum +grandnchild) <= NLEAFTEST)
            {
              if(childSum == 0)
              {
                //Fist modification
                newChildIdx   = j;
                newChildStart = grandchild;
              }
              
              childSum        += grandnchild;
              node_bodies[j]  = (uint2){0,0};  
              removed++;
            }//end if childSUm+grandnchild <= NLEAF
            else if(childSum > 0)
            {
              //Sum is too large but a node is being modified, finish that and make a new one 
              node_bodies[newChildIdx] = (uint2) {newChildStart | (childlevel << BITLEVELS), newChildStart+childSum};
              removed--;
              
              childSum      = grandnchild;
              newChildStart = grandchild; 
              newChildIdx   = j;
              node_bodies[j]  = (uint2){0,0};  
              removed++;
            }
          }
          else
          {
            //This is a non-leaf node
            if(childSum > 0)
            {
              //Sum is too large but a node is being modified, finish that and 
              node_bodies[newChildIdx] = (uint2) {newChildStart | (childlevel << BITLEVELS), newChildStart+childSum};
              removed--;
              
              childSum = 0;              
              newChildStart = -1; newChildIdx = -1;
            }
          }
          
          if((j+1) == nchild+child)
          {
            //This was the last child finish the partial node
            if(childSum > 0)
            { 
              //Sum is too large but a node is being modified, finish that and 
              node_bodies[newChildIdx] = (uint2) {newChildStart | (childlevel << BITLEVELS), newChildStart+childSum};
              removed--;
              childSum = 0;
              newChildStart = -1; newChildIdx = -1;
            }
          }//end if last child
          
          nextLevel.push_back(j);
          
        }

        //Increase the nodeCount, since this node will be part of the tree-structure
        nodeCount++;        
      } //end for curLevel.size
      
      //Put next level stack into current level and continue
      curLevel.clear();
      
      fprintf(stderr ,"Removed on curLevel: %d \n", removed);

      curLevel.assign(nextLevel.begin(), nextLevel.end());
      nextLevel.clear();
    }//end while
    
    particles    = particleCount;
    nodes        = nodeCount;  
    
    fprintf(stderr, "COUNT found: %d particles and %d nodes  can remove: %d \n", particles, nodes, removed);   
    nRemovedNodes = removed;
}


void compress_tree(uint* n_children, uint2* node_bodies, uint* validList, uint* removeOffsets, int n_nodes)
{
  //Walk over all the nodes and remove nodes with 0 children
  
  int storeIdx = 0;
  int removed  = 0;
  
  for(int i=0; i < n_nodes; i++)
  {
    uint child    =   n_children[i] & 0x0FFFFFFF;         //TODO make this name/define?
    uint nchild   = ((n_children[i] & 0xF0000000) >> 28); //TODO make this name/define?          
    
    //Only modify non-leaf nodes
    bool leaf =  (validList[i] >> 31);  
    
    if(leaf)
    {
      uint2 node_body = node_bodies[i];
      
      if(node_body.y == 0)
      {
        //Remove this node
        removed++;       
      }      
      else
      {
        //Keep this node
        uint temp  = (validList[i] >> 31);  
        uint valid = storeIdx| (uint)(temp << 31);   //Distinguish leaves and nodes
        
        n_children[storeIdx]     = n_children[i];        
        node_bodies[storeIdx]    = node_bodies[i];        
//         validList[storeIdx]      = validList[i];
        validList[storeIdx]      = valid;
//         fprintf(stderr, "Leaf: %d goes to: %d || %d %d %d %d \n", i, storeIdx,  n_children[i], node_bodies[i].x, node_bodies[i].y, (validList[i] >> 31) );          
        
        storeIdx++;
            
      }
    }
    else
    {    
      //Walk over the children to count the number that have to be removed      
      int newChildCount = 0;
      for(unsigned int j=child; j < child+nchild; j++)
      {
        bool childleaf =  (validList[j] >> 31); 
        
//         fprintf(stderr, "Child node: %d details: %d %d %d \n", j), 
        
        if(childleaf)
        {
            uint2 child_body = node_bodies[j];      
            if(child_body.y != 0)
              newChildCount++;
        }
        else
        {
          newChildCount++;
        }        
      }
      
//       fprintf(stderr, "Node %d had: (%d, %d)  but becomes: (%d, %d) \n", i, child, nchild, child-removeOffsets[i], newChildCount);
      
      child -= removeOffsets[i]; 
      
      uint temp  = (validList[i] >> 31);  
      uint valid = storeIdx| (uint)(temp << 31);   //Distinguish leaves and nodes   
      
      child = child | (newChildCount << 28);
      n_children[storeIdx]     = child;
      node_bodies[storeIdx]    = node_bodies[i];        
//       validList[storeIdx]      = validList[i];   
      validList[storeIdx]      = valid;
      
//        fprintf(stderr, "Node: %d goes to: %d \n", i, storeIdx);      
      
      storeIdx++;
    }    
  }//end for loop
  
  
  fprintf(stderr, "Number of removed nodes: %d %d %d\n", removed, storeIdx, n_nodes);
  
}


/*
void create_local_essential_tree_count_modify(uint* n_children, uint2* node_bodies, uint* validList,
                                              int &particles, int &nodes)
{
    //Walk the tree as is done on device, level by level
    vector<int> curLevel;
    vector<int> nextLevel;
    
    vector<uint>  newTreeChildren;
    vector<uint2> newTreeNodeBodies;
    vector<uint>  newValidList;
    
    vecotr<int> processed(nodes);
    invalidList.assign(0,0, nodes.size());
   
    int particleCount   = 0;
    int nodeCount       = 0;    
    int removed         = 0;
    
    //Add the top node
    curLevel.push_back(0);
    nodeCount++;
    
    //Start the tree-walk       
    while(curLevel.size() > 0)
    {      
      for(unsigned int i=0; i < curLevel.size(); i++)
      {
        int node = curLevel[i];
        //Read node data        
        bool leaf =  (validList[node] >> 31);
        
        //If the node / leaf is processed we do not add it to the final tree-structure
        if(invalidList[node] == 1)
          continue;
        
        int child, nchild;
        if(!leaf)
        {
          //Node
          child    =   n_children[node] & 0x0FFFFFFF;         //TODO make this name/define?
          nchild   = ((n_children[node] & 0xF0000000) >> 28); //TODO make this name/define?
        }
        else
        {
          //Leaf
          uint2 bij = node_bodies[node];          
          child     = bij.x & ILEVELMASK; 
          nchild    = bij.y - (bij.x & ILEVELMASK);    
        }
        
        bool split = true;

        //if split & node add children to next lvl stack
        if(split && !leaf)
        { 
          //Check if we can compress the child nodes
          int childSum      = 0;          
          int newChildStart = -1;
          int nCount        = 0;
          int newChildIdx   = -1;
                             
          for(int j=child; j < child+nchild; j++)
          {
            //Get info about the child
            bool childIsLeaf = (validList[j] >> 31);

            if(childIsLeaf)
            {
              //Leaf
              uint2 grandchildinfo = node_bodies[j];          
              int   grandchild     = grandchildinfo.x & ILEVELMASK;                          //the first body in the leaf
              int   grandnchild    = grandchildinfo.y - (grandchildinfo.x & ILEVELMASK);     //number of bodies in the leaf 
              
              //Check if we can merge this child with the previous child
              if((childSum + grandnchild) <= NLEAF)
              {
                //We can!
                if(newChildStart < 0)
                {
                  newChildStart = grandchild;   //First time we set it otherwise no need                
                  invalidList[j] = 1;
                  newChildIdx = j;
                }
                childSum += grandnchild;              
                
                invalidList[j] = 1;
                nCount++;
              }
              else
              {
                //Sum becomes too large, write the new node if there is any
                if(newChildStart >= 0)
                {
                  //We started a new node, now store it and remove any previous ones
                  fprintf(stderr, "New merged node goes from: %d and has %d children (Nodes: %d sum: %d)\n", newChildStart, childSum, nCount, removed);
                  
                  
                  
                  
                  removed += nCount-1;
                  
                  invalidList[newChildIdx] = 0;
                  
                  
                  newChildIdx   = -1;
                  childSum      = 0;
                  newChildStart = -1;
                  nCount        = 0;
                  
                  if(grandnchild <= NLEAF)
                  {
                    childSum      = grandnchild;
                    newChildStart = grandchild; 
                    
                    invalidList[j] = 1;
                    
                    nCount++;
                  }                  
                } //if newChildStart >= 0               
              }//end else if (childSum + grandnchild) <= NLEAF)             
            } //end childIsleaf 
            else
            {
              //Check if we merged any previous nodes if so finish them otherwise
              if(newChildStart >= 0)
              {
                //We started a new node, now store it and remove any previous ones
                fprintf(stderr, "New (!leaf) merged node goes from: %d and has %d children (Nodes: %d sum: %d)\n", newChildStart, childSum, nCount, removed);
                removed += nCount-1;
                childSum      = 0;
                newChildStart = -1;
                nCount = 0;
              } //if newChildStart >= 0                   
            }     
            
            //Only do a push_back if the child is not a leaf, otherwise there is the 
            //possibility it has been replaced
            if(!childIsLeaf)
              nextLevel.push_back(j);         
            
            
          }//end for children
          
          if(newChildStart >= 0)
          {
            //We started a new node, now store it and remove any previous ones
            fprintf(stderr, "New (done) merged node goes from: %d and has %d children (Nodes: %d sum: %d)\n", newChildStart, childSum, nCount, removed);
            removed += nCount-1;
            
            childSum      = 0;
            newChildStart = -1;
            nCount = 0;
          } //if newChildStart >= 0                
          
          
//           fprintf(stderr, "%d has %d childs and %d grandchilds\n", node, nchild, childSum);
          
        }
        
        //if split & leaf add particles to particle list
        if(split && leaf)
        { 
          for(int j=child; j < child+nchild; j++)
          {
            particleCount++;
          }
        }

        //Increase the nodeCount, since this node will be part of the tree-structure
        nodeCount++;        
      } //end for curLevel.size
      
      //Put next level stack into current level and continue
      curLevel.clear();

      curLevel.assign(nextLevel.begin(), nextLevel.end());
      nextLevel.clear();
    }//end while
    
    particles    = particleCount;
    nodes        = nodeCount;  
    
    fprintf(stderr, "Count found: %d particles and %d nodes  can remove: %d \n", particles, nodes, removed);   
}
*/

void walkCustomTreeOnHost(tree_structure &tree, my_dev::dev_mem<uint>  &validList)
{
  tree.n_children.d2h();
  validList.d2h();  
  tree.node_bodies.d2h();

  vector<uint> removeOffsets(tree.n_nodes);

//   int nchildren = 0;  
//   fprintf(stderr, "digraph before{ \n");
//   create_dot(&tree.n_children[0],  &tree.node_bodies[0],&validList[0]);
//   fprintf(stderr, "}\n");
  
  
//   for(int i=start; i < end; i++)
  {
//     create_local_essential_tree_host(0, &tree.bodies_Ppos[0], &tree.boxCenterInfo[0],
//   create_local_essential_tree_host_childr(0, &tree.bodies_Ppos[0], &tree.boxCenterInfo[0],
//                                    &tree.boxSizeInfo[0],
//                                    &tree.multipole[0], essentialList, 
//                                    iter, mass, nchildren);
  }
/*
  fprintf(stderr,"Walk done (own tree)!!! selected particles: %ld  Mss sum: %f\tIter: %d  \tChildren: %d\n", 
          essentialList.size(), mass, iter, nchildren);*/

  int particles = 0, nodes = 0, nRemovedNodes = 0;
  create_local_essential_tree_count_modify(&tree.n_children[0],  &tree.node_bodies[0],
                                           &validList[0], particles,nodes, &removeOffsets[0], nRemovedNodes);
                                           

  //Compress the tree
  compress_tree(&tree.n_children[0],  &tree.node_bodies[0], &validList[0], &removeOffsets[0],  tree.n_nodes);
  
  tree.n_nodes = tree.n_nodes - nRemovedNodes;
  tree.n_children.h2d();
  tree.node_bodies.h2d();
  validList.h2d();
  
 /*                                          
  fprintf(stderr, "digraph after{ \n");
  create_dot(&tree.n_children[0],  &tree.node_bodies[0],&validList[0]);
  fprintf(stderr, "}\n");
*/
 
return; 

}


void octree::allocateParticleMemory(tree_structure &tree)
{
  //Allocates the memory to hold the particles data
  //and the arrays that have the same size as there are
  //particles. Eg valid arrays used in tree construction
  int n_bodies = tree.n;
 
  if(nProcs > 1)                //10% extra space, only in parallel when
    n_bodies = n_bodies*1.1;    //number of particles can fluctuate
  
  //Particle properties
  tree.bodies_pos.cmalloc(n_bodies+1, true);   //+1 to set end pos, host mapped? TODO not needed right since we use Ppos
  tree.bodies_key.cmalloc(n_bodies+1, false);   //+1 to set end key
  tree.bodies_ids.cmalloc(n_bodies+1, false);   //+1 to set end key
  
  tree.bodies_Ppos.cmalloc(n_bodies+1, true);   //Memory to store predicted positions, host mapped
  tree.bodies_Pvel.cmalloc(n_bodies+1, true);   //Memory to store predicted velocities, host mapped
  
  tree.bodies_vel.cmalloc(n_bodies, false);     
  tree.bodies_acc0.ccalloc(n_bodies, false);    //ccalloc -> init to 0
  tree.bodies_acc1.ccalloc(n_bodies, false);    //ccalloc -> init to 0
  tree.bodies_time.ccalloc(n_bodies, false);    //ccalloc -> init to 0

  tree.oriParticleOrder.cmalloc(n_bodies, false);	//To desort the bodies tree later on
  
  //iteration properties / information
  tree.activePartlist.ccalloc(n_bodies+1, false);
  tree.ngb.ccalloc(n_bodies, false);  
  tree.interactions.cmalloc(n_bodies, false);
  
  tree.body2group_list.cmalloc(n_bodies, false);
  
  tree.level_list.cmalloc(MAXLEVELS);  
  
  //Added 4096 to give some space for memory allignment
  tree.generalBuffer1.cmalloc(3*n_bodies*4 + 4096, false);
  
//   cmalloc_copy(dev_mem &src_buffer, int offset, int n)
//   tree.body2group_list.cmalloc_copy(tree.generalBuffer1, 0, n_bodies);
  
//   tree.interactions.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
//                                     tree.generalBuffer1.get_flags(), 
//                                     tree.generalBuffer1.get_devMem(),
//                                     &tree.generalBuffer1[0], 0, n_bodies);
//                                     
//                                     
                                      
//     void cmalloc_copy(bool pinned, bool flags, CUdeviceptr cudaMem, 
//                       void* ParentHost_ptr, int offset, int n)  
  
  //Interactions needs; uint2   n 
  //ngb          needs:  uint   n
  //activePartList needs: uint  n  (to be used at same moment as interactions and ngb)
  
  //So we need at least (uint4+uint+uint)*n is uint4*n of memory
  //to be used at the same time, and to be used at the same time as 
  //multipole data


  #if 0
  Vanuit sort bodies

  Keys worden berekend voor de sort bodies, moet het dan nog wel
  in de build.cpp ? Roepen we die niet altijd in combinatie aan?
  

  Could combine it with multiple moments
  none of these are used at the same time
  if we sort, we compute new multiple moments
 
  Have to increase the size of multiple moments
  
  TODO in gpu_iterate
    int blockSize = NBLOCK_REDUCE ;
    my_dev::dev_mem<uint>  nactive(devContext, blockSize);  

  int blockSize = NBLOCK_REDUCE ;
  my_dev::dev_mem<double2>  energy(devContext, blockSize);    

  #endif

//   tree.multipole.cresize(3*tree.n, false);
  
//   tree.multipole.cmalloc(3*tree.n, true); //host alloced
  
  
  //Tree properties, tree size is not known at forehand so
  //allocate worst possible outcome  
  n_bodies = n_bodies / 1;
  tree.node_key.cmalloc(n_bodies, false);
  tree.n_children.cmalloc(n_bodies, false);
  tree.node_bodies.cmalloc(n_bodies, false);  
  //General memory buffers
    
  //Set the context for the memory
  this->tnext.setContext(devContext);  
  this->tnext.ccalloc(NBLOCK_REDUCE,false);  
  
  this->nactive.setContext(devContext);  
  this->nactive.ccalloc(NBLOCK_REDUCE,false);    

  this->devMemRMIN.setContext(devContext);
  this->devMemRMIN.cmalloc(NBLOCK_BOUNDARY, false);  

  this->devMemRMAX.setContext(devContext);
  this->devMemRMAX.cmalloc(NBLOCK_BOUNDARY, false);  
  
  this->devMemCounts.setContext(devContext);
  this->devMemCounts.cmalloc(NBLOCK_PREFIX, false);  

  this->devMemCountsx.setContext(devContext);
  this->devMemCountsx.cmalloc(NBLOCK_PREFIX, false);    
  
  if(mpiGetNProcs() > 1)
  {
//    int remoteSize = (n_bodies*0.1) +  (n_bodies*0.1); //TODO some more realistic number
    int remoteSize = (n_bodies*1); //TODO some more realistic number
    this->remoteTree.fullRemoteTest.cmalloc(remoteSize, true);    
  }
  
}


void octree::reallocateParticleMemory(tree_structure &tree)
{
  //Realloc the memory to hold the particles data
  //and the arrays that have the same size as there are
  //particles. Eg valid arrays used in tree construction
  int n_bodies = tree.n;
  
  bool reduce = false;  //Set this to true to limit memory usage by only allocating what
                        //is required. If its false, then memory is not reduced and a larger
                        //buffer is kept
  
  //Particle properties
  tree.bodies_pos.cresize(n_bodies+1, reduce);   //+1 to set boundary condition
  tree.bodies_key.cresize(n_bodies+1, reduce);   //+1 to set boundary condition
  tree.bodies_ids.cresize(n_bodies+1, reduce);   //

  tree.bodies_Ppos.cresize(n_bodies+1, reduce);   //Memory to store predicted positions
  tree.bodies_Pvel.cresize(n_bodies+1, reduce);   //Memory to store predicted velocities  
  
  tree.bodies_vel.cresize (n_bodies, reduce);     
  tree.bodies_acc0.cresize(n_bodies, reduce);    //ccalloc -> init to 0
  tree.bodies_acc1.cresize(n_bodies, reduce);    //ccalloc -> init to 0
  tree.bodies_time.cresize(n_bodies, reduce);    //ccalloc -> init to 0
  
  tree.oriParticleOrder.cresize(n_bodies, reduce);	//To desort the bodies tree later on
  //iteration properties / information
  tree.activePartlist.cresize(n_bodies+1, reduce);
  tree.ngb.cresize(n_bodies, reduce);  
  tree.interactions.cresize(n_bodies, reduce);
  
  tree.body2group_list.cresize(n_bodies, reduce);  

  
  //Tree properties, tree size is not known at forehand so
  //allocate worst possible outcome  
  n_bodies = n_bodies / 3;
  tree.node_key.cresize(n_bodies, reduce);
  tree.n_children.cresize(n_bodies, reduce);
  tree.node_bodies.cresize(n_bodies, reduce);  
}

void octree::allocateTreePropMemory(tree_structure &tree)
{ 
  int n_nodes = tree.n_nodes;

  //Allocate memory
  if(tree.groupCenterInfo.get_size() > 0)
  {
    //Resize, so we dont alloc if we already have mem alloced
    tree.multipole.cresize(3*n_nodes,     false);
    
    tree.boxSizeInfo.cresize(n_nodes,     false);     //host alloced
    tree.groupSizeInfo.cresize(n_nodes,   false);     

    tree.boxCenterInfo.cresize(n_nodes,   false); //host alloced
    tree.groupCenterInfo.cresize(n_nodes, false);    
  }
  else
  {    
    n_nodes = n_nodes * 1.1; 
    tree.multipole.cmalloc(3*n_nodes, true); //host alloced
        
    tree.boxSizeInfo.cmalloc(n_nodes, true);     //host alloced
    tree.groupSizeInfo.cmalloc(n_nodes, false);     

    tree.boxCenterInfo.cmalloc(n_nodes, true); //host alloced
    tree.groupCenterInfo.cmalloc(n_nodes,false);
  }
}



void octree::build (tree_structure &tree) {

  int level      = 0;
  int validCount = 0;
  int offset     = 0;

//    more process0_profLog |  awk -F "," '{ {if ($3 > $4) SUM+=$3; else SUM+=$4; }} END {print SUM}'


  /******** load kernels **********/

  /******** create memory buffers **********/
  
//   my_dev::dev_mem<uint>  validList(devContext, tree.n*2+10);
//   my_dev::dev_mem<uint>  compactList(devContext, tree.n*2+10);
//   my_dev::dev_mem<uint>  validList(devContext, tree.n*2);
//   my_dev::dev_mem<uint>  compactList(devContext, tree.n*2);

  my_dev::dev_mem<uint>  validList(devContext);
  my_dev::dev_mem<uint>  compactList(devContext);
  
  validList.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[0], 0,
                                    tree.n*2, getAllignmentOffset(0));
                                    
  compactList.cmalloc_copy(tree.generalBuffer1.get_pinned(), 
                                    tree.generalBuffer1.get_flags(), 
                                    tree.generalBuffer1.get_devMem(),
                                    &tree.generalBuffer1[tree.n*2], tree.n*2,
                                    tree.n*2, getAllignmentOffset(tree.n*2));
                                                                        

  
  /******** set kernels parameters **********/
  

  build_key_list.set_arg<cl_mem>(0,   tree.bodies_key.p());
//   build_key_list.set_arg<cl_mem>(1,   tree.bodies_pos.p());
  build_key_list.set_arg<cl_mem>(1,   tree.bodies_Ppos.p());
  build_key_list.set_arg<int>(2, &tree.n);
  build_key_list.set_arg<real4>(3,   &tree.corner);
  
  build_key_list.setWork(tree.n, 128);
  printf("build_key_list: "); build_key_list.printWorkSize();
  
  build_valid_list.set_arg<int>(0, &tree.n);
  build_valid_list.set_arg<int>(1, &level);
  build_valid_list.set_arg<cl_mem>(2,  tree.bodies_key.p());
  build_valid_list.set_arg<cl_mem>(3,  validList.p());  
  
  build_valid_list.setWork(tree.n, 128);
  printf("build_valid_list: "); build_valid_list.printWorkSize();
  
  validList.zeroMem(); 

  build_nodes.set_arg<int>(0, &level);
  build_nodes.set_arg<int>(1, &validCount);
  build_nodes.set_arg<int>(2, &offset);
  build_nodes.set_arg<cl_mem>(3,  compactList.p());
  build_nodes.set_arg<cl_mem>(4,  tree.bodies_key.p());
  build_nodes.set_arg<cl_mem>(5,  tree.node_key.p());
  build_nodes.set_arg<cl_mem>(6,  tree.n_children.p());
  build_nodes.set_arg<cl_mem>(7,  tree.node_bodies.p());

  link_tree.set_arg<int>(0, &offset);
  link_tree.set_arg<cl_mem>(1, tree.n_children.p());
  link_tree.set_arg<cl_mem>(2, tree.node_bodies.p());
//   link_tree.set_arg<cl_mem>(3,  tree.bodies_pos.p());
  link_tree.set_arg<cl_mem>(3,  tree.bodies_Ppos.p());
  link_tree.set_arg<real4>(4,   &tree.corner);
  link_tree.set_arg<cl_mem>(5,  tree.level_list.p());
//   link_tree.set_arg<cl_mem>(6,  id_list.p()); 
  link_tree.set_arg<cl_mem>(6,  validList.p()); 
  link_tree.set_arg<cl_mem>(7,  tree.node_key.p());
  link_tree.set_arg<cl_mem>(8, tree.bodies_key.p());
  link_tree.set_arg<int>(9,  &level);



  /********** build  list of keys ********/
  
  build_key_list.execute();  
  
  
//   tree.bodies_key.d2h();
//   for(int i=0; i < 5; i++)
//   {
//     printf("%d\t%d\t%d\t%d\n", i, tree.bodies_key[i].x, 
//            tree.bodies_key[i].y, tree.bodies_key[i].z);
//   }
//   for(int i=tree.n-5; i < tree.n; i++)
//   {
//     printf("%d\t%d\t%d\t%d\n", i, tree.bodies_key[i].x, 
//            tree.bodies_key[i].y, tree.bodies_key[i].z);
//   }
//     

  /******  build the levels *********/
  
//   my_dev::dev_mem<uint>  testValidList(devContext, tree.n);
  
  int nodeSum = 0;
  for (level = 0; level < MAXLEVELS; level++) {
    // mark bodies to be combined into nodes
    build_valid_list.set_arg<int>(1, &level);
    build_valid_list.execute();
      
    //gpuCompact    
    gpuCompact(devContext, validList, compactList, 
               tree.n*2, &validCount);
                 
    nodeSum += validCount / 2;
    printf("ValidCount (%d): %d \tSum: %d Offset: %d\n", mpiGetRank(), validCount, nodeSum, offset);
    
    validCount /= 2;     
                  
    if (validCount == 0) break;                 
      
    // asssemble nodes           
    build_nodes.setWork(validCount, 128);

    build_nodes.set_arg<int>(0, &level);
    build_nodes.set_arg<int>(1, &validCount);
    build_nodes.set_arg<int>(2, &offset);
    
//     build_nodes.set_arg<cl_mem>(8,  testValidList.p()); Special incase we build dense nodes
    
    build_nodes.execute();
                 
    tree.level_list[level] = (uint2){offset, offset + validCount};
    offset += validCount;
    
    
//       exit(0);
    ////////TEST Make special dense tree  nodes!
    
    //Manually insert leaf nodes, that will become the child
    //nodes of the here detected trees, this is possible to 
    //do on the GPU but for testing lets do it serial on the cpu
    
    #if 0 //Modifies the tree on node level
    int tempNodeTest = 0;
    gpuCompact(devContext, testValidList, compactList, validCount, &tempNodeTest);
               
   printf("TempNodeTest: %d \n",tempNodeTest);
    compactList.d2h();
    tree.node_bodies.d2h();
    tree.bodies_key.d2h();
    
    for(int i=0 ; i < tempNodeTest; i++)
   // int i;if(0)
    {
      uint2 bij = tree.node_bodies[compactList[i]];
      int child     = bij.x & ILEVELMASK; 
      int nchild    = bij.y - (bij.x & ILEVELMASK);         
              
      int newNodes = nchild / NLEAF;
      if((nchild % NLEAF) > 0) newNodes++; 

//       printf("TEST: Node: %d  -> %d (%d, %d )  -> Will be %d  \n", compactList[i], nchild, child, bij.y, newNodes);      
              
      for(int j=0; j < newNodes; j++)
      {
        int start = child + j*NLEAF;
        int end   = min(start+NLEAF, child+nchild);
       
        //Set the particles invalid
        for (int i = start; i < end; i++)
        {
          if((i == start) && (i == end-1))
          {
            test_key[i] = (uint2){0xFFFFFFF4,0xFFFFFFFF}; //sets the key to FF to indicate the body is used
            tree.bodies_key[i] = (uint2){0xFFFFFFF4,0xFFFFFFFF}; //sets the key to FF to indicate the body is used            
          }
          else if(i == start)
          {
//             printf("Border: %d \n", i);
            test_key[i] = (uint2){0xFFFFFFF1,0xFFFFFFFF}; //sets the key to FF to indicate the body is used
            tree.bodies_key[i] = (uint2){0xFFFFFFF1,0xFFFFFFFF}; //sets the key to FF to indicate the body is used            
          }
          else if(i == end-1)
          {
//             printf("Border3: %d \n", i);
            test_key[i] = (uint2){0xFFFFFFF3,0xFFFFFFFF}; //sets the key to FF to indicate the body is used
            tree.bodies_key[i] = (uint2){0xFFFFFFF3,0xFFFFFFFF}; //sets the key to FF to indicate the body is used            
          }    
          else
          {
            test_key[i] = (uint2){0xFFFFFFF2,0xFFFFFFFF}; //sets the key to FF to indicate the body is used          
            tree.bodies_key[i] = (uint2){0xFFFFFFF2,0xFFFFFFFF}; //sets the key to FF to indicate the body is used   
          }
        }
      } //end for new nodes
      
      test_key.h2d();
      tree.bodies_key.h2d();

      
    } //end for to be fixed nodes
    #endif
//       if(level > 10)
//         exit(0);
  } //end for lvl
  

  //Put the last level + 1 index to 0,0 
  //so we dont need an extra if statement in the linking phase
  tree.level_list[level] = (uint2){0, 0};
  tree.level_list.h2d();
    
  int n_nodes  = offset;
  tree.n_nodes = n_nodes;
  
 
  /***** Link the tree ******/
  
  link_tree.set_arg<int>(0, &offset);   //Offset=number of nodes
  link_tree.set_arg<int>(9, &level);   //level=highest number of levels
  
  printf("Max level : %d \n", level);
  if(level >= MAXLEVELS)
  {
    cerr << "The tree has become too deep, the program will exit. \n";
    cerr << "Consider the removal of far away particles to prevent a too large box. \n";
    exit(0);
  }
  

  
  //The maximum number of levels that can be used is MAXLEVEl 
  //if max level is larger than that the program will exit
  
  link_tree.setWork(n_nodes, 128);
  printf("Link_tree: "); link_tree.printWorkSize();
  
  tree.n_levels = level-1;

  for(int i=0; i < level; i++)
    printf("%d\t%d\t%d\n", i, tree.level_list[i].x, tree.level_list[i].y);
 
  //Link the tree      
  link_tree.execute();
  
/*  
  
  tree.n_children.d2h();
  tree.node_key.d2h();
  for(int i=0; i < 10; i++)
  {
      printf("Children info: %d \t %X \t %X\t%X\t%X \n",i, tree.n_children[i],
             tree.node_key[i].x, tree.node_key[i].y, tree.node_key[i].z);
  }
  
//   exit(0);
  */
  
//   if(mpiGetRank() == 0)
  {
  //Test code, merges leaf nodes
//    walkCustomTreeOnHost(tree, validList);
  
  }
//   MPI_Barrier(MPI_COMM_WORLD);
//   exit(0);
  
  

  //After executing link_tree, the id_list contains for each node
  //the ID of its parent.
  //Valid_list contains for each node if its a leaf (valid) or a normal
  //node -> non_valid
  //Execute a split on the validList to get seperate id lists 
  //for the leafs and nodes 
    
  tree.leafNodeIdx.cmalloc(tree.n_nodes , false);
  
  //Split the leaf ids and non-leaf node ids
  gpuSplit(devContext, validList, tree.leafNodeIdx, 
                tree.n_nodes, &tree.n_leafs);     
                 
  printf("Total nodes: %d N_leafs: %d  non-leafs: %d \n", tree.n_nodes, tree.n_leafs, tree.n_nodes - tree.n_leafs);
  

/*  tree.leafPart2Body.zeroMem();
  
  expand_leaflist.set_arg<int>(0,    &tree.n_leafs);
  expand_leaflist.set_arg<cl_mem>(1, tree.leafNodeIdx.p());
  expand_leaflist.set_arg<cl_mem>(2, tree.node_bodies.p());
  expand_leaflist.set_arg<cl_mem>(3, tree.leafPart2Body.p());

  expand_leaflist.setWork(-1, NLEAF, tree.n_leafs); //n_leaf blocks with NLEAF threads each 
  expand_leaflist.printWorkSize();
  expand_leaflist.execute();
//   Dit uitwerken op de GPU, zorgen dat we een grote lijst krijgen die
//   de uitgewerkte cross referenties van leafs naar childs bevat. En die dan
//   opsplitsen voor de groups
  tree.leafPart2Body.d2h();*/
//   for(int i=0; i < tree.n; i++)
//     fprintf(stderr, "%d  \t-> %d \n", i, tree.leafPart2Body[i]);


  build_level_list.set_arg<int>(0, &tree.n_nodes);
  build_level_list.set_arg<int>(1, &tree.n_leafs);
  build_level_list.set_arg<cl_mem>(2, tree.leafNodeIdx.p());
  build_level_list.set_arg<cl_mem>(3, tree.node_bodies.p());
  build_level_list.set_arg<cl_mem>(4, validList.p());
  
  build_level_list.setWork(tree.n_nodes-tree.n_leafs, 128);
  
  validList.zeroMem();
  

  //Build the level list based on the leafIdx list
  //required for easy access in the compute node properties
  build_level_list.execute();  

  tree.node_level_list.cmalloc(level*2 , false);

  int levelThing;
  
  gpuCompact(devContext, validList, tree.node_level_list, 
             2*(tree.n_nodes-tree.n_leafs), &levelThing);             
  
  tree.node_level_list.d2h();
  
  //We only care about end positions, so compress the list:
  int j=0;
  for(int i=0; i < levelThing; i+=2, j++)
    tree.node_level_list[j] = tree.node_level_list[i];
  
  tree.node_level_list[j] =tree.node_level_list[levelThing-1]+1; //Add 1 to make it the end position
  levelThing = j+1;
  tree.node_level_list.h2d();
  
  printf("Finished level list \n");
  
  for(int i=0; i < levelThing; i++)
  {
    printf("node_level_list: %d \t%d\n", i, tree.node_level_list[i]);
  }

/* 
  //Determine if this node is a group...
  //Can only do this after the linking is complete
  //group if nchild <= NCRIT AND parent_nchild > NCRIT
  //Should do one more kernel launch over all nodes
  //and determine if the node is a group or not
  //Question is how the child can look to the parent 
  //     else if ((bj - bi) <= NCRIT)
  //       id_list[id] = id | (1 << 31);   //Group?
  define_groups.set_arg<int>(0, &offset);   //Offset=number of nodes
  define_groups.set_arg<cl_mem>(1, id_list.p());
  define_groups.set_arg<cl_mem>(2, tree.node_bodies.p());
  define_groups.set_arg<cl_mem>(3, validList.p());  
  define_groups.setWork(n_nodes, 128);  
  define_groups.execute();
  
  //Now compact validList to get the list of group ids
  tree.group_list.cmalloc(tree.n_nodes , false);

  gpuCompact(devContext, validList, tree.group_list, 
              tree.n_nodes, &tree.n_groups);
              
  printf("Number of groups according compact: %d \n", tree.n_groups);
*/

/* Peano Hilbert test! */
//   double t1=get_time();
//   fprintf(stderr, "Start creating and sorting PH curve \n");
//   tree.bodies_Ppos.d2h();
// //   sort_local_data_ph(&tree.bodies_Ppos[0], tree.n, tree.corner, &tree.peanoOrder[0]);
// 
//   tree.peanoOrder.h2d();
//   fprintf(stderr, "Done creating and sorting PH curve: Took %g \n", get_time()-t1);
//   
//   
  //Now on the device! woo
/*
  
  build_phkey_list.set_arg<cl_mem>(0,   tree.bodies_key.p());
  build_phkey_list.set_arg<cl_mem>(1,   tree.bodies_Ppos.p());
  build_phkey_list.set_arg<int>(2, &tree.n);
  build_phkey_list.set_arg<real4>(3,   &tree.corner);  
  build_phkey_list.set_arg<cl_mem>(4,   tree.peanoOrder.p());  
  
  build_phkey_list.setWork(tree.n, 256);
  printf("build_phkey_list: "); build_phkey_list.printWorkSize();
    
  build_phkey_list.execute();
  

tree.bodies_key.d2h();
tree.bodies_ids.d2h();
  for(int i=0; i < 10; i++)
{
    fprintf(stderr, "%d  \t (dev)key: %X %X  \t-->\t %d\n", i, tree.bodies_key[i].x, tree.bodies_key[i].y, tree.bodies_ids[i]);
}*/

/*
for(int i=0; i < tree.n; i++)
{
    fprintf(stderr, "%d  \t BODYTEST: %d \n", i, tree.bodies_ids[i]);
}

   exit(0);*/
// Hier gebleven,
// kijken of het mergen van leafs nu wel efficient wordt

  
  
/* End Peano Hilbert test! */

  //Compute the box size, the max length of one of the sides of the rectangle
  real size     = fmax(fabs(rMaxLocalTree.z - rMinLocalTree.z), 
                       fmax(fabs(rMaxLocalTree.y - rMinLocalTree.y),
                            fabs(rMaxLocalTree.x - rMinLocalTree.x)));
  real dist     = ((rMaxLocalTree.z - rMinLocalTree.z) * (rMaxLocalTree.z - rMinLocalTree.z) + 
                   (rMaxLocalTree.y - rMinLocalTree.y) * (rMaxLocalTree.y - rMinLocalTree.y) +
                   (rMaxLocalTree.x - rMinLocalTree.x) * (rMaxLocalTree.x - rMinLocalTree.x));      
                   
  float maxDist = sqrt(dist) / 10;
  maxDist *= maxDist; //Square since we dont do sqrt on device
                       
  fprintf(stderr,"Box max size: %f en max dist: %f \t %f en %f  \n", size, dist, sqrt(dist), maxDist);
  
  //maxDist = 50;
  
  validList.zeroMem();
  //The newest group creation method!
  define_groups.set_arg<int>(0, &tree.n);  
  define_groups.set_arg<cl_mem>(1, validList.p());    
  define_groups.set_arg<cl_mem>(2, tree.bodies_Ppos.p());
  define_groups.set_arg<float>(3, &maxDist);     
  define_groups.setWork(tree.n, 128);  
  define_groups.execute();
  
  //gpuCompact    
  gpuCompact(devContext, validList, compactList, tree.n*2, &validCount);
  
// if(mpiGetRank() == 1)
// {
//   compactList.d2h();
//   for(int i= 0; i < validCount; i++)
//   {
//     if(i > 1090)
//     fprintf(stderr,"Valid: %d \t %d \n", i, compactList[i]);
//   }
// }

  printf("Found number of groups: %d \n", validCount/2);

  tree.n_groups = validCount/2;
  //Now compact validList to get the list of group ids
  tree.group_list_test.cmalloc(tree.n_groups , false);  
  
  store_groups.set_arg<int>(0, &tree.n);  
  store_groups.set_arg<int>(1, &tree.n_groups);  
  store_groups.set_arg<cl_mem>(2, compactList.p());    
  store_groups.set_arg<cl_mem>(3, tree.body2group_list.p());     
  store_groups.set_arg<cl_mem>(4, tree.group_list_test.p());     
  store_groups.setWork(-1, NCRIT, tree.n_groups);  
  
  store_groups.printWorkSize();
  store_groups.execute();  

  /*
  tree.group_list_test.d2h();
  int histogram[NCRIT+1];
  for(int i=0; i < NCRIT+1; i++) histogram[i] = 0;
  for(int i=0; i < tree.n_groups; i++)
  {
    int grpSize = tree.group_list_test[i].y - tree.group_list_test[i].x;
    histogram[grpSize]++;    
  }
  
  int sum = 0;
  for(int i=0; i < NCRIT+1; i++)
  {
    sum += histogram[i];
    printf("%d\t%d\n", i,  histogram[i]);
  }
  printf("sum: %d \n", sum);
//   exit(0);
*/

  /*
  //The new group creation method!
  
  //Compute the number of groups based on the number of particles  
  int n_groups = tree.n / NCRIT;
  if((tree.n % NCRIT) > 0) n_groups++;

//   fprintf(stderr, "Number of particles: %d , number of groups: %d remaining parts: %d \n",
//           tree.n, n_groups, tree.n % NCRIT);
          
  //TODO make this a proper array
  tree.group_list_test.cmalloc(n_groups, false);          
  
  define_groups.set_arg<int>(0, &tree.n);  
  define_groups.set_arg<int>(1, &n_groups);  
  define_groups.set_arg<cl_mem>(2, tree.body2group_list.p());    
  define_groups.set_arg<cl_mem>(3, tree.group_list_test.p());     
  define_groups.setWork(-1, NCRIT, n_groups);  
  
  define_groups.printWorkSize();
  define_groups.execute();
  
  tree.n_groups = n_groups;
  
  
//   group_list.d2h();
//   
//   for(int i=0; i < n_groups+4; i++)
//   {
//     fprintf(stderr, "Grp: %d Start: %d End: %d \n", i, group_list[i].x, group_list[i].y);
//   }
//   
//   tree.body2group_list.d2h();
//   for(int i=0; i < tree.n; i++)
//   {
//     fprintf(stderr, "Particle: %d Grp: %d \n", i, tree.body2group_list[i]);
//   }  
//  
//   exit(0);
  
  
  
  */
  
  
  
  //Group statistics
/*  tree.group_list.d2h();
  tree.node_bodies.d2h();
  int histogram[65];
  for(int i=0; i < 65; i++) histogram[i] = 0;
  for(int i=0; i < tree.n_groups; i++)
  {
      uint2 bij  = tree.node_bodies[tree.group_list[i]];  
      int grpSize    =  bij.y - (bij.x & ILEVELMASK);
      histogram[grpSize]++;    
  }
  
//   int sum = 0;
//   for(int i=0; i < 65; i++)
//   {
//     sum += histogram[i];
//     printf("%d\t%d\n", i,  histogram[i]);
//   }
//   printf("sum: %d \n", sum);
//   exit(0);*/
  
  //Memory allocation for the valid group lists
  if(tree.active_group_list.get_size() > 0)
  {
    tree.active_group_list.cresize(tree.n_groups, false);
    tree.activeGrpList.cresize(tree.n_groups, false);     
    
  }
  else
  {
    tree.active_group_list.cmalloc(tree.n_groups, false);
    tree.activeGrpList.cmalloc(tree.n_groups, false);     
  }
 // tree.group_data.cmalloc(tree.n_groups, false);
  
  //Alloc memory for the body2group list
  //tree.body2group_list.cmalloc(tree.n, false);

  //assign a grp ID to each particle
  
/*  
  build_body2group_list.set_arg<int>(0, &tree.n_groups);  
  build_body2group_list.set_arg<cl_mem>(1, tree.group_list.p());
  build_body2group_list.set_arg<cl_mem>(2, tree.node_bodies.p());
  build_body2group_list.set_arg<cl_mem>(3, tree.body2group_list.p());
  
  //Start as many blocks as there are groups, so each particle in a 
  //group gets that group id block
//   build_body2group_list.setWork(-1, 64, tree.n_groups);
  build_body2group_list.setWork(-1, NCRIT, tree.n_groups);

  //Build the level list based on the leafIdx list
  //required for easy access in the compute node properties
  build_body2group_list.execute();
 */

  printf("Tree built complete!\n");


  /*************************/

}
