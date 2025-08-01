#include <iostream>
// #include "particle.h"
// #include "boundary.h"
#include "node.h"
#include "wtime.h"

std::vector<particle> node::ptcl;
std::vector<node>     node::node_heap;
std::vector<std::pair<node*, node*> > node::pair_list;

void perftest(int nbody, float scale)
{
  std::vector<particle> &ptcl = node::ptcl;
  double t0 = wtime();
  for(int i=0; i<nbody; i++){
    ptcl[i].keygen(scale);
  }
  double t1 = wtime();
  std::sort(ptcl.begin(), ptcl.end(), cmp_particle_key());
  double t2 = wtime();

  node::node_heap.push_back(node());
  node &root = node::node_heap[0];
  for(int i=0; i<nbody; i++){
    root.push_particle(i, 60);
  }
  double t3 = wtime();
  root.make_boundary();
  double t4 = wtime();
#ifdef SLOW
#pragma omp parallel for
  for(int i=0; i<nbody; i++){
    ptcl[i] << root;
  }
#else
#  ifdef _OPENMP
  std::vector<node *> group_list;
  root.find_group_node(2000, group_list);
#  pragma omp parallel for schedule(dynamic)
  for(int i=0; i<(int)group_list.size(); i++){
    *group_list[i] << root;
  }
#  else
  root << root;
#  endif
#endif
  double t5 = wtime();

  std::cerr << "key gen: " << t1-t0 << " sec" << std::endl;
  std::cerr << "sort:    " << t2-t1 << " sec" << std::endl;
  std::cerr << "tree:    " << t3-t2 << " sec" << std::endl;
  std::cerr << "boundary:" << t4-t3 << " sec" << std::endl;
  std::cerr << "search  :" << t5-t4 << " sec" << std::endl;
  std::cerr << "per N:   " << 1.e6*(t5-t4)/nbody << " usec" << std::endl;
  std::cerr << "nbody = " << nbody 
    << ", nnode = " << node::node_heap.size() 
    << std::endl;
  std::cerr << std::endl;
  node::clear();
}

int main(){
  std::vector<particle> &ptcl = node::ptcl;
  std::vector<int> nnb;
  int idum, nbody;
#ifdef BERCZIK
  double fdum;
  std::cin >> idum >> nbody >> fdum;
#else
  std::cin >> idum >> nbody;
#endif
  node::allocate(nbody, nbody);
  nnb.reserve(nbody);
  typedef boundary<float>  Boundary;
  Boundary bound;
  for(int i=0; i<nbody; i++){
    particle::vec pos;
    float h;
    int _nnb;
#ifdef BERCZIK
    particle::vec vel;
    float mass;
    std::cin >> idum >> mass >> pos >> vel;
    assert(i == idum);
    h = 2. * pow(nbody, -1./3.);
#else
    std::cin >> pos >> h >> _nnb;
    h *= 2.f;
#endif
    ptcl.push_back(particle(i, pos, h));
    nnb.push_back(_nnb);
    // ptcl[i].nnb = nnb;
    bound.merge(Boundary(pos, 0.0f));
  }
  std::cerr << bound.min << std::endl;
  std::cerr << bound.max << std::endl;
  fvec3 vmax = maxeach(bound.min.abseach(), bound.max.abseach());
  float rmax = std::max(vmax.x, std::max(vmax.x, vmax.y));
  float scale = 1.f;
  while(scale < rmax) scale *= 2.f;

#ifdef BERCZIK
  for(int run=0; run<5; run++){
    for(int i=0; i<nbody; i++){
      ptcl[i].nnb_gath = 0;
      ptcl[i].nnb_scat = 0;
    }
    perftest(nbody, scale);
    for(int i=0; i<nbody; i++){
      ptcl[i].h *= pow(30./ptcl[i].nnb_gath, 1./3.);
    }
  }
#else
  perftest(nbody, scale);
#endif
  // std::cerr << root.bound_inner.min << std::endl;
  // std::cerr << root.bound_inner.max << std::endl;
  // root.dump_tree(0, ptcl, node_heap, std::cout);
  /*
     for(int i=0; i<nbody; i++){
     std::cout << ptcl[i].id << " " << ptcl[i].pos <<  " " << ptcl[i].key.val << std::endl;
     }
     */
#ifdef BERCZIK
  std::sort(ptcl.begin(), ptcl.end(), cmp_particle_id());
  for(int i=0; i<nbody; i++){
    std::cout << ptcl[i].nnb_gath << " " 
      << ptcl[i].nnb_scat << std::endl;
  }
#else
  std::sort(ptcl.begin(), ptcl.end(), cmp_particle_id());
  for(int i=0; i<nbody; i++){
    if(ptcl[i].nnb_gath != nnb[ptcl[i].id]){
      std::cout << ptcl[i].id << " " 
        << nnb[ptcl[i].id]  << " "
        << ptcl[i].nnb_gath << " " 
        << ptcl[i].nnb_scat << std::endl;
    }
  }
#endif
  return 0;
}
