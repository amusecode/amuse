#include "kdtree.h"

inline double rand(const double min, const double max)
{
  return min + drand48()*(max - min);
}

int range_search_n2(std::vector<int> &id_list, const dvec3 &pos, const double h, const std::vector<kdTree::Particle> &ptcl)
{
  const int n = ptcl.size();
  for (int i = 0; i < n; i++)
    if ((ptcl[i].pos - pos).norm2() < h*h)
      id_list.push_back(i);

  return id_list.size();
}

int main(int argc, char * argv[])
{

  int n = 10323;
  std::vector<kdTree::Particle> ptcl;
  for (int i = 0; i < n; i++)
    ptcl.push_back(kdTree::Particle(i, dvec3(rand(-1,1), rand(-1,1), rand(-1,1))));

#if 0
  kdTree kd(ptcl);
#else
  kdTree kd;
  for (int i = 0; i < n; i++)
    kd.push_ptcl(ptcl[i]);
  kd.build();
#endif

  const double h = 0.2;
  const dvec3 pos(0.0,0.0,0.0);

  std::vector<int> list, list_n2;
  fprintf(stderr, "n_ngb= %d n2= %d\n",
      kd.range_search(list,    pos, h), 
      range_search_n2(list_n2, pos, h, ptcl)
      );

  std::sort(list.begin(), list.end());
  std::sort(list_n2.begin(), list_n2.end());
  assert(list.size() == list_n2.size());

  for (size_t i = 0; i < list.size(); i++)
    assert(list[i] == list_n2[i]);



}
