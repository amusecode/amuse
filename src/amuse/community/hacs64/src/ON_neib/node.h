#include <iostream>
#include <omp.h>
#include "particle.h"
#include "boundary.h"

struct node{
  static const int NLEAF = 64;
  static std::vector<particle> ptcl;
  static std::vector<node> node_heap;
  static std::vector<std::pair<node*, node*> > pair_list;
  static void allocate(int nptcl, int nnode){
    ptcl.reserve(nptcl);
    node_heap.reserve(nnode);
  }
  static void clear(){
    ptcl.clear();
    node_heap.clear();
  }
  typedef boundary<float> Boundary;

  int np;     // number of particle;
  int depth;
  int pfirst; // first particle
  int cfirst; // first child
  Boundary bound_inner;
  Boundary bound_outer;

  node()           : np(0), depth(     0), pfirst(-1), cfirst(-1) {}
  node(int _depth) : np(0), depth(_depth), pfirst(-1), cfirst(-1) {}

  bool is_leaf() const{
    return np < NLEAF;
  }
  void push_particle(
      int paddr, 
      int rshift){
    assert(rshift >= 0);
    if(!is_leaf()){ // assign recursively
      int ic = ptcl[paddr].octkey(rshift);
      node &child = node_heap[cfirst + ic];
      child.push_particle(paddr, rshift-3);
    }else{
      if(pfirst == -1){
        assert(np==0);
        pfirst = paddr;
      }
    }
    np++;
    if(np == NLEAF){ // shi's just become a mother
      assert(pfirst >= 0);
      cfirst = node_heap.size();
#if 0
      for(int ic=0; ic<8; ic++){
        node_heap.push_back(node(1+depth));
      }
#else
      size_t new_size = node_heap.size() + 8;
      assert(node_heap.capacity() >= new_size);
      node_heap.resize(new_size, node(1+depth));
#endif
      for(int addr = pfirst; addr < pfirst+np; addr++){
        int ic = ptcl[addr].octkey(rshift);
        node &child = node_heap[cfirst + ic];
        child.push_particle(addr, rshift-3);
      }
    }
  }
  void dump_tree(
      int level,
      std::ostream &ofs = std::cout) const{
    if(is_leaf()){
      for(int ip=0; ip<np; ip++){
        const particle &p = ptcl[ip+pfirst];
        for(int i=0; i<level; i++) ofs << " ";
        ofs << p.pos << std::endl;
      }
      ofs << std::endl;
    }else{
      for(int i=0; i<level; i++) ofs << ">";
      ofs << std::endl;
      for(int ic=0; ic<8; ic++){
        const node &child = node_heap[cfirst + ic];
        child.dump_tree(level+1, ofs);
      }
    }
  }
  void make_boundary(){
    if(is_leaf()){
      for(int ip=0; ip<np; ip++){
        const particle &p = ptcl[ip+pfirst];
        bound_inner.merge(Boundary(p.pos));
        bound_outer.merge(Boundary(p.pos, p.h));
      }
    }else{
      for(int ic=0; ic<8; ic++){
        node &child = node_heap[cfirst + ic];
        if(child.np > 0){
          child.make_boundary();
          bound_inner.merge(child.bound_inner);
          bound_outer.merge(child.bound_outer);
        }
      }
    }
  }
  static void find_neib_beween_leaves(const node &ileaf, const node &jleaf){
    assert(0);
#if 1 
    for(int i=0; i<ileaf.np; i++){
      particle &ip = ptcl[i+ileaf.pfirst];
      Boundary ibound(ip.pos, ip.h);
      if(not_overlapped(ibound, jleaf.bound_inner)) continue;
      float h2 = ip.h * ip.h;
      for(int j=0; j<jleaf.np; j++){
        particle &jp = ptcl[j+jleaf.pfirst];
        if((jp.pos - ip.pos).norm2() < h2){
#if 1
          ip.nnb_gath++;
          ip.ngb_list.push_back(jp.id);
#endif

//          jp.nnb_scat++;
        }
      }
    }
#else
    typedef float v4sf __attribute__ ((vector_size(16)));
    for(int i=0; i<ileaf.np; i+=4){
      int ni = std::min(4, ileaf.np-i);
      particle &ip0 = ptcl[0+i+ileaf.pfirst];
      particle &ip1 = ptcl[1+i+ileaf.pfirst];
      particle &ip2 = ptcl[2+i+ileaf.pfirst];
      particle &ip3 = ptcl[3+i+ileaf.pfirst];
      if(not_overlapped(Boundary(ip0.pos, ip0.h), jleaf.bound_inner)
          && not_overlapped(Boundary(ip1.pos, ip1.h), jleaf.bound_inner)
          && not_overlapped(Boundary(ip2.pos, ip2.h), jleaf.bound_inner)
          && not_overlapped(Boundary(ip3.pos, ip3.h), jleaf.bound_inner)
        ) continue;
      v4sf xyzh0 = *(v4sf *)&ip0.pos;
      v4sf xyzh1 = *(v4sf *)&ip1.pos;
      v4sf xyzh2 = *(v4sf *)&ip2.pos;
      v4sf xyzh3 = *(v4sf *)&ip3.pos;
      v4sf tmp0 = __builtin_ia32_unpcklps(xyzh0, xyzh2); // y2|y0|x2|x0
      v4sf tmp1 = __builtin_ia32_unpcklps(xyzh1, xyzh3); // y3|y1|x3|x1
      v4sf tmp2 = __builtin_ia32_unpckhps(xyzh0, xyzh2); // h2|h0|z2|z0
      v4sf tmp3 = __builtin_ia32_unpckhps(xyzh1, xyzh3); // h3|h1|z3|z1
      v4sf xi = __builtin_ia32_unpcklps(tmp0, tmp1);
      v4sf yi = __builtin_ia32_unpckhps(tmp0, tmp1);
      v4sf zi = __builtin_ia32_unpcklps(tmp2, tmp3);
      v4sf hi = __builtin_ia32_unpckhps(tmp2, tmp3);
      v4sf h2 = hi * hi;
      int nnb0=0, nnb1=0, nnb2=0, nnb3=0;
      for(int j=0; j<jleaf.np; j++){
        int nnbj = 0;
        particle &jp = ptcl[j+jleaf.pfirst];
        v4sf vj = *(v4sf *)&jp.pos;
        v4sf xj = __builtin_ia32_shufps(vj, vj, 0x00);
        v4sf yj = __builtin_ia32_shufps(vj, vj, 0x55);
        v4sf zj = __builtin_ia32_shufps(vj, vj, 0xaa);
        v4sf dx = xj - xi;
        v4sf dy = yj - yi;
        v4sf dz = zj - zi;
        v4sf r2 = dx*dx + dy*dy + dz*dz;
        int flags = __builtin_ia32_movmskps(
            (v4sf)__builtin_ia32_cmpltps(r2, h2));
        if(flags & 1) nnb0++, nnbj++;
        if(flags & 2) nnb1++, nnbj++;
        if(flags & 4) nnb2++, nnbj++;
        if(flags & 8) nnb3++, nnbj++;
        jp.nnb_scat += nnbj;
      }
      if(ni > 0) ip0.nnb_gath += nnb0;
      if(ni > 1) ip1.nnb_gath += nnb1;
      if(ni > 2) ip2.nnb_gath += nnb2;
      if(ni > 3) ip3.nnb_gath += nnb3;
    }
#endif
  }
  friend void operator << (node &inode, node &jnode){
    if(overlapped(inode.bound_outer, jnode.bound_inner)){
      bool itravel = false;
      bool jtravel = false;
      if(inode.is_leaf()){
        if(jnode.is_leaf()){
          find_neib_beween_leaves(inode, jnode);
          return;
        }else{
          jtravel = true;
        }
      }else{
        if(jnode.is_leaf()){
          itravel = true;
        }else{
          if(inode.depth < jnode.depth){
            itravel = true;
          }else{
            jtravel = true;
          }
        }
      }
      if(itravel){
        for(int i=0; i<8; i++){
          node &ichild = node_heap[i+inode.cfirst];
          if(ichild.np == 0) continue;
          ichild << jnode;
        }
        return;
      }
      if(jtravel){
        for(int j=0; j<8; j++){
          node &jchild = node_heap[j+jnode.cfirst];
          if(jchild.np == 0) continue;
          inode << jchild;
        }
        return;
      }
    }
  }
  void find_group_node(
      int ncrit,
      std::vector<node *> &group_list){
    if(np < ncrit){
      group_list.push_back(this);
    }else{
      for(int ic=0; ic<8; ic++){
        node &child = node_heap[cfirst + ic];
        child.find_group_node(ncrit, group_list);
      }
    }
  }
  friend void operator << (particle &ip, node &jnode){
    Boundary bi(ip.pos, ip.h);
    if(overlapped(bi, jnode.bound_inner)){
      if(jnode.is_leaf()){
        float h2 = ip.h * ip.h;
        for(int j=0; j<jnode.np; j++){
          particle &jp = ptcl[j+jnode.pfirst];
          if((jp.pos - ip.pos).norm2() < h2){
#if 1
            ip.nnb_gath++;
            ip.ngb_list.push_back(jp.id);
#endif
//            jp.nnb_scat++;
          }
        }
      }else{
        for(int j=0; j<8; j++){
          node &jchild = node_heap[j+jnode.cfirst];
          if(jchild.np == 0) continue;
          ip << jchild;
        }
      }
    }
  }
};
