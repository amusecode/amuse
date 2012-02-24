#ifndef _KDTREE_H_
#define _KDTREE_H_

#include <vector>
#include <stack>
#include "localassert.h"
#include <algorithm>
#include "vector3.h"

class kdTree
{
  public:
    struct Particle
    {
      int id;
      dvec3 pos;
      Particle() {}
      Particle(const int _id, const dvec3 _pos) : id(_id), pos(_pos) {}
      template<int SPLIT>
        bool operator() (const Particle &lhs, const Particle &rhs)
        {
          if      (SPLIT == 0) return lhs.pos.x < rhs.pos.x;
          else if (SPLIT == 1) return lhs.pos.y < rhs.pos.y;
          else if (SPLIT == 2) return lhs.pos.z < rhs.pos.z;
          else assert(false);
          return false;
        }
    };

    struct cmp_x
    {
      bool operator() (const Particle &lhs, const Particle &rhs)
      {
        return lhs.pos.x < rhs.pos.x;
      }
    };
    
    struct cmp_y
    {
      bool operator() (const Particle &lhs, const Particle &rhs)
      {
        return lhs.pos.y < rhs.pos.y;
      }
    };
    
    struct cmp_z
    {
      bool operator() (const Particle &lhs, const Particle &rhs)
      {
        return lhs.pos.z < rhs.pos.z;
      }
    };


  protected:
    struct Node
    {
      int split;
      Particle ptcl;
    };

    std::vector<Node> tree;
    std::vector<Particle> ptcl;
    int max_nodes, max_depth; 

    void build_left_balanced_tree(
        const int n_node,
        const int nleft,
        const int n,
        const int depth)
    {
      assert(n >= 0);
      if (n == 0) return;

      max_nodes = std::max(n_node, max_nodes);
      max_depth = std::max(max_depth, depth);

      const int  split = depth%3;

      switch(split)
      {
        case 0:
          std::sort(ptcl.begin()+nleft, ptcl.begin()+nleft+n, cmp_x());
          break;
        case 1:
          std::sort(ptcl.begin()+nleft, ptcl.begin()+nleft+n, cmp_y());
          break;
        case 2:
          std::sort(ptcl.begin()+nleft, ptcl.begin()+nleft+n, cmp_z());
          break;
        default:
          assert(false);
      } 

      int m = 1;
      while (m <= n) 
        m = (m << (1));
      m = (m >> (1));

      const int r = n - (m - 1);
      int lt, rt;
      if (r <= m/2) 
      {
        lt = (m-2)/2 + r;
        rt = (m-2)/2;
      } 
      else 
      {
        lt = (m-2)/2 + m/2;
        rt = (m-2)/2 - m/2 + r;
      }

      const int median  = lt;
      const Particle pm = ptcl[nleft + median];
      const int n_left  = median;
      const int n_right = n - median - 1;
      
      Node &node = tree[n_node];  
      node.split = split;
      node.ptcl  = pm;

      build_left_balanced_tree(2*n_node,   nleft,          n_left,  depth+1);
      build_left_balanced_tree(2*n_node+1, nleft+median+1, n_right, depth+1);
    }

  public:

    kdTree() : max_nodes(0), max_depth(0) {}
    kdTree(const std::vector<Particle> _ptcl) : ptcl(_ptcl) {
      build();
    }
    void push_ptcl(const Particle &p) 
    {
      ptcl.push_back(p);
    }
    const Particle& get_ptcl(const int id) const {return ptcl[id];}

    void build()
    {
      max_depth = max_nodes = 0;
      assert(!ptcl.empty());
      tree.resize(ptcl.size()*2);
      build_left_balanced_tree(1, 0, ptcl.size(), 0);
      fprintf(stderr, "<tree info: max_depth= %d  max_nodes= %d> ... \n", max_depth, max_nodes);
    }

    int range_search(std::vector<int> &id_list, const dvec3 &pos, const double h) const
    {
      std::stack<int> stack;

      stack.push(1);
      while(!stack.empty())
      {
        const int node_id = stack.top();
        stack.pop();
       
        const Node &node = tree[node_id];
        const int   next = node_id << 1;
        
        if (next <= max_nodes)
        {
          switch(node.split)
          {
            case 0:
              if (pos.x - h < node.ptcl.pos.x) stack.push(next);
              if (next + 1 <= max_nodes && (pos.x + h > node.ptcl.pos.x)) stack.push(next+1);
              break;
            case 1:
              if (pos.y - h < node.ptcl.pos.y) stack.push(next);
              if (next + 1 <= max_nodes && (pos.y + h > node.ptcl.pos.y)) stack.push(next+1);
              break;
            case 2:
              if (pos.z - h < node.ptcl.pos.z) stack.push(next);
              if (next + 1 <= max_nodes && (pos.z + h > node.ptcl.pos.z)) stack.push(next+1);
              break;
            default:
              assert(false);
          }
        }

        if ((node.ptcl.pos - pos).norm2() < h*h)
        {
                id_list.push_back(node.ptcl.id);
        }
      }

      return id_list.size();
    }

};
#endif 
