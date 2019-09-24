"""
basic graph class and algorithmes

UnionFind and MinimumSpanningTree taken from PADS:

  a library of Python Algorithms and Data Structures
  implemented by David Eppstein of the University of California, Irvine.

  The current version of PADS may be found at
  <http://www.ics.uci.edu/~eppstein/PADS/>, as individual files or as a
  git repository that may be copied by the command line
  
  git clone http://www.ics.uci.edu/~eppstein/PADS/.git
  
  PADS is hereby placed in the public domain; you may use the code in PADS
  for any purpose whatsoever.  We make no guarantee of quality,
  completeness, correctness, persistence or consistency of APIs, or support.

"""

class UnionFind(object):
    """Union-find data structure.

    Each unionFind instance X maintains a family of disjoint sets of
    hashable objects, supporting the following two methods:

    - X[item] returns a name for the set containing the given item.
      Each set is named by an arbitrarily-chosen one of its members; as
      long as the set remains unchanged it will keep the same name. If
      the item is not yet part of a set in X, a new singleton set is
      created for it.

    - X.union(item1, item2, ...) merges the sets containing each item
      into a single larger set.  If any item is not yet part of a set
      in X, it is added to X as one of the members of the merged set.
    """

    def __init__(self):
        """Create a new empty union-find structure."""
        self.weights = {}
        self.parents = {}

    def __getitem__(self, object):
        """Find and return the name of the set containing the object."""

        # check for previously unknown object
        if object not in self.parents:
            self.parents[object] = object
            self.weights[object] = 1
            return object

        # find path of objects leading to the root
        path = [object]
        root = self.parents[object]
        while root != path[-1]:
            path.append(root)
            root = self.parents[root]

        # compress the path and return
        for ancestor in path:
            self.parents[ancestor] = root
        return root
        
    def __iter__(self):
        """Iterate through all items ever found or unioned by this structure."""
        return iter(self.parents)

    def union(self, *objects):
        """Find the sets containing the objects and merge them all."""
        roots = [self[x] for x in objects]
        heaviest = max([(self.weights[r],r) for r in roots], key = lambda x: x[0])[1]
        for r in roots:
            if r != heaviest:
                self.weights[heaviest] += self.weights[r]
                self.parents[r] = heaviest

    def sets(self):
      sets={}
      for v in self.parents:
        sets.setdefault(self[v],set()).add(v)
      return list(sets.values())

                

class Graph(dict):
    def add_edge(self, n1, n2, w):
        if callable(w): w=w(n1,n2)
        self.setdefault(n1, {}).update({n2: w})
        self.setdefault(n2, {}).update({n1: w})

    def remove_edge(self, n1, n2):
        self[n1].pop(n2)
        self[n2].pop(n1)

    def add_node(self,n):
        self.setdefault(n, {})
    
    def all_edges(self):
        return [(self[u][v],u,v) for u in self for v in self[u]]


def MinimumSpanningTree(G):
    """
    Return the minimum spanning tree of an undirected graph G.
    G should be represented in such a way that G[u][v] gives the
    length of edge u,v, and G[u][v] should always equal G[v][u].
    The tree is returned as a list of edges.
    """
    
    # Kruskal's algorithm: sort edges by weight, and add them one at a time.
    # We use Kruskal's algorithm, first because it is very simple to
    # implement once UnionFind exists, and second, because the only slow
    # part (the sort) is sped up by being built in to Python.
    subtrees = UnionFind()
    tree = []
    edges = [(G[u][v],u,v) for u in G for v in G[u]]
    edges.sort(key=lambda x:x[0])
    for W,u,v in edges:
        if subtrees[u] != subtrees[v]:
            tree.append((W,u,v))
            subtrees.union(u,v)
    return tree        

def MinimumSpanningTreeFromEdges(edges):
    """
    Return the minimum spanning tree of an undirected graph G.
    This version runs directly from an edgelist. An edge is a triple
    (w,u,v), such that u,v are nodes, w is the length of the edge.
    The tree is returned as a list of edges.
    """
    
    # Kruskal's algorithm: sort edges by weight, and add them one at a time.
    # We use Kruskal's algorithm, first because it is very simple to
    # implement once UnionFind exists, and second, because the only slow
    # part (the sort) is sped up by being built in to Python.
    subtrees = UnionFind()
    tree = []
    edges.sort(key=lambda x:x[0])
    for W,u,v in edges:
        if subtrees[u] != subtrees[v]:
            tree.append((W,u,v))
            subtrees.union(u,v)
    return tree        


def ConnectedComponents(G):
    """ 
    Return the connected components of a graph. G should be 
    represented in such a way that G[u] gives the edges from u in a way 
    that and if v in G[u] than u in G[v]. the connected components are 
    returned as sets of nodes. 
    """
    u=UnionFind() 
    for v in G:
      nset=set(G[v])
      nset.add(v)
      u.union(*nset)
    return u.sets()

def ConnectedComponentsFromEdges(edges):
    """ 
    Return the connected components of a graph from a list of egdes. 
    the connected components are returned as sets of nodes. note this does
    not find singletons. 
    """
    u=UnionFind() 
    for e in edges:
      u.union(e[1],e[2])
    return u.sets()


if __name__=="__main__":
  graph = Graph()

  graph.add_edge(0, 1, 1.0)
  graph.add_edge(1, 2, 1.0)
  graph.add_edge(2, 0, 1.0)
  
  graph.add_edge(3, 4, 1.0)
  graph.add_edge(4, 5, 1.0)
  graph.add_edge(5, 3, 1.0)
  
  print(graph[0])
  
  first, second = ConnectedComponents(graph)
  print(first)
  print(second)
  
  print(MinimumSpanningTree(graph))
  
