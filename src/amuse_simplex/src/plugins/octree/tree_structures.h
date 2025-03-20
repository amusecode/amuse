#ifndef TREE_STRUCTURES_H
#define TREE_STRUCTURES_H

#ifdef INDIE
  #include "structs.h"
#else
  #include "Structs.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>
#include <float.h>

#include "octree.h"







using namespace std;

/* class Vertex{ */

/*  public: */

/*   double x; */
/*   double y; */
/*   double z; */


/* }; */

class vertex_tree{

  public:

    //! empty constructor
    vertex_tree();                 

    //! destructor
    ~vertex_tree();
 
    //! initialises empty octree with size tree_size
    void init_octree( unsigned int& tree_size, double& aborder_box );

    //! delete the octree
    void delete_octree();

    //! Insert a vertex in the tree
    void insert_vertex( const float& x, const float& y, const float& z, Vertex& tmp );
  
    //! Search the tree for vertices in a certain box and put the in a vector
    vector<unsigned long long int> find_vertices_in_box( const double& x_min, const double& x_max, const double& y_min, 
					       const double& y_max, const double& z_min, const double& z_max );

    //! Search the tree for vertices in a certain box and put the in a vector
    vector<unsigned int> find_vertices_in_box_2( const double& x_min, const double& x_max, const double& y_min, 
						 const double& y_max, const double& z_min, const double& z_max );

    //! Set the size of the border around the domain
    void set_border_box( const double& _border_box_ ) { border_box = _border_box_; }
    //! Return the size of the octree
    const double& get_border_box() const{ return border_box; }

    //! Set the size of the octree
    void set_size( const unsigned int& _size_ ) { size = _size_; }
    //! Return the size of the octree
    const unsigned int& get_size() const{ return size; }

    //! Set the number of nodes in the octree
    void set_nodes();
    //! Return the node of the octree
    const unsigned int& get_nodes();



 private:
 
    double border_box;                       //!< minimum coordinate of the box
    unsigned int nodes;                      //!< number of nodes in the octree
    unsigned int size;                       //!< size of the octree
    Octree< vector<Vertex> > *octree;        //!< octree containing vertices

};




#endif
