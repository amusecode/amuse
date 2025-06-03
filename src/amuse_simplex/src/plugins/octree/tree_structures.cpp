#include "tree_structures.h"
#include "time.h"

using namespace std;


//constructor
vertex_tree::vertex_tree(){

  size = 0;
  nodes = 0;

}

//destructor
vertex_tree::~vertex_tree(){

  delete octree;

}

// Initialise octree with size tree_size
void vertex_tree::init_octree( unsigned int& tree_size, double& aborder_box ){

   octree = new Octree< vector<Vertex> >( tree_size );
   size = tree_size;
   border_box = aborder_box;
  
}

//delete the octree
void vertex_tree::delete_octree(){

  //only delete the octree when size is zero
  if(size){
    delete octree;
  }


}

//Insert vertex in the octree
void vertex_tree::insert_vertex( const float& x, const float& y, const float& z, Vertex& tmp ){

  //convert to coordinates between 0 and 1
  double temp_x = ( x + border_box )/( 1.0 + 2*border_box );
  double temp_y = ( y + border_box )/( 1.0 + 2*border_box );
  double temp_z = ( z + border_box )/( 1.0 + 2*border_box );

  //convert to tree coordinates
  unsigned int u = (unsigned) floor( temp_x*get_size() );
  unsigned int v = (unsigned) floor( temp_y*get_size() );
  unsigned int w = (unsigned) floor( temp_z*get_size() );

  //make sure positions fit into tree
  if( u >= get_size() ){
    u = get_size() - 1;
  }

  if( v >= get_size() ){
    v = get_size() - 1;
  }

  if( w >= get_size() ){
    w = get_size() - 1;
  }

  (*octree)(u,v,w).push_back( tmp );

}

void vertex_tree::set_nodes(){

  nodes = (*octree).nodes();
  
}

const unsigned int& vertex_tree::get_nodes(){

  set_nodes();
  return nodes;

}

//get the id's of vertices that are inside box
vector<unsigned long long int> vertex_tree::find_vertices_in_box( const double& x_min, const double& x_max, const double& y_min, 
							const double& y_max, const double& z_min, const double& z_max ){
  

  unsigned int tree_size = size;
         
  //convert to coordinates between 0 and 1
  double temp_x_min = ( x_min + border_box )/( 1.0 + 2*border_box );
  double temp_y_min = ( y_min + border_box )/( 1.0 + 2*border_box );
  double temp_z_min = ( z_min + border_box )/( 1.0 + 2*border_box );

  double temp_x_max = ( x_max + border_box )/( 1.0 + 2*border_box );
  double temp_y_max = ( y_max + border_box )/( 1.0 + 2*border_box );
  double temp_z_max = ( z_max + border_box )/( 1.0 + 2*border_box );

  //convert to tree coordinates
  int u_min = (int) floor( temp_x_min*tree_size ) - 1;
  int v_min = (int) floor( temp_y_min*tree_size ) - 1;
  int w_min = (int) floor( temp_z_min*tree_size ) - 1;

  unsigned int u_max = (unsigned) floor( temp_x_max*tree_size ) + 1;
  unsigned int v_max = (unsigned) floor( temp_y_max*tree_size ) + 1;
  unsigned int w_max = (unsigned) floor( temp_z_max*tree_size ) + 1;

  //make sure positions fit into tree
  if( u_min < 0 ){
    u_min = 0;
  }

  if( v_min < 0 ){
    v_min = 0;
  }

  if( w_min < 0 ){
    w_min = 0;
  }

  if( u_max > get_size() ){
    u_max = get_size();
  }

  if( v_max > get_size() ){
    v_max = get_size();
  }

  if( w_max > get_size() ){
    w_max = get_size();
  }

  //vector to store the ids of vertices in the box
  vector<unsigned long long int> ret_vector;
  vector<Vertex>::iterator it;

  //loop through all octree cells looking for the vertices that belong in the selection box

  for( unsigned int k=w_min; k<w_max; k++ ){
    for( unsigned int j=v_min; j<v_max; j++){
      for( unsigned int i=u_min; i<u_max; i++){

	for( it=(*octree)(i,j,k).begin(); it!=(*octree)(i,j,k).end(); it++ ){

	  if( it->get_x() >= x_min && it->get_x() < x_max && 
	      it->get_y() >= y_min && it->get_y() < y_max &&
	      it->get_z() >= z_min && it->get_z() < z_max ){

 	    ret_vector.push_back( it->get_vertex_id() );

	  }



	}


      }
    }
  }
  
  return ret_vector;
  
}

//get the id's of vertices that are inside box
vector<unsigned int> vertex_tree::find_vertices_in_box_2( const double& x_min, const double& x_max, const double& y_min, 
							  const double& y_max, const double& z_min, const double& z_max ){
  

  unsigned int tree_size = size;
         
  //convert to coordinates between 0 and 1
  double temp_x_min = ( x_min + border_box )/( 1.0 + 2*border_box );
  double temp_y_min = ( y_min + border_box )/( 1.0 + 2*border_box );
  double temp_z_min = ( z_min + border_box )/( 1.0 + 2*border_box );

  double temp_x_max = ( x_max + border_box )/( 1.0 + 2*border_box );
  double temp_y_max = ( y_max + border_box )/( 1.0 + 2*border_box );
  double temp_z_max = ( z_max + border_box )/( 1.0 + 2*border_box );

  //convert to tree coordinates
  int u_min = (int) floor( temp_x_min*tree_size ) - 1;
  int v_min = (int) floor( temp_y_min*tree_size ) - 1;
  int w_min = (int) floor( temp_z_min*tree_size ) - 1;

  unsigned int u_max = (unsigned) floor( temp_x_max*tree_size ) + 1;
  unsigned int v_max = (unsigned) floor( temp_y_max*tree_size ) + 1;
  unsigned int w_max = (unsigned) floor( temp_z_max*tree_size ) + 1;

  //make sure positions fit into tree
  if( u_min < 0 ){
    u_min = 0;
  }

  if( v_min < 0 ){
    v_min = 0;
  }

  if( w_min < 0 ){
    w_min = 0;
  }

  if( u_max > get_size() ){
    u_max = get_size();
  }

  if( v_max > get_size() ){
    v_max = get_size();
  }

  if( w_max > get_size() ){
    w_max = get_size();
  }

  double coords[3];
  double half_thick[3];

  coords[0] = 0.5*( x_max - x_min ) + x_min;
  coords[1] = 0.5*( y_max - y_min ) + y_min;
  coords[2] = 0.5*( z_max - z_min ) + z_min;
  half_thick[0] = coords[0] - x_min;
  half_thick[1] = coords[1] - y_min;
  half_thick[2] = coords[2] - z_min;

  int crd;
  bool in_box;
  vector<unsigned int> ret_vector;
  vector<Vertex>::iterator it;
  double vert_coord[3];

   //loop through all octree cells looking for the vertices that belong in the selection box

  for( unsigned int k=w_min; k<w_max; k++ ){
    for( unsigned int j=v_min; j<v_max; j++){
      for( unsigned int i=u_min; i<u_max; i++){

	for( it=(*octree)(i,j,k).begin(); it!=(*octree)(i,j,k).end(); it++ ){

	  vert_coord[0] = it->get_x();  // get the vertex position
	  vert_coord[1] = it->get_y();  // get the vertex position
	  vert_coord[2] = it->get_z();  // get the vertex position

	  // see it is in the box...
 	  in_box = true;

	  for(crd=0; crd<3; crd++){
	    if (fabs(coords[crd]-vert_coord[crd]) > half_thick[crd]){
	      in_box=false;
	      break;
	    }
	  }

 	  if (in_box){
 	    ret_vector.push_back( it->get_vertex_id() );
 	  }


	}


      }
    }
  }
  
  return ret_vector;
  
}

