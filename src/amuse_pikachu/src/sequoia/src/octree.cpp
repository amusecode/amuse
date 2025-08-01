#include "octree.h"

/*********************************/
/*********************************/
/*********************************/

void octree::set_src_directory(string src_dir) {                                                                                                                                 
    this->src_directory = (char*)src_dir.c_str();                                                                                                                                
}   


int octree::getAllignmentOffset(int n)
{
  const int allignBoundary = 128*sizeof(uint); //Fermi,128 bytes 
  
  int offset = 0;
  //Compute the number of bytes  
  offset = n*sizeof(uint); 
  //Compute number of 256 byte blocks  
  offset = (offset / allignBoundary) + (((offset % allignBoundary) > 0) ? 1 : 0); 
  //Compute the number of bytes padded / offset 
  offset = (offset * allignBoundary) - n*sizeof(uint); 
  //Back to the actual number of elements
  offset = offset / sizeof(uint);   
  
  return offset;
}

int octree::getTextureAllignmentOffset(int n, int size)
{
    const int texBoundary = TEXTURE_BOUNDARY; //Fermi
  
    int textOffset = 0;
    //Compute the number of bytes  
    textOffset = n*size; 
    //Compute number of texBoundary byte blocks  
    textOffset = (textOffset / texBoundary) + (((textOffset % texBoundary) > 0) ? 1 : 0); 
    //Compute the number of bytes padded / offset 
    textOffset = (textOffset * texBoundary) - n*size; 
    //Back to the actual number of elements
    textOffset = textOffset / size; 
    
    return textOffset;
}   

/*********************************/
/*********************************/
/*********************************/

