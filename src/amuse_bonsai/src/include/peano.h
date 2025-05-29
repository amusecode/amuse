#ifndef _PEANO_H_
#define _PEANO_H_

#define BITS_PER_DIMENSION 20

namespace peano_hilbert {

  typedef unsigned long peanokey;
  
  struct peano_struct {
    peanokey key;
    int      idx;
  };
  
  
  
  static int quadrants[24][2][2][2] = {
    /* rotx=0, roty=0-3 */
    {{{0, 7}, {1, 6}}, {{3, 4}, {2, 5}}},
    {{{7, 4}, {6, 5}}, {{0, 3}, {1, 2}}},
    {{{4, 3}, {5, 2}}, {{7, 0}, {6, 1}}},
    {{{3, 0}, {2, 1}}, {{4, 7}, {5, 6}}},
    /* rotx=1, roty=0-3 */
    {{{1, 0}, {6, 7}}, {{2, 3}, {5, 4}}},
    {{{0, 3}, {7, 4}}, {{1, 2}, {6, 5}}},
    {{{3, 2}, {4, 5}}, {{0, 1}, {7, 6}}},
    {{{2, 1}, {5, 6}}, {{3, 0}, {4, 7}}},
    /* rotx=2, roty=0-3 */
    {{{6, 1}, {7, 0}}, {{5, 2}, {4, 3}}},
    {{{1, 2}, {0, 3}}, {{6, 5}, {7, 4}}},
    {{{2, 5}, {3, 4}}, {{1, 6}, {0, 7}}},
    {{{5, 6}, {4, 7}}, {{2, 1}, {3, 0}}},
    /* rotx=3, roty=0-3 */
    {{{7, 6}, {0, 1}}, {{4, 5}, {3, 2}}},
    {{{6, 5}, {1, 2}}, {{7, 4}, {0, 3}}},
    {{{5, 4}, {2, 3}}, {{6, 7}, {1, 0}}},
    {{{4, 7}, {3, 0}}, {{5, 6}, {2, 1}}},
    /* rotx=4, roty=0-3 */
    {{{6, 7}, {5, 4}}, {{1, 0}, {2, 3}}},
    {{{7, 0}, {4, 3}}, {{6, 1}, {5, 2}}},
    {{{0, 1}, {3, 2}}, {{7, 6}, {4, 5}}},
    {{{1, 6}, {2, 5}}, {{0, 7}, {3, 4}}},
    /* rotx=5, roty=0-3 */
    {{{2, 3}, {1, 0}}, {{5, 4}, {6, 7}}},
    {{{3, 4}, {0, 7}}, {{2, 5}, {1, 6}}},
    {{{4, 5}, {7, 6}}, {{3, 2}, {0, 1}}},
    {{{5, 2}, {6, 1}}, {{4, 3}, {7, 0}}}
  };
  
  static int rotxmap_table[24] = { 4, 5, 6, 7, 8, 9, 10, 11,
				   12, 13, 14, 15, 0, 1, 2, 3, 17, 18, 19, 16, 23, 20, 21, 22
  };

  static int rotymap_table[24] = { 1, 2, 3, 0, 16, 17, 18, 19,
				   11, 8, 9, 10, 22, 23, 20, 21, 14, 15, 12, 13, 4, 5, 6, 7
  };
  
  static int rotx_table[8]  = { 3, 0, 0, 2, 2, 0, 0, 1 };
  static int roty_table[8]  = { 0, 1, 1, 2, 2, 3, 3, 0 };
  
  static int sense_table[8] = { -1, -1, -1, +1, +1, -1, -1, -1 };
  
  // static int flag_quadrants_inverse = 1;                                                                                                                                                                                                           
  // static char quadrants_inverse_x[24][8];                                                                                                                                                                                                          
  // static char quadrants_inverse_y[24][8];                                                                                                                                                                                                          
  // static char quadrants_inverse_z[24][8];                                                                                                                                                                                                          
  
  /*! This function computes a Peano-Hilbert key for an integer triplet (x,y,z),                                                                                                                                                                      
   *  with x,y,z in the range between 0 and 2^bits-1.                                                                                                                                                                                                 
   */
  
  inline peanokey peano_hilbert_key(const int x, const int y, const int z, const int bits) {
    int i, quad, bitx, bity, bitz;
    int mask, rotation, rotx, roty, sense;
    peanokey key;
    
    mask = 1 << (bits - 1);
    key = 0;
    rotation = 0;
    sense = 1;
    
    for(i = 0; i < bits; i++, mask >>= 1)
      {
	bitx = (x & mask) ? 1 : 0;
	bity = (y & mask) ? 1 : 0;
	bitz = (z & mask) ? 1 : 0;
	
	quad = quadrants[rotation][bitx][bity][bitz];
	
	key <<= 3;
	key += (sense == 1) ? (quad) : (7 - quad);
	
	rotx = rotx_table[quad];
	roty = roty_table[quad];
	sense *= sense_table[quad];
	
	while(rotx > 0)
	  {
	    rotation = rotxmap_table[rotation];
	    rotx--;
	  }
	
	while(roty > 0)
	  {
	    rotation = rotymap_table[rotation];
	    roty--;
	  }
      }
    
    return key;
  }
  
  
  inline peanokey peano_hilbert_key_new(int x, int y, int z, const int bits) {
    int i, bitx, bity, bitz;
    int mask;
    peanokey key;
    
    //0= 000, 1=001, 2=011, 3=010, 4=110, 5=111, 6=101, 7=100
    //000=0=0, 001=1=1, 011=3=2, 010=2=3, 110=6=4, 111=7=5, 101=5=6, 100=4=7
    const int C[8] = {0, 1, 3, 2, 7, 6, 4, 5};
    
    int temp;
    
    mask = 1 << (bits - 1);
    key  = 0;
    

    for(i = 0; i < bits; i++, mask >>= 1)
      {
        bitx = (x & mask) ? 1 : 0;
        bity = (y & mask) ? 1 : 0;
        bitz = (z & mask) ? 1 : 0;
        
        int index = (bitx << 2) + (bity << 1) + (bitz);
        
        key  = (key << 3) + C[index];
        
        if(index == 0 || index == 4)  //Switch X and Z
        {
          temp = z; z = x; x = temp;
        }
        else if(index == 1) //switch Y and Z
        {
          temp = z; z = y; y = temp;
        }
        else if(index == 3 || index == 7) //switch X and Z, inverse X and Z
        {
          temp = x^(-1);
          x = z^(-1);
          z = temp;
        }
        else if(index == 6) //switch Y and Z, inverse Y and Z
        {
          temp  = y^(-1);
          y     = z^(-1);
          z     = temp;
        }
        //Index 2 and 5 dont do anything
      }
    
    return key;
  }  

 inline peanokey peano_hilbert_key_new_fix(int x, int y, int z, const int bits) {
    int i, bitx, bity, bitz;
    int mask;
    peanokey key;
    
    //0= 000, 1=001, 2=011, 3=010, 4=110, 5=111, 6=101, 7=100
    //000=0=0, 001=1=1, 011=3=2, 010=2=3, 110=6=4, 111=7=5, 101=5=6, 100=4=7
    const int C[8] = {0, 1, 3, 2, 7, 6, 4, 5};
    
    int temp;
    
    mask = 1 << (bits - 1);
    key  = 0;
    

    for(i = 0; i < bits; i++, mask >>= 1)
      {
        bitx = (x & mask) ? 1 : 0;
        bity = (y & mask) ? 1 : 0;
        bitz = (z & mask) ? 1 : 0;
        
        int index = (bitx << 2) + (bity << 1) + (bitz);
        
        key  = (key << 3) + C[index];
        
        if(index == 0 || index == 4)  //Switch X and Z
        {
          temp = z; z = x; x = temp;
        }
        else if(index == 1) //switch Y and X
        {
          temp = x; x = y; y = temp;
        }
        else if(index == 3 || index == 7) //switch X and Z, inverse X and Z
        {
          temp = x^(-1);
          x = z^(-1);
          z = temp;
        }
        else if(index == 6) //switch Y and X, inverse Y and X
        {
          temp  = y^(-1);
          y     = x^(-1);
          x     = temp;
        }
        //Index 2 and 5 dont do anything
      }
    
    return key;
  } 
  
 inline peanokey peano_hilbert_key_new_thesis(int x, int y, int z, const int bits) {
    int i,xi, yi, zi;
    int mask;
    peanokey key;
    
    //0= 000, 1=001, 2=011, 3=010, 4=110, 5=111, 6=101, 7=100
    //000=0=0, 001=1=1, 011=3=2, 010=2=3, 110=6=4, 111=7=5, 101=5=6, 100=4=7
    const int C[8] = {0, 1, 7, 6, 3, 2, 4, 5};
    
    int temp;
    
    mask = 1 << (bits - 1);
    key  = 0;
    

//     for(i = 0; i < bits; i++, mask >>= 1)
//       {
//         xi = ((x) >> i) & 1;
//         yi = ((y) >> i) & 1;
//         zi = ((z) >> i) & 1;
    for(i = 0; i < bits; i++, mask >>= 1)
      {
        xi = (x & mask) ? 1 : 0;
        yi = (y & mask) ? 1 : 0;
        zi = (z & mask) ? 1 : 0;        
               
        if(xi == 0 && yi == 0 && zi == 0)
        {
          temp = z; z = y; y = temp;
        }
        else  if(xi == 0 && yi == 0 && zi == 1)
        {
          temp = x; x = y; y = temp;
        }
        else  if(xi == 1 && yi == 0 && zi == 1)
        {
          temp = x; x = y; y = temp;
        }
        else  if(xi == 1 && yi == 0 && zi == 0)
        {
         x = (x) ^ (-1);
         z = (z) ^ (-1);
        }
        else  if(xi == 1 && yi == 1 && zi == 0)
        {
         x = (x) ^ (-1);
         z = (z) ^ (-1);
        }
        else  if(xi == 1 && yi == 1 && zi == 1)
        {
         temp = (x) ^ (-1);         
         x = (y) ^ (-1);
         y = temp;
        }
        else  if(xi == 0 && yi == 1 && zi == 1)
        {
         temp = (x) ^ (-1);         
         x = (y) ^ (-1);
         y = temp;
        }
        else
        {
         temp = (z) ^ (-1);         
         z = (y) ^ (-1);
         y = temp;          
        }
        
        int index = (xi << 2) + (yi << 1) + zi;
        key = (key << 3) + C[index];
      }
    
    return key;
  }   
  

inline peanokey peano_hilbert_key_new_thesis2(int x, int y, int z, const int bits) {
    int i,xi, yi, zi;
    int mask;
    peanokey key;
    
    //0= 000, 1=001, 2=011, 3=010, 4=110, 5=111, 6=101, 7=100
    //000=0=0, 001=1=1, 011=3=2, 010=2=3, 110=6=4, 111=7=5, 101=5=6, 100=4=7
    const int C[8] = {0, 1, 7, 6, 3, 2, 4, 5};
    
    int temp;
    
    mask = 1 << (bits - 1);
    key  = 0;
    
    
//     fprintf(stderr, "TEST MASK: %d\t%X\t%d  \n", mask, mask, bits);
//     for(i = 0; i < bits; i++, mask >>= 1)
//       {
//         xi = ((x) >> i) & 1;
//         yi = ((y) >> i) & 1;
//         zi = ((z) >> i) & 1;
    for(i = 0; i < bits; i++, mask >>= 1)
    {
        xi = (x & mask) ? 1 : 0;
        yi = (y & mask) ? 1 : 0;
        zi = (z & mask) ? 1 : 0;        

        int index = (xi << 2) + (yi << 1) + zi;
        
//         fprintf(stderr, "TEST MASK2: %d\t%X\t%d  \n", mask, mask, bits);    
//         fprintf(stderr, "TEST MASK3: %d %d %d \t %d %d %d (%X %X %X) \t i: %d  \n", xi, yi, zi, x, y, z,x, y, z, i);    

        
        if(index == 0)
        {
          temp = z; z = y; y = temp;
        }
        else  if(index == 1 || index == 5)
        {
          temp = x; x = y; y = temp;
        }
        else  if(index == 4 || index == 6)
        {
         x = (x) ^ (-1);
         z = (z) ^ (-1);
        }
        else  if(index == 7 || index == 3)
        {
         temp = (x) ^ (-1);         
         x = (y) ^ (-1);
         y = temp;
        }
        else
        {
         temp = (z) ^ (-1);         
         z = (y) ^ (-1);
         y = temp;          
        }
        

        key = (key << 3) + C[index];
      }
//     exit(0);
    return key;
  }   
  
};


#endif // _PEANO_H_
