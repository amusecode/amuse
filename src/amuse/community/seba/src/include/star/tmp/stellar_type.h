       //=======================================================// _\|/_
      //  __  _____           ___                    ___       //   /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     // _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //   /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       // _\|/_
//=======================================================//   /|\ ~

/*
 *  stellar_type.h: 
 *
 *.............................................................................
 *    version 1:  Jan 1994   Simon F. Portegies Zwart
 *    version 2:  Aug 1996   Simon Portegies Zwart
 *
 *.............................................................................
 *     This file includes:
 *  1) definition of stellar types
 *
 *
 *.............................................................................
 */
#ifndef     _STELLAR_TYPE
#   define  _STELLAR_TYPE

enum stellar_type {NAS = 0, Main_Sequence = 1, Hertzsprung_Gap,
                   Sub_Giant, Horizontal_Branch, Giant, Super_Giant,
                   Helium_Star, Wolf_Rayet, White_Dwarf, Thorn_Zytkow,
                   Neutron_Star, Black_Hole, Brown_Dwarf, Disintegrated,
                   Double, no_of_stellar_type};

#endif
