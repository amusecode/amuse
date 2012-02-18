
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/// @file hydrobase.h  Underlying class for hydro systems, on same level as node.
//
//  version 1:  Jan 1993   Piet Hut
//  version 2:
//
//  This file includes:
//  1) definition of class hydrobase

#ifndef  STARLAB_HYDROBASE_H
#  define  STARLAB_HYDROBASE_H

#include  "starlab_vector.h"
#include  "story.h"

/// \a hydrobase: The underlying class for all hydrodynamical systems.

class  hydrobase
{
    protected:

	story * hydro_story;

    public:

        hydrobase()  {hydro_story = mk_story_chapter();}

	virtual ~hydrobase()  {delete hydro_story;}
	    
	story * get_hydro_story()               {return hydro_story;}

	void  set_hydro_story(story * ss)       {hydro_story = ss;}

	virtual ostream & print_hydro_story(ostream&);
	virtual istream & scan_hydro_story(istream&);
};

typedef  hydrobase *(*hbpfp)();

inline  hydrobase * new_hydrobase()    {return  new hydrobase;}

#endif
 
