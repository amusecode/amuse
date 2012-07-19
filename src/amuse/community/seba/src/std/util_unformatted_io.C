
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

// util_unformatted_io.C:  functions for unformatted I/O.  May require lowering
//			   the optimization level with some compilers...

#include "starlab_vector.h"
#include "util_io.h"
#include "story.h"
#include <ctype.h>

#undef isalnum	    // hacks for Irix 6.5 <ctype.h> backward compatibility
#undef isspace

void write_unformatted_real( ostream & s, real v )
{
#if WORDS_BIGENDIAN
    s.write( (char *)&v, 8 );
#else
    unsigned long long lv = *(unsigned long long *)&v;
    lv = (lv>>32) | (lv<<32);
    lv = (lv&0x0000FFFF0000FFFFLL)<<16
       | (lv>>16)&0x0000FFFF0000FFFFLL;
    lv = (lv&0x00FF00FF00FF00FFLL)<<8
       | (lv>>8)&0x00FF00FF00FF00FFLL;
    s.write( (char *)&lv, 8 );
#endif
}

void write_unformatted32_real( ostream & s, real v )
{
    float f = v;
#if WORDS_BIGENDIAN
    s.write( (char *)&f, 4 );
#else
    unsigned int l = (*(unsigned int *)&f)>>16 | (*(unsigned int *)&f)<<16;
    l = (l&0x00FF00FF)<<8
      | (l>>8)&0x00FF00FF;
    s.write( (char *)&l, 4 );
#endif
}

void write_unformatted_vector( ostream & s, vec & v )
{
    write_unformatted_real( s, v[0] );
    write_unformatted_real( s, v[1] );
    write_unformatted_real( s, v[2] );
}

void write_unformatted32_vector( ostream & s, vec & v )
{
    write_unformatted32_real( s, v[0] );
    write_unformatted32_real( s, v[1] );
    write_unformatted32_real( s, v[2] );
}

real read_unformatted_real( istream & s )
{
#if WORDS_BIGENDIAN
    real r;
    s.read( (char *)&r, 8 );
    return r;
#else

    // Note from Steve (12/04): the old version of this, in which only
    // a single lv was used in place of all these, now seems to fail.
    // Using lv1, ..., lv4 seems to fix the problem.  Not clear why it
    // occurred -- presumably it afflicts similar operations in other
    // functions in this file.
	
    unsigned long long lv1;
    s.read( (char *)&lv1, 8 );
    unsigned long long lv2 = (lv1>>32) | (lv1<<32);
    unsigned long long lv3 = (lv2&0x0000FFFF0000FFFFLL)<<16
				| (lv2>>16)&0x0000FFFF0000FFFFLL;
    unsigned long long lv4 = (lv3&0x00FF00FF00FF00FFLL)<<8
				| (lv3>>8)&0x00FF00FF00FF00FFLL;
    return *(real *)&lv4;
#endif
}

real read_unformatted32_real( istream & s )
{
#if WORDS_BIGENDIAN
    float f;
    s.read( (char *)&f, 4 );
    return f;
#else
    unsigned int iv;
    s.read( (char *)&iv, 4 );
    iv = (iv>>16) | (iv<<16);
    iv = (iv&0x00FF00FF)<<8
	| (iv>>8)&0x00FF00FF;
    return *(float *)&iv;
#endif
}

void read_unformatted_vector( istream & s, vec & v )
{
    v[0] = read_unformatted_real( s );
    v[1] = read_unformatted_real( s );
    v[2] = read_unformatted_real( s );
}

void read_unformatted32_vector( istream & s, vec & v )
{
    v[0] = read_unformatted32_real( s );
    v[1] = read_unformatted32_real( s );
    v[2] = read_unformatted32_real( s );
}
