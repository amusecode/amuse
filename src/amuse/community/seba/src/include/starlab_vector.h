
       //=======================================================//    _\|/_
      //  __  _____           ___                    ___       //      /|\ ~
     //  /      |      ^     |   \  |         ^     |   \     //          _\|/_
    //   \__    |     / \    |___/  |        / \    |___/    //            /|\ ~
   //       \   |    /___\   |  \   |       /___\   |   \   // _\|/_
  //     ___/   |   /     \  |   \  |____  /     \  |___/  //   /|\ ~
 //                                                       //            _\|/_
//=======================================================//              /|\ ~

/// @file starlab_vector.h  3D vector operations.
//
//  This file includes:
//  1) definition of class vec

#ifndef  STARLAB_vector_H
#  define  STARLAB_vector_H

#include "stdinc.h"

/// \a vec: A class for 3-dimensional vectors.

///    version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino<br>
///    version 2:  Jun 2003   Steve McMillan: vector now vec to avoid
///                           conflict with upstart STL vector type.

class vec
{
    private:

	/// 3-D coordinates.

        real element[3];

    public:

	/// Default: initialize to (0,0,0).  Note that vec(1) is (1,1,1).

        vec(real c = 0)
	    {element[0] = element[1] = element[2] = c;}

	/// Initialize to (x,y,z).

	vec(real x, real y, real z)
	    {element[0] = x; element[1] = y; element[2] = z;}

	/// Element access.

	/// []: Allow array-style access.  The return type is declared as a
	/// reference (&), so that it can be used on the left-hand side of an
	/// asignment, as well as on the right-hand side, i.e.v[1] = 3.14  and
	/// x = v[2]  are both allowed and work as expected.

        real & operator [] (int i)       {return element[i];}

	/// Print to cout.  Note that cout << vec alwo works.

	inline void print() {cout << element[0] << " " << element[1] << " "
				  << element[2] << "\n";}

	/// Unary -

        vec operator - ()
	    {return vec(-element[0], -element[1], -element[2]);}

	/// Dot product.

        real operator * (const vec& b)
	    {return element[0]*b.element[0]
		  + element[1]*b.element[1]
		  + element[2]*b.element[2];}

	/// Outer product.

        vec operator ^ (const vec &b)
	    {return vec(element[1]*b.element[2] - element[2]*b.element[1],
			element[2]*b.element[0] - element[0]*b.element[2],
			element[0]*b.element[1] - element[1]*b.element[0]);}

	/// Vector +

        vec operator + (const vec &b)
	    {return vec(element[0]+b.element[0],
			element[1]+b.element[1],
			element[2]+b.element[2]);}
	/// Vector -

        vec operator - (const vec &b)
	    {return vec(element[0]-b.element[0],
			element[1]-b.element[1],
			element[2]-b.element[2]);}

        friend vec operator + (real, const vec & );
        friend vec operator + (const vec &, real);

	/// Scalar premultiply.

        friend vec operator * (real, const vec & );

	/// Scalar postmultiply.

        friend vec operator * (const vec &, real);

	/// Scalar divide.

        friend vec operator / (const vec &, real);

//	Bool ==, !=: apparently can't do bit-by-bit compare because the
//	member data are private.  Hence these need to be friends in order
//	to apply these operators.  Almost never used, however...

	/// Logical comparison ==

	friend bool operator == (const vec &, const vec &);

	/// Logical comparison !=

	friend bool operator != (const vec &, const vec &);

//      Post operations (again...) -- the GCC3 december 2002 sage:
//      this is a some interesting way to resolve the ambiguity of
//      resolving  v*s vs. s*v and the confusion of having a v*v with
//      a non-explicit vec constructor from a scalar.....

	/// Real postmultiply.

	vec operator * (const real b)
	  {return vec(element[0]*b, element[1]*b, element[2]*b); }

	/// Real postadd.

	vec operator + (const real b)
	  {return vec(element[0]+b, element[1]+b, element[2]+b); }

	/// Vector +=

        vec& operator += (const vec& b)
	    {element[0] += b.element[0];       
	     element[1] += b.element[1];
	     element[2] += b.element[2];
	     return *this;}

	/// Vector -=

	vec& operator -= (const vec& b)
	    {element[0] -= b.element[0];
	     element[1] -= b.element[1];
	     element[2] -= b.element[2];
	     return *this;}

	/// Vector *=

	vec& operator *= (const real b)
	    {element[0] *= b;
	     element[1] *= b;
	     element[2] *= b;
	     return *this;}

	/// Vector /=

	vec& operator /= (const real b)
	    {element[0] /= b;
	     element[1] /= b;
	     element[2] /= b;
	     return *this;}

	/// Input

        friend ostream & operator << (ostream & , const vec & );

	/// Output.

	friend istream & operator >> (istream & , vec & );
};

inline  ostream & operator << (ostream & s, const vec & v)
	    {return s << v.element[0] << "  " << v.element[1]
		      << "  " << v.element[2];}

inline  istream & operator >> (istream & s, vec & v)
	    {s >> v.element[0] >> v.element[1] >> v.element[2];
	     return s;}

inline  real square(vec v) {return v*v;}
inline  real abs(vec v)    {return sqrt(v*v);}

// Another measure of vector magnitude; less work than abs():

inline  real abs1(vec v)   {return abs(v[0]) + abs(v[1]) + abs(v[2]);}

inline  vec operator + (real b, const vec & v)
	    {return vec(b+v.element[0],
			b+v.element[1],
			b+v.element[2]);}

inline  vec operator + (const vec & v, real b)
	    {return vec(b+v.element[0],
			b+v.element[1],
			b+v.element[2]);}

inline  vec operator * (real b, const vec & v)
	    {return vec(b*v.element[0],
			b*v.element[1],
			b*v.element[2]);}

inline  vec operator * (const vec & v, real b)
	    {return vec(b*v.element[0],
			b*v.element[1],
			b*v.element[2]);}

inline  vec operator / (const vec & v, real b)
	    {return vec(v.element[0]/b,
			v.element[1]/b,
			v.element[2]/b);}

inline  bool operator == (const vec& u, const vec& v)
	    {return (u.element[0] == v.element[0]
		     && u.element[1] == v.element[1]
		     && u.element[2] == v.element[2]);}

inline  bool operator != (const vec& u, const vec& v)
	    {return !(u==v);}

#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------\\ ~
//  |  the end of:  |         /|\         |  inc/starlab_vector.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================\\ ~
 
