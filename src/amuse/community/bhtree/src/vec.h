/*
 *  vector.h: 3D vector operations include file
 *.............................................................................
 *    version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class vector
 *.............................................................................
 */

// This is slightly modified version of vector header file
// taken from STARLAB
// -- J. Makino

#ifndef  STARLAB_VECTOR_H
#  define  STARLAB_VECTOR_H

#include "stdinc.h"

/*-----------------------------------------------------------------------------
 *  vec  --  a class for 3-dimensional vectors
 *-----------------------------------------------------------------------------
 */

const int ndim = 3;

class vec
{
private:
    
    real element[3];
    
public:
    
    //	Default: initialize to zero.
    
    vec(real c = 0)
    {element[0] = element[1] = element[2] = c;}
    
    vec(real x, real y, real z)
    {element[0] = x; element[1] = y; element[2] = z;}
    
    //  []: the return type is declared as a reference (&), so that it can be used
    //  on the left-hand side of an asignment, as well as on the right-hand side,
    //  i.e.  v[1] = 3.14  and  x = v[2]  are both allowed and work as expected.
    
    real & operator [] (int i)       {return element[i];}
    
    inline void print() {cout << element[0] << " " << element[1] << " "
			      << element[2] << "\n";}

    vec  readjust(){
	return vec(*this);
    }
	
//	Unary -
	
        vec operator - ()
	    {return vec(-element[0], -element[1], -element[2]);}
	
	//	Dot product.
	
        real operator * (const vec& b)
	    {return element[0]*b.element[0]
		  + element[1]*b.element[1]
		  + element[2]*b.element[2];}

//	Outer product.

        vec operator ^ (const vec &b)
	    {return vec(element[1]*b.element[2] - element[2]*b.element[1],
			   element[2]*b.element[0] - element[0]*b.element[2],
			   element[0]*b.element[1] - element[1]*b.element[0]);}

//	Vec +, -

        vec operator + (const vec &b)
	    {return vec(element[0]+b.element[0],
			   element[1]+b.element[1],
			   element[2]+b.element[2]);}
        vec operator - (const vec &b)
	    {return vec(element[0]-b.element[0],
		 	   element[1]-b.element[1],
			   element[2]-b.element[2]);}

        friend vec operator + (real, const vec & );
        friend vec operator + (const vec &, real);

//	Scalar *, /

        friend vec operator * (real, const vec & );
        friend vec operator * (const vec &, real);
        friend vec operator / (const vec &, real);

//	Vec +=, -=, *=, /=

        vec& operator += (const vec& b)
	    {element[0] += b.element[0];       
	     element[1] += b.element[1];
	     element[2] += b.element[2];
	     return *this;}

	vec& operator -= (const vec& b)
	    {element[0] -= b.element[0];
	     element[1] -= b.element[1];
	     element[2] -= b.element[2];
	     return *this;}

	vec& operator *= (const real b)
	    {element[0] *= b;
	     element[1] *= b;
	     element[2] *= b;
	     return *this;}

	vec& operator /= (const real b)
	    {element[0] /= b;
	     element[1] /= b;
	     element[2] /= b;
	     return *this;}

//      Input / Output

        friend ostream & operator << (ostream & , const vec & );

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

//inline real readjust_r(real x) {
//  return x;
//}

//vec  readjust(vec v) {
//  return vec(v);
//}

#endif
