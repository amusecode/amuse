/*
 *  vec.h: 3D vec operations include file
 *.............................................................................
 *    version 1:  Dec 1992   Piet Hut, Steve McMillan, Jun Makino
 *    version 2:
 *.............................................................................
 *     This file includes:
 *  1) definition of class vec
 *.............................................................................
 */

// This is slightly modified version of vec header file
// taken from STARLAB
// -- J. Makino

#ifndef  STARLAB_VEC_H
#  define  STARLAB_VEC_H

typedef double real;
//#include "stdinc.h"

#define ONE_THIRD (0.3333333333333333333333333333L)
#define ONE_SIXTH (0.1666666666666666666666666666L)
#define PR(x)  cerr << #x << " = " << x << " "
#define PRC(x) cerr << #x << " = " << x << ",  "
#define PRL(x) cerr << #x << " = " << x << "\n"

#ifdef PERIODIC
inline real readjust_r(real x)
{
    return fmod(x+2.5, 1.0)-0.5;
}
#else
inline real readjust_r(real x)
{
    return x;
}
#endif

/*-----------------------------------------------------------------------------
 *  vec  --  a class for 3-dimensional vecs
 *-----------------------------------------------------------------------------
 */

#define vec my3dvec

const int ndim = 3;

class vec
{
private:
    
    real element[3];
    
public:
    
    //	Default: initialize to zero.
    
    inline vec(real c = 0)
    {element[0] = element[1] = element[2] = c;}
    
    inline vec(real x, real y, real z)
    {element[0] = x; element[1] = y; element[2] = z;}
    
    inline vec(real x[3])
    {element[0] = x[0]; element[1] = x[1]; element[2] = x[2];}

    //  []: the return type is declared as a reference (&), so that it can be used
    //  on the left-hand side of an asignment, as well as on the right-hand side,
    //  i.e.  v[1] = 3.14  and  x = v[2]  are both allowed and work as expected.
    
    inline real & operator [] (int i)       {return element[i];}
    
    inline void print() {cout << element[0] << " " << element[1] << " "
			      << element[2] << "\n";}
#ifdef PERIODIC
    // PERIODIC basic reajustment menber function
    vec readjust(){
	return vec(readjust_r(element[0]),
		      readjust_r(element[1]),
		      readjust_r(element[2]));
    }
#else
    inline vec  readjust(){
	return vec(*this);
    }
#endif    
	
	
	
//	Unary -
	
        inline vec operator - ()
	    {return vec(-element[0], -element[1], -element[2]);}
	
	//	Dot product.
	
        inline real operator * (const vec& b)
	    {return element[0]*b.element[0]
		  + element[1]*b.element[1]
		  + element[2]*b.element[2];}

//	Outer product.

    inline vec operator ^ (const vec &b)
	    {return vec(element[1]*b.element[2] - element[2]*b.element[1],
			   element[2]*b.element[0] - element[0]*b.element[2],
			   element[0]*b.element[1] - element[1]*b.element[0]);}

//	Vec +, -

        inline vec operator + (const vec &b)
	    {return vec(element[0]+b.element[0],
			   element[1]+b.element[1],
			   element[2]+b.element[2]);}
        inline vec operator - (const vec &b)
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

        inline vec& operator += (const vec& b)
	    {element[0] += b.element[0];       
	     element[1] += b.element[1];
	     element[2] += b.element[2];
	     return *this;}

	inline vec& operator -= (const vec& b)
	    {element[0] -= b.element[0];
	     element[1] -= b.element[1];
	     element[2] -= b.element[2];
	     return *this;}

	inline vec& operator *= (const real b)
	    {element[0] *= b;
	     element[1] *= b;
	     element[2] *= b;
	     return *this;}

	inline vec& operator /= (const real b)
	    {register real binv = 1.0/b;
             element[0] *= binv;
	     element[1] *= binv;
	     element[2] *= binv;
	     return *this;}

//      Input / Output

        friend ostream & operator << (ostream & , const vec & );

	friend istream & operator >> (istream & , vec & );

//      Comare  08/09/02 M.Fujii

	friend class compare_vec;
};


class compare_vec{
  public:
    const int dim;
    compare_vec(int n): dim(n){}

    bool operator()(const vec&a, const vec&b) const {
        return b.element[dim] > a.element[dim];
    }
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


#endif

//=======================================================================//
//  +---------------+        _\|/_        +------------------------------//
//  |  the end of:  |         /|\         |  inc/vec.h
//  +---------------+                     +------------------------------//
//========================= STARLAB =====================================//
 
