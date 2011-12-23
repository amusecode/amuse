#include <cassert>
#include <vector>
#include "vector3.h"
#include "morton_key.h"

//#define USE_SSE
struct particle{
	typedef vector3<float> vec;
	typedef std::vector<particle>::iterator it;
  std::vector<int> ngb_list;

	morton_key <vec, float>  key; //  8   8
	int id;         //  4  12
	short nnb_gath; //  2  14
	short nnb_scat;  //  2  16
	vec pos;        // 12  28
	float h;        //  4  32

  particle() {}
	particle(int _id, vec &_pos, float _h) : 
		id(_id), nnb_gath(0), nnb_scat(0), pos(_pos), h(_h) {
#ifdef USE_SSE
			assert(sizeof(particle) == 32);
#endif
		}
	void keygen(float size){
		key = morton_key<vec, float> (pos, size);
	}
	int octkey(int rshift){
		return 7 & (key.val >> rshift);
	}
/* An optimization for the STL sort */
#ifdef USE_SSE
	const particle operator = (const particle &rhs){
		typedef float v4sf __attribute__ ((vector_size(16)));
		v4sf *lp =(v4sf *)this;
		v4sf *rp =(v4sf *)(&rhs);
		lp[0] = rp[0];
		lp[1] = rp[1];
		return *this;
	}
#endif
};

struct cmp_particle_key{
	bool operator () (const particle &a, const particle &b){
		return a.key.val < b.key.val;
	}
};
struct cmp_particle_id{
	bool operator () (const particle &a, const particle &b){
		return a.id < b.id;
	}
};

/* An optimization for the STL sort */
#ifdef USE_SSE
namespace std{
	template <> 
	inline void iter_swap <particle::it, particle::it> (particle::it a, particle::it b){
		typedef float v4sf __attribute__ ((vector_size(16)));
		v4sf *ap =(v4sf *)&(*a);
		v4sf *bp =(v4sf *)&(*b);
		v4sf tmpa0 = ap[0];
		v4sf tmpa1 = ap[1];
		v4sf tmpb0 = bp[0];
		v4sf tmpb1 = bp[1];
		ap[0] = tmpb0;
		ap[1] = tmpb1;
		bp[0] = tmpa0;
		bp[1] = tmpa1;
	}
}
#endif
