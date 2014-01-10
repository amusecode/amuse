#include "vector3.h"

#ifndef HUGE
#define HUGE HUGE_VAL
#endif
template <typename REAL>
struct boundary{
	typedef vector3<REAL> vec;
	vec min, max;
	boundary() : min(HUGE), max(-HUGE) {}
	boundary(const vec &_min, const vec &_max) : min(_min), max(_max) {}
	boundary(const vec &pos, const REAL &h = 0.0) : min(pos - vec(h)), max(pos + vec(h)) {}

	static const boundary merge(const boundary &a, const boundary &b){
		return boundary(mineach(a.min, b.min), maxeach(a.max, b.max));
	}
	void merge(const boundary &b){
		*this = merge(*this, b);
	}
	friend bool not_overlapped(const boundary &a, const boundary &b){
		return (a.max.x < b.min.x) || (b.max.x < a.min.x)
		    || (a.max.y < b.min.y) || (b.max.y < a.min.y)
		    || (a.max.z < b.min.z) || (b.max.z < a.min.z);
	}
	friend bool overlapped(const boundary &a, const boundary &b){
		return !not_overlapped(a, b);
	}
	const vec center() const {
		return REAL(0.5) * (max + min);
	}
	const vec hlen() const {
		return REAL(0.5) * (max - min);
	}
	const REAL separation2_from(const vec &pos) const{
		vec dr = center() - pos;
		dr = dr.abseach() - hlen();
		dr = vec::maxeach(dr, vec(0.0));
		return dr.norm2();
	}
};
#if 0
template <>
struct boundary<float>{
	typedef vector3<float> vec;
	typedef float v4sf __attribute__ ((vector_size(16)));

	vec min; 
	float p0;
	vec max;
	float p1;

	boundary() : min(HUGE), p0(0.f), max(-HUGE), p1(0.f) {}
	boundary(const vec &_min, const vec &_max) : min(_min), p0(0.f), max(_max), p1(0.f) {}
	boundary(const vec &pos, float h = 0.f) : 
		min(pos - vec(h)), p0(0.f), max(pos + vec(h)), p1(0.f) {}
	boundary(v4sf _min, v4sf _max){
		*(v4sf *)&min = _min;
		*(v4sf *)&max = _max;
	}

	static const boundary merge(const boundary &a, const boundary &b){
		return boundary(
				__builtin_ia32_minps(*(v4sf *)&a.min, *(v4sf *)&b.min),
				__builtin_ia32_maxps(*(v4sf *)&a.max, *(v4sf *)&b.max));
	}
	void merge(const boundary &b){
		*this = merge(*this, b);
	}
	friend bool not_overlapped(const boundary &a, const boundary &b){
		return __builtin_ia32_movmskps(
				(v4sf)(__builtin_ia32_cmpltps(
						*(v4sf *)&a.max, *(v4sf *)&b.min)))
		   ||  __builtin_ia32_movmskps(
				(v4sf)(__builtin_ia32_cmpltps(
						*(v4sf *)&b.max, *(v4sf *)&a.min)));
	}
	friend bool overlapped(const boundary &a, const boundary &b){
		return !not_overlapped(a, b);
	}
};
#endif
