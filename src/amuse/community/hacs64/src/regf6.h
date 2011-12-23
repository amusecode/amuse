#ifndef __REGF6_H__
#define __REGF6_H__

#include "NGBlist.h"
#include "Timer.h"
#include "vector3.h"

namespace regf6
{

	struct Particle
	{
#if 1
		double mass, time, h2;
		dvec3 pos, vel, acc, jrk, snp, crk;
#else
    double time;
		float mass, h2;
		dvec3 pos;
    fvec3 vel, acc, jrk, snp, crk;
#endif
		Particle() {}
		Particle(
				const double _mass,
				const double _time,
				const double _h2,
				const dvec3 &_pos,
				const dvec3 &_vel,
				const dvec3 &_acc,
				const dvec3 &_jrk,
				const dvec3 &_snp,
				const dvec3 &_crk) : 
			mass(_mass), time(_time), h2(_h2),
			pos(_pos), vel(_vel), acc(_acc), jrk(_jrk), snp(_snp), crk(_crk) {}
	};


	struct Force
	{
		float h2;
		dvec3 acc;
    fvec3 jrk, snp;
	};

#if  0
	enum {
		NGBMIN  = 48, 
		NGBMEAN = 64,
		NGBMAX  = 96,	
	};
#endif

#if 1
	enum {
		NGBMIN  = 16,
		NGBMEAN = 32,
		NGBMAX  = 48,	
	};
#endif

	struct regf
	{
		public:

			regf(const int ni_max = 0, const double h2max = 0.0);
			regf(const int ni_max, const double h2max, const double dt_tick);
			~regf();

			int resize(const int ni);
			int set_ti(const double ti);
			int set_jp(const int iaddr, const Particle &pi);
			int set_jp(const std::vector<int>&, const std::vector<Particle>&);
			int set_list(const int iaddr, const NGBlist &ilist);
			int set_list(const std::vector<int>&, const std::vector<NGBlist>&);
			int get_list(const int iaddr, NGBlist &list);
			int get_list(const std::vector<int>&, std::vector<NGBlist>&);
			int force_first(const std::vector<int>&, std::vector<Force>&, const double eps2);
			int force_last();
			int potential_first(std::vector<double> &gpot, const double eps2);
			int potential_last();
	};

}

#endif // __REGF6_H__


