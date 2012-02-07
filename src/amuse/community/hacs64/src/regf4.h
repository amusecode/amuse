#ifndef __REGF4_H__
#define __REGF4_H__

#include "NGBlist.h"
#include "Timer.h"
#include "vector3.h"

namespace regf4
{

	struct Particle
	{
		double mass, time, h2;
		dvec3 pos, vel, acc, jrk;
		Particle() {}
		Particle(
				const double _mass,
				const double _time,
				const double _h2,
				const dvec3 &_pos,
				const dvec3 &_vel,
				const dvec3 &_acc,
				const dvec3 &_jrk) :
			mass(_mass), time(_time), h2(_h2),
			pos(_pos), vel(_vel), acc(_acc), jrk(_jrk) {}
	};


	struct Force
	{
		double h2;
		dvec3 acc, jrk;
	};

#if 1
	enum {
		NGBMIN  = 48,
		NGBMEAN = 64,
		NGBMAX  = 96,
	};
#endif

#if 0
	enum {
		NGBMIN  = 16,
		NGBMEAN = 32,
		NGBMAX  = 48,
	};
#endif

#if 0
  enum {
		NGBMIN  = 8,
		NGBMEAN = 16,
		NGBMAX  = 24,
	};
#endif


	struct regf
	{
		regf(const int _ni_max = 0, const double _h2max = 0.0);
		regf(const int _ni_max, const double _h2max, const double dt_tick);
		~regf();

		int resize(const int ni);
		int set_ti(const double ti);
    int commit_changes();
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

#endif // __REGF4_H__


