#ifndef __IRRF6_H__
#define __IRRF6_H__

#include "NGBlist.h"
#include "Timer.h"
#include "vector3.h"

namespace irrf6
{

	struct Particle
	{
		double mass, time;
		dvec3 pos, vel, acc, jrk, snp, crk;
		Particle() {}
		Particle(
				const double _mass,
				const double _time,
				const dvec3 &_pos,
				const dvec3 &_vel,
				const dvec3 &_acc,
				const dvec3 &_jrk,
				const dvec3 &_snp,
				const dvec3 &_crk) : 
			mass(_mass), time(_time), 
			pos(_pos), vel(_vel), acc(_acc), jrk(_jrk), snp(_snp), crk(_crk) {}
	};

	struct Force
	{
		dvec3 acc, jrk, snp;
	};

	struct Interpolate
	{
		dvec3 crk, pop, d5a;
	};

	struct irrf
	{
		public:

			irrf(const int _ni_max = 0);
			~irrf();
			
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
			int interpolate_first(std::vector<int>&,	std::vector<Interpolate>&,	const double);
			int interpolate_last();
	};
}

#endif // __IRRF6_H__
