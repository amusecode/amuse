#include "regf6.h"

namespace regf6
{
	int ni_max;
  float dt_tick;
	double h2max;
	double t_predictor, t_interaction;
	unsigned long long n_interaction;

	double t_global;
	std::vector<Particle > ptcl;
	std::vector<NGBlist  > list;
	struct Predictor
	{
		float mass;
		dvec3 pos;
    fvec3 vel, acc;
		Predictor() {}
		Predictor(const Particle &p, const double ti)
		{
#if 0
			const double dt  = ti - p.time;
			const double dt2 = dt*(1.0/2.0);
			const double dt3 = dt*(1.0/3.0);
			const double dt4 = dt*(1.0/4.0);
			const double dt5 = dt*(1.0/5.0);
			pos  = p.pos + dt*(p.vel + dt2*(p.acc + dt3*(p.jrk + dt4*(p.snp + dt5*p.crk))));
			vel  = p.vel + dt*(p.acc + dt2*(p.jrk + dt3*(p.snp + dt4* p.crk)));
			acc  = p.acc + dt*(p.jrk + dt2*(p.snp + dt3* p.crk));
#else
			const float dt  = dt_tick*((unsigned int)(ti/(double)dt_tick) - (unsigned int)(p.time/(double)dt_tick));
			const float dt2 = dt*(1.0/2.0);
			const float dt3 = dt*(1.0/3.0);
			const float dt4 = dt*(1.0/4.0);
			const float dt5 = dt*(1.0/5.0);
      pos = p.pos;
      vel = p.vel;
      acc = p.acc;
      fvec3 jrk = p.jrk;
      fvec3 snp = p.snp;
      fvec3 crk = p.crk;
			pos  = pos + dt*(vel + dt2*(acc + dt3*(jrk + dt4*(snp + dt5*crk))));
			vel  = vel + dt*(acc + dt2*(jrk + dt3*(snp + dt4* crk)));
			acc  = acc + dt*(jrk + dt2*(snp + dt3* crk));
#endif
			mass = p.mass;
		}
	};
	std::vector<Predictor> pred;

	regf::regf(const int _ni_max, const double _h2max)
	{
		ni_max = _ni_max;
		h2max  = _h2max;
		ptcl.resize(ni_max);
		pred.resize(ni_max);
		list.resize(ni_max);

		t_predictor = t_interaction = 0.0;
		n_interaction = 0;
	}
	regf::regf(const int _ni_max, const double _h2max, const double _dt_tick)
	{
		ni_max = _ni_max;
		h2max  = _h2max;
    dt_tick = _dt_tick;
		ptcl.resize(ni_max);
		pred.resize(ni_max);
		list.resize(ni_max);

		t_predictor = t_interaction = 0.0;
		n_interaction = 0;
	}
	regf::~regf() {}

	int regf::resize(const int ni)
	{
		assert(ni < ni_max);
		ptcl.resize(ni);
		pred.resize(ni);
		list.resize(ni);
		return 0;
	}

	int regf::set_ti(const double ti)
	{
		t_global = ti;
		return 0;
	}

	int regf::set_jp(const int iaddr, const Particle &pi)
	{
		assert(iaddr < ni_max);
		ptcl[iaddr] = pi;
		return 0;
	}
	int regf::set_jp(const std::vector<int> &ilist, const std::vector<Particle> &ptcl_list)
	{
		for (int i = 0; i < (const int)ilist.size(); i++)
		{
			assert(ilist[i] < ni_max);
			ptcl[ilist[i]] = ptcl_list[i];
		}
		return 0;
	}

	int regf::set_list(const int iaddr, const NGBlist &ngb)
	{
		assert(iaddr < ni_max);
		list[iaddr] = ngb;
		return 0;
	}
	int regf::set_list(const std::vector<int> &ilist, const std::vector<NGBlist> &ngb_list)
	{
		for (int i = 0; i < (const int)ilist.size(); i++)
		{
			assert(ilist[i] < ni_max);
			list[ilist[i]] = ngb_list[i];
		}
		return 0;
	}

	int regf::get_list(const int iaddr, NGBlist &ngb) 
	{
		assert(iaddr < ni_max);
		ngb = list[iaddr];
		return 0;
	}
	int regf::get_list(const std::vector<int>&ilist, std::vector<NGBlist> &ngb_list)
	{
		ngb_list.resize(ilist.size());
		for (int i = 0; i < (const int)ilist.size(); i++)
		{
			assert(ilist[i] < ni_max);
			ngb_list[i] = list[ilist[i]];
		}
		return 0;
	}

	std::vector<int  >  iptcl_list;
	std::vector<Force> *force_result;
	double eps2_force;

	int regf::force_first(const std::vector<int> &ilist, std::vector<Force> &force, const double eps2_in)
	{
		iptcl_list = ilist;
		force_result = &force;
		eps2_force = eps2_in;

		// predict_all

		const int nmax = ptcl.size();
		const int jbeg = 0;
		const int jend = nmax;
		assert(jbeg >= 0);
		assert(jend <= nmax);
		const double t0 = get_wtime();
#pragma omp parallel for
		for (int j = jbeg; j < jend; j++)
			pred[j] = Predictor(ptcl[j], t_global);
		const double t1 = get_wtime();
		t_predictor += t1 - t0;

		return 0;
	}

	int regf::force_last()
	{
		const std::vector<int> &ilist = iptcl_list;
		std::vector<Force> &force = *force_result;
		const double eps2 = eps2_force;

		const int ni = ilist.size();
		const int nj = pred.size();
		force.resize(ni);

		const double t0 = get_wtime();
#pragma omp parallel 
		for (int ix = 0; ix < ni; ix++)
		{
			const int i = ilist[ix];
			const Predictor &pi = pred[i];
			const double dt_reg = t_global - ptcl[i].time;

			const int ngbi = list[i].size();
			double h2 = ptcl[i].h2;
#if 0
			if (ngbi > NGBMAX || ngbi < NGBMIN) 
				h2 *= std::pow(NGBMEAN*1.0/(ngbi + 1), 2.0/3.0); 
			h2 = std::min(h2, h2max);
#endif

			Force &fi = force[ix];

			const int jiter_max = 10;
			int jiter = 0;
			for (jiter = 0; jiter < jiter_max; jiter++) 
			{
				fi.acc = fi.jrk = fi.snp = dvec3(0.0);	
				list[i].clear();

				for (int j = 0; j < nj; j++) 
				{
					const Predictor &pj = pred[j];
					const dvec3 dr = pj.pos - pi.pos;
					const fvec3 dv = pj.vel - pi.vel;
					const fvec3 da = pj.acc - pi.acc;
					double r2 = dr.norm2();

//					if (std::min(r2, std::min((dr+dv*dt_reg).norm2(),(dr + dv*dt_reg + da*dt_reg*dt_reg*0.5).norm2())) < h2) 
					if (r2 < h2)
					{
						if ((list[i].size() < NGBlist::NGB_MAX) && (i != j))
							list[i].push_back(j);
						if (list[i].size() >= NGBlist::NGB_MAX)
							break;
						continue;
					}
          if (r2 == 0.0) continue;
					assert(r2 > 0.0);
					r2 += eps2;	
					const float rv = dr * dv;
					const float v2 = dv * dv;
					const float ra = dr * da;

					const double rinv1  = 1.0/std::sqrt(r2);
					const double rinv2  = rinv1*rinv1;	
					const double mrinv3 = pj.mass * (rinv1 * rinv2);

					const float alpha  = rv * rinv2;
					const float alphalpha = alpha * alpha;
					const float beta   = ((v2 + ra) * rinv2 + alphalpha);

					const dvec3 Aij = mrinv3*dr;
					const dvec3 Jij = mrinv3*dv - (3.0*alpha)*Aij;
					const dvec3 Sij = mrinv3*da - (6.0*alpha)*Jij - (3.0*beta)*Aij;
					fi.acc += Aij;
					fi.jrk += Jij;
					fi.snp += Sij;
				}
        
        double fac = (std::pow(NGBMEAN*1.0/(list[i].size()+1), 2.0/3.0) + 1)*0.5;
        if (fac > 1.0) fac = std::min(fac, 1.25);
        else           fac = std::max(fac, 1.0/1.25);
        h2 *= fac;

				if (list[i].size() >= NGBlist::NGB_MAX) 
				{
#if 0
					h2 *= std::pow(NGBMEAN*1.0/(list.size() + 1), 2.0/3.0);
					h2  = std::min(h2max, h2);
#endif
					fprintf(stderr, " ** WARNING **  new_ngbi >= NGBBUFF, i= %d jiter= %d < jiter_max= %d\n",
							i, jiter, jiter_max);
				}
				else
					break;
			} 
			assert(jiter < jiter_max);
			fi.h2 = h2;
		}
		n_interaction += ni*ptcl.size();
		const double t1 = get_wtime();
		t_interaction += t1 - t0;

		return 0;
	}

	std::vector<double> *gpot_result;
	double eps2_pot;

	int regf::potential_first(std::vector<double> &gpot, const double eps2_pot_in)
	{
		gpot_result = &gpot;
		eps2_pot    = eps2_pot_in;
		return 0;
	}

	int regf::potential_last()
	{
		std::vector<double> &gpot = *gpot_result;
		const double eps2 = eps2_pot;

		const int nbody = ptcl.size();
		gpot.resize(nbody);
#pragma omp parallel for
		for(int i = 0; i < nbody; i++)
		{
			double p = 0.0;
			for(int j = 0; j < nbody; j++)
			{
				if (j == i) continue;
				const dvec3 dr = ptcl[j].pos - ptcl[i].pos;
				const double rinv = 1.0 / std::sqrt(eps2 + dr*dr);
				p += ptcl[j].mass * rinv;
			}
			gpot[i] = -p;
		}
		return 0;
	}

}
