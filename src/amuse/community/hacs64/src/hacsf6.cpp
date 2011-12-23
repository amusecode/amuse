#include "regf6.h"
#include "irrf6.h"

namespace hacsf6
{
	typedef regf6::Particle Particle;
	struct Predictor
	{
		double mass;
		dvec3 pos, vel, acc;
		Predictor() {}
		Predictor(const Particle &p, const double ti)
		{
			const double dt  = ti - p.time;
			const double dt2 = dt*(1.0/2.0);
			const double dt3 = dt*(1.0/3.0);
			const double dt4 = dt*(1.0/4.0);
			const double dt5 = dt*(1.0/5.0);
			pos  = p.pos + dt*(p.vel + dt2*(p.acc + dt3*(p.jrk + dt4*(p.snp + dt5*p.crk))));
			vel  = p.vel + dt*(p.acc + dt2*(p.jrk + dt3*(p.snp + dt4* p.crk)));
			acc  = p.acc + dt*(p.jrk + dt2*(p.snp + dt3* p.crk));
			mass = p.mass;
		}
	};


	bool predict_all = false;

	double t_global;
	int ni_max;
	double h2max;
	std::vector<Particle > ptcl;
	std::vector<Predictor> pred;
	std::vector<NGBlist  > list;

	double t_predictor, t_reg_force, t_irr_force, t_irr_interp;
	unsigned long long  n_reg_force, n_irr_force, n_irr_interp;

	int open(const int _ni_max, const double _h2max) 
	{
		ni_max = _ni_max; h2max = _h2max;

		ptcl.resize(ni_max);
		pred.resize(ni_max);
		list.resize(ni_max);

		t_predictor = t_reg_force = t_irr_force = t_irr_interp = 0.0;
		n_reg_force = n_irr_force = n_irr_interp = 0;
		return 0;
	}

	int close()
	{
		return 0;
	};

	int resize(const int ni)
	{
		assert(ni < ni_max);
		ptcl.resize(ni);
		pred.resize(ni);
		list.resize(ni);
		return 0;
	}

	int set_ti(const double ti)
	{
		t_global = ti;
		predict_all = false;
		return 0;
	}

	int set_jp(const int iaddr, const Particle &pi)
	{
		assert(iaddr < ni_max);
		ptcl[iaddr] = pi;
		predict_all = false;
		return 0;
	}
	int set_jp(const std::vector<int> &ilist, const std::vector<Particle> &ptcl_list)
	{
		for (int i = 0; i < (const int)ilist.size(); i++)
		{
			assert(ilist[i] < ni_max);
			ptcl[ilist[i]] = ptcl_list[i];
		}
		predict_all = false;
		return 0;
	}

	int set_list(const int iaddr, const NGBlist &ngb)
	{
		assert(iaddr < ni_max);
		list[iaddr] = ngb;
		return 0;
	}
	int set_list(const std::vector<int> &ilist, const std::vector<NGBlist> &ngb_list)
	{
		for (int i = 0; i < (const int)ilist.size(); i++)
		{
			assert(ilist[i] < ni_max);
			list[ilist[i]] = ngb_list[i];
		}
		return 0;
	}

	int get_list(const int iaddr, NGBlist &ngb) 
	{
		assert(iaddr < ni_max);
		ngb = list[iaddr];
		return 0;
	}
	int get_list(const std::vector<int>&ilist, std::vector<NGBlist> &ngb_list)
	{
		ngb_list.resize(ilist.size());
		for (int i = 0; i < (const int)ilist.size(); i++)
		{
			assert(ilist[i] < ni_max);
			ngb_list[i] = list[ilist[i]];
		}
		return 0;
	}

	std::vector<int  >   ilist_freg;
	std::vector<regf6::Force> *freg;
	double                eps2_freg;

	std::vector<int  >   ilist_firr;
	std::vector<irrf6::Force> *firr;
	double                eps2_firr;
	
	std::vector<int  >         ilist_finterp;
	std::vector<irrf6::Interpolate> *finterp;
	double                      eps2_finterp;

	
	int freg_first(const std::vector<int> &ilist, std::vector<regf6::Force> &force, const double eps2)
	{
		ilist_freg = ilist;
		freg       = &force;
		eps2_freg  = eps2;

		if (predict_all) return 0;

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

		predict_all = true;

		return 0;
	}
	
	int firr_first(const std::vector<int> &ilist, std::vector<irrf6::Force> &force, const double eps2)
	{
		ilist_firr = ilist;
		firr       = &force;
		eps2_firr  = eps2;

		assert(predict_all);
		return 0;
	}


	int interpolate_first(const std::vector<int> &ilist, std::vector<irrf6::Interpolate> &force, const double eps2)
	{
		ilist_finterp = ilist;
		finterp       = &force;
		eps2_finterp  = eps2;

		assert(predict_all);
		return 0;
	}

	//////////
	
	int freg_last()
	{
		const std::vector<int> &ilist    = ilist_freg;
		std::vector<regf6::Force> &force =      *freg;
		const double eps2                =  eps2_freg;

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
			if (ngbi > regf6::NGBMAX || ngbi < regf6::NGBMIN) 
				h2 *= std::pow(regf6::NGBMEAN*1.0/(ngbi + 1), 2.0/3.0); 
			h2 = std::min(h2, h2max);

			regf6::Force &fi = force[ix];

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
					const dvec3 dv = pj.vel - pi.vel;
					const dvec3 da = pj.acc - pi.acc;
					double r2 = dr.norm2();

					if (std::min(r2, std::min((dr+dv*dt_reg).norm2(),(dr + dv*dt_reg + da*dt_reg*dt_reg*0.5).norm2())) < h2) 
					{
						if ((list[i].size() < NGBlist::NGB_MAX) && (i != j))
							list[i].push_back(j);
						if (list[i].size() >= NGBlist::NGB_MAX)
							break;
						continue;
					}
					assert(r2 > 0.0);
					r2 += eps2;	
					const double rv = dr * dv;
					const double v2 = dv * dv;
					const double ra = dr * da;

					const double rinv1  = 1.0/std::sqrt(r2);
					const double rinv2  = rinv1*rinv1;	
					const double mrinv3 = pj.mass * (rinv1 * rinv2);

					const double alpha  = rv * rinv2;
					const double alphalpha = alpha * alpha;
					const double beta   = ((v2 + ra) * rinv2 + alphalpha);

					const dvec3 Aij = mrinv3*dr;
					const dvec3 Jij = mrinv3*dv - (3.0*alpha)*Aij;
					const dvec3 Sij = mrinv3*da - (6.0*alpha)*Jij - (3.0*beta)*Aij;
					fi.acc += Aij;
					fi.jrk += Jij;
					fi.snp += Sij;
				}

				if (list[i].size() >= NGBlist::NGB_MAX) 
				{
					h2 *= std::pow(regf6::NGBMEAN*1.0/(list.size() + 1), 2.0/3.0);
					h2  = std::min(h2max, h2);
					fprintf(stderr, " ** WARNING **  new_ngbi >= NGBBUFF, i= %d jiter= %d < jiter_max= %d\n",
							i, jiter, jiter_max);
				}
				else
					break;
			} 
			assert(jiter < jiter_max);
			fi.h2 = h2;
		}
		n_reg_force += ni*ptcl.size();
		const double t1 = get_wtime();
		t_reg_force += t1 - t0;

		return 0;
	}
	
	std::vector<double> *gpot_result;
	double eps2_pot;

	int potential_first(std::vector<double> &gpot, const double eps2_pot_in)
	{
		gpot_result = &gpot;
		eps2_pot    = eps2_pot_in;
		return 0;
	}

	int potential_last()
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

	int firr_last()	
	{
		const std::vector<int> &ilist    = ilist_firr;
		std::vector<irrf6::Force> &force =      *firr;
		const double eps2                =  eps2_firr;

		const int ni = ilist.size();
		force.resize(ni);

		const double t0 = get_wtime();
		unsigned long long nint = 0;
#pragma omp parallel for reduction(+: nint)
		for (int ix = 0; ix < ni; ix++)
		{
			irrf6::Force &fi = force[ix];
			const int i = ilist[ix];
			const int nj = list[i].size();
			const Predictor &pi = pred[i];

			fi.acc = fi.jrk = fi.snp = 0.0;

			for (int jx = 0; jx < nj; jx++) 
			{
				const int j = list[i][jx];

				const Predictor &pj = pred[j];

				const dvec3 dr = pj.pos - pi.pos;
				const dvec3 dv = pj.vel - pi.vel;
				const dvec3 da = pj.acc - pi.acc;

				double r2 = dr * dr;
				assert(r2 > 0.0);
				r2 += eps2;
				const double rv = dr * dv;
				const double v2 = dv * dv;
				const double ra = dr * da;

				const double rinv1  = 1.0/std::sqrt(r2);
				const double rinv2  = rinv1*rinv1;	
				const double mrinv3 = pj.mass * (rinv1 * rinv2);


				const double alpha  = rv * rinv2;
				const double alphalpha = alpha * alpha;
				const double beta   = ((v2 + ra) * rinv2 + alphalpha);

				const dvec3 acc = mrinv3*dr;
				const dvec3 jrk = mrinv3*dv - (3.0*alpha)*acc;
				const dvec3 snp = mrinv3*da - (6.0*alpha)*jrk - (3.0*beta)*acc;

				fi.acc += acc;
				fi.jrk += jrk;
				fi.snp += snp;
			}

			nint += list[i].size();
		}
		n_irr_force += nint;
		const double t1 = get_wtime();
		t_irr_force+= t1 - t0;
		return 0;
	}

	int interpolate_last()
	{
		const std::vector<int> &ilist                = ilist_finterp;
		std::vector<irrf6::Interpolate> &interpolate =      *finterp;
		const double eps2                            =  eps2_finterp;

		const int ni = ilist.size();
		interpolate.resize(ni);

		const double t0 = get_wtime();
		unsigned long long nint = 0;
#pragma omp parallel for reduction(+: nint)
		for (int ix = 0; ix < ni; ix++)
		{
			const int i = ilist[ix];
			irrf6::Interpolate &interp = interpolate[ix];
			const int nj = list[i].size();

			interp.crk = interp.pop = interp.d5a = dvec3(0.0);

			const Predictor &pi = pred[i];
			const Particle  &ip = ptcl[i];

			for (int jx = 0; jx < nj; jx++) 
			{
				const int j = list[i][jx];

				const Predictor &pj = pred[j];
				const Particle  &jp = ptcl[j];
				const dvec3 dr = pj.pos - pi.pos;
				const dvec3 dv = pj.vel - pi.vel;
				const dvec3 da = pj.acc - pi.acc;
				const dvec3 dj = jp.jrk - ip.jrk;
				const dvec3 ds = jp.snp - ip.snp;
				const dvec3 dc = jp.crk - ip.crk;

				double r2 = dr * dr;
				assert(r2 > 0.0);
				r2 += eps2;

				const double rv = dr * dv;
				const double v2 = dv * dv;
				const double ra = dr * da;
				const double va = dv * da;
				const double rj = dr * dj;

				const double rinv1  = 1.0/std::sqrt(r2);
				const double rinv2  = rinv1*rinv1;	
				const double mrinv3 = pj.mass * (rinv1 * rinv2);


				const double alpha  = rv * rinv2;
				const double alphalpha = alpha * alpha;
				const double beta   = ((v2 + ra) * rinv2 + alphalpha);
				const double gamma  = ((3.0*va + rj) * rinv2 + alpha * (3.0*beta - 4.0*alphalpha));

				const dvec3 acc = mrinv3*dr;
				const dvec3 jrk = mrinv3*dv - (3.0*alpha)*acc;
				const dvec3 snp = mrinv3*da - (6.0*alpha)*jrk - (3.0*beta)*acc;
				const dvec3 crk = mrinv3*dj - (9.0*alpha)*snp - (9.0*beta)*jrk - (3.0*gamma)*acc;

				const double a2 = da*da;
				const double vj = dv*dj;
				const double aj = da*dj;
				const double rs = dr*ds;
				const double vs = dr*ds;
				const double rc = dr*dc;

				const double i0 = -rinv2;
				const double i1 = 2.0*rv;
				const double i2 = v2 + ra;
				const double i3 = va + rj/3.0;
				const double i4 = a2/4.0 + vj/3.0  + rs/12.0;
				const double i5 = aj/6.0 + vs/12.0 + rc/60.0;

				const double w0 = mrinv3;
				const double w1 = i0/1.0 * (1.5*i1*w0);
				const double w2 = i0/2.0 * (3.0*i2*w0 + 2.5*i1*w1);
				const double w3 = i0/3.0 * (4.5*i3*w0 + 4.0*i2*w1 + 3.5*i1*w2);
				const double w4 = i0/4.0 * (6.0*i4*w0 + 5.5*i3*w1 + 5.0*i2*w2 + 4.5*i1*w3);
				const double w5 = i0/5.0 * (7.5*i5*w0 + 7.0*i4*w1 + 6.5*i3*w2 + 6.0*i2*w3 + 5.5*i1*w4);

				interp.crk += crk;
				interp.pop += ds*w0 + 4.0*dj*w1 + 12.0*da*w2 + 24.0*(dv*w3 + dr*w4);
				interp.d5a += dc*w0 + 5.0*ds*w1 + 20.0*dj*w2 + 60.0*da*w3 + 120.0*(dv*w4 + dr*w5);
			}

			nint += list[i].size();
		}
		n_irr_interp += nint;
		const double t1 = get_wtime();
		t_irr_interp += t1 - t0;
		return 0;
	}


}

namespace regf6
{
	regf::regf(const int ni_max, const double h2max)
	{
		assert(hacsf6::open(ni_max, h2max) == 0);
	}
	regf::~regf() {
		assert(hacsf6::close() == 0);
	};

	int regf::resize(const int ni)
	{
		return hacsf6::resize(ni);
	}
	int regf::set_ti(const double ti)
	{
		return hacsf6::set_ti(ti);
	}
	int regf::set_jp(const int iaddr, const Particle &pi)
	{
		return hacsf6::set_jp(iaddr, pi);
	}
	int regf::set_jp(const std::vector<int> &ilist, const std::vector<Particle> &ptcl_list)
	{
		return hacsf6::set_jp(ilist, ptcl_list);
	}

	int regf::set_list(const int iaddr, const NGBlist &ngb)
	{
		return hacsf6::set_list(iaddr, ngb);
	}
	int regf::set_list(const std::vector<int> &ilist, const std::vector<NGBlist> &ngb_list)
	{
		return hacsf6::set_list(ilist, ngb_list);
	}

	int regf::get_list(const int iaddr, NGBlist &ngb) 
	{
		return hacsf6::get_list(iaddr, ngb);
	}
	int regf::get_list(const std::vector<int>&ilist, std::vector<NGBlist> &ngb_list)
	{
		return hacsf6::get_list(ilist, ngb_list);
	}


	int regf::force_first(const std::vector<int> &ilist, std::vector<Force> &force, const double eps2)
	{
		return hacsf6::freg_first(ilist, force, eps2);
	}

	int regf::force_last()
	{
		return hacsf6::freg_last();
	}

	std::vector<double> *gpot_result;
	double eps2_pot;

	int regf::potential_first(std::vector<double> &gpot, const double eps2)
	{
		return hacsf6::potential_first(gpot, eps2);
	}

	int regf::potential_last()
	{
		return hacsf6::potential_last();
	}
};

namespace irrf6
{

	irrf::irrf(const int _ni_max) {}
	irrf::~irrf() {}

	int irrf::resize(const int ni) {return 0;}
	int irrf::set_ti(const double ti) {return 0;}
	int irrf::set_jp(const int iaddr, const Particle &pi) {return 0;}
	int irrf::set_jp(const std::vector<int> &ilist, const std::vector<Particle> &ptcl_list) {return 0;}
	int irrf::set_list(const int iaddr, const NGBlist &ngb) {return 0;}
	int irrf::set_list(const std::vector<int> &ilist, const std::vector<NGBlist> &ngb_list) {return 0;}
	int irrf::get_list(const int iaddr, NGBlist &ngb) {return 0;}
	int irrf::get_list(const std::vector<int>&ilist, std::vector<NGBlist> &ngb_list) {return 0;}

	int irrf::force_first(const std::vector<int> &ilist, std::vector<Force> &force, const double eps2)
	{
		return hacsf6::firr_first(ilist, force, eps2);
	}
	int irrf::force_last() 
	{
		return hacsf6::firr_last();
	}
	int irrf::interpolate_first(std::vector<int> &ilist,	std::vector<Interpolate> &interpolate,	const double eps2)
	{
		return hacsf6::interpolate_first(ilist, interpolate, eps2);
	}
	int irrf::interpolate_last()
	{
		return hacsf6::interpolate_last();
	}
}

