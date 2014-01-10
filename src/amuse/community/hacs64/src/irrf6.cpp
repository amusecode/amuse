#include "irrf6.h"
#include <algorithm>

#ifndef HUGE
#define HUGE HUGE_VAL
#endif

namespace irrf6
{
	int ni_max;
	double t_predictor, t_interaction, t_interpolate;
	unsigned long long n_interaction, n_interpolate;

	double t_global;
	struct Predictor
	{
		double mass;
		dvec3 pos, vel, acc;
    dvec3 jrk, snp, crk;
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
			jrk  = p.jrk + dt*(p.snp + dt2* p.crk);
#if 0
			snp  = p.snp + dt* p.crk;
			crk  = p.crk;
#endif
			mass = p.mass;
		}
	};
	std::vector<Predictor> pred;
	std::vector<NGBlist  > list;
	std::vector<Particle > ptcl;
	std::vector<bool>      predicted_hash;
	std::vector<int>       predicted_list;

	irrf::irrf(const int _ni_max) 
	{
		ni_max = _ni_max;
		ptcl.resize(ni_max);
		pred.resize(ni_max);
		list.resize(ni_max);

		t_predictor = t_interaction = 0.0;
		n_interaction = 0;
	}
	irrf::~irrf() {}

	int irrf::resize(const int ni) 
	{
		assert(ni <= ni_max);
		ptcl.resize(ni);
		pred.resize(ni);
		list.resize(ni);
		return 0;
	}

	int irrf::set_ti(const double ti) 
	{
		t_global = ti;
		return 0;
	}

	int irrf::set_jp(const int iaddr, const Particle &pi)
	{
		assert(iaddr < ni_max);
		ptcl[iaddr] = pi;
		return 0;
	}
	int irrf::set_jp(const std::vector<int> &ilist, const std::vector<Particle> &ptcl_list)
	{
		for (int i = 0; i < (const int)ilist.size(); i++)
		{
			assert(ilist[i] < ni_max);
			ptcl[ilist[i]] = ptcl_list[i];
		}
		return 0;
	}

  int irrf::commit_changes()
  {
    return 0;
  }

	int irrf::set_list(const int iaddr, const NGBlist &ngb)
	{
		assert(iaddr < ni_max);
		list[iaddr] = ngb;
		return 0;
	}
	int irrf::set_list(const std::vector<int> &ilist, const std::vector<NGBlist> &ngb_list)
	{
		for (int i = 0; i < (const int)ilist.size(); i++)
		{
			assert(ilist[i] < ni_max);
			list[ilist[i]] = ngb_list[i];
		}
		return 0;
	}

	int irrf::get_list(const int iaddr, NGBlist &ngb) 
	{
		assert(iaddr < ni_max);
		ngb = list[iaddr];
		return 0;
	}
	int irrf::get_list(const std::vector<int>&ilist, std::vector<NGBlist> &ngb_list)
	{
		ngb_list.resize(ilist.size());
		for (int i = 0; i < (const int)ilist.size(); i++)
		{
			assert(ilist[i] < ni_max);
			ngb_list[i] = list[ilist[i]];
		}
		return 0;
	}

	std::vector<int  >  iptcl_force_list;
	std::vector<Force> *force_result;
	double eps2_force;

	int irrf::force_first(const std::vector<int> &ilist, std::vector<Force> &force, const double eps2_in)
	{
		iptcl_force_list = ilist;
		force_result = &force;
		eps2_force = eps2_in;

		const int nmax = ptcl.size();
		const int jbeg = 0;
		const int jend = nmax;
		assert(jbeg >= 0);
		assert(jend <= nmax);
		const double t0 = get_wtime();
    if (ilist.size() > 	0.01*nmax)
    {
#pragma omp parallel for
      for (int j = jbeg; j < jend; j++)
        pred[j] = Predictor(ptcl[j], t_global);
    }
    else
    {
      static std::vector<bool> flag(nmax, false);
      static std::vector<int> idx;
      const int ni = ilist.size();
      for (int ix = 0; ix < ni; ix++)
      {
        const int i = ilist[ix];
        if (!flag[i])
        {
          idx.push_back(i);
          flag[i] = true;
        }

        const int nj = list[i].size();
        for (int jx = 0; jx < nj; jx++) 
        {
          const int j = list[i][jx];
          if (!flag[j])
          {
            idx.push_back(j);
            flag[j] = true;
          }
        }
      }
//      std::sort(idx.begin(), idx.end());

      const int nidx = idx.size();
      for (int j = 0; j < nidx; j++)
      {
        const int i = idx[j];
        flag[i] = false;
        pred[i] = Predictor(ptcl[i], t_global);
      }
      idx.clear();
    }
		const double t1 = get_wtime();
		t_predictor += t1 - t0;
		return 0;
	}

	int irrf::force_last() 
	{
		const std::vector<int> &ilist = iptcl_force_list;
		std::vector<Force> &force = *force_result;
		const double eps2 = eps2_force;

		const int ni = ilist.size();
		force.resize(ni);
		const double t0 = get_wtime();
		unsigned long long nint = 0;
#pragma omp parallel for reduction(+: nint)
		for (int ix = 0; ix < ni; ix++)
		{
			Force &fi = force[ix];
      fi.jr2 = HUGE;
			const int i = ilist[ix];
			const int nj = list[i].size();
			const Predictor &pi = pred[i];

			fi.acc = fi.jrk = fi.snp = 0.0;
			
			for (int jx = 0; jx < nj; jx++) 
			{
				const int j = list[i][jx];

				const Predictor &pj = pred[j];
#ifdef __GNUC__
        __builtin_prefetch(&pred[list[i][jx+1]]);
#endif

        const dvec3 dr = pj.pos - pi.pos;
        const dvec3 dv = pj.vel - pi.vel;
        const dvec3 da = pj.acc - pi.acc;

        double r2 = dr * dr;
        assert(r2 > 0.0);
        if (r2 < fi.jr2)
        {
          fi.jr2 = r2;
          fi.jnb = j;
        }
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
    n_interaction += nint;
    const double t1 = get_wtime();
    t_interaction += t1 - t0;
    return 0;
  }

  std::vector<int> iptcl_interp_list;
  std::vector<Interpolate> *interpolate_result;
  double eps2_interpolate;

  int irrf::interpolate_first(std::vector<int> &ilist,	std::vector<Interpolate> &interpolate,	const double eps2)
  {
    iptcl_interp_list = ilist;
    interpolate_result = &interpolate;
    eps2_interpolate = eps2;
    return 0;
  }


  int irrf::interpolate_last()
  {
    const std::vector<int> &ilist = iptcl_interp_list;
    std::vector<Interpolate> &interpolate = *interpolate_result;
#if 0
    std::vector<Force> &force = *force_result;
#endif
    const double eps2 = eps2_interpolate;

    const int ni = ilist.size();
    interpolate.resize(ni);

    const double t0 = get_wtime();
    unsigned long long nint = 0;
#pragma omp parallel for reduction(+: nint)
    for (int ix = 0; ix < ni; ix++)
    {
      const int i = ilist[ix];
      Interpolate &interp = interpolate[ix];
      const int nj = list[i].size();

#if 0
      Force &fi = force[ix];
      fi.acc = fi.jrk = fi.snp = 0.0;
#endif
      interp.crk = interp.pop = interp.d5a = dvec3(0.0);

      const Predictor &pi = pred[i];

      for (int jx = 0; jx < nj; jx++) 
      {
        const int j = list[i][jx];

        const Predictor &pj = pred[j];
        const dvec3 dr = pj.pos - pi.pos;
        const dvec3 dv = pj.vel - pi.vel;
        const dvec3 da = pj.acc - pi.acc;
        const dvec3 dj = pj.jrk - pi.jrk;
        const dvec3 ds = pj.snp - pi.snp;
        const dvec3 dc = pj.crk - pi.crk;

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

#if 0
        fi.acc += acc;
        fi.jrk += jrk;
        fi.snp += snp;
#endif

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
    n_interpolate += nint;
    const double t1 = get_wtime();
    t_interpolate += t1 - t0;
    return 0;
  }
}
