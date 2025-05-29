#ifndef __HACS6_FORCE__
#define __HACS6_FORCE__

namespace hacs6
{
	struct Force 
	{
		dvec3 acc, jrk, snp, crk;
		Force() {};
    Force(const double val) : acc(val), jrk(val), snp(val), crk(val) {}
		~Force() {};

		Force(const dvec3 &_acc, const dvec3 &_jrk, const dvec3 &_snp, const dvec3 &_crk) :
		 	acc(_acc), jrk(_jrk), snp(_snp), crk(_crk) {}

    Force(const hacs4::Force &f) : acc(f.acc), jrk(f.jrk), snp(0.0), crk(0.0) {}

		Force operator+(const Force &rhs) const 
		{
			return Force(acc + rhs.acc, jrk + rhs.jrk, snp + rhs.snp, crk + rhs.crk);
		}
		Force operator-(const Force &rhs) const 
		{
			return Force(acc - rhs.acc, jrk - rhs.jrk, snp - rhs.snp, crk - rhs.crk);
		}
		friend Force operator << (const Force &f, const double h)
		{
			const double h2 = h *h * (1.0/2.0);
			const double h3 = h2*h * (1.0/3.0);
			return Force(f.acc - h*f.jrk + h2*f.snp - h3*f.crk, f.jrk - h*f.snp + h2*f.crk, f.snp - h*f.crk, f.crk);
		}
		friend Force operator >> (const Force &f, const double h)
		{
			const double h2 = h *h * (1.0/2.0);
			const double h3 = h2*h * (1.0/3.0);
			return Force(f.acc + h*f.jrk + h2*f.snp + h3*f.crk, f.jrk + h*f.snp + h2*f.crk, f.snp + h*f.crk, f.crk);
		}
	};

	struct Interpolate 
	{
		dvec3 crk, pop, d5a;
		Interpolate() {};
    Interpolate(const double val) : crk(val), pop(val), d5a(val) {}
		~Interpolate() {};

		Interpolate(const dvec3 &_crk, const dvec3 &_pop, const dvec3 &_d5a) : crk(_crk), pop(_pop), d5a(_d5a) {}
		Interpolate(const Force &f0, const Force &f1, const double dt) 
		{
			const double h    = 0.5*dt;
			const double h2   = h*h;
			const double hinv = 2.0/dt;
			const double hinv2 = hinv *hinv;
			const double hinv3 = hinv *hinv2;
			const double hinv4 = hinv2*hinv2;
			const double hinv5 = hinv4*hinv;

			const dvec3 Am = (f1.acc - f0.acc);
			const dvec3 Ap = (f1.acc + f0.acc);
			const dvec3 Jm = (f1.jrk - f0.jrk)*h;
			const dvec3 Jp = (f1.jrk + f0.jrk)*h;
			const dvec3 Sm = (f1.snp - f0.snp)*h2;
			const dvec3 Sp = (f1.snp + f0.snp)*h2;

			// timestep
			crk = (  6.0/ 8.0)*hinv3 * (-5.0*Am + 5.0*Jp - Sm);
			pop = ( 24.0/16.0)*hinv4 * (        Sp -     Jm);
			d5a = (120.0/16.0)*hinv5 * ( 3.0*Am - 3.0*Jp + Sm);
    }

    Interpolate friend operator << (const Interpolate &ip, const double h) 
    {
      const double h2 = h*h * (1.0/2.0);
      return Interpolate(ip.crk - h*ip.pop + h2*ip.d5a, ip.pop - h*ip.d5a, ip.d5a);
    }
    Interpolate friend operator >> (const Interpolate &ip, const double h) 
    {
      const double h2 = h*h * (1.0/2.0);
      return Interpolate(ip.crk + h*ip.pop + h2*ip.d5a, ip.pop + h*ip.d5a, ip.d5a);
    }

  };

  inline double aarseth_step(
      const dvec3 &acc, const dvec3 &jrk,
      const dvec3 &snp, const dvec3 &crk, 
      const dvec3 &pop, const dvec3 &d5a,
      const double eta) 
  {
    const double s0 = acc.norm2();
    const double s1 = jrk.norm2();
    const double s2 = snp.norm2();
    const double s3 = crk.norm2();
    const double s4 = pop.norm2();
    const double s5 = d5a.norm2();

    const double h = std::sqrt(s0*s2) + s1;
    const double l = std::sqrt(s3*s5) + s4;
    assert(l > 0.0);
    return eta*std::pow(h/l, 1.0/6.0);
  }
}

#endif 
