#ifndef __HACS4_FORCE__
#define __HACS4_FORCE__

namespace hacs4
{
	struct Force 
	{
		dvec3 acc, jrk;
		Force() {};
		~Force() {};

    Force(const double val) : acc(val), jrk(val) {}
		Force(const dvec3 &_acc, const dvec3 &_jrk) : acc(_acc), jrk(_jrk) {}

		Force operator+(const Force &rhs) const {
			return Force(acc + rhs.acc, jrk + rhs.jrk);
		}
		Force operator-(const Force &rhs) const {
			return Force(acc - rhs.acc, jrk - rhs.jrk);
		}
		friend Force operator << (const Force &f, const double h){
			return Force(f.acc - h*f.jrk, f.jrk);
		}
		friend Force operator >> (const Force &f, const double h){
			return Force(f.acc + h*f.jrk, f.jrk);
		}
	}; 

	struct Interpolate 
	{
		dvec3 snp, crk;
		Interpolate() {};
    Interpolate(const double val) : snp(val), crk(val) {}
		~Interpolate() {};

		Interpolate(const dvec3 &_snp, const dvec3 &_crk) : snp(_snp), crk(_crk) {}
		Interpolate(const Force &f0, const Force &f1, const double dt) {
			const double h    = 0.5*dt;
			const double hinv = 2.0/dt;

			const dvec3 Am = (f1.acc - f0.acc);
			const dvec3 Ap = (f1.acc + f0.acc);
			const dvec3 Jm = (f1.jrk - f0.jrk)*h;
			const dvec3 Jp = (f1.jrk + f0.jrk)*h;

			// timestep
			snp = (0.5 * hinv*hinv     ) *  Jm;
			crk = (1.5 * hinv*hinv*hinv) * (Jp - Am);

		}

		Interpolate friend operator << (const Interpolate &ip, const double h) {
			return Interpolate(ip.snp - h*ip.crk, ip.crk);
		}
		Interpolate friend operator >> (const Interpolate &ip, const double h) {
			return Interpolate(ip.snp + h*ip.crk, ip.crk);
		}
	};
  
  inline double aarseth_step(
      const dvec3 &acc, const dvec3 &jrk,
      const dvec3 &snp, const dvec3 &crk, 
      const double eta) 
  {
    const double s0 = acc.norm2();
    const double s1 = jrk.norm2();
    const double s2 = snp.norm2();
    const double s3 = crk.norm2();
    const double u = std::sqrt(s0*s2) + s1;
    const double l = std::sqrt(s1*s3) + s2;
    if(l>0.0) {eta*std::sqrt(u/l);}
    else {return 0.0;}
  }
}

#endif //
