#ifndef __ENERGY_H__
#define __ENERGY_H__

struct Energy
{
	double ke, pe, e;
	dvec3  psum, lsum;
	Energy(const Particle::Vector &ptcl, const double eps2){
		const int n = ptcl.size();
		std::vector<double> pot(n);
		calc_pot(ptcl, pot, eps2);
		ke = pe = 0.0;
		psum = lsum = dvec3(0.0);
		for(int i=0; i<n; i++){
			const Particle &p = ptcl[i];
			ke += 0.5 * p.mass * p.vel.norm2();
			pe += 0.5 * p.mass * pot[i];
			psum += p.mass * ptcl[i].vel;
			lsum += p.mass * (p.pos  % p.vel);
		}
		e = ke + pe;
	}
	void print(FILE *fp = stderr, const char *prefix = "##") const {
		fprintf(fp, "%s ke= %g pe= %g  e= %g\n", prefix, ke, pe, e);
	}
	void print_mom(FILE *fp = stderr, const char *prefix = "##") const {
		fprintf(fp, "%s P = [%g %g %g], L = [%g %g %g]\n", 
				prefix, psum.x, psum.y, psum.z, lsum.x, lsum.y, lsum.z);
	}

	Energy(const double _ke, const double _pe, const double _e,
			const dvec3 &_psum, const dvec3 &_lsum) : 
		ke(_ke), pe(_pe), e(_e), psum(_psum), lsum(_lsum) {}
	Energy operator- (const Energy &rhs) const {
		return Energy(ke-rhs.ke, pe-rhs.pe, e-rhs.e, psum-rhs.psum, lsum-rhs.lsum);
	}
};

#endif

