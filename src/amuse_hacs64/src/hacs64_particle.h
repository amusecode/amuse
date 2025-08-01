#ifndef __HACS64_PARTICLE__
#define __HACS64_PARTICLE__

namespace hacs64
{
  struct Particle
  {
    enum action {MASS, RADIUS, POS, VEL, ALL};
    typedef std::vector<Particle> Vector;
    int id;
    double mass;
    dvec3  pos, vel;
    float  h2;
    float  radius;
    hacs6::Force firr;
    hacs4::Force freg;
    hacs6::Force ftot;
    int    irr_rung; 
    int    reg_rung; 
    double pot;
    double t_reg;
    double t_irr;
    int    jnb;
    float  jr2;

    Particle() {}
    ~Particle() {}      
    Particle(
        const double _mass, 
        const dvec3 &_pos, 
        const dvec3 &_vel)
      : id(-1), mass(_mass), pos(_pos), vel(_vel), radius(0.0), firr(0.0), freg(0.0), ftot(0.0), irr_rung(-1), reg_rung(-1) 
    {
      t_reg = t_irr = 0.0;
    }
    Particle(
        const double _mass, 
        const double _radius,
        const dvec3 &_pos, 
        const dvec3 &_vel,
        const int   _id)
      : id(_id), mass(_mass), pos(_pos), vel(_vel), radius(_radius), firr(0.0), freg(0.0), ftot(0.0), irr_rung(-1), reg_rung(-1) 
    {
      t_reg = t_irr = 0.0;
    }

    double correct_irr(
        const hacs6::Force &firr_new, 
        const double dt_irr, 
        const double dt_reg,
        const double eta_irr) 
    {
      const hacs6::Force firr_old = firr;
      const hacs6::Interpolate ip(firr_old, firr_new, dt_irr);

      const hacs6::Interpolate ip0 = ip << (0.5*dt_irr);
      const hacs6::Interpolate ip1 = ip >> (0.5*dt_irr);

      const double dt  = dt_irr;
      const double dt2 = dt *dt * (1.0/2.0);
      const double dt3 = dt2*dt * (1.0/3.0);
      const double dt4 = dt3*dt * (1.0/4.0);
      const double dt5 = dt4*dt * (1.0/5.0);
      const double dt6 = dt5*dt * (1.0/6.0);
      const double dt7 = dt6*dt * (1.0/7.0);

      pos +=      vel*dt + ftot.acc*dt2 + ftot.jrk*dt3 + ftot.snp*dt4 + (ftot.crk - firr.crk  + ip0.crk)*dt5 + ip0.pop*dt6 + ip0.d5a*dt7;
      vel += ftot.acc*dt + ftot.jrk*dt2 + ftot.snp*dt3 + (ftot.crk - firr.crk + ip0.crk)*dt4  + ip0.pop *dt5 + ip0.d5a*dt6;

      firr     = firr_new;
      firr.crk = ip1.crk;

      ftot = firr + (freg >> dt_reg);		

      const double dt_irr_new = 
        hacs6::aarseth_step(ftot.acc, firr.jrk, firr.snp, ip1.crk, ip1.pop, ip1.d5a, eta_irr);
      return dt_irr_new;
    }

    void correct_irr_zero(const double dt_irr)   /* this corrector is called when n_ngb = 0 */
    {
      const double dt  = dt_irr;
      const double dt2 = dt *dt * (1.0/2.0);
      const double dt3 = dt2*dt * (1.0/3.0);
      const double dt4 = dt3*dt * (1.0/4.0);
      const double dt5 = dt4*dt * (1.0/5.0);

      pos +=      vel*dt + ftot.acc*dt2 + ftot.jrk*dt3 + ftot.snp*dt4 + ftot.crk*dt5;
      vel += ftot.acc*dt + ftot.jrk*dt2 + ftot.snp*dt3 + ftot.crk*dt4;
    }

    double correct_reg(
        const hacs4::Force &freg_new, 
        const hacs6::Force &firr_new, 
        const double dt_reg,
        const double eta_reg)
    {
      const hacs6::Force freg_old = (firr_new - firr) + freg_new;
      const hacs4::Force f0(freg    .acc, freg    .jrk);  /* f0 with old ngb */
      const hacs4::Force f1(freg_old.acc, freg_old.jrk);  /* f1 with old ngb */
      const hacs4::Interpolate ip(f0, f1, dt_reg);
      const hacs4::Interpolate ip0 = ip << (0.5*dt_reg);
      const hacs4::Interpolate ip1 = ip >> (0.5*dt_reg);
			
      const double dt  = dt_reg;
			const double dt2 = dt *dt * (1.0/2.0);
			const double dt3 = dt2*dt * (1.0/3.0);
			const double dt4 = dt3*dt * (1.0/4.0);
			const double dt5 = dt4*dt * (1.0/5.0);

			pos += dt4*ip0.snp + dt5*ip0.crk;
			vel += dt3*ip0.snp + dt4*ip0.crk;
			
      freg = hacs4::Force(freg_new.acc, freg_new.jrk);
			firr = firr_new;

			ftot = firr + freg;
			
      const double dt_reg_new = 
				hacs4::aarseth_step(freg.acc, freg.jrk, ip1.snp, ip1.crk, eta_reg);
      
      if(dt_reg_new > 0.0) {return dt_reg_new;}
      else {return  dt_reg;}
    }	
  }; // 64 words * sizeof(int)
	
  struct Predictor
	{
    typedef std::vector<Predictor> Vector;
		double mass;
		dvec3 pos, vel, acc;
    dvec3 jrk, snp, crk;
		Predictor() {}
		Predictor(const Particle &p, const double dt)
		{
			const double dt2 = dt*(1.0/2.0);
			const double dt3 = dt*(1.0/3.0);
			const double dt4 = dt*(1.0/4.0);
			const double dt5 = dt*(1.0/5.0);
      const dvec3 _acc = p.ftot.acc;
      const dvec3 _jrk = p.ftot.jrk;
      const dvec3 _snp = p.ftot.snp;
      const dvec3 _crk = p.ftot.crk;
			pos  = p.pos + dt*(p.vel + dt2*(_acc + dt3*(_jrk + dt4*(_snp + dt5*_crk))));
			vel  = p.vel + dt*(_acc + dt2*(_jrk + dt3*(_snp + dt4* _crk)));
			acc  = _acc + dt*(_jrk + dt2*(_snp + dt3* _crk));
			jrk  = _jrk + dt*(_snp + dt2* _crk);
			mass = p.mass;
		}
	};

}


#endif
