#ifndef __HACS6_H__
#define __HACS6_H__

#include <map>
#include "irrf6.h"
#include "regf4.h"
#include "Scheduler.h"

#include "node.h"
#include "hacs4_force.h"
#include "hacs6_force.h"
#include "hacs64_particle.h"


namespace hacs64
{
  struct Nbody
  {
    int local_n;

    Particle::Vector ptcl;
    Particle::Vector ptcl2add;
    Predictor::Vector predicor;
    std::vector<int> ptcl2remove;
    std::vector< std::pair<int, Particle::action> >  ptcl2modify;
    std::map<int, int>         index2id_map;
    std::vector<irrf6::Interpolate> ip_irr;
    std::vector<int>           irr_list, reg_list;
    std::vector<NGBlist>       ngb_list;
    unsigned long long cyclical_idx;
    int nmax;
    double dtmax;

    std::vector<regf4::Force>  freg;
    std::vector<irrf6::Force>  firr;

    irrf6::irrf *irr_ptr;
    regf4::regf *reg_ptr;


    Scheduler scheduler;

    long long nsteps_reg, nsteps_irr, nsteps;
    long long ninter_reg, ninter_irr;
    double eta_reg, eta_irr;
    double t_global, dt_global;
    double eps2, h2max;
    int iteration;
    bool is_synched;

    double dt_reg, dt_irr, dt_all, dt_ngb;
    double dt_corr_irr, dt_corr_reg, dt_corr_mix;

    Nbody() : irr_ptr(NULL), reg_ptr(NULL)
    {
      nsteps_reg = nsteps_irr = nsteps = 0;
      ninter_reg = ninter_irr = 0;
      eta_reg = 0.14;
      eta_irr = 0.8;
      t_global = dt_global = 0;
      eps2 = 0;
      h2max = 0.5;
      iteration = 0;
      nmax = -1;
      dtmax = -1.0;
      cyclical_idx = 0;
      
      dt_all = dt_reg = dt_irr = dt_ngb = 0;
      dt_corr_irr = dt_corr_reg = dt_corr_mix = 0;
      is_synched = false;
    }
    ~Nbody() 
    {
      if (irr_ptr != NULL) delete irr_ptr;
      if (reg_ptr != NULL) delete reg_ptr;
    }
    bool is_sane()
    {
      assert(irr_ptr != NULL);
      assert(reg_ptr != NULL);
      assert((int)ptcl.size() <= nmax);
      assert(scheduler.dtmax == dtmax);
      return true;
    }

    void iterate()
    {

#if 0
      assert(scheduler.get_nreg() == local_n);
      assert(scheduler.get_nirr() == local_n);
#endif

      dt_global = scheduler.pull_active_list(irr_list, reg_list);
      t_global += dt_global;
      
      iterate_active();
    }

    void iterate_active()
    {
      double t0;
      const double t0_all = get_wtime();

      irrf6::irrf &irr = *irr_ptr;
      regf4::regf &reg = *reg_ptr;

      iteration++;

      t0 = get_wtime();
      reg.set_ti(t_global);
      reg.force_first(reg_list, freg, eps2);
      dt_reg += get_wtime() - t0;

      t0 = get_wtime();
      irr.set_ti(t_global);
      irr.force_first(irr_list, firr, eps2);
      irr.force_last();
      dt_irr += get_wtime() - t0;

      t0 = get_wtime();
      const int nirr = irr_list.size();
      int nsteps_irr_th = 0, nsteps_th = 0, ninter_irr_th = 0;
#pragma omp parallel for reduction(+:nsteps_irr_th,nsteps_th, ninter_irr_th)
      for (int ix = 0; ix < nirr; ix++) 
      {
        const int i = irr_list[ix];
        if (ptcl[i].t_irr == t_global)
        {
          ptcl[i].irr_rung = -1-ptcl[i].irr_rung;
          continue;
        }

        if (ngb_list[i].size() > 0) 
        {
          const hacs6::Force firr_i(firr[ix].acc, firr[ix].jrk, firr[ix].snp, dvec3(0.0));
          const double dt_irr = ptcl[i].correct_irr(
              firr_i,
              scheduler.get_dt_corr(ptcl[i].irr_rung),
              scheduler.get_dt_pred(ptcl[i].reg_rung), 
              eta_irr);
          ptcl[i].irr_rung = -1-scheduler.get_rung(dt_irr);

          nsteps_irr_th++;
          nsteps_th++;
          ninter_irr_th += ngb_list[i].size();
        } 
        else 
        { // no neighbour
          ptcl[i].correct_irr_zero(scheduler.get_dt_corr(ptcl[i].irr_rung));
          ptcl[i].irr_rung = -1-scheduler.get_rung(ptcl[i].reg_rung); 
        }
      }

      nsteps_irr += nsteps_irr_th;
      nsteps     += nsteps_th;
      ninter_irr += ninter_irr_th;

      dt_corr_irr += get_wtime() - t0;

      t0 = get_wtime();
      reg.force_last();
      dt_reg += get_wtime() - t0;

      t0 = get_wtime();
      const int nreg = reg_list.size();
      {
        NGBlist list_tmp;
        for (int ix = 0; ix < nreg; ix++)
        {
          if (ptcl[reg_list[ix]].t_reg == t_global)
            continue;

          reg.get_list(reg_list[ix],   list_tmp);
          irr.set_list(reg_list[ix],   list_tmp);
          ngb_list    [reg_list[ix]] = list_tmp;
        }
      }
      dt_ngb += get_wtime() - t0;

      t0 = get_wtime();
      irr.force_first      (reg_list,   firr, eps2);
      irr.interpolate_first(reg_list, ip_irr, eps2);
//      irr.force_last();
      irr.interpolate_last();
       dt_irr += get_wtime() - t0;

       is_synched = false;
       if (nreg == local_n)
       {
         int ngb_min = local_n;
         int ngb_max = 0;
         int ngb_mean = 0;
         std::vector<double> r2(local_n);
         for (int i = 0; i < local_n; i++)
         {
           ngb_min = std::min(ngb_min, (int)ngb_list[i].size());
           ngb_max = std::max(ngb_max, (int)ngb_list[i].size());
           ngb_mean += ngb_list[i].size();
           r2[i] = ptcl[i].pos.norm2();
         }
         fprintf(stderr, " *** ngb_min= %d  ngb_max= %d  ngb_mean= %g  *** r: ",
             ngb_min, ngb_max, 1.0*ngb_mean/local_n);
#if 0
         std::sort(r2.begin(), r2.end());
         for (int i = 0; i < 10; i++)
           fprintf(stderr, "%g ", std::sqrt(r2[local_n-1-i]));
#endif
         fprintf(stderr, "\n");
         is_synched = true;
       }

       t0 = get_wtime();
       int nsteps_reg_th = 0; 
       nsteps_th = 0; 
       int ninter_reg_th = 0;
#pragma omp parallel for reduction(+:nsteps_reg_th,nsteps_th, ninter_reg_th)
       for (int ix = 0; ix < nreg; ix++) 
       {
         const int i = reg_list[ix];
         if (ptcl[i].t_reg == t_global) 
         {
           ptcl[i].reg_rung = -1-ptcl[i].reg_rung;
           continue;
         }

         ptcl[i].h2  = freg[ix].h2;

         const hacs4::Force freg_i(freg[ix].acc, freg[ix].jrk);
         const hacs6::Force firr_i(firr[ix].acc, firr[ix].jrk, firr[ix].snp, ip_irr[ix].crk);

         const double dt_reg = ptcl[i].correct_reg(
             freg_i, firr_i,
             scheduler.get_dt_corr(ptcl[i].irr_rung), 
             scheduler.get_dt_corr(ptcl[i].reg_rung), 
             eta_reg);
         ptcl[i].reg_rung = -1-scheduler.get_rung(dt_reg);

         nsteps_reg_th++;
         nsteps_th++;
         ninter_reg_th += local_n; // - ngb_list[i].size();
       }

       nsteps_reg += nsteps_reg_th;
       nsteps     += nsteps_th;
       ninter_reg += ninter_reg_th;

       dt_corr_reg += get_wtime() - t0;

       t0 = get_wtime();
#pragma omp parallel for
       for (int ix = 0; ix < nreg; ix++)
       {
         const int i = reg_list[ix];
         if (ptcl[i].t_reg == t_global) 
         {
           if (ngb_list[i].size() != 0) 
           {
             if (-1-ptcl[i].reg_rung > -1-ptcl[i].irr_rung)
               ptcl[i].irr_rung = -1-((-1-ptcl[i].reg_rung)+1);
           }
           else
             ptcl[i].irr_rung = ptcl[i].reg_rung;
           continue;
         }


         if (ngb_list[i].size() != 0) 
         {
           const double dt_irr =	
             hacs6::aarseth_step(
                 ptcl[i].ftot.acc,
                 ptcl[i].firr.jrk, 
                 ptcl[i].firr.snp,
                 ip_irr[ix].crk, 
                 ip_irr[ix].pop, 
                 ip_irr[ix].d5a, 
                 eta_irr);


           assert(ptcl[i].irr_rung < 0);
           assert(ptcl[i].reg_rung < 0);
           ptcl[i].irr_rung = -1-scheduler.get_rung(dt_irr);
           if (-1-ptcl[i].reg_rung > -1-ptcl[i].irr_rung)
             ptcl[i].irr_rung = -1-((-1-ptcl[i].reg_rung)+1);
         }
         else
           ptcl[i].irr_rung = ptcl[i].reg_rung;
       }

       for (int ix = 0; ix < nreg; ix++) 
       {
         const int i = reg_list[ix];
         assert(ptcl[i].irr_rung < 0);
         assert(ptcl[i].reg_rung < 0);

         ptcl[i].irr_rung = -1-ptcl[i].irr_rung;
         ptcl[i].reg_rung = -1-ptcl[i].reg_rung;

         ptcl[i].irr_rung = scheduler.push_particle(i, ptcl[i].irr_rung, false);
         ptcl[i].reg_rung = scheduler.push_particle(i, ptcl[i].reg_rung, true);

         const Particle &pi = ptcl[i];
         reg.set_jp(i, regf4::Particle(
               pi.mass, t_global, pi.h2,
               pi.pos, pi.vel,
               pi.ftot.acc, pi.ftot.jrk));
         irr.set_jp(i, irrf6::Particle(
               pi.mass, t_global,
               pi.pos, pi.vel,
               pi.ftot.acc, pi.ftot.jrk, pi.ftot.snp, pi.ftot.crk));
         ptcl[i].t_reg = t_global;
       }

       for (int ix = 0; ix < nirr; ix++)
       {
         const int i = irr_list[ix];
         if (ptcl[i].irr_rung < 0)
         {
           ptcl[i].irr_rung = -1-ptcl[i].irr_rung;
           ptcl[i].irr_rung = scheduler.push_particle(i, ptcl[i].irr_rung, false);

           const Particle &pi = ptcl[i];
           reg.set_jp(i, regf4::Particle(
                 pi.mass, t_global, pi.h2,
                 pi.pos, pi.vel,
                 pi.ftot.acc, pi.ftot.jrk));
           irr.set_jp(i, irrf6::Particle(
                 pi.mass, t_global,
                 pi.pos, pi.vel,
                 pi.ftot.acc, pi.ftot.jrk, pi.ftot.snp, pi.ftot.crk));
         }
         ptcl[i].t_irr = t_global;
       }
       dt_corr_mix += get_wtime() - t0;

       const double t1_all = get_wtime();
       dt_all += t1_all - t0_all;
    }

    struct cmp_dist_ngb {
      bool operator() (const std::pair<int, float> &lhs, const std::pair<int, float> &rhs){
        return lhs.second < rhs.second;
      }
    };
    
    void initialize(
        const int   nbodies,
        const double mass[],
        const dvec3  pos [],
        const dvec3  vel [],
        const double eps2,
        const double hmax,
        const double eta_irr,
        const double eta_reg,
        const double dt_max)
    {
      local_n = nbodies;


      this->eps2    = eps2;
      this->h2max   = hmax*hmax;
      this->eta_irr = eta_irr;
      this->eta_reg = eta_reg;

      scheduler     = Scheduler(dt_max);

      irr_ptr = new irrf6::irrf(local_n);
      reg_ptr = new regf4::regf(local_n, h2max, scheduler.dt_tick);
      


      irrf6::irrf &irr = *irr_ptr;
      regf4::regf &reg = *reg_ptr;
      
      ptcl.    resize(local_n);
      ngb_list.resize(local_n);
      irr.resize(local_n);
      reg.resize(local_n);

      nsteps_reg = nsteps_irr = nsteps = 0;
      ninter_reg = ninter_irr = 0;
      t_global = dt_global = 0.0;
      iteration = 0;

      for (int i = 0; i < nbodies; i++) 
        ptcl[i] = Particle(mass[i], pos[i], vel[i]);

      std::vector<int> ilist(nbodies);
      for (int i = 0; i < nbodies; i++) ilist[i] = i;

#if 0
      for (int i = 0 ; i < nbodies; i++) 
      {
#if 0 
        std::vector< std::pair<int, float> > ingb_list(nbodies);

        for (int j = 0; j < nbodies; j++) 
          ingb_list[j] = std::pair<int, float>(j, (ptcl[j].pos - ptcl[i].pos).norm2());

        std::sort(ingb_list.begin(), ingb_list.end(), cmp_dist_ngb());
        ptcl[i].h2 = ingb_list[regf6::NGBMEAN].second;

        ngb_list[i].clear();
        for (int j = 0; j < regf6::NGBMEAN; j++)
          if (i != ingb_list[j].first)
            ngb_list[i].push_back(ingb_list[j].first);
#else
        ptcl[i].h2 = 0.0;
        ngb_list[i].clear();
#endif
      }
#else
      {
        fprintf(stderr,  "-- searching for neighbours --\n");

        // Initialize

        const int nbody = nbodies;

        node::allocate(nbody, nbody);
        std::vector<int> num_neib(nbody);
        typedef boundary<float>  Boundary;

        // Initialize
        Boundary bound;
        for(int i=0; i<nbody; i++)
        {
          particle::vec ipos(ptcl[i].pos);
          float h = 0.0f;;
          node::ptcl.push_back(particle(i, ipos, h));
          bound.merge(Boundary(ipos, 0.0f));
        }
        std::cerr << bound.max << std::endl;
        fvec3 vmax = maxeach(bound.min.abseach(), bound.max.abseach());
        float rmax = std::max(vmax.x, std::max(vmax.x, vmax.y));
        float scale = 1.f;
        while(scale < rmax) scale *= 2.f;

        // Setup tree
        //
        for(int i=0; i<nbody; i++)
          node::ptcl[i].keygen(scale);
        std::sort(node::ptcl.begin(), node::ptcl.end(), cmp_particle_key());

        node::node_heap.push_back(node());
        node &root = node::node_heap[0];
        for(int i=0; i<nbody; i++)
          root.push_particle(i, 60);
        root.make_boundary();

        /**** apply range-search ****/

        for (int i = 0; i < nbodies; i++)
        {
          particle::vec ipos(ptcl[i].pos);
          particle pi(i, ipos, 0.1);
          int iter = 0;
          while(1)
          {
            iter++;
            assert(iter < 1000);

            pi.ngb_list.clear();
            pi.nnb_gath = 0;
            pi << root;
            const int nngb = pi.ngb_list.size();
            assert(nngb == pi.nnb_gath);

            assert(nngb > 0);
            if (nngb > regf4::NGBMAX)
            {
              pi.h *= 0.75; //std::max(std::pow(regf4::NGBMEAN*1.0/(nngb + 1), 1.0/3.0), 0.75);
            }
            else if (nngb < regf4::NGBMIN)
            {
              pi.h *= 1.25; //std::min(std::pow(regf4::NGBMEAN*1.0/(nngb + 1), 1.0/3.0), 1.25);
            }
            else
              break;
          }

          const int nngb = pi.ngb_list.size();
          ngb_list[i].clear();
          for (int j = 0; j < nngb; j++)
            if (i != pi.ngb_list[j])
              ngb_list[i].push_back(pi.ngb_list[j]);
          ptcl[i].h2 = pi.h*pi.h;
          if (i == 0)
          {
            int nj = 0;
            for (int j = 0; j < nbodies; j++) 
              nj += (ptcl[j].pos - ptcl[i].pos).norm2() < ptcl[i].h2 ? 1 : 0;
            fprintf(stderr, "i= %d  nngb= %d nj= %d h= %g \n",
                i, nngb, nj, pi.h);
            //            assert(0);
          }
        }
        fprintf(stderr, " -- done -- \n");
      }
#endif

      for (int i = 0; i < nbodies; i++)
      {
        const Particle &pi = ptcl[i];
        irr.set_jp(i, irrf6::Particle(
              pi.mass, t_global,
              pi.pos, pi.vel,
              dvec3(0.0), dvec3(0.0), dvec3(0.0), dvec3(0.0)));
        reg.set_jp(i, regf4::Particle(
              pi.mass, t_global, pi.h2,
              pi.pos, pi.vel,
              dvec3(0.0), dvec3(0.0)));
        irr.set_list(i, ngb_list[i]);
        reg.set_list(i, ngb_list[i]);
      }

      reg.set_ti(t_global);
      irr.set_ti(t_global);

      reg.force_first(ilist, freg, eps2);
      reg.force_last();
      for (int i = 0; i < nbodies; i++)
      {
        reg.get_list(i, ngb_list[i]);
        irr.set_list(i, ngb_list[i]);
      }

      irr.force_first(ilist, firr, eps2);
      irr.force_last();

#pragma omp parallel for
      for (int i = 0; i < local_n; i++) 
      {
        ptcl[i].h2   = freg[i].h2;
        ptcl[i].freg = hacs4::Force(freg[i].acc, freg[i].jrk);
        ptcl[i].firr = hacs6::Force(firr[i].acc, firr[i].jrk, firr[i].snp, dvec3(0.0));
        ptcl[i].ftot = ptcl[i].firr + ptcl[i].freg;
      }

      {
        for (int i = 0; i < nbodies; i++)
        {
          const Particle &pi = ptcl[i];
          irr.set_jp(i, irrf6::Particle(
                pi.mass, t_global,
                pi.pos, pi.vel,
                pi.ftot.acc, pi.ftot.jrk, pi.ftot.snp, pi.ftot.crk));
          reg.set_jp(i, regf4::Particle(
                pi.mass, t_global, pi.h2,
                pi.pos, pi.vel,
                pi.ftot.acc, pi.ftot.jrk));
        }

        irr.force_first(ilist, firr, eps2);
        irr.force_last();

        for (int i = 0; i < local_n; i++) 
        {
          ptcl[i].h2   = freg[i].h2;
          ptcl[i].freg = hacs4::Force(freg[i].acc, freg[i].jrk);
          ptcl[i].firr = hacs6::Force(firr[i].acc, firr[i].jrk, firr[i].snp, dvec3(0.0));
          ptcl[i].ftot = ptcl[i].firr + ptcl[i].freg;
        }

        for (int i = 0; i < local_n; i++) 
        {
          const double fac = 0.03;
          double dt_irr, dt_reg;
          {
            const dvec3 acc = ptcl[i].freg.acc;
            const dvec3 jrk = ptcl[i].freg.jrk;
            const double s0 = acc.norm2();
            const double s1 = jrk.norm2();

            const double u = s1;
            const double l = s0;
            assert(l > 0.0);
            dt_irr = dt_reg = 0.001*eta_reg*std::sqrt(u/l);
          }
          if (ngb_list[i].size() > 0)
          {
            const dvec3 acc = ptcl[i].ftot.acc;
            const dvec3 jrk = ptcl[i].firr.jrk;
            const dvec3 snp = ptcl[i].firr.snp;
            const dvec3 crk = ptcl[i].firr.crk;
            const double s0 = acc.norm2();
            const double s1 = jrk.norm2();
            const double s2 = snp.norm2();
            const double s3 = crk.norm2();

            const double u = std::sqrt(s0*s2) + s1;
            const double l = std::sqrt(s1*s3) + s2;
            assert(l > 0.0);
            dt_irr = dt_reg = fac*eta_irr*std::sqrt(u/l);
          }

          ptcl[i].reg_rung = scheduler.push_particle(i, dt_reg, true);
          if (ngb_list[i].size() == 0)
            ptcl[i].irr_rung = scheduler.push_particle(i, dt_reg, false);
          else 
          {
            const int irr_rung = scheduler.get_rung(dt_irr);
            if (ptcl[i].reg_rung > irr_rung)
              ptcl[i].irr_rung = scheduler.push_particle(i, ptcl[i].reg_rung+1, false);
            else
              ptcl[i].irr_rung = scheduler.push_particle(i, irr_rung, false);
          }
        }
      }
#if 0
      scheduler.debug_dump();
#endif
    }

    double potential()
    {
      std::vector<double> gpot(local_n);
      reg_ptr->potential_first(gpot, eps2);
      reg_ptr->potential_last();
      double gpot_tot = 0.0;
      for (int i = 0;  i < local_n; i++)
      {
        ptcl[i].pot = gpot[i];
        gpot_tot   += ptcl[i].mass*gpot[i]*0.5;
      }
      return gpot_tot;
    }

    void write_dumbp(const char *filename)
    {

      FILE *fout; 
      if (!(fout = fopen(filename, "w"))) {
        std::cerr << "Cannot open file " << filename << std::endl;
        exit(-1);
      }

      const int ni = ptcl.size();
      for (int i = 0 ;i < ni; i++)
      {
        fprintf(fout, " %d  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  %.15g  \n",
            i,
            ptcl[i].mass,
            ptcl[i].pos.x,
            ptcl[i].pos.y,
            ptcl[i].pos.z,
            ptcl[i].vel.x,
            ptcl[i].vel.y,
            ptcl[i].vel.z);
      }

      fclose(fout);

    }

    /************ AMUSE INTERFACE ***************/

    void commit_parameters();
    void recommit_parameters();

    void commit_particles();
    void recommit_particles();

    void evolve_model(const double t_next)
    {
      while (t_global < t_next)
        iterate();
    };

    void __synchmodel();

    double get_total_mass()
    {
      const int ni = ptcl.size();
      double mass = 0;
      for (int i = 0; i < ni; i++)
        mass += ptcl[i].mass;
      return mass;
    }

    double get_epot()
    {
      return  potential();
    }

    double get_ekin()
    {
      const int ni = ptcl.size();
      double ekin = 0;
      for (int i = 0; i < ni; i++)
        ekin += ptcl[i].mass*0.5*ptcl[i].vel.norm2();
      return ekin;
    }

    dvec3 get_com_pos()
    {
      const int ni = ptcl.size();
      double mass = 0.0;
      dvec3  mpos(0.0);
      for (int i = 0; i < ni; i++)
      {
        mass += ptcl[i].mass;
        mpos += ptcl[i].mass*ptcl[i].pos;
      }
      return mpos*(1.0/mass);
    }

    dvec3 get_com_vel()
    {
      const int ni = ptcl.size();
      double mass = 0;
      dvec3  mvel(0.0);
      for (int i = 0; i < ni; i++)
      {
        mass += ptcl[i].mass;
        mvel += ptcl[i].mass*ptcl[i].vel;
      }
      return mvel*(1.0/mass);
    }

    double get_total_radius()
    {
      const int ni = ptcl.size();
      double r2max = 0.0;
      const dvec3 cpos = get_com_pos();
      for (int i = 0; i < ni; i++)
        r2max = std::max(r2max, (ptcl[i].pos - cpos).norm2());
      return std::sqrt(r2max);
    }



  };

  struct Energy
  {
    double ke, pe, e;
    dvec3  psum, lsum;
    Energy(const Particle::Vector &ptcl, const double eps2){
      const int n = ptcl.size();
      ke = pe = 0.0;
      psum = lsum = dvec3(0.0);
      for(int i=0; i<n; i++){
        const Particle &p = ptcl[i];
        ke += 0.5 * p.mass * p.vel.norm2();
        pe += 0.5 * p.mass * p.pot;
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




};

#endif // __HACS6_H__
