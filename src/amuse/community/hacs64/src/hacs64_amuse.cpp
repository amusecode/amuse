#include "hacs64.h"

namespace hacs64
{
  void Nbody::commit_parameters()
  {
    fprintf(stderr, " hacs64::commit_parameters() -- \n");
    assert(irr_ptr == NULL);
    assert(reg_ptr == NULL);

    fprintf(stderr, "   nmax= %d \n", nmax);
    fprintf(stderr, "   dtmax= %g \n", dtmax);
    fprintf(stderr, "   n_ngb= %d \n", regf4::NGBMEAN);
    fprintf(stderr, "   h2max= %g \n", h2max);
    fprintf(stderr, "   eta_irr= %g \n", eta_irr);
    fprintf(stderr, "   eta_reg= %g \n", eta_reg);
    fprintf(stderr, "   eps2= %g \n", eps2);

    assert(nmax   > 0  );
    assert(dtmax  > 0.0); 
    assert(h2max  > 0.0);


    scheduler     = Scheduler(dtmax);
    irr_ptr = new irrf6::irrf(nmax);
    reg_ptr = new regf4::regf(nmax, h2max, scheduler.dt_tick);
  }
    
  void Nbody::recommit_parameters()
  {
    fprintf(stderr, " hacs64::recommit_parameters() -- \n");
  }
    
  void Nbody::commit_particles()
  {
    fprintf(stderr, " hacs64::commit_particles() -- \n");

    assert((int)ptcl2add.size() <= nmax);
    assert(ptcl.empty());
    assert(ptcl2remove.empty());
    assert(ptcl2modify.empty());
    assert(index2id_map.empty());

    assert(is_sane());

    irrf6::irrf &irr = *irr_ptr;
    regf4::regf &reg = *reg_ptr;


    const int nbodies = ptcl2add.size();

    reg.resize(nbodies);
    irr.resize(nbodies);

    local_n = nbodies;

    ptcl.    resize(nbodies);
    ngb_list.resize(nbodies);

    for (int i = 0; i < nbodies; i++)
    {
      ptcl[i] = ptcl2add[i];
      index2id_map[ptcl[i].id] = i;
    }
    ptcl2add.clear();

    std::vector<int> ilist(nbodies);
    for (int i = 0; i < nbodies; i++) ilist[i] = i;

    {
      /*********** search neighbours ********/

      fprintf(stderr,  " hacs64::commit_particles() -- initialising neighbour lists --\n");

#if 0
      for (int i = 0 ; i < nbodies; i++) 
      {
#if 0          /* O(N^2) approach */
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

        /* using octree */

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
      }
#endif

      fprintf(stderr,  " hacs64::commit_particles() -- done initialising neighbour lists--\n");
    }

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

    is_synched = true;
  }
    
  void Nbody::recommit_particles()
  {
    fprintf(stderr, " hacs64::commit_particles() -- \n");
    assert(false);
  }
    
  void Nbody::__synchmodel()
  {
    if (!is_synched)
    {
      const int ni = ptcl.size();

      irr_list.clear();
      for (int i = 0; i < ni; i++)
        irr_list.push_back(i);
      reg_list = irr_list;


      for (int i = 0; i < ni; i++)
      {
        ptcl[i].reg_rung = t_global > ptcl[i].t_reg ? scheduler.get_rung(t_global - ptcl[i].t_reg) : ptcl[i].reg_rung;
        ptcl[i].irr_rung = t_global > ptcl[i].t_irr ? scheduler.get_rung(t_global - ptcl[i].t_irr) : ptcl[i].irr_rung;
      }

      scheduler.flush();
      iteration--;
      iterate_active();

      assert(is_synched);
    }
  }


}
