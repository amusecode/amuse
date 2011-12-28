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



    const int nbodies = ptcl2add.size();


    local_n = nbodies;

    ptcl.    resize(nbodies);

    for (int i = 0; i < nbodies; i++)
    {
      ptcl[i] = ptcl2add[i];
      assert(index2id_map.find(ptcl[i].id) == index2id_map.end());
      ptcl[i].t_reg = ptcl[i].t_irr = t_global;
      index2id_map[ptcl[i].id] = i;
    }
    ptcl2add.clear();
    
    local_n = ptcl.size();
    index2id_map.clear();
    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      assert(index2id_map.find(it->id) == index2id_map.end());
      index2id_map[it->id] = it - ptcl.begin();
    }

    init_model();
  }
    
  void Nbody::recommit_particles()
  {
    fprintf(stderr, " hacs64::recommit_particles() -- \n");
    /***** first we remove particles ******/

#if 1
    {
      std::sort(ptcl2remove.begin(), ptcl2remove.end());
      std::vector<int> unique_id2rm;
      unique_id2rm.push_back(ptcl2remove[0]);
      for (std::vector<int>::iterator it = ptcl2remove.begin() + 1; it != ptcl2remove.end(); it++)
        if (*it != unique_id2rm.back())
          unique_id2rm.push_back(*it);
      if (ptcl2remove.size() != unique_id2rm.size())
          fprintf(stderr, "  warrning: it appears that some particles were asked to be removed twice ...\n");

      for (std::vector<int>::iterator it = unique_id2rm.begin(); it != unique_id2rm.end(); it++)
        ptcl.erase(ptcl.begin() + *it);
    }

    {
      const int n2add = ptcl2add.size();
      assert((int)ptcl.size() + n2add <= nmax);
      for (int i = 0; i < n2add; i++)
      {
        ptcl.push_back(ptcl2add[i]);
        ptcl.back().t_reg = ptcl.back().t_irr = t_global;
      }
    }
#endif
    ptcl2add.clear();
    ptcl2remove.clear();

    fprintf(stderr, " old_n= %d :: new_n= %d \n", local_n, (int)ptcl.size());
    fprintf(stderr, " cyclical_idx= %u\n", cyclical_idx);
    local_n = ptcl.size();
    index2id_map.clear();
    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      assert(index2id_map.find(it->id) == index2id_map.end());
      index2id_map[it->id] = it - ptcl.begin();
    }


    for (int i = 0; i < local_n; i++)
    {
      assert(index2id_map.find(ptcl[i].id) != index2id_map.end());
      assert(index2id_map[ptcl[i].id] == i);
    }


    init_model();

  }

  void Nbody::__predictmodel()
  {
    fprintf(stderr, " hacs64::__predictmodel() \n");
    const int ni = ptcl.size();
    predictor.resize(ni);
    for (int i = 0; i < ni; i++)
      predictor[i] = Predictor(ptcl[i], t_global - ptcl[i].t_irr);
  }

  void Nbody::__synchmodel()
  {
    fprintf(stderr, " hacs64::__synchmodel() \n");
    if (!is_synched)
    {
      fprintf(stderr, "  force sync \n");
      const int ni = ptcl.size();

#if 0

      irr_list.clear();
      for (int i = 0; i < ni; i++)
        irr_list.push_back(i);
      reg_list = irr_list;


#if 0
      for (int i = 0; i < ni; i++)
      {
        ptcl[i].reg_rung = t_global > ptcl[i].t_reg ? scheduler.get_rung(t_global - ptcl[i].t_reg) : ptcl[i].reg_rung;
        ptcl[i].irr_rung = t_global > ptcl[i].t_irr ? scheduler.get_rung(t_global - ptcl[i].t_irr) : ptcl[i].irr_rung;
        assert(scheduler.get_dt_pred(ptcl[i].reg_rung) == t_global - ptcl[i].t_reg);
        assert(scheduler.get_dt_corr(ptcl[i].irr_rung) == t_global - ptcl[i].t_irr);
      }
#endif

      scheduler.flush();
      iteration--;
      iterate_active();

#if 0
      for (int i = 0; i < ni; i++)
      {
        assert(ptcl[i].t_reg == t_global);
        assert(ptcl[i].t_irr == t_global);
      }
#endif
#endif

#if 1
      assert(!ptcl_pre.empty());
      for (Particle::Vector::iterator it = ptcl_pre.begin(); it != ptcl_pre.end(); it++)
      {
        assert(index2id_map.find(it->id) != index2id_map.end());
        const int i = index2id_map[it->id];
        assert(i >= 0);
        assert(i < ni);
        
        const Particle &pi = *it;
        assert(ptcl[i].id == it->id);
        ptcl[i] = pi;
        
        reg_ptr->set_jp(i, regf4::Particle(
              pi.mass, t_global - dt_global, pi.h2,
              pi.pos, pi.vel,
              pi.ftot.acc, pi.ftot.jrk));
        irr_ptr->set_jp(i, irrf6::Particle(
              pi.mass, t_global - dt_global,
              pi.pos, pi.vel,
              pi.ftot.acc, pi.ftot.jrk, pi.ftot.snp, pi.ftot.crk));
      }

      assert(ngb_list_id_pre.size() == ngb_list_pre.size());
      const int nngb = ngb_list_id_pre.size();
      for (int ix = 0 ; ix < nngb; ix++)
      {
        const int i = ngb_list_id_pre[ix];
        ngb_list[i] = ngb_list_pre   [ix];
        irr_ptr->set_list(i, ngb_list[i]);
      }


      reg_ptr->commit_changes();
      irr_ptr->commit_changes();

      irr_list.clear();
      for (int i = 0; i < ni; i++)
        irr_list.push_back(i);
      reg_list = irr_list;
      
      scheduler.flush();
      iteration--;
      iterate_active();
      
#endif

      assert(is_synched);
    }
  }
    


}
