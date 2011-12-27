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
      index2id_map[ptcl[i].id] = i;
    }
    ptcl2add.clear();

    init_model();
  }
    
  void Nbody::recommit_particles()
  {
    fprintf(stderr, " hacs64::recommit_particles() -- \n");
    /***** first we remove particles ******/

    {
      std::sort(ptcl2remove.begin(), ptcl2remove.end());
      std::vector<int> unique_id2rm;
      unique_id2rm.push_back(ptcl2remove[0]);
      for (std::vector<int>::iterator it = ptcl2remove.begin() + 1; it != ptcl2remove.end(); it++)
        if (*it != unique_id2rm.back())
          unique_id2rm.push_back(*it);
      if (ptcl2remove.size() != unique_id2rm.size())
          fprintf(stderr, "  warrning: it appears that some particles were asked to be removed twice ...\n");

      ptcl2remove.clear();
      for (std::vector<int>::iterator it = unique_id2rm.begin(); it != unique_id2rm.end(); it++)
        ptcl.erase(ptcl.begin() + *it);
    }

    {
      const int n2add = ptcl2add.size();
      assert((int)ptcl.size() + n2add <= nmax);
      for (int i = 0; i < n2add; i++)
        ptcl.push_back(ptcl2add[i]);
      ptcl2add.clear();
    }

    index2id_map.clear();
    for (Particle::Vector::iterator it = ptcl.begin(); it != ptcl.end(); it++)
    {
      assert(index2id_map.find(it->id) == index2id_map.end());
      index2id_map[it->id] = it - ptcl.begin();
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
