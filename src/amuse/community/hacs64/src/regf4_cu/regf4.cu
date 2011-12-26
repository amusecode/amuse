#include <algorithm>
#include <map>
#include "regf4.h"
#include <cassert>
#include "Timer.h"


#include "cuda_pointer.h"
#include "cutil_inline.h"

namespace regf4
{

#include "dev_regf4.cu"

#define NJBLOCK  32
#define NJBLOCK2 32
#define NJTHREAD 128
#define NPIPE    512

  enum {NGB_PER_BLOCK = NGBlist::NGB_MAX};

  /************************************************/

  int ni_tot, ni_tot_round, ni_max;
  double h2max;
  double t_predictor, t_interaction;
  unsigned long long n_interaction;

  /************************************************/

  std::vector<NGBlist  > list;

  std::vector<Force > *force_result;
  std::vector<double> *gpot_result;
  float eps2;

  /***************** data vectors BEGIN *******************/

  float dt_tick;                                   // basic time unit
  unsigned int tsys;                               // system time in basic time units (32 bit mantissa)

  cuVector<dev_particle,  false> ptcl_list;        // particle list
  //cuVector<int,           false> active_list;      // list of active particles addr
  std::vector<int> active_list;
  cuVector<dev_predictor, false> predictor_list;   // predictor list
  cuVector<float,         false> gpot_list;        // list of potential

#if 0
  std::vector<dev_particle> jp_ptcl_vector;
  std::vector<int         > jp_addr_vector;
#else
  std::map<int, dev_particle> jp_ptcl_map;
#endif

  cuVector<dev_particle,  false> jp_ptcl_list;          // list of jp_particles
  cuVector<int         ,  false> jp_addr_list;          // list of jp_particles' addresses

  cuVector<unsigned int, false> ngb_list, ngb_reduced_list;
  cuVector<dev_force,    false> force_list, force_reduced_list;

  cuVector<int2, false> ngb_offset_list;

  /***************** data vectors  END  *******************/

  regf::regf(const int _ni_max, const double _h2max) {}

  regf::regf(const int _ni_max, const double _h2max, const double _dt_tick)
  {
    ni_max  = ((_ni_max - 1)/(NJTHREAD*NJBLOCK) + 1)*NJBLOCK*NJTHREAD;

    ni_tot       = _ni_max;
    ni_tot_round =  ni_max;
    fprintf(stderr, " ni_tot= %d   ni_tot_round= %d \n", ni_tot, ni_tot_round);
    h2max   = _h2max;
    dt_tick = (float)_dt_tick;

    assert(dt_tick > 0.0);

    list.resize(ni_max);

    /////////////////

    int ndevice;
    assert(cudaGetDeviceCount(&ndevice) == 0);
    fprintf(stderr, " regf::regf -  %d CUDA devices found \n", ndevice);
    fprintf(stderr, " regf::regf - using last (%d) device \n", ndevice-1);

    assert(cudaSetDevice(ndevice-1) == cudaSuccess);

    //////////////

    ptcl_list     .allocate(ni_max);
    //    active_list   .allocate(ni_max);
    predictor_list.allocate(ni_max);
    gpot_list     .allocate(ni_max);

    jp_addr_list  .allocate(ni_max);
    jp_ptcl_list  .allocate(ni_max);

    for (int i = 0; i < ni_max; i++)
    {
      ptcl_list[i].pos  = dcuvec3(0.0f);
      ptcl_list[i].vel  = fcuvec3(0.0f);
      ptcl_list[i].acc  = fcuvec3(0.0f);
      ptcl_list[i].jrk  = fcuvec3(0.0f);
      ptcl_list[i].mass = 0.0;
      ptcl_list[i].h2   = 0.0;
      ptcl_list[i].time = 0;
      ptcl_list[i].id   = i;
    }
    ptcl_list.h2d();

    ngb_list        .allocate(NPIPE*NJBLOCK*NGB_PER_BLOCK);
    ngb_reduced_list.allocate(NPIPE*NGBlist::NGB_MAX);

    force_list        .allocate(NPIPE*NJBLOCK2);
    force_reduced_list.allocate(NPIPE);

    ngb_offset_list.allocate(NPIPE*NJBLOCK);

    //////////////

    t_predictor = t_interaction = 0.0;
    n_interaction = 0;

    /////////////////

    CUDA_SAFE_CALL(cudaMemcpyToSymbol("regf4::DT_TICK", &dt_tick, sizeof(float)));
  }

  regf::~regf() 
  {
    ptcl_list.cufree();
    //    active_list.cufree();
    predictor_list.cufree();
    gpot_list.cufree();

    jp_addr_list.cufree();
    jp_ptcl_list.cufree();

    ngb_list.cufree();
    ngb_reduced_list.cufree();

    force_list.cufree();
    force_reduced_list.cufree();

    ngb_offset_list.cufree();
  }

  int regf::resize(const int _ni)
  {
    assert(_ni <= ni_max);
    ni_tot        = _ni;
    ni_tot_round  = ((_ni - 1)/(NJTHREAD*NJBLOCK) + 1)*NJBLOCK*NJTHREAD;
    ptcl_list.d2h();
    for (int i = ni_tot; i < ni_max; i++)
    {
      ptcl_list[i].pos  = dcuvec3(0.0f);
      ptcl_list[i].vel  = fcuvec3(0.0f);
      ptcl_list[i].acc  = fcuvec3(0.0f);
      ptcl_list[i].jrk  = fcuvec3(0.0f);
      ptcl_list[i].mass = 0.0;
      ptcl_list[i].h2   = 0.0;
      ptcl_list[i].time = 0;
      ptcl_list[i].id   = -1;
    }
    ptcl_list.h2d();
    return 0;
  }

  int regf::set_ti(const double ti)
  {
    tsys = (unsigned int)(ti/(double)dt_tick);
    return 0;
  }

  /********** jp-particles ********/

  int regf::set_jp(const int iaddr, const Particle &pi)
  {
    assert(iaddr < ni_tot);
#if 0
    jp_ptcl_vector.push_back(pi);
    jp_addr_vector.push_back(iaddr);
#else
    jp_ptcl_map[iaddr] = pi;
#endif

    return 0;
  }
  int regf::set_jp(const std::vector<int> &ilist, const std::vector<Particle> &ptcl_list)
  {
    for (size_t i = 0; i < ilist.size(); i++)
    {
      assert(ilist[i] < ni_tot);
#if 0
      jp_ptcl_vector.push_back(ptcl_list[i]);
      jp_addr_vector.push_back(ilist[i]);
#else
      jp_ptcl_map[ilist[i]] = ptcl_list[i];
#endif
    }
    return 0;
  }

  /********** NGB lists ********/

  int regf::set_list(const int iaddr, const NGBlist &ngb)
  {
    assert(iaddr < ni_tot);
    list[iaddr] = ngb;
    return 0;
  }
  int regf::set_list(const std::vector<int> &ilist, const std::vector<NGBlist> &ngb_list)
  {
    for (size_t i = 0; i < ilist.size(); i++)
    {
      assert(ilist[i] < ni_tot);
      list[ilist[i]] = ngb_list[i];
    }
    return 0;
  }

  int regf::get_list(const int iaddr, NGBlist &ngb) 
  {
    assert(iaddr < ni_tot);
    ngb = list[iaddr];
    return 0;
  }
  int regf::get_list(const std::vector<int>&ilist, std::vector<NGBlist> &ngb_list)
  {
    ngb_list.resize(ilist.size());
    for (size_t i = 0; i < ilist.size(); i++)
    {
      assert(ilist[i] < ni_tot);
      ngb_list[i] = list[ilist[i]];
    }
    return 0;
  }

  /************************************/

  void copy_jp_to_device()
  {
    if (jp_ptcl_map.empty()) return;

    /********** copy new/updated jp-particles from the host to the device ****/

    const int nj = jp_ptcl_map.size();
    assert(nj <= ni_tot);

    std::vector<int> iaddr(nj);
    std::vector<dev_particle> iptcl(nj);
    int cnt = 0;
    for (std::map<int, dev_particle>::iterator it = jp_ptcl_map.begin(); it != jp_ptcl_map.end(); it++)
    {
      iaddr[cnt] = it->first;
      iptcl[cnt] = it->second;
      cnt++;
    }

    for (int j = 0; j < nj; j++)
    {
      iptcl[j].id   = iaddr[j];
      iptcl[j].iPad = j;
    }

    jp_ptcl_list.copy(&iptcl[0], nj);
    jp_ptcl_list.h2d(nj);

    jp_addr_list.copy(&iaddr[0], nj);
    jp_addr_list.h2d(nj);

    const int nthreads = 256;
    const int nblocks  = (nj-1)/nthreads + 1;
    const dim3 block(nthreads, 1, 1);
    const dim3 grid (nblocks,  1, 1);
    dev_move_particles<<<grid, block>>>(nj, jp_addr_list, jp_ptcl_list, ptcl_list);

    jp_ptcl_map.clear();
  }

  /************************************/

  int regf::force_first(const std::vector<int> &ilist, std::vector<Force> &force, const double eps2_in)
  {
    force_result     = &force;
    eps2             = eps2_in;

    active_list = ilist;

    if (ilist.empty()) return 0;


    copy_jp_to_device();

    /******** execute predictor kernel ********/

    const int nthreads = 128;
    const int nblocks  = (ni_tot_round - 1)/nthreads + 1;
    const dim3 block(nthreads, 1, 1);
    const dim3 grid (nblocks,  1, 1);
    dev_predict_ptcl<<<grid, block>>>(ni_tot_round, tsys, ptcl_list, predictor_list, gpot_list);

    return 0;
  }

  /************************************/

  int regf::force_last()
  {
    const double t0 = get_wtime();

    std::vector<Force> &force = *force_result;

    if (active_list.empty()) return 0;
    force.resize(active_list.size());

    CUDA_SAFE_CALL(cudaMemcpyToSymbol("regf4::EPS2", &eps2, sizeof(float)));

    cuVector<int, false> dev_active_list;
    dev_active_list.allocate(active_list.size());

    std::vector<int> active_idx(active_list.size());
    for (size_t i = 0; i < active_list.size(); i++)
      active_idx[i] = i;


    
    int nev_tot = 0, nev_tot1 = 0;

    while(!active_list.empty())
    {
      std::vector<int  > failed_idx;
      std::vector<float> new_h2_list;

      const int ni_active = active_list.size();
      dev_active_list.init(&active_list[0], ni_active);
      dev_active_list.h2d(ni_active);
      active_list.clear();

      for (int ix = 0; ix < ni_active; ix += NPIPE)
      {
        const int ni = std::min(NPIPE, ni_active - ix);
        nev_tot  += NPIPE;
        nev_tot1 += ni;

        {
          const int nthreads     = NJTHREAD;
          const int niblocks     = (ni-1)/nthreads + 1;
          const int njblocks     = NJBLOCK;
          const int nj_per_block = (ni_tot_round-1)/njblocks + 1;
          const dim3 block(nthreads, 1, 1);
          const dim3 grid (niblocks, njblocks, 1);
#if 0
          cudaThreadSetCacheConfig(cudaFuncCachePreferL1);
#endif
          dev_regf<NJTHREAD, NJBLOCK, NJBLOCK2, NGB_PER_BLOCK><<<grid, block>>>(
              ni,
              nj_per_block,
              dev_active_list+ix,
              predictor_list,
              gpot_list,  // dt_list
              force_list,
              ngb_list);
#if 0
          cudaThreadSynchronize();
          cutilCheckMsg("dev_regf failed");
#endif
        }

        {
          const int nthreads = NJBLOCK2 > 64 ? NJBLOCK2 : 64; 
          const dim3 block(nthreads, 1, 1);
          const dim3 grid (NPIPE,    1, 1);
          dev_reduce_regf<nthreads, NJBLOCK, NJBLOCK2><<<grid, block>>>(
              force_list,
              ngb_offset_list,
              force_reduced_list);
          force_reduced_list.d2h(ni);

          dev_reduce_ngb<nthreads, NJBLOCK, NGB_PER_BLOCK, NGBlist::NGB_MAX><<<grid, block>>>(
              ngb_offset_list,
              ngb_list,
              ngb_reduced_list);
          ngb_reduced_list.d2h(ni*NGBlist::NGB_MAX);
        }


        for (int ii = 0; ii < ni; ii++)
        {
          Force &fi = force[active_idx[ix+ii]];
          fi.acc = dvec3(force_reduced_list[ii].accx.dbl, force_reduced_list[ii].accy.dbl, force_reduced_list[ii].accz.dbl);
          fi.jrk = dvec3(force_reduced_list[ii].jrk.x, force_reduced_list[ii].jrk.y, force_reduced_list[ii].jrk.z);
          fi.h2  = force_reduced_list[ii].h2;
          const int nj = force_reduced_list[ii].nngb;
          double fac = (std::pow(NGBMEAN*1.0/(nj+1), 2.0/3.0) + 1)*0.5;
          if (fac > 1.0) fac = std::min(fac, 1.25);
          else           fac = std::max(fac, 1.0/1.25);
          fi.h2 *= fac;
#if 0
          fi.h2  = std::min(fi.h2, h2max);
#endif
          const int addr = dev_active_list[ix+ii];
          list[addr].clear();
          if (nj >= NGBlist::NGB_MAX) 
//          if (nj > NGBMAX || nj < NGBMIN) 
          {
            fprintf(stderr, " ** WARNING **  new_ngbi >= NGBBUFF, addr= %d  nj= %d\n",
                addr, nj);
            failed_idx .push_back(active_idx[ix+ii]);
            active_list.push_back(addr);
            new_h2_list.push_back(fi.h2);
          }
          else
            for (int jj = 0; jj < nj; jj++)
              list[addr].push_back(ngb_reduced_list[jj+ii*NGBlist::NGB_MAX]);
        }
      }
      if (!new_h2_list.empty())
      {
        predictor_list.d2h(ni_tot);
        for (size_t i = 0; i < active_list.size(); i++)
          predictor_list[active_list[i]].h2 = new_h2_list[i];
        predictor_list.h2d(ni_tot);
        new_h2_list.clear();
      }
      active_idx = failed_idx;
      assert(active_idx.size() == active_list.size());
      failed_idx.clear();
    }
#if 0
    const double t1 = get_wtime();
    const double dt = t1 -t0;
    fprintf(stderr, "done in  %g sec  %g/%g GFLOP/s  [ %d %d/%d ]\n",
        dt, 
        (double)ni_tot_round*(double)nev_tot1*60.0/1.0e9/dt,
        (double)ni_tot_round*(double)nev_tot*60.0/1.0e9/dt,
        ni_tot_round, nev_tot1, nev_tot);
#endif

    return 0;
  }

  /************************************/

  /**********************/

  int regf::potential_first(std::vector<double> &gpot, const double eps2_pot_in)
  {
    gpot_result = &gpot;
    eps2        = eps2_pot_in;

    copy_jp_to_device();

    return 0;
  }

  int regf::potential_last()
  {

    CUDA_SAFE_CALL(cudaMemcpyToSymbol("regf4::EPS2", &eps2, sizeof(float)));

    std::vector<double> &gpot = *gpot_result;

    const int nthreads = 256;
    const int nblocks  = (ni_tot - 1)/nthreads + 1;
    const int nj       = nblocks * nthreads;
    const dim3 block(nthreads, 1, 1);
    const dim3 grid (nblocks,  1, 1);
    dev_compute_potential<256><<<grid, block>>>(ni_tot, nj, ptcl_list, gpot_list);
    cudaThreadSynchronize();
    cutilCheckMsg("dev_compute_potential failed");

    gpot_list.d2h(ni_tot);

    gpot.resize(ni_tot); 
    for (int i = 0; i < ni_tot; i++)
      gpot[i] = gpot_list[i];
    return 0;
  }


  __host__ dev_particle::dev_particle(const regf4::Particle &p)
  {
    assert(sizeof(dev_particle) == 32*sizeof(int));
    pos  = dcuvec3(p.pos.x, p.pos.y, p.pos.z);
    vel  = fcuvec3(p.vel.x, p.vel.y, p.vel.z);
    acc  = fcuvec3(p.acc.x, p.acc.y, p.acc.z);
    jrk  = fcuvec3(p.jrk.x, p.jrk.y, p.jrk.z);
    mass = p.mass;
    h2   = p.h2;
    time = (unsigned int)(p.time/(double)regf4::dt_tick);
  }
};

