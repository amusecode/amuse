/* ################################################################################## */
/* ###                                                                            ### */
/* ###                                 Gadgetmp2                                  ### */
/* ###                                                                            ### */
/* ###   Original: Gadget2 in the version used in Amuse                           ### */
/* ###   Author: Gadget2 and Amuse contributors                                   ### */
/* ###                                                                            ### */
/* ###   Modified: July 2020                                                      ### */
/* ###   Author: Thomas Schano                                                    ### */
/* ###                                                                            ### */
/* ###   Changes are intended to enable precise calculations in                   ### */
/* ###   non periodic small domain simulations in which comoving parts            ### */
/* ###   are simulated in std precision                                           ### */
/* ###                                                                            ### */
/* ################################################################################## */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef NOMPI
#include <mpi.h>
#endif

//#include "allvars.hpp"
#include "proto.hpp"

/*! \file allocate.c
 *  \brief routines for allocating particle and tree storage
 */

/*! Allocates a number of small buffers and arrays, the largest one being
 *  the communication buffer. The communication buffer itself is mapped
 *  onto various tables used in the different parts of the force
 *  algorithms. We further allocate space for the top-level tree nodes, and
 *  auxiliary arrays for the domain decomposition algorithm.
 */


void gadgetmp2::allocate_commbuffers(void)
{
    size_t bytes, i;
    if (CommBuffer == nullptr)
        if(!(CommBuffer = malloc(bytes = All.BufferSize * 1024 * 1024)))
        {
            printf("failed to allocate memory for `CommBuffer' (%g MB).\n", bytes / (1024.0 * 1024.0));
            endrun(2);
        };
    if (Exportflag == nullptr) Exportflag = new char[NTask];
    if (DomainStartList == nullptr) DomainStartList = new int[NTask];
    if (DomainEndList == nullptr) DomainEndList = new int[NTask];

    if (TopNodes == nullptr) TopNodes = new topnode_data[MAXTOPNODES];

    if (DomainWork == nullptr){
        DomainWork = new my_float[MAXTOPNODES];
    }else{
        for(i=0;i<MAXTOPNODES; i++)
            DomainWork[i].set_prec(my_float::get_default_prec());
    }
    if (DomainCount == nullptr) DomainCount = new int[MAXTOPNODES];
    if (DomainCountSph == nullptr) DomainCountSph = new int[MAXTOPNODES];
    if (DomainTask == nullptr) DomainTask = new int[MAXTOPNODES];
    if (DomainNodeIndex == nullptr) DomainNodeIndex = new int[MAXTOPNODES];
    if (DomainTreeNodeLen == nullptr){
        DomainTreeNodeLen = new my_float[MAXTOPNODES];
    }else{
        for(i=0;i<MAXTOPNODES; i++)
            DomainTreeNodeLen[i].set_prec(my_float::get_default_prec());
    }
    size_t all_reduce_size = All_Reduce_buff::gen_size();
    if (DomainHmax != nullptr) free(DomainHmax);
    DomainHmax = (All_Reduce_buff*)malloc(MAXTOPNODES*all_reduce_size);
    if (DomainMoment != nullptr) free(DomainMoment);
    DomainMoment = (DomainNODE*)malloc(MAXTOPNODES*DomainNODE::gen_size());

    size_t grav_size = gravdata_in::gen_size();
    All.BunchSizeForce = (All.BufferSize * 1024 * 1024) / (sizeof(struct gravdata_index) + 2 * grav_size);

    if(All.BunchSizeForce & 1)
    {
        All.BunchSizeForce -= 1;	/* make sure that All.BunchSizeForce is an even number --> 8-byte alignment for 64bit processors */
    };
    GravDataIndexTable = (struct gravdata_index *) CommBuffer;
    GravDataIn = (gravdata_in *) (GravDataIndexTable + All.BunchSizeForce);
    GravDataGet =(gravdata_in *)((size_t)GravDataIn + All.BunchSizeForce*grav_size);
    GravDataOut = GravDataIn;	/* this will overwrite the GravDataIn-Table */
    GravDataResult = GravDataGet;	/* this will overwrite the GravDataGet-Table */

    size_t dens_size_in = densdata_in::gen_size();
    size_t dens_size_out = densdata_out::gen_size();
    All.BunchSizeDensity =
    (All.BufferSize * 1024 * 1024) / (2 * dens_size_in + 2 * dens_size_out);

    DensDataIn = (densdata_in *) CommBuffer;
    //DensDataGet = DensDataIn + All.BunchSizeDensity;
    DensDataGet = (densdata_in *)((size_t)DensDataIn + All.BunchSizeDensity*dens_size_in);
    //DensDataResult = (struct densdata_out *) (DensDataGet + All.BunchSizeDensity);
    DensDataResult = (densdata_out *)((size_t)DensDataGet + All.BunchSizeDensity*dens_size_in);
    //DensDataPartialResult = DensDataResult + All.BunchSizeDensity;
    DensDataPartialResult = (densdata_out *)((size_t)DensDataResult + All.BunchSizeDensity*dens_size_out);


    size_t hydro_size_in = hydrodata_in::gen_size();
    size_t hydro_size_out = hydrodata_out::gen_size();
    All.BunchSizeHydro = (All.BufferSize * 1024 * 1024) / (2 * hydro_size_in + 2 * hydro_size_out);

    HydroDataIn = (hydrodata_in *) CommBuffer;
    //HydroDataGet = HydroDataIn + All.BunchSizeHydro;
    HydroDataGet = (hydrodata_in *)((size_t)HydroDataIn + All.BunchSizeHydro*hydro_size_in);
    //HydroDataResult = (struct hydrodata_out *) (HydroDataGet + All.BunchSizeHydro);
    HydroDataResult = (hydrodata_out *)((size_t)HydroDataGet + All.BunchSizeHydro*hydro_size_in);
    //HydroDataPartialResult = HydroDataResult + All.BunchSizeHydro;
    HydroDataPartialResult = (hydrodata_out *)((size_t)HydroDataResult + All.BunchSizeHydro*hydro_size_out);

    size_t particle_size=particle_data_buff::gen_size();
    size_t sphparticle_size=sph_particle_data_buff::gen_size();
    All.BunchSizeDomain = (All.BufferSize * 1024 * 1024) / (2*(particle_size + sphparticle_size) + sizeof(peanokey));

    if(All.BunchSizeDomain & 1)
    {
        All.BunchSizeDomain -= 1;	/* make sure that All.BunchSizeDomain is even --> 8-byte alignment of DomainKeyBuf for 64bit processors */
    };
    DomainPartBuf_s = (particle_data_buff *) CommBuffer;
    DomainPartBuf_r = (particle_data_buff *)((size_t)DomainPartBuf_s + All.BunchSizeDomain*particle_size);
    //DomainSphBuf = (struct sph_particle_data *) (DomainPartBuf + All.BunchSizeDomain);
    DomainSphBuf_s = (sph_particle_data_buff *)((size_t)DomainPartBuf_r + All.BunchSizeDomain*particle_size);
    DomainSphBuf_r = (sph_particle_data_buff *)((size_t)DomainSphBuf_s + All.BunchSizeDomain*sphparticle_size);
    //DomainKeyBuf = (peanokey *) (DomainSphBuf + All.BunchSizeDomain);
    DomainKeyBuf = (peanokey *)((size_t)DomainSphBuf_r + All.BunchSizeDomain*sphparticle_size);

    #ifdef TIMESTEP_LIMITER
    size_t time_size = timedata_in::gen_size();
    All.BunchSizeTime = (All.BufferSize * 1024 * 1024) / (2 * time_size);
    TimeDataIn = (timedata_in *) CommBuffer;
    //  TimeDataGet = TimeDataIn + All.BunchSizeTime;
    TimeDataGet =(timedata_in *)((size_t)TimeDataIn + All.BunchSizeTime*time_size);
    #endif


    if((All.BufferSize * 1024 * 1024) < all_reduce_size * NTask)
        exit(0);

    all_reduce_buff = (All_Reduce_buff *) CommBuffer;

    if(ThisTask == 0)
    {
        printf("\nAllocated %d MByte communication buffer per processor.\n\n", All.BufferSize);
        printf("Communication buffer has room for %d particles in gravity computation\n", All.BunchSizeForce);
        printf("Communication buffer has room for %d particles in density computation\n", All.BunchSizeDensity);
        printf("Communication buffer has room for %d particles in hydro computation\n", All.BunchSizeHydro);
        printf("Communication buffer has room for %d particles in domain decomposition\n", All.BunchSizeDomain);
        printf("\n");
    }
//    DEBUG << "sph" <<  hydro_size_in << " sppppp " << All.BunchSizeHydro<<"\n";
//    DEBUG.flush();
}



/*! This routine allocates memory for particle storage, both the
 *  collisionless and the SPH particles.
 */
void gadgetmp2::allocate_memory(void)
{
    size_t bytes;
    long long bytes_tot = 0;

    if(All.MaxPart > 0)
    {
        if(!(P = new particle_data[All.MaxPart]))
        {
            bytes = All.MaxPart * particle_data::get_size_();
            printf("failed to allocate memory for `P' (%g MB).\n", bytes / (1024.0 * 1024.0));
            endrun(1);
        }
        bytes = All.MaxPart * particle_data::get_size_();
        bytes_tot += bytes;

        if(ThisTask == 0)
            printf("\nAllocated %g MByte for particle storage. %ld\n\n", bytes_tot / (1024.0 * 1024.0), sizeof(struct particle_data));
    }

    if(All.MaxPartSph > 0)
    {
        bytes_tot = 0;

        if(!(SphP = new sph_particle_data[All.MaxPartSph]))
        {
            bytes = All.MaxPartSph * sph_particle_data::get_size_();
            printf("failed to allocate memory for `SphP' (%g MB) %ld.\n", bytes / (1024.0 * 1024.0), sizeof(struct sph_particle_data));
            endrun(1);
        }
        bytes = All.MaxPartSph * sph_particle_data::get_size_();
        bytes_tot += bytes;

        if(ThisTask == 0)
            printf("Allocated %g MByte for storage of SPH data. %ld\n\n", bytes_tot / (1024.0 * 1024.0), sizeof(struct sph_particle_data));
    }
}




/*! This routine frees the memory for the particle storage.  Note: We don't
 *  actually bother to call it in the code...  When the program terminats,
 *  the memory will be automatically freed by the operating system.
 */
void gadgetmp2::free_memory(void)
{
    if(All.MaxPartSph > 0)
    {    //free(SphP);
        delete[] SphP;
        SphP=nullptr;
    }

    if(All.MaxPart > 0)
    {    //free(P);
        delete[] P;
        P=nullptr;
    }
}


