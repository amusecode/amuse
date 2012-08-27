#ifndef GRAVITY_PROTOTYPES_H
#define GRAVITY_PROTOTYPES_H 
#include "../copyright.h"
/*============================================================================*/
/*! \file prototypes.h
 *  \brief Prototypes for all public functions in the /src/gravity directory */
/*============================================================================*/
#include <stdio.h>
#include <stdarg.h>
#include "../athena.h"
#include "../defs.h"

#include "../config.h"

/* bvals_grav.c  */
#ifdef SELF_GRAVITY
void bvals_grav_init(MeshS *pM);
void bvals_grav_fun(DomainS *pD, enum BCDirection dir, VGFun_t prob_bc);
void bvals_grav(DomainS *pDomain);
#endif

/* selfg.c  */
#ifdef SELF_GRAVITY
VDFun_t selfg_init(MeshS *pM);
void selfg_flux_correction(GridS *pG);
#endif /* SELF_GRAVITY */

/* selfg_multigrid.c  */
#ifdef SELF_GRAVITY
void selfg_multig_1d(DomainS *pD);
void selfg_multig_2d(DomainS *pD);
void selfg_multig_3d(DomainS *pD);
void selfg_multig_3d_init(MeshS *pM);
#endif /* SELF_GRAVITY */

/* selfg_fft.c  */
#ifdef SELF_GRAVITY
#if defined(FFT_ENABLED) && defined(SELF_GRAVITY_USING_FFT)
void selfg_fft_1d(DomainS *pD);
void selfg_fft_2d(DomainS *pD);
void selfg_fft_3d(DomainS *pD);
void selfg_fft_2d_init(MeshS *pM);
void selfg_fft_3d_init(MeshS *pM);
#endif /* FFT_ENABLED */
#if defined(FFT_ENABLED) && defined(SELF_GRAVITY_USING_FFT_OBC)
void selfg_fft_obc_3d(DomainS *pD);
void selfg_fft_obc_3d_init(MeshS *pM);
#endif /* FFT_ENABLED SELF_GRAVITY_USING_FFT_OBC */



#endif /* SELF_GRAVITY */

#endif /* GRAVITY_PROTOTYPES_H */
