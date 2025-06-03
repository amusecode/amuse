#include "../copyright.h"

// #define PM_CIC
#define PM_TSC

/*============================================================================*/
/*! \file selfg_fft_obc.c
 *  \brief Contains functions to solve Poisson's equation for self-gravity in
 *   3D using FFTs, using OPEN BCs in all three directions 
 *
 *
 *   The function uses FFTW3.x, and for MPI parallel use Steve Plimpton's
 *   block decomposition routines added by N. Lemaster to /athena/fftsrc.
 *   This means to use these fns the code must be
 *   - (1) configured with --with-gravity=fft_obc --enable-fft
 *   - (2) compiled with links to FFTW libraries (may need to edit Makeoptions)
 *
 *   For PERIODIC BCs, use selfg_fft() functions.
 *
 * CONTAINS PUBLIC FUNCTIONS:
 *   - selfg_fft_obc_2d() - 2D Poisson solver using FFTs
 *   - selfg_fft_obc_3d() - 3D Poisson solver using FFTs
 *   - selfg_fft_obc_2d_init() - initializes FFT plans for 2D
 *   - selfg_fft_obc_3d_init() - initializes FFT plans for 3D
 *============================================================================*/

#include <math.h>
#include <float.h>
#include "../defs.h"
#include "../athena.h"
#include "../globals.h"
#include "prototypes.h"
#include "../prototypes.h"

#ifdef SELF_GRAVITY_USING_FFT_OBC

#ifndef FFT_ENABLED
#error self gravity with FFT requires configure --enable-fft
#endif /* FFT_ENABLED */

/* plans for forward and backward FFTs; work space for FFTW */
static struct ath_2d_fft_plan *fplan2d, *bplan2d;
static struct ath_3d_fft_plan *fplan3d, *bplan3d;
static ath_fft_data *work=NULL;


#ifdef STATIC_MESH_REFINEMENT
#error self gravity with FFT not yet implemented to work with SMR
#endif


/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_obc_2d(DomainS *pD)
 *  \brief Only works for uniform grid, periodic boundary conditions
 */
void selfg_fft_obc_2d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i,ioff,ip, is = pG->is, ie = pG->ie;
  int j,joff,jp, js = pG->js, je = pG->je;
  int ks = pG->ks;
  Real pcoeff,offset,kx,ky;
  int Nx1=pD->Nx[0], Nx2=pD->Nx[1];
  int hNx1=Nx1/2, hNx2=Nx2/2;
  Real dkx=2.0*PI/(double)(Nx1), dky=2.0*PI/(double)(Nx2);
#ifdef RADIATION
  RadParticleS *pRP=NULL;
  Real x1,x2,x3;
  Real tmp,W1[3],W2[3],W3[3],d;
#endif

  /* Copy current potential into old and zero-out */
  for (j=js-nghost; j<=je+nghost; j++) {
    for (i=is-nghost; i<=ie+nghost; i++) {
      pG->Phi_old[ks][j][i] = pG->Phi[ks][j][i];
      pG->Phi[ks][j][i] = 0.0;
    }
  }

  /* Loop over the index offsets (0 for even, 1 for odd) */
  for (joff=0; joff<=1; joff++) {
    for (ioff=0; ioff<=1; ioff++) {

      /* STEP 1: Forward FFT of 4\piG*(d-d0)*(\Delta z) */
      /* Add gas density into work array 0. */
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0] = pG->U[ks][j][i].d;
        }
      }

#ifdef RADIATION
      /* Interpolate particle density onto the mesh using CIC or TSC shape
       * functions as described in Hockney & Eastwood.  NOTE:  There are more
       * efficient ways to implement this, but since the radiation particle lists
       * will typically be small, there is no real need. */
      pRP = pG->radparticles;
      while (pRP) {
        /* Calculate indices and spatial coordinates of the lower nearest
         * neighbor grid cell.
         * x1 <= pRP->x1 < x1+pG->dx1, x2 <= pRP->x2 < x2+pG->dx2 */
        cc_ijk(pG,pRP->x1,pRP->x2,pRP->x3,&ip,&jp,&ks);
        cc_pos(pG,ip,jp,ks,&x1,&x2,&x3);

#ifdef PM_CIC
        /* Compute CIC weights.   */
        if (pRP->x1 > x1)
          W1[1] = (pRP->x1 - x1)/pG->dx1;  W1[0] = 1.0 - W1[1];
        else
          W1[0] = (x1 - pRP->x1)/pG->dx1;  W1[1] = 1.0 - W1[0];
        if (pRP->x2 > x2)
          W2[1] = (pRP->x2 - x2)/pG->dx2;  W2[0] = 1.0 - W2[1];
        else
          W2[0] = (x2 - pRP->x2)/pG->dx2;  W2[1] = 1.0 - W2[0];

        /* Add weighted particle density to work array 0. */
        d = pRP->m/(pG->dx1*pG->dx2);
        for (j=0; j<=1; j++) {
          for (i=0; i<=1; i++) {
            work[F2DI(ip+i-is,jp+j-js,pG->Nx[0],pG->Nx[1])][0] += d*W1[i]*W2[j];
          }
        }
#elif defined(PM_TSC)
        /* Compute TSC weights.   */
        tmp = 0.5 + (pRP->x1 - x1)/pG->dx1;
        W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2] - tmp + 0.5;  W1[1] = 1.0 - W1[0] - W1[2];
        tmp = 0.5 + (pRP->x2 - x2)/pG->dx2;
        W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2] - tmp + 0.5;  W2[1] = 1.0 - W2[0] - W2[2];

        /* Add weighted particle density to work array 0. */
        d = pRP->m/(pG->dx1*pG->dx2);
        for (j=0; j<=2; j++) {
          for (i=0; i<=2; i++) {
            work[F2DI(ip-1+i-is,jp-1+j-js,pG->Nx[0],pG->Nx[1])][0] += d*W1[i]*W2[j];
          }
        }
#endif

        pRP = pRP->next;
      }
#endif

      /* Copy work array 0 into work array 1, then multiply by complex offsets.
       * To compute offsets, note that indices relative to whole Domain are
       * needed. */
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1] =
            work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0];

          offset = 0.5*((i-is+pG->Disp[0])*ioff*dkx
                      + (j-js+pG->Disp[1])*joff*dky);
          work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0] *=  cos(offset);
          work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1] *= -sin(offset);
        }
      }

      /* Forward FFT */
      ath_2d_fft(fplan2d, work);

      /* STEP 2:  Compute potential in Fourier space.  Multiple loops are used
       * to avoid divide by zero at i=is,j=js, and to avoid if statement in
       * loop.  To compute kx,ky, note that indices relative to whole Domain
       * are needed. */
      if ((pG->Disp[0])==0 && (pG->Disp[1])==0 && ioff==0 && joff==0) {
        work[F2DI(0,0,pG->Nx[0],pG->Nx[1])][0] = 0.0;
        work[F2DI(0,0,pG->Nx[0],pG->Nx[1])][1] = 0.0;
      } else {
        ip = (pG->Disp[0] + hNx1) % Nx1 - hNx1;
        jp = (pG->Disp[1] + hNx2) % Nx2 - hNx2;
        kx = (ip + 0.5*ioff)*dkx/pG->dx1;
        ky = (jp + 0.5*joff)*dky/pG->dx2;
        pcoeff = -0.5/sqrt(SQR(kx)+SQR(ky));
        work[F2DI(0,0,pG->Nx[0],pG->Nx[1])][0] *= pcoeff;
        work[F2DI(0,0,pG->Nx[0],pG->Nx[1])][1] *= pcoeff;
      }

      for (j=js+1; j<=je; j++) {
        ip = (         pG->Disp[0] + hNx1) % Nx1 - hNx1;
        jp = ((j-js) + pG->Disp[1] + hNx2) % Nx2 - hNx2;
        kx = (ip + 0.5*ioff)*dkx/pG->dx1;
        ky = (jp + 0.5*joff)*dky/pG->dx2;
        pcoeff = -0.5/sqrt(SQR(kx)+SQR(ky));
        work[F2DI(0,j-js,pG->Nx[0],pG->Nx[1])][0] *= pcoeff;
        work[F2DI(0,j-js,pG->Nx[0],pG->Nx[1])][1] *= pcoeff;
      }

      for (i=is+1; i<=ie; i++) {
        for (j=js; j<=je; j++) {
          ip = ((i-is) + pG->Disp[0] + hNx1) % Nx1 - hNx1;
          jp = ((j-js) + pG->Disp[1] + hNx2) % Nx2 - hNx2;
          kx = (ip + 0.5*ioff)*dkx/pG->dx1;
          ky = (jp + 0.5*joff)*dky/pG->dx2;
          pcoeff = -0.5/sqrt(SQR(kx)+SQR(ky));
          work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0] *= pcoeff;
          work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1] *= pcoeff;
        }
      }

      /* STEP 3:  Backward FFT and set potential in real space */
      ath_2d_fft(bplan2d, work);

      /* Multiply by complex offsets and add real part to Phi.  To compute
       * offsets, note that indices relative to whole Domain are needed. */
      for (j=js; j<=je; j++) {
        for (i=is; i<=ie; i++) {
          offset = 0.5*((i-is+pG->Disp[0])*ioff*dkx
                      + (j-js+pG->Disp[1])*joff*dky);
          pG->Phi[ks][j][i] +=
            cos(offset)*work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][0]
          - sin(offset)*work[F2DI(i-is,j-js,pG->Nx[0],pG->Nx[1])][1];
        }
      }
    }
  }

  /* Finally, normalize the transforms. NOTE: There is a 4.0 here because
   * formally the transforms are performed over the extended domain of size
   * (2Nx)(2Ny). */
  pcoeff = four_pi_G/(4.0*bplan2d->gcnt);
  for (j=js; j<=je; j++) {
    for (i=is; i<=ie; i++) {
      pG->Phi[ks][j][i] *= pcoeff;
    }
  }

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_obc_3d(DomainS *pD)
 *  \brief Only works for uniform grid, periodic boundary conditions
 */
void selfg_fft_obc_3d(DomainS *pD)
{
  GridS *pG = (pD->Grid);
  int i,ioff,ip, is = pG->is, ie = pG->ie;
  int j,joff,jp, js = pG->js, je = pG->je;
  int k,koff,kp, ks = pG->ks, ke = pG->ke;
  Real pcoeff,offset;
  int Nx1=pD->Nx[0], Nx2=pD->Nx[1], Nx3=pD->Nx[2];
  int hNx1=Nx1/2, hNx2=Nx2/2, hNx3=Nx3/2;
  Real idx1sq=1.0/SQR(pG->dx1),idx2sq=1.0/SQR(pG->dx2),idx3sq=1.0/SQR(pG->dx3);
  Real dkx=2.0*PI/(double)(Nx1),dky=2.0*PI/(double)(Nx2),dkz=2.0*PI/(double)(Nx3);
#ifdef RADIATION
  RadParticleS *pRP=NULL;
  Real x1,x2,x3;
  Real tmp,W1[3],W2[3],W3[3],d;
#endif

  /* Copy current potential into old and zero-out */
  for (k=ks-nghost; k<=ke+nghost; k++){
    for (j=js-nghost; j<=je+nghost; j++){
      for (i=is-nghost; i<=ie+nghost; i++){
        pG->Phi_old[k][j][i] = pG->Phi[k][j][i];
        pG->Phi[k][j][i] = 0.0;
      }
    }
  }

  /* Loop over the index offsets (0 for even, 1 for odd) */
  for (koff=0; koff<=1; koff++) {
  for (joff=0; joff<=1; joff++) {
    for (ioff=0; ioff<=1; ioff++) {

      /* STEP 1: Forward FFT of 4\piG*(d-d0) */
      /* Add gas density into work array 0. */
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = pG->U[k][j][i].d;
          }
        }
      }

#ifdef RADIATION
      /* Interpolate particle density onto the mesh using CIC or TSC shape
       * functions as described in Hockney & Eastwood.  NOTE:  There are more
       * efficient ways to implement this, but since the radiation particle lists
       * will typically be small, there is no real need. */
      pRP = pG->radparticles;
      while (pRP) {
        /* Calculate indices and spatial coordinates of the lower nearest
         * neighbor grid cell.
         * x1 <= pRP->x1 < x1+pG->dx1, x2 <= pRP->x2 < x2+pG->dx2 */
        cc_ijk(pG,pRP->x1,pRP->x2,pRP->x3,&ip,&jp,&kp);
        cc_pos(pG,ip,jp,kp,&x1,&x2,&x3);

#ifdef PM_CIC
        /* Compute CIC weights.   */
        if (pRP->x1 > x1)
          W1[1] = (pRP->x1 - x1)/pG->dx1;  W1[0] = 1.0 - W1[1];
        else
          W1[0] = (x1 - pRP->x1)/pG->dx1;  W1[1] = 1.0 - W1[0];
        if (pRP->x2 > x2)
          W2[1] = (pRP->x2 - x2)/pG->dx2;  W2[0] = 1.0 - W2[1];
        else
          W2[0] = (x2 - pRP->x2)/pG->dx2;  W2[1] = 1.0 - W2[0];
        if (pRP->x3 > x3)
          W3[1] = (pRP->x3 - x3)/pG->dx3;  W3[0] = 1.0 - W3[1];
        else
          W3[0] = (x3 - pRP->x3)/pG->dx3;  W3[1] = 1.0 - W3[0];

        /* Add weighted particle density to work array 0. */
        d = pRP->m/(pG->dx1*pG->dx2*pG->dx3);
        for (k=0; k<=1; k++) {
          for (j=0; j<=1; j++) {
            for (i=0; i<=1; i++) {
              work[F3DI(ip+i-is,jp+j-js,kp+k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] += d*W1[i]*W2[j]*W3[k];
            }
          }
        }

#elif defined(PM_TSC)
        /* Compute TSC weights.   */
        tmp = 0.5 + (pRP->x1 - x1)/pG->dx1;
        W1[2] = 0.5*SQR(tmp);  W1[0] = W1[2] - tmp + 0.5;  W1[1] = 1.0 - W1[0] - W1[2];
        tmp = 0.5 + (pRP->x2 - x2)/pG->dx2;
        W2[2] = 0.5*SQR(tmp);  W2[0] = W2[2] - tmp + 0.5;  W2[1] = 1.0 - W2[0] - W2[2];
        tmp = 0.5 + (pRP->x3 - x3)/pG->dx3;
        W3[2] = 0.5*SQR(tmp);  W3[0] = W3[2] - tmp + 0.5;  W3[1] = 1.0 - W3[0] - W3[2];

        /* Add weighted particle density to work array 0. */
        d = pRP->m/(pG->dx1*pG->dx2*pG->dx3);
        for (k=0; k<=2; k++) {
          for (j=0; j<=2; j++) {
            for (i=0; i<=2; i++) {
              work[F3DI(ip-1+i-is,jp-1+j-js,kp-1+k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] += d*W1[i]*W2[j]*W3[k];
            }
          }
        }
#endif

        pRP = pRP->next;
      }
#endif

      /* Copy work array 0 into work array 1, then multiply by complex offsets.
       * To compute offsets, note that indices relative to whole Domain are
       * needed. */
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] =
              work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0];

            offset = 0.5*((i-is+pG->Disp[0])*ioff*dkx
                        + (j-js+pG->Disp[1])*joff*dky
                        + (k-ks+pG->Disp[2])*koff*dkz);
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *=  cos(offset);
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= -sin(offset);
          }
        }
      }

      /* Forward FFT */
      ath_3d_fft(fplan3d, work);

      /* STEP 2:  Compute potential in Fourier space.  Multiple loops are used
       * to avoid divide by zero at i=is,j=js,k=ks, and to avoid if statement in
       * loop.  To compute kx,ky,kz, note that indices relative to whole Domain
       * are needed. */
      if ((pG->Disp[0])==0 && (pG->Disp[1])==0 && (pG->Disp[2])==0
           && ioff==0 && joff==0 && koff==0) {
        work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] = 0.0;
        work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] = 0.0;
      } else {
        pcoeff = -0.5/((1.0-cos((pG->Disp[0] + 0.5*ioff)*dkx))*idx1sq +
                       (1.0-cos((pG->Disp[1] + 0.5*joff)*dky))*idx2sq +
                       (1.0-cos((pG->Disp[2] + 0.5*koff)*dkz))*idx3sq);
        work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
        work[F3DI(0,0,0,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
      }

      for (k=ks+1; k<=ke; k++) {
        pcoeff = -0.5/((1.0-cos((         pG->Disp[0] + 0.5*ioff)*dkx))*idx1sq +
                       (1.0-cos((         pG->Disp[1] + 0.5*joff)*dky))*idx2sq +
                       (1.0-cos(((k-ks) + pG->Disp[2] + 0.5*koff)*dkz))*idx3sq);
        work[F3DI(0,0,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
        work[F3DI(0,0,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
      }

      for (j=js+1; j<=je; j++) {
        for (k=ks; k<=ke; k++) {
          pcoeff = -0.5/((1.0-cos((         pG->Disp[0] + 0.5*ioff)*dkx))*idx1sq +
                         (1.0-cos(((j-js) + pG->Disp[1] + 0.5*joff)*dky))*idx2sq +
                         (1.0-cos(((k-ks) + pG->Disp[2] + 0.5*koff)*dkz))*idx3sq);
          work[F3DI(0,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
          work[F3DI(0,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
        }
      }

      for (i=is+1; i<=ie; i++) {
        for (j=js; j<=je; j++) {
          for (k=ks; k<=ke; k++) {
            pcoeff = -0.5/((1.0-cos(((i-is) + pG->Disp[0] + 0.5*ioff)*dkx))*idx1sq +
                           (1.0-cos(((j-js) + pG->Disp[1] + 0.5*joff)*dky))*idx2sq +
                           (1.0-cos(((k-ks) + pG->Disp[2] + 0.5*koff)*dkz))*idx3sq);
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0] *= pcoeff;
            work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1] *= pcoeff;
          }
        }
      }


      /* STEP 3:  Backward FFT and set potential in real space */
      ath_3d_fft(bplan3d, work);

      /* Multiply by complex offsets and add real part to Phi.  To compute
       * offsets, note that indices relative to whole Domain are needed. */
      for (k=ks; k<=ke; k++) {
        for (j=js; j<=je; j++) {
          for (i=is; i<=ie; i++) {
            offset = 0.5*((i-is+pG->Disp[0])*ioff*dkx
                        + (j-js+pG->Disp[1])*joff*dky
                        + (k-ks+pG->Disp[2])*koff*dkz);
            pG->Phi[k][j][i] +=
              cos(offset)*work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][0]
            - sin(offset)*work[F3DI(i-is,j-js,k-ks,pG->Nx[0],pG->Nx[1],pG->Nx[2])][1];
          }
        }
      }
    }
  }}

  /* Finally, normalize the transforms. NOTE: There is an 8.0 here because
   * formally the transforms are performed over the extended domain of size
   * (2Nx)(2Ny)(2Nz). */
  pcoeff = four_pi_G/(8.0*bplan3d->gcnt);
  for (k=ks; k<=ke; k++){
    for (j=js; j<=je; j++){
      for (i=is; i<=ie; i++){
        pG->Phi[k][j][i] *= pcoeff;
      }
    }
  }

  return;
}


/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_2d_init(MeshS *pM)
 *  \brief Initializes plans for forward/backward FFTs, and allocates memory 
 *   needed by FFTW.
 */
void selfg_fft_2d_init(MeshS *pM)
{
  DomainS *pD;
  int nl,nd;
  for (nl=0; nl<(pM->NLevels); nl++) {
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        fplan2d = ath_2d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
        bplan2d = ath_2d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
        work = ath_2d_fft_malloc(fplan2d);
      }
    }
  }
}


/*----------------------------------------------------------------------------*/
/*! \fn void selfg_fft_3d_init(MeshS *pM)
 *  \brief Initializes plans for forward/backward FFTs, and allocates memory 
 *   needed by FFTW.
 */
void selfg_fft_3d_init(MeshS *pM)
{
  DomainS *pD;
  int nl,nd;
  for (nl=0; nl<(pM->NLevels); nl++) {
    for (nd=0; nd<(pM->DomainsPerLevel[nl]); nd++) {
      if (pM->Domain[nl][nd].Grid != NULL) {
        pD = (DomainS*)&(pM->Domain[nl][nd]);
        fplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_FORWARD);
        bplan3d = ath_3d_fft_quick_plan(pD, NULL, ATH_FFT_BACKWARD);
        work = ath_3d_fft_malloc(fplan3d);
      }
    }
  }
}

#endif /* SELF_GRAVITY_USING_FFT_OBC */
