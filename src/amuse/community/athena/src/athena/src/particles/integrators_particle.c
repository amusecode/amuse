#include "../copyright.h"
/*===========================================================================*/
/*! \file integrators_particle.c
 *  \brief Provide three kinds of particle integrators.
 *
 * PURPOSE: provide three kinds of particle integrators, namely, 2nd order
 *   explicit, 2nd order semi-implicit and 2nd order fully implicit.
 * 
 * CONTAINS PUBLIC FUNCTIONS:
 * - Integrate_Particles();
 * - int_par_exp   ()
 * - int_par_semimp()
 * - int_par_fulimp()
 * - feedback_predictor()
 * - feedback_corrector()
 *
 * PRIVATE FUNCTION PROTOTYPES:
 * - Delete_Ghost()   - delete ghost particles
 * - JudgeCrossing()  - judge if the particle cross the grid boundary
 * - Get_Drag()       - calculate the drag force
 * - Get_Force()      - calculate forces other than the drag
 * - Get_Term()       - calculate the termination particle velocity
 * - Get_ForceDiff()  - calculate the force difference between particle and gas
 * 
 * History:
 * - Written by Xuening Bai, Mar.2009					      */
/*============================================================================*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../defs.h"
#include "../athena.h"
#include "../prototypes.h"
#include "prototypes.h"
#include "particle.h"
#include "../globals.h"

#ifdef PARTICLES         /* endif at the end of the file */

/*==============================================================================
 * PRIVATE FUNCTION PROTOTYPES:
 *   Delete_Ghost()   - delete ghost particles
 *   JudgeCrossing()  - judge if the particle cross the grid boundary
 *   Get_Drag()       - calculate the drag force
 *   Get_Force()      - calculate forces other than the drag
 *   Get_Term()       - calculate the termination particle velocity
 *   Get_ForceDiff()  - calculate the force difference between particle and gas
 *============================================================================*/
void   Delete_Ghost(Grid *pG);
void   JudgeCrossing(Grid *pG, Real x1, Real x2, Real x3, Grain *gr);
Vector Get_Drag(Grid *pG, int type, Real x1, Real x2, Real x3,
                Real v1, Real v2, Real v3, Vector cell1, Real *tstop1);
Vector Get_Force(Grid *pG, Real x1, Real x2, Real x3,
                           Real v1, Real v2, Real v3);
Vector Get_Term(Grid *pG, int type, Real x1, Real x2, Real x3, Vector cell1,
                                                      Real *tstop);
Vector Get_ForceDiff(Grid *pG, Real x1, Real x2, Real x3,
                               Real v1, Real v2, Real v3);

/*=========================== PUBLIC FUNCTIONS ===============================*/
/*----------------------------------------------------------------------------*/

/*---------------------------------- Main Integrator -------------------------*/
/*! \fn void Integrate_Particles(Grid *pG, Domain *pD)
 *  \brief Main particle integrator.
 *
 * Input: Grid which is already evolved in half a time step. Paricles unevolved.
 * Output: particle updated for one full time step; feedback array for corrector
 *         updated.
 * Note: This routine allows has the flexibility to freely choose the particle
 *       integrator.
 * Should use fully implicit integrator for tightly coupoled particles.
 * Otherwise the semi-implicit integrator performs better.
 */
void Integrate_Particles(Grid *pG, Domain *pD)
{
  Grain *curG, *curP, mygr;     /* pointer of the current working position */
  long p;                       /* particle index */
  Real dv1, dv2, dv3, ts, t1;   /* amount of velocity update, stopping time */
  Vector cell1;                 /* one over dx1, dx2, dx3 */
  Vector vterm;                 /* termination velocity */

/* Initialization */
#ifdef FEEDBACK
  feedback_clear(pG);   /* clean the feedback array */
#endif /* FEEDBACK */

  curP = &(mygr);       /* temperory particle */

  /* cell1 is a shortcut expressions as well as dimension indicator */
  if (pG->Nx1 > 1)  cell1.x1 = 1.0/pG->dx1;  else cell1.x1 = 0.0;
  if (pG->Nx2 > 1)  cell1.x2 = 1.0/pG->dx2;  else cell1.x2 = 0.0;
  if (pG->Nx3 > 1)  cell1.x3 = 1.0/pG->dx3;  else cell1.x3 = 0.0;

  /* delete all ghost particles */
  Delete_Ghost(pG);

  p = 0;
  while (p<pG->nparticle)
  {/* loop over all particles */
    curG = &(pG->particle[p]);

/* Step 1: Calculate velocity update */
    switch(pG->grproperty[curG->property].integrator)
    {
      case 1: /* 2nd order explicit integrator */
        int_par_exp(pG,curG,cell1, &dv1,&dv2,&dv3, &ts);
        break;

      case 2: /* 2nd order semi-implicit integrator */
        int_par_semimp(pG,curG,cell1, &dv1,&dv2,&dv3, &ts);
        break;

      case 3: /* 2nd order fully implicit integrator */
        int_par_fulimp(pG,curG,cell1, &dv1,&dv2,&dv3, &ts);
        break;

      case 4: /* specific integrator in the strong coupling limit */
        int_par_spec(pG,curG,cell1, &dv1,&dv2,&dv3, &ts);
        break;

      default:
        ath_error("[integrate_particle]: unknown integrator type!");
    }

/* Step 2: particle update to curP */

    /* velocity update */
    curP->v1 = curG->v1 + dv1;
    curP->v2 = curG->v2 + dv2;
    curP->v3 = curG->v3 + dv3;

    /* position update */
    if (pG->Nx1 > 1)
      curP->x1 = curG->x1 + 0.5*pG->dt*(curG->v1 + curP->v1);
    else /* do not move if this dimension collapses */
      curP->x1 = curG->x1;

    if (pG->Nx2 > 1)
      curP->x2 = curG->x2 + 0.5*pG->dt*(curG->v2 + curP->v2);
    else /* do not move if this dimension collapses */
      curP->x2 = curG->x2;

   if (pG->Nx3 > 1)
      curP->x3 = curG->x3 + 0.5*pG->dt*(curG->v3 + curP->v3);
    else /* do not move if this dimension collapses */
      curP->x3 = curG->x3;

#ifdef FARGO
    /* shift = -qshear * Omega_0 * x * dt */
    pG->parsub[p].shift = -0.5*qshear*Omega_0*(curG->x1+curP->x1)*pG->dt;
#endif

    /* special treatment for integrator #4 */
    if (pG->grproperty[curG->property].integrator == 4)
    {
       vterm = Get_Term(pG,curG->property,curP->x1,curP->x2,curP->x3,cell1,&t1);
       curP->v1 = vterm.x1;     curP->v2 = vterm.x2;     curP->v2 = vterm.x2;

       dv1=curP->v1-curG->v1;   dv2=curP->v2-curG->v2;   dv3=curP->v3-curG->v3;
    }

/* Step 3: calculate feedback force to the gas */
#ifdef FEEDBACK
    feedback_corrector(pG, curG, curP, cell1, dv1, dv2, dv3, ts);
#endif /* FEEDBACK */

/* Step 4: Final update of the particle */
    /* update particle status (crossing boundary or not) */
    JudgeCrossing(pG, curP->x1, curP->x2, curP->x3, curG);

    /* update the particle */
    curG->x1 = curP->x1;
    curG->x2 = curP->x2;
    curG->x3 = curP->x3;
    curG->v1 = curP->v1;
    curG->v2 = curP->v2;
    curG->v3 = curP->v3;
    p++;

  } /* end of the for loop */

/* output the status */
  ath_pout(0, "In processor %d, there are %ld particles.\n",
                             pG->my_id, pG->nparticle);

  return;
}

/* ------------ 2nd order fully implicit particle integrator -----------------*/
/*! \fn void int_par_fulimp(Grid *pG, Grain *curG, Vector cell1, 
 *                            Real *dv1, Real *dv2, Real *dv3, Real *ts)
 *  \brief 2nd order fully implicit particle integrator
 *
 * Input: 
 *   grid pointer (pG), grain pointer (curG), cell size indicator (cell1)
 * Output:
 *   dv1,dv2,dv3: velocity update
 */
void int_par_fulimp(Grid *pG, Grain *curG, Vector cell1, 
                              Real *dv1, Real *dv2, Real *dv3, Real *ts)
{
  Real x1n, x2n, x3n;	/* first order new position at half a time step */
  Vector fd, fr;	/* drag force and other forces */
  Vector fc, fp, ft;	/* force at current & predicted position, total force */
  Real ts11, ts12;	/* 1/stopping time */
  Real b0,A,B,C,D,Det1;	/* matrix elements and determinant */
#ifdef SHEARING_BOX
  Real oh, oh2;		/* Omega_0*dt and its square */
#endif

/* step 1: predict of the particle position after one time step */
  if (pG->Nx1 > 1)  x1n = curG->x1+curG->v1*pG->dt;
  else x1n = curG->x1;
  if (pG->Nx2 > 1)  x2n = curG->x2+curG->v2*pG->dt;
  else x2n = curG->x2;
  if (pG->Nx3 > 1)  x3n = curG->x3+curG->v3*pG->dt;
  else x3n = curG->x3;

#ifdef SHEARING_BOX
#ifndef FARGO
  if (pG->Nx3 > 1) x2n -= 0.5*qshear*curG->v1*SQR(pG->dt); /* advection part */
#endif
#endif

/* step 2: calculate the force at current position */
  fd = Get_Drag(pG, curG->property, curG->x1, curG->x2, curG->x3,
                                    curG->v1, curG->v2, curG->v3, cell1, &ts11);

  fr = Get_Force(pG, curG->x1, curG->x2, curG->x3,
                     curG->v1, curG->v2, curG->v3);

  fc.x1 = fd.x1+fr.x1;
  fc.x2 = fd.x2+fr.x2;
  fc.x3 = fd.x3+fr.x3;

/* step 3: calculate the force at the predicted positoin */
  fd = Get_Drag(pG, curG->property, x1n, x2n, x3n,
                                    curG->v1, curG->v2, curG->v3, cell1, &ts12);

  fr = Get_Force(pG, x1n, x2n, x3n, curG->v1, curG->v2, curG->v3);

  fp.x1 = fd.x1+fr.x1;
  fp.x2 = fd.x2+fr.x2;
  fp.x3 = fd.x3+fr.x3;

/* step 4: calculate the velocity update */
  /* shortcut expressions */
  b0 = 1.0+pG->dt*ts11;

  /* Total force */
  ft.x1 = 0.5*(fc.x1+b0*fp.x1);
  ft.x2 = 0.5*(fc.x2+b0*fp.x2);
  ft.x3 = 0.5*(fc.x3+b0*fp.x3);

#ifdef SHEARING_BOX
  oh = Omega_0*pG->dt;
  if (pG->Nx3 > 1) {/* 3D shearing sheet (x1,x2,x3)=(X,Y,Z) */
    ft.x1 += -oh*fp.x2;
  #ifdef FARGO
    ft.x2 += 0.5*(2.0-qshear)*oh*fp.x1;
  #else
    ft.x2 += oh*fp.x1;
  #endif
  } else {         /* 2D shearing sheet (x1,x2,x3)=(X,Z,Y) */
    ft.x1 += -oh*fp.x3;
  #ifdef FARGO
    ft.x3 += 0.5*(2.0-qshear)*oh*fp.x1;
  #else
    ft.x3 += oh*fp.x1;
  #endif
  }
#endif /* SHEARING_BOX */

  /* calculate the inverse matrix elements */
  D = 1.0+0.5*pG->dt*(ts11 + ts12 + pG->dt*ts11*ts12);
#ifdef SHEARING_BOX
  oh2 = SQR(oh);
  B = oh * (-2.0-(ts11+ts12)*pG->dt);
#ifdef FARGO
  A = D - (2.0-qshear)*oh2;
  C = 0.5*(qshear-2.0)*B;
#else /* FARGO */
  A = D - 2.0*oh2;
  C = -B;
#endif /* FARGO */
  Det1 = 1.0/(SQR(A)-B*C);
  if (pG->Nx3>1) {
    *dv1 = pG->dt*Det1*(ft.x1*A-ft.x2*B);
    *dv2 = pG->dt*Det1*(-ft.x1*C+ft.x2*A);
    *dv3 = pG->dt*ft.x3/D;
  } else {
    *dv1 = pG->dt*Det1*(ft.x1*A-ft.x3*B);
    *dv3 = pG->dt*Det1*(-ft.x1*C+ft.x3*A);
    *dv2 = pG->dt*ft.x2/D;
  }
#else /* SHEARING_BOX */
  D = 1.0/D;
  *dv1 = pG->dt*ft.x1*D;
  *dv2 = pG->dt*ft.x2*D;
  *dv3 = pG->dt*ft.x3*D;
#endif /* SHEARING_BOX */

  *ts = 0.5/ts11+0.5/ts12;

  return;
}


/*--------------- 2nd order semi-implicit particle integrator ----------------*/
/*! \fn void int_par_semimp(Grid *pG, Grain *curG, Vector cell1, 
 *                            Real *dv1, Real *dv2, Real *dv3, Real *ts)
 *  \brief 2nd order semi-implicit particle integrator 
 *
 * Input: 
 *   grid pointer (pG), grain pointer (curG), cell size indicator (cell1)
 * Output:
 *   dv1,dv2,dv3: velocity update
 */
void int_par_semimp(Grid *pG, Grain *curG, Vector cell1, 
                              Real *dv1, Real *dv2, Real *dv3, Real *ts)
{
  Vector fd, fr, ft;	/* drag force and other forces, total force */
  Real ts1, b, b2;	/* other shortcut expressions */
  Real x1n, x2n, x3n;	/* first order new position at half a time step */
#ifdef SHEARING_BOX
  Real b1, oh;		/* Omega_0*h */
#endif

/* step 1: predict of the particle position after half a time step */
  if (pG->Nx1 > 1)  x1n = curG->x1+0.5*curG->v1*pG->dt;
  else x1n = curG->x1;
  if (pG->Nx2 > 1)  x2n = curG->x2+0.5*curG->v2*pG->dt;
  else x2n = curG->x2;
  if (pG->Nx3 > 1)  x3n = curG->x3+0.5*curG->v3*pG->dt;
  else x3n = curG->x3;

#ifdef SHEARING_BOX
#ifndef FARGO
  if (pG->Nx3 > 1) x2n -= 0.125*qshear*curG->v1*SQR(pG->dt);/* advection part */
#endif
#endif

/* Step 2: interpolation to get fluid density, velocity and the sound speed at\  * predicted position
 */
  fd = Get_Drag(pG, curG->property, x1n, x2n, x3n,
                                    curG->v1, curG->v2, curG->v3, cell1, &ts1);

  fr = Get_Force(pG, x1n, x2n, x3n, curG->v1, curG->v2, curG->v3);

  ft.x1 = fd.x1+fr.x1;
  ft.x2 = fd.x2+fr.x2;
  ft.x3 = fd.x3+fr.x3;

/* step 3: calculate velocity update */

  /* shortcut expressions */
  b = pG->dt*ts1+2.0;
#ifdef SHEARING_BOX
  oh = Omega_0*pG->dt;
#ifdef FARGO
  b1 = 1.0/(SQR(b)+2.0*(2.0-qshear)*SQR(oh));
#else
  b1 = 1.0/(SQR(b)+4.0*SQR(oh));
#endif /* FARGO */
  b2 = b*b1;
#else
  b2 = 1.0/b;
#endif /* SHEARING BOX */

    /* velocity evolution */
#ifdef SHEARING_BOX
  if (pG->Nx3>1)
  {/* 3D shearing sheet (x1,x2,x3)=(X,Y,Z) */
    *dv1 = pG->dt*2.0*b2*ft.x1 + pG->dt*4.0*oh*b1*ft.x2;
  #ifdef FARGO
    *dv2 = pG->dt*2.0*b2*ft.x2 - 2.0*(2.0-qshear)*pG->dt*oh*b1*ft.x1;
  #else
    *dv2 = pG->dt*2.0*b2*ft.x2 - 4.0*pG->dt*oh*b1*ft.x1;
  #endif /* FARGO */
    *dv3 = pG->dt*2.0*ft.x3/b;
  }
  else
  {/* 2D shearing sheet (x1,x2,x3)=(X,Z,Y) */
    *dv1 = pG->dt*2.0*b2*ft.x1 + pG->dt*4.0*oh*b1*ft.x3;
    *dv2 = pG->dt*2.0*ft.x2/b;
  #ifdef FARGO
    *dv3 = pG->dt*2.0*b2*ft.x3 - 2.0*(2.0-qshear)*pG->dt*oh*b1*ft.x1;
  #else
    *dv3 = pG->dt*2.0*b2*ft.x3 - 4.0*pG->dt*oh*b1*ft.x1;
  #endif
  }
#else
  *dv1 = pG->dt*2.0*b2*ft.x1;
  *dv2 = pG->dt*2.0*b2*ft.x2;
  *dv3 = pG->dt*2.0*b2*ft.x3;
#endif /* SHEARING_BOX */

  *ts = 1.0/ts1;

  return;
}


/*------------------- 2nd order explicit particle integrator -----------------*/
/*! \fn void int_par_exp(Grid *pG, Grain *curG, Vector cell1,
 *                         Real *dv1, Real *dv2, Real *dv3, Real *ts)
 *  \brief 2nd order explicit particle integrator 
 *
 * Input: 
 *   grid pointer (pG), grain pointer (curG), cell size indicator (cell1)
 * Output:
 *   dv1,dv2,dv3: velocity update
 */
void int_par_exp(Grid *pG, Grain *curG, Vector cell1,
                           Real *dv1, Real *dv2, Real *dv3, Real *ts)
{
  Vector fd, fr, ft;	/* drag force and other forces, total force */
  Real ts1;		/* 1/stopping time */
  Real x1n, x2n, x3n;	/* first order new position at half a time step */
  Real v1n, v2n, v3n;	/* first order new velocity at half a time step */

/* step 1: predict of the particle position after half a time step */
  if (pG->Nx1 > 1)
    x1n = curG->x1+0.5*curG->v1*pG->dt;
  else x1n = curG->x1;
  if (pG->Nx2 > 1)
    x2n = curG->x2+0.5*curG->v2*pG->dt;
  else x2n = curG->x2;
  if (pG->Nx3 > 1)
    x3n = curG->x3+0.5*curG->v3*pG->dt;
  else x3n = curG->x3;

#ifdef SHEARING_BOX
#ifndef FARGO
  if (pG->Nx3 > 1) x2n -= 0.125*qshear*curG->v1*SQR(pG->dt);/* advection part */
#endif
#endif

/* step 2: calculate the force at current position */
  fd = Get_Drag(pG, curG->property, curG->x1, curG->x2, curG->x3,
                                    curG->v1, curG->v2, curG->v3, cell1, &ts1);

  fr = Get_Force(pG, curG->x1, curG->x2, curG->x3,
                     curG->v1, curG->v2, curG->v3);

  ft.x1 = fd.x1+fr.x1;
  ft.x2 = fd.x2+fr.x2;
  ft.x3 = fd.x3+fr.x3;

  v1n = curG->v1 + 0.5*ft.x1*pG->dt;
  v2n = curG->v2 + 0.5*ft.x2*pG->dt;
  v3n = curG->v3 + 0.5*ft.x3*pG->dt;

/* step 3: calculate the force at the predicted positoin */
  fd = Get_Drag(pG, curG->property, x1n, x2n, x3n, v1n, v2n, v3n, cell1, &ts1);

  fr = Get_Force(pG, x1n, x2n, x3n, v1n, v2n, v3n);

  ft.x1 = fd.x1+fr.x1;
  ft.x2 = fd.x2+fr.x2;
  ft.x3 = fd.x3+fr.x3;

/* step 4: calculate velocity update */
  *dv1 = ft.x1*pG->dt;
  *dv2 = ft.x2*pG->dt;
  *dv3 = ft.x3*pG->dt;

  *ts = 1.0/ts1;

  return;
}

/*------------------- 2nd order specific particle integrator -----------------*/
/*! \fn void int_par_spec(Grid *pG, Grain *curG, Vector cell1,
 *                          Real *dv1, Real *dv2, Real *dv3, Real *ts)
 *  \brief 2nd order specific particle integrator;
 *  This integrator works ONLY in the strong coupling regime (t_stop<h)
 *
 * Input:
 *   grid pointer (pG), grain pointer (curG), cell size indicator (cell1)
 * Output:
 *   dv1,dv2,dv3: velocity update
 */
void int_par_spec(Grid *pG, Grain *curG, Vector cell1,
                            Real *dv1, Real *dv2, Real *dv3, Real *ts)
{
  Vector vterm;         /* termination velocity */
  Real x1n, x2n, x3n;   /* predicted position at half a time step */

/* step 1: predict of the particle position after half a time step */
  if (pG->Nx1 > 1)
    x1n = curG->x1+0.5*curG->v1*pG->dt;
  else x1n = curG->x1;
  if (pG->Nx2 > 1)
    x2n = curG->x2+0.5*curG->v2*pG->dt;
  else x2n = curG->x2;
  if (pG->Nx3 > 1)
    x3n = curG->x3+0.5*curG->v3*pG->dt;
  else x3n = curG->x3;

#ifdef SHEARING_BOX
#ifndef FARGO
  if (pG->Nx3 > 1) x2n -= 0.125*qshear*curG->v1*SQR(pG->dt);/* advection part */
#endif
#endif

/* step 2: get gas termination velocity */
  vterm = Get_Term(pG, curG->property, x1n, x2n, x3n, cell1, ts);

/* step 3: calculate the velocity difference */
  *dv1 = 2.0*(vterm.x1 - curG->v1);
  *dv2 = 2.0*(vterm.x2 - curG->v2);
  *dv3 = 2.0*(vterm.x3 - curG->v3);

  return;
}

#ifdef FEEDBACK

/*! \fn void feedback_predictor(Grid* pG)
 *  \brief Calculate the feedback of the drag force from the particle to the gas
 *
 * Serves for the predictor step. It deals with all the particles.
 * Input: pG: grid with particles
 * Output: pG: the array of drag forces exerted by the particle is updated
*/
void feedback_predictor(Grid* pG)
{
  int is,js,ks,i,j,k;
  long p;                   /* particle index */
  Real weight[3][3][3];     /* weight function */
  Real rho, cs, tstop;      /* density, sound speed, stopping time */
  Real u1, u2, u3;
  Real vd1, vd2, vd3, vd;   /* velocity difference between particle and gas */
  Real f1, f2, f3;          /* feedback force */
  Real m, ts1h;             /* grain mass, 0.5*dt/tstop */
  Vector cell1;             /* one over dx1, dx2, dx3 */
  Vector fb;                /* drag force, fluid velocity */
#ifndef BAROTROPIC
  Real Elosspar;            /* energy dissipation rate due to drag */
#endif
  Real stiffness;           /* stiffness parameter of feedback */
  Grain *cur;               /* pointer of the current working position */

  /* initialization */
  get_gasinfo(pG);		/* calculate gas information */

  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++) {
        /* clean the feedback array */
        pG->Coup[k][j][i].fb1 = 0.0;
        pG->Coup[k][j][i].fb2 = 0.0;
        pG->Coup[k][j][i].fb3 = 0.0;
#ifndef BAROTROPIC
        pG->Coup[k][j][i].Eloss = 0.0;
#endif
        pG->Coup[k][j][i].FBstiff = 0.0;
      }

  /* convenient expressions */
  if (pG->Nx1 > 1)  cell1.x1 = 1.0/pG->dx1;
  else              cell1.x1 = 0.0;

  if (pG->Nx2 > 1)  cell1.x2 = 1.0/pG->dx2;
  else              cell1.x2 = 0.0;

  if (pG->Nx3 > 1)  cell1.x3 = 1.0/pG->dx3;
  else              cell1.x3 = 0.0;

  /* loop over all particles to calculate the drag force */
  for (p=0; p<pG->nparticle; p++)
  {/* loop over all particle */
    cur = &(pG->particle[p]);

    /* interpolation to get fluid density and velocity */
    getweight(pG, cur->x1, cur->x2, cur->x3, cell1, weight, &is, &js, &ks);
    if (getvalues(pG, weight, is, js, ks,
                              &rho, &u1, &u2, &u3, &cs, &stiffness) == 0)
    { /* particle is in the grid */

      /* apply gas velocity shift due to pressure gradient */
      gasvshift(cur->x1, cur->x2, cur->x3, &u1, &u2, &u3);
      /* velocity difference */
      vd1 = u1-cur->v1;
      vd2 = u2-cur->v2;
      vd3 = u3-cur->v3;
      vd = sqrt(vd1*vd1 + vd2*vd2 + vd3*vd3);

      /* calculate particle stopping time */
      tstop = get_ts(pG, cur->property, rho, cs, vd);
      ts1h = 0.5*pG->dt/tstop;

      /* Drag force density */
      m = pG->grproperty[cur->property].m;
      fb.x1 = m * vd1 * ts1h;
      fb.x2 = m * vd2 * ts1h;
      fb.x3 = m * vd3 * ts1h;

      /* calculate feedback stiffness */
       stiffness = 2.0*m*ts1h;

#ifndef BAROTROPIC
      Elosspar = fb.x1*vd1 + fb.x2*vd2 + fb.x3*vd3;
      /* distribute the drag force (density) to the grid */
      distrFB_pred(pG, weight, is, js, ks, fb, stiffness, Elosspar);
#else
      /* distribute the drag force (density) to the grid */
      distrFB_pred(pG, weight, is, js, ks, fb, stiffness);
#endif
    }
  }/* end of the for loop */

/* normalize stiffness and correct for feedback */
  for (k=klp; k<=kup; k++)
    for (j=jlp; j<=jup; j++)
      for (i=ilp; i<=iup; i++)
      {
        pG->Coup[k][j][i].FBstiff /= pG->U[k][j][i].d;

//        stiffness = 1.0/MAX(1.0, pG->Coup[k][j][i].FBstiff);

//        pG->Coup[k][j][i].fb1 *= stiffness;
//        pG->Coup[k][j][i].fb2 *= stiffness;
//        pG->Coup[k][j][i].fb3 *= stiffness;
      }

  return;
}

/*----------------------------------------------------------------------------*/
/*! \fn void feedback_corrector(Grid *pG, Grain *gri, Grain *grf, Vector cell1,
 *                                Real dv1, Real dv2, Real dv3, Real ts)
 *  \brief  Calculate the feedback of the drag force from the particle 
 *	    to the gas.
 *
 * Serves for the corrector step. It deals with one particle at a time.
 * Input: pG: grid with particles; gri,grf: initial and final particles;
 *        dv: velocity difference between gri and grf.
 * Output: pG: the array of drag forces exerted by the particle is updated
*/
void feedback_corrector(Grid *pG, Grain *gri, Grain *grf, Vector cell1,
                                  Real dv1, Real dv2, Real dv3, Real ts)
{
  int is, js, ks;
  Real x1, x2, x3, v1, v2, v3;
  Real mgr;
  Real weight[3][3][3];
  Real Elosspar;                        /* particle energy dissipation */
  Vector fb;

  mgr = pG->grproperty[gri->property].m;
  x1 = 0.5*(gri->x1+grf->x1);
  x2 = 0.5*(gri->x2+grf->x2);
  x3 = 0.5*(gri->x3+grf->x3);
  v1 = 0.5*(gri->v1+grf->v1);
  v2 = 0.5*(gri->v2+grf->v2);
  v3 = 0.5*(gri->v3+grf->v3);

  /* Force other than drag force */
  fb = Get_Force(pG, x1, x2, x3, v1, v2, v3);

  fb.x1 = dv1 - pG->dt*fb.x1;
  fb.x2 = dv2 - pG->dt*fb.x2;
  fb.x3 = dv3 - pG->dt*fb.x3;

  /* energy dissipation */
  Elosspar = mgr*(SQR(fb.x1)+SQR(fb.x2)+SQR(fb.x3))*ts;

  /* Drag force density */
  fb.x1 = mgr*fb.x1;
  fb.x2 = mgr*fb.x2;
  fb.x3 = mgr*fb.x3;

  /* distribute the drag force (density) to the grid */
  getweight(pG, x1, x2, x3, cell1, weight, &is, &js, &ks);
  distrFB_corr(pG, weight, is, js, ks, fb, Elosspar);

  return;

}

#endif /* FEEDBACK */


/*=========================== PRIVATE FUNCTIONS ==============================*/
/*----------------------------------------------------------------------------*/
/*! \fn void Delete_Ghost(Grid *pG)
 *  \brief Delete ghost particles */
void Delete_Ghost(Grid *pG)
{
  long p;
  Grain *gr;

  p = 0;
  while (p<pG->nparticle)
  {/* loop over all particles */
    gr = &(pG->particle[p]);

    if (gr->pos == 0)
    {/* gr is a ghost particle */
      pG->nparticle -= 1;
      pG->grproperty[gr->property].num -= 1;
      pG->particle[p] = pG->particle[pG->nparticle];
    }
    else
      p++;
  }

  return;
}

/*--------------------------------------------------------------------------- */
/*! \fn void JudgeCrossing(Grid *pG, Real x1, Real x2, Real x3, Grain *gr)
 *  \brief Judge if the particle is a crossing particle */
void JudgeCrossing(Grid *pG, Real x1, Real x2, Real x3, Grain *gr)
{
#ifndef FARGO
    /* if it crosses the grid boundary, mark it as a crossing out particle */
    if ((x1>=x1upar) || (x1<x1lpar) || (x2>=x2upar) || (x2<x2lpar) ||
                                       (x3>=x3upar) || (x3<x3lpar) )
      gr->pos = 10;
#else
    /* FARGO will naturally return the "crossing out" particles in the x2
        direction to the grid
     */
    if ((x1>=x1upar) || (x1<x1lpar) || (x3>=x3upar) || (x3<x3lpar))
      gr->pos = 10;
    else if ((pG->Nx3==1) && ((x2>=x2upar) || (x2<x2lpar))) /* 2D problem */
      gr->pos = 10;
#endif
    return;
}

/*--------------------------------------------------------------------------- */
/*! \fn Vector Get_Drag(Grid *pG, int type, Real x1, Real x2, Real x3,
 *              Real v1, Real v2, Real v3, Vector cell1, Real *tstop1)
 *  \brief Calculate the drag force to the particles 
 *
 * Input:
 *   pG: grid;	type: particle type;	cell1: 1/dx1,1/dx2,1/dx3;
 *   x1,x2,x3,v1,v2,v3: particle position and velocity;
 * Output:
 *   tstop1: 1/stopping time;
 * Return:
 *   drag force;
 */
Vector Get_Drag(Grid *pG, int type, Real x1, Real x2, Real x3,
                Real v1, Real v2, Real v3, Vector cell1, Real *tstop1)
{
  int is, js, ks;
  Real rho, u1, u2, u3, cs;
  Real vd1, vd2, vd3, vd, tstop, ts1;
#ifdef FEEDBACK
  Real stiffness;
#endif
  Real weight[3][3][3];		/* weight function */
  Vector fd;

  /* interpolation to get fluid density, velocity and the sound speed */
  getweight(pG, x1, x2, x3, cell1, weight, &is, &js, &ks);

#ifndef FEEDBACK
  if (getvalues(pG, weight, is, js, ks, &rho, &u1, &u2, &u3, &cs) == 0)
#else
  if (getvalues(pG, weight, is, js, ks, &rho,&u1,&u2,&u3,&cs, &stiffness) == 0)
#endif
  { /* particle in the grid */

    /* apply possible gas velocity shift (e.g., for fake gas velocity field) */
    gasvshift(x1, x2, x3, &u1, &u2, &u3);

    /* velocity difference */
    vd1 = v1-u1;
    vd2 = v2-u2;
    vd3 = v3-u3;
    vd = sqrt(SQR(vd1) + SQR(vd2) + SQR(vd3)); /* dimension independent */

    /* particle stopping time */
    tstop = get_ts(pG, type, rho, cs, vd);
#ifdef FEEDBACK
//    tstop *= MAX(1.0,stiffness);
#endif
    ts1 = 1.0/tstop;
  }
  else
  { /* particle out of the grid, free motion, with warning sign */
    vd1 = 0.0;	vd2 = 0.0;	vd3 = 0.0;	ts1 = 0.0;
    ath_perr(0, "Particle move out of grid %d with position (%f,%f,%f)!\n",
                                           pG->my_id,x1,x2,x3); /* warning! */
  }

  *tstop1 = ts1;

  /* Drag force */
  fd.x1 = -ts1*vd1;
  fd.x2 = -ts1*vd2;
  fd.x3 = -ts1*vd3;

  return fd;
}

/*--------------------------------------------------------------------------- */
/*! \fn Vector Get_Force(Grid *pG, Real x1, Real x2, Real x3,
 *                         Real v1, Real v2, Real v3)
 *  \brief Calculate the forces to the particle other than the gas drag
 *
 * Input:
 *   pG: grid;
 *   x1,x2,x3,v1,v2,v3: particle position and velocity;
 * Return:
 *   forces;
 */
Vector Get_Force(Grid *pG, Real x1, Real x2, Real x3,
                           Real v1, Real v2, Real v3)
{
  Vector ft;

  ft.x1 = ft.x2 = ft.x3 = 0.0;

/* User defined forces
 * Should be independent of velocity, or the integrators must be modified
 * Can also change the velocity to account for velocity difference.
 */
  Userforce_particle(&ft, x1, x2, x3, v1, v2, v3);

#ifdef SHEARING_BOX
  Real omg2 = SQR(Omega_0);

  if (pG->Nx3 > 1)
  {/* 3D shearing sheet (x1,x2,x3)=(X,Y,Z) */
  #ifdef FARGO
    ft.x1 += 2.0*v2*Omega_0;
    ft.x2 += (qshear-2.0)*v1*Omega_0;
  #else
    ft.x1 += 2.0*(qshear*omg2*x1 + v2*Omega_0);
    ft.x2 += -2.0*v1*Omega_0;
  #endif /* FARGO */
  }
  else
  { /* 2D shearing sheet (x1,x2,x3)=(X,Z,Y) */
  #ifdef FARGO
    ft.x1 += 2.0*v3*Omega_0;
    ft.x3 += (qshear-2.0)*v1*Omega_0;
  #else
    ft.x1 += 2.0*(qshear*omg2*x1 + v3*Omega_0);
    ft.x3 += -2.0*v1*Omega_0;
  #endif /* FARGO */
  }
#endif /* SHEARING_BOX */

  return ft;
}

/*! \fn Vector Get_Term(Grid *pG, int type, Real x1, Real x2, Real x3, 
 *                      Vector cell1, Real *tstop)
 *  \brief Calculate the termination velocity of strongly coupled particles
 *
 * Used for the special integrator
 * Force difference include pressure gradient and momentum feedback
 * Input:
 *   pG: grid;  type: particle type;  x1,x2,x3: particle position;
 * Return:
 *   termination velocity, and the stopping time.
 */
Vector Get_Term(Grid *pG, int type, Real x1, Real x2, Real x3, Vector cell1,
                                                       Real *tstop)
{
  Vector vterm;             /* termination velocity */
  Vector ft;                /* force difference */
  Real rho, u1, u2, u3, cs; /* gas velocity */
  int is, js, ks;
#ifdef FEEDBACK
  Real stiffness;
#endif
  Real weight[3][3][3];     /* weight function */

  /* interpolation to get fluid density, velocity and the sound speed */
  getweight(pG, x1, x2, x3, cell1, weight, &is, &js, &ks);

#ifndef FEEDBACK
  if (getvalues(pG, weight, is, js, ks, &rho, &u1, &u2, &u3, &cs) == 0)
#else
  if (getvalues(pG, weight, is, js, ks, &rho,&u1,&u2,&u3,&cs, &stiffness) == 0)
#endif
  { /* position in the grid */

    /* apply possible gas velocity shift (e.g., for fake gas velocity field) */
    gasvshift(x1, x2, x3, &u1, &u2, &u3);

    /* particle stopping time */
    *tstop = get_ts(pG, type, rho, cs, 0.0);

    /* force difference */
    ft = Get_ForceDiff(pG, x1, x2, x3, u1, u2, u3);

    /* termination velocity */
    vterm.x1 = u1 + *tstop*ft.x1;
    vterm.x2 = u2 + *tstop*ft.x2;
    vterm.x3 = u3 + *tstop*ft.x3;
  }
  else
  {
    vterm.x1 = 0.0;    vterm.x2 = 0.0;    vterm.x3 = 0.0;
    ath_perr(0, "[get_term]: Position (%f,%f,%f) is out of grid!\n",
                                           pG->my_id,x1,x2,x3); /* warning! */
  }

  return vterm;
}

/*! \fn Vector Get_ForceDiff(Grid *pG, Real x1, Real x2, Real x3,
 *                             Real v1, Real v2, Real v3)
 *  \brief Calculate the force (density) difference between particles and gas
 *
 * Used ONLY for the special integrator.
 * THIS ROUTINE MUST BE EDITTED BY THE USER!
 * Force differences due to gas pressure gradient and momentum feedback are
 * automatically included. The user must provide other sources.
 *
 * Input:
 *   pG: grid;
 *   x1,x2,x3,v1,v2,v3: particle position and velocity;
 * Return:
 *   forces;
 */
Vector Get_ForceDiff(Grid *pG, Real x1, Real x2, Real x3,
                               Real v1, Real v2, Real v3)
{
  Vector fd;

  fd.x1 = 0.0;    fd.x2 = 0.0;    fd.x3 = 0.0;

  Userforce_particle(&fd, x1, x2, x3, v1, v2, v3);

/*
  fd.x1 += x1;
  fd.x2 += x2;
  fd.x3 += x3;
*/

  return fd;
}

#endif /*PARTICLES*/
