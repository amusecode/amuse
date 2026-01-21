/*!
 * \copyright   This file is part of the public version of the AREPO code.
 * \copyright   Copyright (C) 2009-2019, Max-Planck Institute for Astrophysics
 * \copyright   Developed by Volker Springel (vspringel@MPA-Garching.MPG.DE) and
 *              contributing authors.
 * \copyright   Arepo is free software: you can redistribute it and/or modify
 *              it under the terms of the GNU General Public License as published by
 *              the Free Software Foundation, either version 3 of the License, or
 *              (at your option) any later version.
 *
 *              Arepo is distributed in the hope that it will be useful,
 *              but WITHOUT ANY WARRANTY; without even the implied warranty of
 *              MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *              GNU General Public License for more details.
 *
 *              A copy of the GNU General Public License is available under
 *              LICENSE as part of this program.  See also
 *              <https://www.gnu.org/licenses/>.
 *
 * \file        src/riemann_hllc.c
 * \date        05/2018
 * \brief       Routines for a HLLC Riemann solver.
 * \details     contains functions:
 *                static void hllc_get_fluxes_from_state(struct state *st,
 *                  struct fluxes *flux)
 *                static double get_hllc_star_fluxes(const struct state *st,
 *                  const struct fluxes *flux, struct fluxes *hllc_flux,
 *                  double S_star, double S)
 *                double godunov_flux_3d_hllc(struct state *st_L, struct state
 *                  *st_R, struct state_face *st_face, struct fluxes *flux)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 21.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"

#if defined(RIEMANN_HLLC)

#if defined(RIEMANN_HLLD)
#error option RIEMANN_HLLC is incompatible with option RIEMANN_HLLD.
Only one Riemann solver can be chosen among the above options.If none of them is selected,
    the exact Riemann solver will be used.
#endif /* #if defined(RIEMANN_HLLD) */
    /*! \brief Calculates the flux from a state.
     *
     *  Mass, momentum and energy flux.
     *
     *  \param[in] st State.
     *  \param[out] flux Flux corresponding to the state.
     *
     *  \return void
     */
    static void
    hllc_get_fluxes_from_state(struct state *st, struct fluxes *flux)
{
  flux->mass        = st->rho * st->velx;
  flux->momentum[0] = st->rho * st->velx * st->velx + st->press;
  flux->momentum[1] = st->rho * st->velx * st->vely;
  flux->momentum[2] = st->rho * st->velx * st->velz;

  st->Energy   = st->press / GAMMA_MINUS1 + 0.5 * st->rho * (st->velx * st->velx + st->vely * st->vely + st->velz * st->velz);
  flux->energy = (st->Energy + st->press) * st->velx;
}

/*! \brief Calculates a central flux in HLLC approximation.
 *
 *  \param[in] st State of the Riemann problem (either left or right).
 *  \param[in] flux Flux through face (either left or right).
 *  \param[out] hllc_flux State at the face (determined by this routine).
 *  \param[in] S_star speed of characteristics in central region.
 *  \param[in] S speed of characteristics in outside state (left or right).
 *
 *  \return Central density.
 */
static double get_hllc_star_fluxes(const struct state *st, const struct fluxes *flux, struct fluxes *hllc_flux, double S_star,
                                   double S)
{
  double Q0 = st->rho * (S - st->velx) / (S - S_star);
  double Q1 = Q0 * S_star;
  double Q2 = Q0 * st->vely;
  double Q3 = Q0 * st->velz;
  double Q4 = Q0 * (st->Energy / st->rho + (S_star - st->velx) * (S_star + st->press / (st->rho * (S - st->velx))));

  hllc_flux->mass = flux->mass + S * (Q0 - st->rho);

  hllc_flux->momentum[0] = flux->momentum[0] + S * (Q1 - st->rho * st->velx);

  hllc_flux->momentum[1] = flux->momentum[1] + S * (Q2 - st->rho * st->vely);

  hllc_flux->momentum[2] = flux->momentum[2] + S * (Q3 - st->rho * st->velz);

  hllc_flux->energy = flux->energy + S * (Q4 - st->Energy);

  return Q0;
}

/*! \brief Main routine for the hllc Riemann solver.
 *
 *  Called in finite_volume_solver.c
 *
 *  \param[in] st_L Left state of the Riemann problem.
 *  \param[in] st_R Right state of the Riemann problem.
 *  \param[out] st_face State at face.
 *  \param[out] flux Flux through face.
 *
 *  \return Pressure.
 */
double godunov_flux_3d_hllc(struct state *st_L, struct state *st_R, struct state_face *st_face, struct fluxes *flux)
{
  double S_L, S_R, S_star;
  double Press_star, rho_star;
  double rho_hat, csnd_hat;

  if(st_L->rho > 0 && st_R->rho > 0)
    {
      struct fluxes flux_L, flux_R;

      st_L->csnd = sqrt(GAMMA * st_L->press / st_L->rho);
      st_R->csnd = sqrt(GAMMA * st_R->press / st_R->rho);

      /* first estimate wave speeds */
      S_L = dmin(st_L->velx - st_L->csnd, st_R->velx - st_R->csnd);
      S_R = dmax(st_L->velx + st_L->csnd, st_R->velx + st_R->csnd);

      rho_hat    = 0.5 * (st_L->rho + st_R->rho);
      csnd_hat   = 0.5 * (st_L->csnd + st_R->csnd);
      Press_star = 0.5 * ((st_L->press + st_R->press) + (st_L->velx - st_R->velx) * (rho_hat * csnd_hat));
      S_star     = 0.5 * ((st_L->velx + st_R->velx) + (st_L->press - st_R->press) / (rho_hat * csnd_hat));

      /* compute fluxes for the left and right states */
      hllc_get_fluxes_from_state(st_L, &flux_L);
      hllc_get_fluxes_from_state(st_R, &flux_R);

      if(S_L >= 0.0) /* F_hllc = F_L */
        {
          /* copy the fluxes from the left state */
          flux->mass        = flux_L.mass;
          flux->momentum[0] = flux_L.momentum[0];
          flux->momentum[1] = flux_L.momentum[1];
          flux->momentum[2] = flux_L.momentum[2];
          flux->energy      = flux_L.energy;

          /* set the primitive variables at the face */
          st_face->rho   = st_L->rho;
          st_face->velx  = st_L->velx;
          st_face->vely  = st_L->vely;
          st_face->velz  = st_L->velz;
          st_face->press = st_L->press;
        }
      else if(S_R <= 0.0) /* F_hllc = F_R */
        {
          /* copy the fluxes from the left state */
          flux->mass        = flux_R.mass;
          flux->momentum[0] = flux_R.momentum[0];
          flux->momentum[1] = flux_R.momentum[1];
          flux->momentum[2] = flux_R.momentum[2];
          flux->energy      = flux_R.energy;

          /* set the primitive variables at the face */
          st_face->rho   = st_R->rho;
          st_face->velx  = st_R->velx;
          st_face->vely  = st_R->vely;
          st_face->velz  = st_R->velz;
          st_face->press = st_R->press;
        }
      else if(S_L <= 0.0 && S_star >= 0.0) /* F_hllc = F*_L */
        {
          /* compute star flux */
          rho_star = get_hllc_star_fluxes(st_L, &flux_L, flux, S_star, S_L);

          /* set the primitive variables at the face */
          st_face->rho   = rho_star;
          st_face->velx  = S_star;
          st_face->vely  = st_L->vely;
          st_face->velz  = st_L->velz;
          st_face->press = Press_star;
        }
      else /* F_hllc = F*_R */
        {
          /* compute star flux */
          rho_star = get_hllc_star_fluxes(st_R, &flux_R, flux, S_star, S_R);

          /* set the primitive variables at the face */
          st_face->rho   = rho_star;
          st_face->velx  = S_star;
          st_face->vely  = st_R->vely;
          st_face->velz  = st_R->velz;
          st_face->press = Press_star;
        }
    }
  else
    {
      printf("Left:  st_L->press=%g st_L->rho=%g  st_L->velx=%g\n", st_L->press, st_L->rho, st_L->velx);
      printf("Right: st_R->press=%g st_R->rho=%g  st_R->velx=%g\n", st_R->press, st_R->rho, st_R->velx);
      terminate("density is zero\n");
      return 0;
    }

  return st_face->press;
}

#endif /* #if defined(RIEMANN_HLLC) */
