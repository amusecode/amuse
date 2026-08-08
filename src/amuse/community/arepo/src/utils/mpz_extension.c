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
 * \file        src/utils/mpz_extension.c
 * \date        05/2018
 * \brief       Auxiliary functions to facilitate usage of mpz functions.
 * \details     Integer arithmetic used by Voronoi mesh construction.
 *              contains functions:
 *                void MY_mpz_set_si(mpz_t dest, signed long long int val)
 *                void MY_mpz_mul_si(mpz_t prod, mpz_t mult, signed long long
 *                  int val)
 *                void MY_mpz_sub_ui(mpz_t prod, mpz_t mult,
 *                  unsigned long long int val)
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 20.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include <gmp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../main/allvars.h"
#include "../main/proto.h"

#include "../mesh/voronoi/voronoi.h"

#if USEDBITS > 31

/*! \brief Sets mpz variable from signed long long int.
 *
 *  For Arepo-internal use of mpz.
 *
 *  \param[out] dest Variable to be set.
 *  \param[in] val Value in signed long long int.
 *
 *  \return void
 */
void MY_mpz_set_si(mpz_t dest, signed long long int val)
{
  mpz_t tmp, tmp2;

  unsigned long int lower = (unsigned long int)(val & 0xffffffffL);
  signed long int higher  = (signed long int)(val >> 32);

  mpz_init(tmp);
  mpz_init(tmp2);

  mpz_set_si(tmp, higher);
  mpz_mul_2exp(tmp2, tmp, 32);
  mpz_add_ui(dest, tmp2, lower);

  mpz_clear(tmp2);
  mpz_clear(tmp);
}

/*! \brief Multiplies an mpz type with a signed long long int.
 *
 *  \param[out] pred Result of multiplication.
 *  \param[in] mult Multiplicator (mpz_t).
 *  \param[in] val Multiplicand (signed long long int).
 *
 *  \return void
 */
void MY_mpz_mul_si(mpz_t prod, mpz_t mult, signed long long int val)
{
  mpz_t tmp;

  mpz_init(tmp);

  MY_mpz_set_si(tmp, val);

  mpz_mul(prod, mult, tmp);

  mpz_clear(tmp);
}

/*! \brief Subtracts 'val' from 'mult'.
 *
 *  \param[out] prod Result of subtraction.
 *  \param[in] mult Minuend (mpz_t).
 *  \param[in] val Subtrahend (unsigned long long int).
 *
 *  \return void
 */
void MY_mpz_sub_ui(mpz_t prod, mpz_t mult, unsigned long long int val)
{
  mpz_t tmp;

  mpz_init(tmp);

  MY_mpz_set_si(tmp, val);

  mpz_sub(prod, mult, tmp);

  mpz_clear(tmp);
}

#endif
