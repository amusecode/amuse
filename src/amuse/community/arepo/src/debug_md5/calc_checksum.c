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
 * \file        src/debug_md5/calc_checksum.c
 * \date        05/2018
 * \brief       Functions to calculate an MD5 checksum from a dataset.
 * \details     contains functions:
 *                void calc_memory_checksum(void *base, size_t bytes)
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 24.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#include "../main/allvars.h"
#include "../main/proto.h"

#include "Md5.h"

/*! \brief Calculates a md5 checksum (on all MPI tasks) and prints it.
 *
 *  \param[in] base Pointer to start of data.
 *  \param[in] bytes Number of bytes to be checked.
 *
 *  \return void
 */
void calc_memory_checksum(void *base, size_t bytes)
{
  MD5_CTX sum;
  union
  {
    unsigned char digest[16];
    int val[4];
  } u, uglob;

  MD5Init(&sum);
  MD5UpdateLong(&sum, base, bytes);
  MD5Final(&sum);

  int i;

  for(i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob.val, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      printf("Step=%d  MD5=", All.NumCurrentTiStep);
      for(i = 0; i < 16; i++)
        printf("%02x", uglob.digest[i]);
      printf("\n");
    }
}

#ifdef RESTART_DEBUG
/*! \brief Calculates md5 checksums of main data structures of a restart file.
 *
 *  \return void
 */
void log_restart_debug(void)
{
  MD5_CTX sum;
  union
  {
    unsigned char digest[16];
    int val[4];
  } u, uglob_P, uglob_SphP;
  int i;

  MD5Init(&sum);
  MD5UpdateLong(&sum, (void *)P, NumPart * sizeof(struct particle_data));
  MD5Final(&sum);

  for(i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob_P.val, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  MD5Init(&sum);
  MD5UpdateLong(&sum, (void *)SphP, NumGas * sizeof(struct sph_particle_data));
  MD5Final(&sum);

  for(i = 0; i < 16; i++)
    u.digest[i] = sum.digest[i];

  MPI_Allreduce(u.val, uglob_SphP.val, 4, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask == 0)
    {
      fprintf(FdRestartTest, "\n");
      fprintf(FdRestartTest, "Step=%8d  P[]        ", All.NumCurrentTiStep);
      for(i = 0; i < 16; i++)
        fprintf(FdRestartTest, "%02x", uglob_P.digest[i]);
      fprintf(FdRestartTest, "\n");
      fprintf(FdRestartTest, "               SphP[]     ");
      for(i = 0; i < 16; i++)
        fprintf(FdRestartTest, "%02x", uglob_SphP.digest[i]);
      fprintf(FdRestartTest, "\n");
      fflush(FdRestartTest);
    }
}
#endif
