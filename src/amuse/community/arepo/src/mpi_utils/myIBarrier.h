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
 * \file        src/mpi_utils/myIBarrier.h
 * \date        05/2018
 * \brief       Header for myIBarrier functions.
 * \details
 *
 *
 * \par Major modifications and contributions:
 *
 * - DD.MM.YYYY Description
 * - 27.05.2018 Prepared file for public release -- Rainer Weinberger
 */

#ifndef MYIBARRIER_H
#define MYIBARRIER_H

#ifdef MYIBARRIER
#define MPI_TAG_IBARRIER 0x666

struct sMyIBarrier
{
  MPI_Comm comm;
  int rank;
  int nTasks;
  int nLevels;
  char *LevelDone;
};

void myIBarrier(MPI_Comm comm, struct sMyIBarrier *barrier);
void myIBarrierTest(struct sMyIBarrier *barrier, int *flag, MPI_Status *unused);
#endif /* #ifdef MYIBARRIER */

#endif /* #ifndef MYIBARRIER_H */
