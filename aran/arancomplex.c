/* LIBARAN - Fast Multipole Method library
 * Copyright (C) 2006-2007 Pierre Gay
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "arancomplex.h"

#ifdef VSG_HAVE_MPI
/**
 * ARAN_MPI_TYPE_GCOMPLEX64:
 *
 * The #MPI_Datatype associated to #gcomplex64
 */

MPI_Datatype aran_gcomplex64_get_mpi_type (void)
{
  static MPI_Datatype gcomplex64_mpi_type = MPI_DATATYPE_NULL;

  if (gcomplex64_mpi_type == MPI_DATATYPE_NULL)
    {
      MPI_Type_contiguous (2, MPI_FLOAT, &gcomplex64_mpi_type);
      MPI_Type_commit (&gcomplex64_mpi_type);
    }

  return gcomplex64_mpi_type;
}

/**
 * ARAN_MPI_TYPE_GCOMPLEX128:
 *
 * The #MPI_Datatype associated to #gcomplex128
 */

MPI_Datatype aran_gcomplex128_get_mpi_type (void)
{
  static MPI_Datatype gcomplex128_mpi_type = MPI_DATATYPE_NULL;

  if (gcomplex128_mpi_type == MPI_DATATYPE_NULL)
    {
      MPI_Type_contiguous (2, MPI_DOUBLE, &gcomplex128_mpi_type);
      MPI_Type_commit (&gcomplex128_mpi_type);
    }

  return gcomplex128_mpi_type;
}
#endif
