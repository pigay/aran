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

#include "aran.h"
#include "aranbinomial.h"

/**
 * AranZeroFunc:
 * @devel: some data.
 *
 * Functions used to set some development to zero in #AranSolver2d and
 * #AranSolver.
 */

/**
 * gcomplex64:
 *
 * C99 64 bytes (single precision) complex numbers.
 */

/**
 * gcomplex128:
 *
 * C99 128 bytes (double precision) complex numbers.
 */

/**
 * G_I:
 *
 * C99 complex number %i (satifies equation: "%i=sqrt(-1)").
 */

/* from aransolver2d.c */
void aran_solver2d_init ();

/* from aransolver3d.c */
void aran_solver3d_init ();

/**
 * aran_init:
 *
 * Call this before you use functions from #Aran library.
 */
void aran_init()
{
  g_type_init ();

  aran_solver2d_init ();
  aran_solver3d_init ();

  /* require some precomputed binomial values */
  aran_binomial_require (20);
}
