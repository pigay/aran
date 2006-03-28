/* LIBARAN - Fast Multipole Method library
 * Copyright (C) 2006 Pierre Gay
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

#include "aranbinomial.h"

#include "aranbinomialbufferd.h"

/**
 * aran_slow_binomial: 
 * @n: a #guint.
 * @p: a #guint.
 *
 * Computes the binomial coefficient (choice of @p in @n) using slow
 * recurrence formula.
 *
 * Returns: C_@p^@n.
 */
gdouble aran_slow_binomial (guint n, guint p)
{
  if (p > n) return 0.;
  if (p == 0) return 1.;

  return aran_slow_binomial (n-1, p) + aran_slow_binomial (n-1, p-1);
}

static AranBinomialBufferd *_binomial = NULL;

static gdouble _binomial_generator (guint l, guint m, AranBinomialBufferd *buf)
{
  if (l == 0) return 1.;
  if ((m == 0) || (m == l)) return 1.;

  return *aran_binomial_bufferd_get_unsafe (buf, l-1, m) +
    *aran_binomial_bufferd_get_unsafe (buf, l-1, m-1);
}

static void _binomial_atexit ()
{
  if (_binomial != NULL)
    aran_binomial_bufferd_free (_binomial);

  _binomial = NULL;

}

/**
 * aran_binomial: 
 * @n: a #guint.
 * @p: a #guint.
 *
 * Computes the binomial coefficient (choice of @p in @n).
 */

/**
 * aran_fast_binomial: 
 * @n: a #guint.
 * @p: a #guint.
 *
 * Computes the binomial coefficient (choice of @p in @n) using fast
 * buffered values.
 *
 * Returns: C_@p^@n.
 */
gdouble aran_fast_binomial (guint n, guint p)
{
  aran_binomial_require (n);

  return *aran_binomial_bufferd_get_unsafe (_binomial, n, p);
}

/**
 * aran_binomial_require: 
 * @max: a #guint.
 *
 * Preallocates binomial coefficients up to n=@max. Calling this for a
 * sufficiently large @max can improve your program efficiency, since
 * successive calls with slowly increasing values of @max can lead to
 * numerous buffer reallocations.
 */
void aran_binomial_require (guint max)
{
  if (_binomial == NULL)
    {
      _binomial = aran_binomial_bufferd_new (_binomial_generator, max);
      g_atexit (_binomial_atexit);
    }
  else
    {
      aran_binomial_bufferd_require (_binomial, max);
    }
}
