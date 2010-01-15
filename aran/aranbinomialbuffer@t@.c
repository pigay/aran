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

#include "aranbinomialbuffer@t@.h"

/**
 * AranBinomialBuffer@t@: 
 *
 * Opaque structure. Contains only private fields.
 */

struct _AranBinomialBuffer@t@ {
  gint l;
  @type@ *buffer;
  @type@ **direct;
  AranBinomialBuffer@t@Generator generator;
};

/**
 * AranBinomialBuffer@t@Generator: 
 * @l: degree of term.
 * @m: order of term.
 * @buf: the #AranBinomialBuffer@t@.
 *
 * Generator function for #AranBinomialBuffer@t@ structure. These functions
 * are user provided to compute #AranBinomialBuffer@t@ terms. Argument @buf
 * is passed to the genrator to allow recurrent formulas.
 *
 * Returns: @buf[@l][@m] value.
 */

/**
 * aran_binomial_buffer@t@_new: 
 * @generator: the buffer's generation function.
 * @l: initial allocation order.
 *
 * Allocates an #AranBinomialBuffer@t@ structure from @generator function.
 *
 * Returns: newly allocated structure.
 */
AranBinomialBuffer@t@ *
aran_binomial_buffer@t@_new (AranBinomialBuffer@t@Generator generator,
                             guint l)
{
  AranBinomialBuffer@t@ *ret;

  g_return_val_if_fail (generator != NULL, NULL);

  ret = g_malloc (sizeof (AranBinomialBuffer@t@));

  ret->l = -1;
  ret->generator = generator;
  ret->buffer = NULL;
  ret->direct = NULL;

  aran_binomial_buffer@t@_require (ret, l);

  return ret;
}

/**
 * aran_binomial_buffer@t@_free: 
 * @buf: an #AranBinomialBuffer@t@.
 *
 * Deallocates @buf and all associated memory.
 */
void aran_binomial_buffer@t@_free (AranBinomialBuffer@t@ *buf)
{
  if (buf != NULL)
    {
      if (buf->buffer != NULL)
        g_free (buf->buffer);

      if (buf->direct != NULL)
        g_free (buf->direct);

      buf->buffer = NULL;
      buf->direct = NULL;

      g_free (buf);
    }
}

/**
 * aran_binomial_buffer@t@_require: 
 * @buf: an #AranBinomialBuffer@t@.
 * @max: a #guint.
 *
 * Preallocates @buf coefficients up to n=@max. 
 */
void aran_binomial_buffer@t@_require (AranBinomialBuffer@t@ *buf,
                                      guint max)
{
  guint l, size;
  gint i, j, oldl;

  g_return_if_fail (buf != NULL);

  if (buf->l > (gint) max) return;

  l = MAX (buf->l, 1);

  while (l < (gint) max) l *= 2;

  size = ((l+1)*(l+2))/2;

  buf->direct = g_realloc (buf->direct, (l+1) * sizeof (@type@ *));
  buf->buffer = g_realloc (buf->buffer, size * sizeof (@type@));

  for (i=0; i<=l; i ++)
    {
      buf->direct[i] = buf->buffer + (i*(i+1))/2;
    }

  oldl = buf->l;
  buf->l = l;

  for (i=oldl+1; i<=l; i ++)
    {
      for (j=0; j<=i; j ++)
        buf->direct[i][j] = buf->generator (i, j, buf);
    }
}

/**
 * aran_binomial_buffer@t@_get: 
 * @buf: an #AranBinomialBuffer@t@.
 * @l: degree.
 * @m: order.
 *
 * Localizes the specified @buf[@l][@m] term.
 *
 * Returns: address of the specified term.
 */
@type@ *aran_binomial_buffer@t@_get (AranBinomialBuffer@t@ *buf,
                                     guint l, guint m)
{
  g_return_val_if_fail (buf != NULL, NULL);

  g_return_val_if_fail (m <= l, NULL);

  aran_binomial_buffer@t@_require (buf, l);

  return aran_binomial_buffer@t@_get_unsafe (buf, l, m);
}

/**
 * aran_binomial_buffer@t@_get_unsafe: 
 * @buf: an #AranBinomialBuffer@t@.
 * @l: degree.
 * @m: order.
 *
 * Localizes the specified @buf[@l][@m] term without checking @l and @m are in
 * the correct intervals. if @l happens to be greater than the maximum
 * allocated degree for @buf or if @m > @l, the result is undefined.
 *
 * Returns: address of the specified term.
 */
@type@ *aran_binomial_buffer@t@_get_unsafe (AranBinomialBuffer@t@ *buf,
                                            guint l, guint m)
{
  return &(buf->direct[l][m]);
}


