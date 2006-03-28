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

#include "arancoefficientbuffer@t@.h"

/**
 * AranCoefficientBuffer@t@: 
 *
 * Opaque structure. Contains only private fields.
 */

struct _AranCoefficientBuffer@t@ {
  guint n;
  @type@ *buffer;
  AranCoefficientBuffer@t@Generator generator;
};

/**
 * AranCoefficientBuffer@t@Generator: 
 * @n: order of term.
 * @buf: the #AranCoefficientBuffer@t@.
 *
 * Generator function for #AranCoefficientBuffer@t@ structure. These functions
 * are user provided to compute #AranCoefficientBuffer@t@ terms. Argument @buf
 * is passed to the genrator to allow recurrent formulas.
 *
 * Returns: @buf[@n] value.
 */

/**
 * aran_coefficient_buffer@t@_new:
 * @generator: generation function.
 * @n: preallocation size.
 *
 * Allocates a new #AranCoefficientBuffer@t@.
 *
 * Returns: newly allocated structure.
 */
AranCoefficientBuffer@t@ *
aran_coefficient_buffer@t@_new (AranCoefficientBuffer@t@Generator generator,
                                guint n)
{
  AranCoefficientBuffer@t@ *ret;

  g_return_val_if_fail (generator != NULL, NULL);

  ret = g_malloc (sizeof (AranCoefficientBuffer@t@));

  ret->n = 0;
  ret->generator = generator;
  ret->buffer = NULL;

  aran_coefficient_buffer@t@_require (ret, n);

  return ret;
}

/**
 * aran_coefficient_buffer@t@_free:
 * @buf: an #AranCoefficientBuffer@t@.
 *
 * Deallocates @buf and all associated memory.
 */
void aran_coefficient_buffer@t@_free (AranCoefficientBuffer@t@ *buf)
{
  g_free (buf->buffer);
  g_free (buf);
}

/**
 * aran_coefficient_buffer@t@_require:
 * @buf: an #AranCoefficientBuffer@t@.
 * @max: new size.
 *
 * Computes missing terms in @buf up to @max.
 */
void aran_coefficient_buffer@t@_require (AranCoefficientBuffer@t@ *buf,
                                         guint max)
{
  guint size, i, oldn;

  g_return_if_fail (buf != NULL);

  if (buf->n >= max) return;

  size = MAX (buf->n+1, 1);

  while (size <= max) size *= 2;

  buf->buffer = g_realloc (buf->buffer, size * sizeof (@type@));

  oldn = buf->n;
  buf->n = size-1;

  for (i=oldn; i<size; i ++)
    {
      buf->buffer[i] = buf->generator (i, buf);
    }

}

/**
 * aran_coefficient_buffer@t@_get:
 * @buf: an #AranCoefficientBuffer@t@.
 * @n: a #guint.
 *
 * Localizes the specified @buf[@n] term.
 *
 * Returns: address of @buf[@n].
 */
@type@ *aran_coefficient_buffer@t@_get (AranCoefficientBuffer@t@ *buf,
                                        guint n)
{
  g_return_val_if_fail (buf != NULL, NULL);

  aran_coefficient_buffer@t@_require (buf, n);

  return aran_coefficient_buffer@t@_get_unsafe (buf, n);
}

/**
 * aran_coefficient_buffer@t@_get_unsafe:
 * @buf: an #AranCoefficientBuffer@t@.
 * @n: a #guint.
 *
 * As aran_coefficient_buffer@t@_get() except this function will not check
 * if @buf is allocated up to @n. Use for performance at your own risk.
 *
 * Returns: address of @buf[@n].
 */
@type@ *aran_coefficient_buffer@t@_get_unsafe (AranCoefficientBuffer@t@ *buf,
                                               guint n)
{
  return buf->buffer + n;
}


