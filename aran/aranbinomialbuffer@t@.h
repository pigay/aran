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

#ifndef __ARAN_BINOMIALBUFFER@T@_H__
#define __ARAN_BINOMIALBUFFER@T@_H__

#include <glib.h>

#include <aran/arancomplex.h>

G_BEGIN_DECLS;

/* typedefs */

typedef struct _AranBinomialBuffer@t@ AranBinomialBuffer@t@;

typedef @type@ (*AranBinomialBuffer@t@Generator) (guint l, guint m,
                                                  AranBinomialBuffer@t@ *buf);


/* functions */

AranBinomialBuffer@t@ *
aran_binomial_buffer@t@_new (AranBinomialBuffer@t@Generator generator,
                             guint l);

void aran_binomial_buffer@t@_free (AranBinomialBuffer@t@ *buf);

void aran_binomial_buffer@t@_require (AranBinomialBuffer@t@ *buf,
                                      guint max);

@type@ *aran_binomial_buffer@t@_get (AranBinomialBuffer@t@ *buf,
                                     guint l, guint m);

@type@ *aran_binomial_buffer@t@_get_unsafe (AranBinomialBuffer@t@ *buf,
                                            guint l, guint m);



G_END_DECLS;

#endif /* __ARAN_BINOMIALBUFFER@T@_H__ */
