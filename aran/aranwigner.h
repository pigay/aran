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

#ifndef __ARAN_WIGNER_H__
#define __ARAN_WIGNER_H__

#include <glib.h>
#include <stdio.h>

#include <aran/arancomplex.h>

typedef struct _AranWigner AranWigner;

AranWigner *aran_wigner_new (gdouble alpha, gdouble beta, gdouble gamma, gint l);
void aran_wigner_free (AranWigner * aw);

void aran_wigner_copy (AranWigner *src, AranWigner *dst);

void aran_wigner_require (AranWigner * aw, guint l);

gcomplex128 *aran_wigner_term (AranWigner * aw, guint l, guint mprime, gint m);

void aran_wigner_write (AranWigner * aw, FILE * file);

#endif /* __ARAN_WIGNER_H__ */
