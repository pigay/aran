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

#ifndef __ARAN_POLY1D_H__
#define __ARAN_POLY1D_H__

#include <glib.h>
#include <stdio.h>

G_BEGIN_DECLS;

typedef struct _AranPoly1d AranPoly1d;

struct _AranPoly1d {
  guint degree;
  gdouble *terms;
};

/* macros */
#define aran_poly1d_term(a,i) ((a)->terms[i])

/* functions */
AranPoly1d *aran_poly1d_new (guint deg);

AranPoly1d *aran_poly1d_new_with_terms (int degree, const double *terms);

void aran_poly1d_free (AranPoly1d *ap1d);

void aran_poly1d_write (AranPoly1d *ap1d, FILE *file);

void aran_poly1d_add (AranPoly1d *one, AranPoly1d *other, AranPoly1d *ret);
void aran_poly1d_sub (AranPoly1d *one, AranPoly1d *other, AranPoly1d *ret);

void aran_poly1d_scalp (AranPoly1d *ap1d, gdouble factor, AranPoly1d *ret);

gdouble aran_poly1d_eval (AranPoly1d *ap1d, gdouble x);

void aran_poly1d_write_key_file (AranPoly1d *ap1d, GKeyFile *kf, gchar *group,
                                 gchar *key);
G_END_DECLS;

#endif /* __ARAN_POLY1D_H__ */
