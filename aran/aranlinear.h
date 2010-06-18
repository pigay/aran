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

#ifndef __ARAN_LINEAR_H__
#define __ARAN_LINEAR_H__

#include <glib.h>
#include <float.h>

G_BEGIN_DECLS;

gdouble aran_vec_norm (unsigned int n, gdouble *vec);

void aran_vec_sub (unsigned int n, gdouble *left, gdouble *right,
                   gdouble *result);

void aran_mat_vec_mult (unsigned int n, gdouble **mat, gdouble *vec,
                        gdouble *result);

int aran_linear_solve (unsigned int n, gdouble **mat, gdouble *vec);

void aran_pivot_init (unsigned int n, int *piv);

int aran_lu_factorize (unsigned int n, gdouble **mat, int *piv);
void aran_lu_updown (unsigned int n, gdouble **mat, int *piv, gdouble *rhs);

gdouble aran_lu_determinant (unsigned int n, gdouble **mat, int *piv);

void aran_mat_mat_mult (unsigned int m, unsigned int n, unsigned int o,
                        gdouble **left, gdouble **right,
                        gdouble **result);

void aran_mat_transpose (unsigned int n, gdouble **mat, gdouble **result);

void aran_svd (unsigned int m, unsigned int n, gdouble **a, gdouble *w,
               gdouble **v);

void aran_svd_solve (unsigned int m, unsigned int n, gdouble **u, gdouble *w,
                     gdouble **v, gdouble *sol, gdouble *rhs);

G_END_DECLS;

#endif /* __ARAN_LINEAR_H__ */
