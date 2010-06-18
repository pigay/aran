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

#ifndef __ARAN_FIT_H__
#define __ARAN_FIT_H__

#include <aranlinear.h>


typedef void (*AranFittingFunc) (gint n, gdouble x, gdouble *f);

void aran_svd_fit (gint ndata, gint nfunc, gdouble *x, gdouble *y,
                   gdouble *sigma, gdouble *a, gdouble **u, gdouble **v,
                   gdouble *w, gdouble *chisq, AranFittingFunc func);

#endif /* __ARAN_FIT_H__ */
