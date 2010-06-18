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

#ifndef __ARAN_POLYNOMIAL_FIT_H__
#define __ARAN_POLYNOMIAL_FIT_H__

#include "aranpoly1d.h"

gdouble aran_poly1d_fit (gint nsamples, gdouble *x, gdouble *f,
                         AranPoly1d *ap1d);

gdouble aran_poly1d_fit_sigma (gint nsamples, gdouble *x, gdouble *f,
                               gdouble *sigma, AranPoly1d *ap1d);

#endif /* __ARAN_POLYNOMIAL_FIT_H__ */


