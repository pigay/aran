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

#ifndef __ARAN_LEGENDRE_H__
#define __ARAN_LEGENDRE_H__

#include <glib.h>

G_BEGIN_DECLS;

/* functions */
gdouble aran_legendre_associated_evaluate_internal (guint l, guint m,
						    gdouble x,
						    gdouble sqrt_1_m_x2);

gdouble aran_legendre_associated_evaluate (guint l, guint m, gdouble x);

void aran_legendre_evaluate_multiple (guint l, gdouble x, gdouble *result);

gdouble aran_legendre_evaluate (guint l, gdouble x);

void aran_legendre_associated_evaluate_multiple_internal (guint l,
							  gdouble x,
							  gdouble sqrt_1_m_x2,
							  gdouble *result);

void aran_legendre_associated_evaluate_multiple (guint l,
						 gdouble x,
						 gdouble *result);

void
aran_legendre_associated_evaluate_special_internal (guint l,
                                                    gdouble x,
                                                    gdouble sqrt_1_m_x2,
                                                    gdouble *legendre,
                                                    gdouble *special_result);

gdouble *aran_legendre_associated_multiple_get_term (gint l, gint m,
						     gdouble *result);

G_END_DECLS;

#endif /* __ARAN_LEGENDRE_H__ */
