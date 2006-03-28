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

#ifndef __ARAN_SPHERICAL_HARMONIC_H__
#define __ARAN_SPHERICAL_HARMONIC_H__

#include <glib.h>

#include <aran/arancomplex.h>

G_BEGIN_DECLS;

/* functions */
gcomplex128 aran_spherical_harmonic_evaluate_internal (guint l, gint m,
						       gdouble cost,
						       gdouble sint,
						       gcomplex128 expmp);

gcomplex128 aran_spherical_harmonic_evaluate (guint l, gint m,
					      gdouble theta, gdouble phi);

void aran_spherical_harmonic_evaluate_multiple_internal (guint l,
							 gdouble cost,
							 gdouble sint,
							 gcomplex128 expp,
							 gcomplex128 *result);

void aran_spherical_harmonic_evaluate_multiple (guint l,
						gdouble theta,
						gdouble phi,
						gcomplex128 *result);

void
aran_spherical_harmonic_pre_gradient_multiple_internal (guint l,
                                                        gdouble cost,
                                                        gdouble sint,
                                                        gcomplex128 expp,
                                                        gcomplex128 *harmonics,
                                                        gcomplex128 *special);

void
aran_spherical_harmonic_pre_gradient_multiple (guint l,
                                               gdouble theta,
                                               gdouble phi,
                                               gcomplex128 *harmonics,
                                               gcomplex128 *special);

gcomplex128 *aran_spherical_harmonic_multiple_get_term (gint l, gint m,
						        gcomplex128 *result);

void aran_spherical_harmonic_require (guint degree);

G_END_DECLS;

#endif /* __ARAN_SPHERICAL_HARMONIC_H__ */
