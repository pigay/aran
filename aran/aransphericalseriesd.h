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

#ifndef __ARAN_SPHERICAL_SERIESD_H__
#define __ARAN_SPHERICAL_SERIESD_H__

#include <glib.h>

#include <vsg/vsgd.h>
#ifdef VSG_HAVE_MPI
#include <vsg/vsgpackedmsg.h>
#endif

#include <aran/arancomplex.h>

#include <aran/aranwigner.h>

G_BEGIN_DECLS;

/* typedefs */

typedef struct _AranSphericalSeriesd AranSphericalSeriesd;

/* functions */

AranSphericalSeriesd *aran_spherical_seriesd_new (guint8 posdeg, guint8 negdeg);

void aran_spherical_seriesd_free (AranSphericalSeriesd *ass);

void aran_spherical_seriesd_copy (const AranSphericalSeriesd *src,
                                  AranSphericalSeriesd *dst);

AranSphericalSeriesd *
aran_spherical_seriesd_clone (const AranSphericalSeriesd *src);

gcomplex128 *aran_spherical_seriesd_get_term (AranSphericalSeriesd *ass,
					      gint i, gint j);

guint8 aran_spherical_seriesd_get_posdeg (const AranSphericalSeriesd *ass);

guint8 aran_spherical_seriesd_get_negdeg (const AranSphericalSeriesd *ass);

void aran_spherical_seriesd_set_zero (AranSphericalSeriesd *ass);

void aran_spherical_seriesd_write (const AranSphericalSeriesd *ass,
				   FILE *file);

#ifdef VSG_HAVE_MPI
void aran_spherical_seriesd_pack (AranSphericalSeriesd *ass, VsgPackedMsg *pm);
void aran_spherical_seriesd_unpack (AranSphericalSeriesd *ass,
                                    VsgPackedMsg *pm);
#endif

gcomplex128
aran_spherical_seriesd_evaluate_internal (const AranSphericalSeriesd *ass,
					  gdouble r,
					  gdouble cost, gdouble sint,
					  gdouble cosp, gdouble sinp);

gcomplex128 aran_spherical_seriesd_evaluate (const AranSphericalSeriesd *ass,
					     const VsgVector3d *x);

void aran_spherical_seriesd_local_gradient_evaluate_internal
(const AranSphericalSeriesd *ass,
 gdouble r,
 gdouble cost, gdouble sint,
 gdouble cosp, gdouble sinp,
 gcomplex128 *dr, gcomplex128 *dt, gcomplex128 *dp);

void aran_spherical_seriesd_local_gradient_evaluate
(const AranSphericalSeriesd *ass, const VsgVector3d *x, VsgVector3d *grad);

void aran_spherical_seriesd_add (AranSphericalSeriesd *one,
                               AranSphericalSeriesd *other,
                               AranSphericalSeriesd *result);

void aran_spherical_seriesd_translate (const AranSphericalSeriesd *src,
				       const VsgVector3d *xsrc,
				       AranSphericalSeriesd *dst,
				       const VsgVector3d *xdst);

void aran_spherical_seriesd_to_local (const AranSphericalSeriesd *src,
				      const VsgVector3d *xsrc,
				      AranSphericalSeriesd *dst,
				      const VsgVector3d *xdst);

void aran_spherical_seriesd_translate_kkylin (const AranSphericalSeriesd *src,
                                              const VsgVector3d *xsrc,
                                              AranSphericalSeriesd *dst,
                                              const VsgVector3d *xdst);

void aran_spherical_seriesd_to_local_kkylin (const AranSphericalSeriesd *src,
                                             const VsgVector3d *xsrc,
                                             AranSphericalSeriesd *dst,
                                             const VsgVector3d *xdst);

void aran_spherical_seriesd_rotate (const AranSphericalSeriesd *src,
                                    gdouble alpha, gdouble beta,
                                    gdouble gamma,
                                    AranSphericalSeriesd *dst);

void aran_spherical_seriesd_rotate_inverse (const AranSphericalSeriesd * src,
                                            gdouble alpha, gdouble beta,
                                            gdouble gamma,
                                            AranSphericalSeriesd * dst);

void aran_spherical_seriesd_translate_rotate (const AranSphericalSeriesd *src,
                                              const VsgVector3d *xsrc,
                                              AranSphericalSeriesd *dst,
                                              const VsgVector3d *xdst);

void aran_spherical_seriesd_to_local_rotate (const AranSphericalSeriesd *src,
                                             const VsgVector3d *xsrc,
                                             AranSphericalSeriesd *dst,
                                             const VsgVector3d *xdst);

G_END_DECLS;

#endif /* __ARAN_SPHERICAL_SERIESD_H__ */
