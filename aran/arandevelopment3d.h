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

#ifndef __ARAN_DEVELOPMENT3D_H__
#define __ARAN_DEVELOPMENT3D_H__

#include <glib-object.h>

#include <vsg/vsgd.h>
#include <aran/aransphericalharmonic.h>
#include <aran/aransphericalseriesd.h>

G_BEGIN_DECLS;

/* macros */
#define ARAN_TYPE_DEVELOPMENT3D (aran_development3d_get_type ())

/* typedefs */

typedef struct _AranDevelopment3d AranDevelopment3d;

struct _AranDevelopment3d
{
  AranSphericalSeriesd *multipole;
  AranSphericalSeriesd *local;
};

/* functions */
GType aran_development3d_get_type ();

AranDevelopment3d *aran_development3d_new (guint8 posdeg, guint8 negdeg);

void aran_development3d_free (AranDevelopment3d *ad);

void aran_development3d_copy (const AranDevelopment3d *src,
                              AranDevelopment3d *dst);

AranDevelopment3d *aran_development3d_clone (AranDevelopment3d *src);

void aran_development3d_set_zero (AranDevelopment3d *ad);

void aran_development3d_write (AranDevelopment3d *ad, FILE *file);

void aran_development3d_m2m (const VsgVector3d *src_center,
			     AranDevelopment3d *src,
			     const VsgVector3d *dst_center,
			     AranDevelopment3d *dst);

void aran_development3d_m2l (const VsgVector3d *src_center,
			     AranDevelopment3d *src,
			     const VsgVector3d *dst_center,
			     AranDevelopment3d *dst);

void aran_development3d_l2l (const VsgVector3d *src_center,
			     AranDevelopment3d *src,
			     const VsgVector3d *dst_center,
			     AranDevelopment3d *dst);

gcomplex128
aran_development3d_multipole_evaluate (const VsgVector3d *devel_center,
                                       AranDevelopment3d *devel,
                                       const VsgVector3d *pos);

gcomplex128 aran_development3d_local_evaluate (const VsgVector3d *devel_center,
					       AranDevelopment3d *devel,
					       const VsgVector3d *pos);

void
aran_development3d_local_gradient_evaluate (const VsgVector3d *devel_center,
                                            AranDevelopment3d *devel,
                                            const VsgVector3d *pos,
                                            VsgVector3d *grad);

void aran_development3d_m2m_kkylin (const VsgVector3d *src_center,
                                    AranDevelopment3d *src,
                                    const VsgVector3d *dst_center,
                                    AranDevelopment3d *dst);

void aran_development3d_m2l_kkylin (const VsgVector3d *src_center,
                                    AranDevelopment3d *src,
                                    const VsgVector3d *dst_center,
                                    AranDevelopment3d *dst);

void aran_development3d_l2l_kkylin (const VsgVector3d *src_center,
                                    AranDevelopment3d *src,
                                    const VsgVector3d *dst_center,
                                    AranDevelopment3d *dst);

void aran_development3d_m2m_rotate (const VsgVector3d *src_center,
                                    AranDevelopment3d *src,
                                    const VsgVector3d *dst_center,
                                    AranDevelopment3d *dst);

void aran_development3d_m2l_rotate (const VsgVector3d *src_center,
                                    AranDevelopment3d *src,
                                    const VsgVector3d *dst_center,
                                    AranDevelopment3d *dst);

void aran_development3d_l2l_rotate (const VsgVector3d *src_center,
                                    AranDevelopment3d *src,
                                    const VsgVector3d *dst_center,
                                    AranDevelopment3d *dst);

G_END_DECLS;

#endif /* __ARAN_DEVELOPMENT3D_H__ */
