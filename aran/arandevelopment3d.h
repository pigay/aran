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

void aran_development3d_p2m (const VsgVector3d *position, const gdouble charge,
                             const VsgPRTree3dNodeInfo *dst_node,
                             AranDevelopment3d *dst);

void aran_development3d_p2l (const VsgVector3d *position, const gdouble charge,
                             const VsgPRTree3dNodeInfo *dst_node,
                             AranDevelopment3d *dst);

void aran_development3d_m2m (const VsgPRTree3dNodeInfo *src_node,
			     AranDevelopment3d *src,
			     const VsgPRTree3dNodeInfo *dst_node,
			     AranDevelopment3d *dst);

void aran_development3d_m2l (const VsgPRTree3dNodeInfo *src_node,
                             AranDevelopment3d *src,
                             const VsgPRTree3dNodeInfo *dst_node,
                             AranDevelopment3d *dst);

void aran_development3d_l2l (const VsgPRTree3dNodeInfo *src_node,
			     AranDevelopment3d *src,
			     const VsgPRTree3dNodeInfo *dst_node,
			     AranDevelopment3d *dst);

gcomplex128
aran_development3d_multipole_evaluate (const VsgPRTree3dNodeInfo *devel_node,
                                       AranDevelopment3d *devel,
                                       const VsgVector3d *pos);

gcomplex128 aran_development3d_local_evaluate (const VsgPRTree3dNodeInfo *devel_node,
					       AranDevelopment3d *devel,
					       const VsgVector3d *pos);

gcomplex128 aran_development3d_l2p (const VsgPRTree3dNodeInfo *devel_node,
                                    AranDevelopment3d *devel,
                                    const VsgVector3d *pos);

void
aran_development3d_multipole_gradient_evaluate (const VsgPRTree3dNodeInfo *devel_node,
                                                AranDevelopment3d *devel,
                                                const VsgVector3d *pos,
                                                VsgVector3d *grad);

void
aran_development3d_local_gradient_evaluate (const VsgPRTree3dNodeInfo *devel_node,
                                            AranDevelopment3d *devel,
                                            const VsgVector3d *pos,
                                            VsgVector3d *grad);

void
aran_development3d_l2pv (const VsgPRTree3dNodeInfo *devel_node,
                         AranDevelopment3d *devel,
                         const VsgVector3d *pos,
                         VsgVector3d *grad);

void aran_development3d_m2m_kkylin (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst);

void aran_development3d_m2l_kkylin (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst);

void aran_development3d_l2l_kkylin (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst);

void aran_development3d_m2m_rotate (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst);

void aran_development3d_m2l_rotate (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst);

void aran_development3d_l2l_rotate (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst);

#ifdef VSG_HAVE_MPI

void aran_development3d_vtable_init (VsgParallelVTable *vtable, guint8 posdeg,
                                     guint8 negdeg);

void aran_development3d_vtable_clear (VsgParallelVTable *vtable);

gpointer aran_development3d_alloc (gboolean resident, AranDevelopment3d *src);

void aran_development3d_destroy (gpointer data, gboolean resident,
                                 gpointer user_data);

void aran_development3d_migrate_pack (AranDevelopment3d *devel,
                                      VsgPackedMsg *pm,
                                      gpointer user_data);

void aran_development3d_migrate_unpack (AranDevelopment3d *devel,
                                        VsgPackedMsg *pm,
                                        gpointer user_data);

void aran_development3d_visit_fw_pack (AranDevelopment3d *devel,
                                       VsgPackedMsg *pm,
                                       gpointer user_data);

void aran_development3d_visit_fw_unpack (AranDevelopment3d *devel,
                                         VsgPackedMsg *pm,
                                         gpointer user_data);

void aran_development3d_visit_fw_reduce (AranDevelopment3d *a,
                                         AranDevelopment3d  *b,
                                         gpointer user_data);

void aran_development3d_visit_bw_pack (AranDevelopment3d *devel,
                                       VsgPackedMsg *pm,
                                       gpointer user_data);


void aran_development3d_visit_bw_unpack (AranDevelopment3d *devel,
                                         VsgPackedMsg *pm,
                                         gpointer user_data);

void aran_development3d_visit_bw_reduce (AranDevelopment3d *a,
                                         AranDevelopment3d  *b,
                                         gpointer user_data);

#endif /* VSG_HAVE_MPI */

G_END_DECLS;

#endif /* __ARAN_DEVELOPMENT3D_H__ */
