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

#ifndef __ARAN_SOLVER3D_H__
#define __ARAN_SOLVER3D_H__

#include <glib.h>

#include <vsg/vsgd.h>

#include <aran/aran.h>
#include <aran/arandevelopment3d.h>

G_BEGIN_DECLS;

/* typedefs */
typedef struct _AranSolver3d AranSolver3d;

/* Particle functions */
typedef void (*AranParticle2ParticleFunc3d) (VsgPoint3 src, VsgPoint3 dst);

typedef void (*AranParticle2MultipoleFunc3d) (VsgPoint3 src,
                                              const VsgPRTree3dNodeInfo *dst_node,
                                              gpointer dst);

typedef void (*AranLocal2ParticleFunc3d) (const VsgPRTree3dNodeInfo *src_node,
                                          gpointer src,
                                          VsgPoint3 dst);

/* Translation functions */
typedef void (*AranMultipole2MultipoleFunc3d) (const VsgPRTree3dNodeInfo *src_node,
                                                gpointer src,
                                                const VsgPRTree3dNodeInfo *dst_node,
                                                gpointer dst);

typedef gboolean (*AranMultipole2LocalFunc3d) (const VsgPRTree3dNodeInfo *src_node,
                                                gpointer src,
                                                const VsgPRTree3dNodeInfo *dst_node,
                                                gpointer dst);

typedef void (*AranLocal2LocalFunc3d) (const VsgPRTree3dNodeInfo *src_node,
                                        gpointer src,
                                        const VsgPRTree3dNodeInfo *dst_node,
                                        gpointer dst);

/* functions */

AranSolver3d *aran_solver3d_new (VsgPRTree3d *prtree,
				 GType devel_type,
				 gpointer devel,
				 AranZeroFunc zero);

void aran_solver3d_free (AranSolver3d *solver);

void aran_solver3d_set_development (AranSolver3d *solver,
				    GType devel_type,
				    gpointer devel,
				    AranZeroFunc zero);

void aran_solver3d_set_functions (AranSolver3d *solver,
				  AranParticle2ParticleFunc3d p2p,
				  AranParticle2MultipoleFunc3d p2m,
				  AranMultipole2MultipoleFunc3d m2m,
				  AranMultipole2LocalFunc3d m2l,
				  AranLocal2LocalFunc3d l2l,
				  AranLocal2ParticleFunc3d l2p);

gdouble aran_solver3d_get_tolerance (AranSolver3d *solver);

void aran_solver3d_set_tolerance (AranSolver3d *solver,
				  gdouble tolerance);

void aran_solver3d_get_bounds (AranSolver3d *solver,
			       VsgVector3d *lbound,
			       VsgVector3d *ubound);

guint aran_solver3d_depth (const AranSolver3d *solver);

guint aran_solver3d_point_count (const AranSolver3d *solver);

void aran_solver3d_insert_point (AranSolver3d *solver,
                                 VsgPoint3 point);

gboolean aran_solver3d_insert_point_local (AranSolver3d *solver,
                                           VsgPoint3 point);

gboolean aran_solver3d_remove_point (AranSolver3d *solver,
                                     VsgPoint3 point);

VsgPoint3 aran_solver3d_find_point (AranSolver3d *solver,
                                    VsgPoint3 selector);

void aran_solver3d_foreach_point (AranSolver3d *solver,
                                  GFunc func,
                                  gpointer user_data);

void aran_solver3d_foreach_point_custom (AranSolver3d *solver,
                                         VsgRegion3 selector,
                                         VsgRegion3Point3LocFunc locfunc,
                                         GFunc func,
                                         gpointer user_data);

void aran_solver3d_set_children_order_hilbert (AranSolver3d *solver);

void aran_solver3d_set_children_order_default (AranSolver3d *solver);

void aran_solver3d_solve (AranSolver3d *solver);

#ifdef VSG_HAVE_MPI

void aran_solver3d_set_parallel (AranSolver3d *solver,
                                 VsgPRTreeParallelConfig *pconfig);

void aran_solver3d_migrate_flush (AranSolver3d *solver);

void aran_solver3d_distribute_nodes (AranSolver3d *solver,
                                     VsgPRTree3dDistributionFunc func,
                                     gpointer user_data);

void aran_solver3d_distribute_contiguous_leaves (AranSolver3d *solver);

#endif /* VSG_HAVE_MPI */

G_END_DECLS;

#endif /* __ARAN_SOLVER3D_H__ */
