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

#ifndef __ARAN_SOLVER2D_H__
#define __ARAN_SOLVER2D_H__

#include <glib.h>

#include <vsg/vsgd.h>

#include <aran/aran.h>
#include <aran/arandevelopment2d.h>

G_BEGIN_DECLS;

/* typedefs */
typedef struct _AranSolver2d AranSolver2d;

/* Particle functions */
typedef void (*AranParticle2ParticleFunc2d) (VsgPoint2 src, VsgPoint2 dst);

typedef void (*AranParticle2MultipoleFunc2d) (VsgPoint2 src,
                                              const VsgPRTree2dNodeInfo *dst_node,
                                              gpointer dst);

typedef void (*AranLocal2ParticleFunc2d) (const VsgPRTree2dNodeInfo *src_node,
                                          gpointer src,
                                          VsgPoint2 dst);

/* Translation functions */
typedef void (*AranMultipole2MultipoleFunc2d) (const VsgPRTree2dNodeInfo *src_node,
                                                gpointer src,
                                                const VsgPRTree2dNodeInfo *dst_node,
                                                gpointer dst);

typedef gboolean (*AranMultipole2LocalFunc2d) (const VsgPRTree2dNodeInfo *src_node,
                                                gpointer src,
                                                const VsgPRTree2dNodeInfo *dst_node,
                                                gpointer dst);

typedef void (*AranLocal2LocalFunc2d) (const VsgPRTree2dNodeInfo *src_node,
                                        gpointer src,
                                        const VsgPRTree2dNodeInfo *dst_node,
                                        gpointer dst);


/* functions */

AranSolver2d *aran_solver2d_new (VsgPRTree2d *prtree,
				 GType devel_type,
				 gpointer devel,
				 AranZeroFunc zero);

void aran_solver2d_free (AranSolver2d *solver);

void aran_solver2d_set_development (AranSolver2d *solver,
				    GType devel_type,
				    gpointer devel,
				    AranZeroFunc zero);

void aran_solver2d_set_functions (AranSolver2d *solver,
				  AranParticle2ParticleFunc2d p2p,
				  AranParticle2MultipoleFunc2d p2m,
				  AranMultipole2MultipoleFunc2d m2m,
				  AranMultipole2LocalFunc2d m2l,
				  AranLocal2LocalFunc2d l2l,
				  AranLocal2ParticleFunc2d l2p);

void aran_solver2d_reinit_stats (AranSolver2d *solver);
void aran_solver2d_get_stats (AranSolver2d *solver, glong *zero_count,
			      glong *p2p_count, glong *p2m_count,
			      glong *m2m_count, glong *m2l_count,
			      glong *l2l_count, glong *l2p_count);

gdouble aran_solver2d_get_tolerance (AranSolver2d *solver);

void aran_solver2d_set_tolerance (AranSolver2d *solver,
				  gdouble tolerance);

void aran_solver2d_get_bounds (AranSolver2d *solver,
			       VsgVector2d *lbound,
			       VsgVector2d *ubound);

guint aran_solver2d_depth (const AranSolver2d *solver);

guint aran_solver2d_point_count (const AranSolver2d *solver);

void aran_solver2d_insert_point (AranSolver2d *solver,
                                 VsgPoint2 point);

gboolean aran_solver2d_insert_point_local (AranSolver2d *solver,
                                           VsgPoint2 point);

gboolean aran_solver2d_remove_point (AranSolver2d *solver,
                                     VsgPoint2 point);

VsgPoint2 aran_solver2d_find_point (AranSolver2d *solver,
                                    VsgPoint2 selector);

void aran_solver2d_foreach_point (AranSolver2d *solver,
                                  GFunc func,
                                  gpointer user_data);

void aran_solver2d_foreach_point_custom (AranSolver2d *solver,
                                         VsgRegion2 selector,
                                         VsgRegion2Point2LocFunc locfunc,
                                         GFunc func,
                                         gpointer user_data);

void aran_solver2d_traverse (AranSolver2d *solver,
			     GTraverseType order,
			     VsgPRTree2dFunc func,
			     gpointer user_data);

void aran_solver2d_set_children_order_hilbert (AranSolver2d *solver);

void aran_solver2d_set_children_order_default (AranSolver2d *solver);


void aran_solver2d_solve (AranSolver2d *solver);

#ifdef VSG_HAVE_MPI

void aran_solver2d_set_parallel (AranSolver2d *solver,
                                 VsgPRTreeParallelConfig *pconfig);

void aran_solver2d_get_parallel (AranSolver2d *solver,
                                 VsgPRTreeParallelConfig *pconfig);

void aran_solver2d_migrate_flush (AranSolver2d *solver);

void aran_solver2d_distribute_nodes (AranSolver2d *solver,
                                     VsgPRTree2dDistributionFunc func,
                                     gpointer user_data);

void aran_solver2d_distribute_contiguous_leaves (AranSolver2d *solver);

#endif /* VSG_HAVE_MPI */

G_END_DECLS;

#endif /* __ARAN_SOLVER2D_H__ */
