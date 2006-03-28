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

#ifndef __ARAN_SOLVER2D_H__
#define __ARAN_SOLVER2D_H__

#include <glib.h>

#include <vsg/vsgd.h>

#include <aran/aran.h>
#include <aran/arandevelopment2d.h>

G_BEGIN_DECLS;

/* typedefs */
typedef struct _AranSolver2d AranSolver2d;

typedef void (*AranParticle2ParticleFunc2d) (VsgPoint2 src, VsgPoint2 dst);

typedef void (*AranParticle2MultipoleFunc2d) (VsgPoint2 src,
                                              const VsgVector2d *dst_center,
                                              gpointer dst);

typedef void (*AranMultipole2MultipoleFunc2d) (const VsgVector2d *src_center,
                                               gpointer src,
                                               const VsgVector2d *dst_center,
                                               gpointer dst);

typedef void (*AranMultipole2LocalFunc2d) (const VsgVector2d *src_center,
                                           gpointer src,
                                           const VsgVector2d *dst_center,
                                           gpointer dst);

typedef void (*AranLocal2LocalFunc2d) (const VsgVector2d *src_center,
                                       gpointer src,
                                       const VsgVector2d *dst_center,
                                       gpointer dst);

typedef void (*AranLocal2ParticleFunc2d) (const VsgVector2d *src_center,
                                          gpointer src,
                                          VsgPoint2 dst);

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

void aran_solver2d_solve (AranSolver2d *solver);

G_END_DECLS;

#endif /* __ARAN_SOLVER2D_H__ */
