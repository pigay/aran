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

#include "aransolver2d.h"

#include <glib-object.h>
#include <vsg/vsgd.h>

#include <string.h>

/**
 * AranSolver2d:
 *
 * Opaque structure. Possesses only private data.
 */
struct _AranSolver2d
{
  VsgPRTree2d *prtree;

  GType devel_type;
  gpointer devel;
  AranZeroFunc zero;

  AranParticle2ParticleFunc2d p2p;

  AranParticle2MultipoleFunc2d p2m;
  AranMultipole2MultipoleFunc2d m2m;
  AranMultipole2LocalFunc2d m2l;
  AranLocal2LocalFunc2d l2l;
  AranLocal2ParticleFunc2d l2p;

};

#define ARAN_SOLVER2D_PREALLOC 4

/**
 * AranParticle2ParticleFunc2d:
 * @src: source particle.
 * @dst: destination particle.
 *
 * Function provided to compute the particle/particle direct interaction.
 */

/**
 * AranParticle2MultipoleFunc2d:
 * @src: source particle.
 * @dst_center: @dst center.
 * @dst: destination development.
 *
 * Function provided to accumulate @src contribution into @dst multipole
 * expansion.
 */

/**
 * AranMultipole2MultipoleFunc2d:
 * @src_center: @src center.
 * @src: source development.
 * @dst_center: @dst center.
 * @dst: destination development.
 *
 * Function provided to translate a multipole expansion from @src to @dst.
 */

/**
 * AranMultipole2LocalFunc2d:
 * @src_center: @src center.
 * @src: source development.
 * @dst_center: @dst center.
 * @dst: destination development.
 *
 * Function provided to translate a multipole expansion from @src to @dst local
 * expansion.
 */

/**
 * AranLocal2LocalFunc2d:
 * @src_center: @src center.
 * @src: source development.
 * @dst_center: @dst center.
 * @dst: destination development.
 *
 * 
 * Function provided to translate a local expansion from @src to @dst.
 */

/**
 * AranLocal2ParticleFunc2d:
 * @src_center: @src center.
 * @src: source development.
 * @dst: destination particle.
 *
 * Function provided to evaluate @src contribuition and accumulate it into
 * particle @dst.
 */

static GMemChunk *aran_solver2d_mem_chunk = 0;
static guint aran_solver2d_instances_count = 0;

/* static functions: */
static void aran_solver2d_finalize ();
static AranSolver2d *_solver2d_alloc ();

static void aran_solver2d_finalize ()
{
  if (aran_solver2d_mem_chunk)
    {
      g_mem_chunk_destroy (aran_solver2d_mem_chunk);
      aran_solver2d_mem_chunk = 0;
    }
}

static AranSolver2d *_solver2d_alloc ()
{
  AranSolver2d *solver;

  if (!aran_solver2d_mem_chunk)
    {
      aran_solver2d_mem_chunk = g_mem_chunk_create (AranSolver2d,
						    ARAN_SOLVER2D_PREALLOC,
						    G_ALLOC_ONLY);
    }

  aran_solver2d_instances_count ++;

  solver = g_chunk_new (AranSolver2d, aran_solver2d_mem_chunk);

  solver->devel = NULL;
  solver->zero = NULL;

  solver->p2p = NULL;

  solver->p2m = NULL;
  solver->m2m = NULL;
  solver->m2l = NULL;
  solver->l2l = NULL;
  solver->l2p = NULL;

  return solver;
}

static void _solver2d_dealloc (AranSolver2d *solver)
{
  g_chunk_free (solver, aran_solver2d_mem_chunk);

  aran_solver2d_instances_count --;

  if (aran_solver2d_instances_count == 0)
    aran_solver2d_finalize ();
}

/* private function */
void aran_solver2d_init ()
{
  static gboolean wasinit = FALSE;

  if (! wasinit)
    {
      wasinit = TRUE;
      g_atexit (aran_solver2d_finalize);
    }
}


static void aran_solver2d_clear_development (AranSolver2d *solver,
					     gpointer devel)
{
  solver->zero (devel);
}
/*----------------------------------------------------*/

static void near_func (const VsgPRTree2dNodeInfo *one_info,
		       const VsgPRTree2dNodeInfo *other_info,
		       AranSolver2d *solver)
{
  GSList *one_list = one_info->point_list;

  if (solver->p2p == NULL) return;

  while (one_list)
    {
      VsgPoint2 one_point = (VsgPoint2) one_list->data;
      GSList *other_list =
	(one_info != other_info)? other_info->point_list : one_list;

      while (other_list)
	{
	  VsgPoint2 other_point = (VsgPoint2) other_list->data;

	  /* Particle to Particle interaction */
	  solver->p2p (one_point, other_point);

	  other_list = other_list->next;
	}

      one_list = one_list->next;
    }
}

static void far_func (const VsgPRTree2dNodeInfo *one_info,
		      const VsgPRTree2dNodeInfo *other_info,
		      AranSolver2d *solver)
{
  /* Multipole to Local transformation */
  if (solver->m2l != NULL)
    {
      gpointer one_dev = one_info->user_data;
      gpointer other_dev = other_info->user_data;

      /* both ways in order to get symmetric exchange */
      solver->m2l (&one_info->center, one_dev,
		   &other_info->center, other_dev);

      solver->m2l (&other_info->center, other_dev,
		   &one_info->center, one_dev);
    }
}


static void clear_func (const VsgPRTree2dNodeInfo *node_info,
                        AranSolver2d *solver)
{
  gpointer node_dev = node_info->user_data;

  aran_solver2d_clear_development (solver, node_dev);
}

static void up_func (const VsgPRTree2dNodeInfo *node_info,
		     AranSolver2d *solver)
{
  gpointer node_dev = node_info->user_data;

  if (node_info->isleaf)
    {
      if (solver->p2m != NULL)
	{
	  GSList *node_list = node_info->point_list;

	  while (node_list)
	    {
	      VsgPoint2 node_point = (VsgPoint2) node_list->data;

	      /* Particle to Multipole gathering */
	      solver->p2m (node_point, &node_info->center, node_dev);

	      node_list = node_list->next;
	    }
	}
    }

  if (solver->m2m != NULL && node_info->point_count != 0 &&
      node_info->father_info)
    {
      /* Multipole to Multipole translation */
      solver->m2m (&node_info->center,
                   node_dev,
                   &(node_info->father_info->center),
                   node_info->father_info->user_data);
    }
}

static void down_func (const VsgPRTree2dNodeInfo *node_info,
		       AranSolver2d *solver)
{
  gpointer node_dev = node_info->user_data;

  if (solver->l2l != NULL && node_info->point_count != 0 &&
      node_info->father_info)
    {
      /* Local to Local translation */
      solver->l2l (&(node_info->father_info->center),
                   node_info->father_info->user_data,
                   &node_info->center,
                   node_dev);
    }

  if ((node_info->isleaf))
    {
      if (solver->l2p != NULL)
	{
	  GSList *node_list = node_info->point_list;

	  while (node_list)
	    {
	      VsgPoint2 node_point = (VsgPoint2) node_list->data;

	      /* Local to Particle distribution */
	      solver->l2p (&node_info->center, node_dev, node_point);

	      node_list = node_list->next;
	    }
	}
    }
}



/* public functions */

/**
 * aran_solver2d_new:
 * @prtree: a #VsgPRTree2d.
 * @devel_type: #GType of development to be used by the solver.
 * @devel: instance of type @devel_type for cloning.
 * @zero: @devel_type zeroing function.
 *
 * Creates a new #AranSolver2d from @prtree and a development type. @prtree
 * is considered as owned by the new #AranSolver2d (ie. #VsgPRTree2d cannot be
 * shared between different #AranSolver2d nor being reused elsewhere).
 *
 * Returns: newly allocated structure.
 */
AranSolver2d *aran_solver2d_new (VsgPRTree2d *prtree,
				 GType devel_type,
				 gpointer devel,
				 AranZeroFunc zero)
{
  VsgVector2d lbound = {-1., -1.};
  VsgVector2d ubound = {1., 1.};
  AranSolver2d *solver;

  solver = _solver2d_alloc ();

  if (prtree == NULL)
    solver->prtree = vsg_prtree2d_new (&lbound, &ubound, NULL);
  else
    solver->prtree = prtree;

  aran_solver2d_set_development (solver, devel_type, devel, zero);

  return solver;
}

/**
 * aran_solver2d_free:
 * @solver: an #AranSolver2d.
 *
 * Deallocates all memory associated with @solver (Even the #VsgPRTree2d is
 * freed).
 */
void aran_solver2d_free (AranSolver2d *solver)
{
  if (solver == NULL) return;

  vsg_prtree2d_free (solver->prtree);

  _solver2d_dealloc (solver);
}

/**
 * aran_solver2d_set_development:
 * @solver: an #AranSolver2d.
 * @devel_type: #GType of development to be used by the solver.
 * @devel: instance of type @devel_type for cloning.
 * @zero: @devel_type zeroing function.
 *
 * Associates @solver with a new development definition.
 */
void aran_solver2d_set_development (AranSolver2d *solver,
				    GType devel_type,
				    gpointer devel,
				    AranZeroFunc zero)
{
  g_return_if_fail (solver != NULL);

  g_return_if_fail ((devel_type == G_TYPE_NONE) ||
		    G_TYPE_IS_BOXED (devel_type));

  if (solver->devel != NULL)
    g_boxed_free (solver->devel_type, solver->devel);

  solver->devel_type = devel_type;
  solver->devel = devel;
  solver->zero = zero;

  if (devel != NULL)
    {
      vsg_prtree2d_set_node_data (solver->prtree, devel_type, devel);
    }
  else
    {
      vsg_prtree2d_set_node_data (solver->prtree,
                                  G_TYPE_NONE,
                                  NULL);
    }
}

/**
 * aran_solver2d_set_functions:
 * @solver: an #AranSolver2d.
 * @p2p: particle 2 particle function.
 * @p2m: particle 2 multipole function.
 * @m2m: multipole 2 multipole function.
 * @m2l: multipole 2 local function.
 * @l2l: local 2 local function.
 * @l2p: local 2 particle function.
 *
 * Associates @solver with a set of FMM functions.
 */
void aran_solver2d_set_functions (AranSolver2d *solver,
				  AranParticle2ParticleFunc2d p2p,
				  AranParticle2MultipoleFunc2d p2m,
				  AranMultipole2MultipoleFunc2d m2m,
				  AranMultipole2LocalFunc2d m2l,
				  AranLocal2LocalFunc2d l2l,
				  AranLocal2ParticleFunc2d l2p)
{
  g_return_if_fail (solver != NULL);

  solver->p2p = p2p;

  solver->p2m = p2m;
  solver->m2m = m2m;
  solver->m2l = m2l;
  solver->l2l = l2l;
  solver->l2p = l2p;
}

/**
 * aran_solver2d_get_tolerance:
 * @solver: an #AranSolver2d.
 *
 * Inquiry of the spatial tolerance used in the associated #VsgPRTree2d.
 *
 * Returns: spatial tolerance.
 */
gdouble aran_solver2d_get_tolerance (AranSolver2d *solver)
{
  g_return_val_if_fail (solver != NULL, -1.);

  return vsg_prtree2d_get_tolerance (solver->prtree);
}

/**
 * aran_solver2d_set_tolerance:
 * @solver: an #AranSolver2d.
 * @tolerance: a positive #gdouble.
 *
 * Sets the spatial tolerance used in the associated #VsgPRTree2d.
 */
void aran_solver2d_set_tolerance (AranSolver2d *solver,
				  gdouble tolerance)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_set_tolerance (solver->prtree, tolerance);
}

/**
 * aran_solver2d_get_bounds:
 * @solver: an #AranSolver2d.
 * @lbound: bbox lower bound.
 * @ubound: bbox upper bound.
 *
 * Gets the bounding box for the associated #VsgPRTree2d.
 */
void aran_solver2d_get_bounds (AranSolver2d *solver,
			       VsgVector2d *lbound,
			       VsgVector2d *ubound)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_get_bounds (solver->prtree, lbound, ubound);
}

/**
 * aran_solver2d_depth:
 * @solver: an #AranSolver2d.
 *
 * Computes the depth of the associated #VsgPRTree2d.
 *
 * Returns: tree depth.
 */
guint aran_solver2d_depth (const AranSolver2d *solver)
{
  g_return_val_if_fail (solver != NULL, 0);

  return vsg_prtree2d_depth (solver->prtree);
}

/**
 * aran_solver2d_point_count:
 * @solver: an #AranSolver2d.
 *
 * Counts the number of points (particles) in the associated #VsgPRTree2d.
 *
 * Returns: point count.
 */
guint aran_solver2d_point_count (const AranSolver2d *solver)
{
  g_return_val_if_fail (solver != NULL, 0);

  return vsg_prtree2d_point_count (solver->prtree);
}

/**
 * aran_solver2d_insert_point:
 * @solver: an #AranSolver2d.
 * @point: a particle.
 *
 * Inserts a particle @point into the associated #VsgPRTree2d.
 */
void aran_solver2d_insert_point (AranSolver2d *solver,
                                 VsgPoint2 point)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_insert_point (solver->prtree, point);
}

/**
 * aran_solver2d_remove_point:
 * @solver: an #AranSolver2d.
 * @point: a particle.
 *
 * Removes a particle from the associated #VsgPRTree2d.
 *
 * Returns: %TRUE on success.
 */
gboolean aran_solver2d_remove_point (AranSolver2d *solver,
                                     VsgPoint2 point)
{
  g_return_val_if_fail (solver != NULL, FALSE);

  return vsg_prtree2d_remove_point (solver->prtree, point);
}

/**
 * aran_solver2d_find_point:
 * @solver: an #AranSolver2d.
 * @selector: a particle.
 *
 * See vsg_prtree2d_find_point() in Vsg API docs.
 *
 * Returns: found #VsgPoint2 or NULL.
 */
VsgPoint2 aran_solver2d_find_point (AranSolver2d *solver,
                                    VsgPoint2 selector)
{
  g_return_val_if_fail (solver != NULL, NULL);

  return vsg_prtree2d_find_point (solver->prtree, selector);
}

/**
 * aran_solver2d_foreach_point:
 * @solver: an #AranSolver2d.
 * @func: foreach function.
 * @user_data: pointer to pass to @func.
 *
 * See vsg_prtree2d_foreach_point() in Vsg API docs.
 */
void aran_solver2d_foreach_point (AranSolver2d *solver,
                                  GFunc func,
                                  gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_foreach_point (solver->prtree, func, user_data);
}

/**
 * aran_solver2d_foreach_point_custom:
 * @solver: an #AranSolver2d.
 * @selector: a #VsgRegion2.
 * @locfunc: localization function for @selector.
 * @func: foreach function.
 * @user_data: pointer to pass to @func.
 *
 * See vsg_prtree2d_foreach_point_custom() in Vsg API docs.
 */
void aran_solver2d_foreach_point_custom (AranSolver2d *solver,
                                         VsgRegion2 selector,
                                         VsgRegion2Point2LocFunc locfunc,
                                         GFunc func,
                                         gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_foreach_point_custom (solver->prtree,
                                     selector, locfunc,
                                     func, user_data);
}

/**
 * aran_solver2d_solve:
 * @solver: an #AranSolver2d.
 *
 * Solves the FMM problem for @solver.
 */
void aran_solver2d_solve (AranSolver2d *solver)
{
  g_return_if_fail (solver != NULL);

  /* clear multipole and local develoments before the big work */
  vsg_prtree2d_traverse (solver->prtree, G_POST_ORDER,
                         (VsgPRTree2dFunc) clear_func,
                         solver);

  /* gather information in Multipole development */
  vsg_prtree2d_traverse (solver->prtree, G_POST_ORDER,
                         (VsgPRTree2dFunc) up_func,
                         solver);

  /* transmit info from Multipole to Local developments */
  vsg_prtree2d_near_far_traversal (solver->prtree,
                                   (VsgPRTree2dInteractionFunc) far_func,
                                   (VsgPRTree2dInteractionFunc) near_func,
                                   solver);

  /* distribute information through Local developments towards particles */
  vsg_prtree2d_traverse (solver->prtree, G_PRE_ORDER,
                         (VsgPRTree2dFunc) down_func,
                         solver);
}
