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

#include "aransolver3d.h"

#include <glib-object.h>
/* #include <vsg/vsgd-inline.h> */

#include <string.h>

/**
 * AranSolver3d:
 *
 * Opaque structure. Possesses only private data.
 */
struct _AranSolver3d
{
  VsgPRTree3d *prtree;

  GType devel_type;
  gpointer devel;
  AranZeroFunc zero;

  AranParticle2ParticleFunc3d p2p;

  AranParticle2MultipoleFunc3d p2m;
  AranMultipole2MultipoleFunc3d m2m;
  AranMultipole2LocalFunc3d m2l;
  AranLocal2LocalFunc3d l2l;
  AranLocal2ParticleFunc3d l2p;

  glong zero_counter;
  glong p2p_counter;
  glong p2m_counter;
  glong m2m_counter;
  glong m2l_counter;
  glong l2l_counter;
  glong l2p_counter;
};

#define ARAN_SOLVER3D_PREALLOC 4

/**
 * AranParticle2ParticleFunc3d:
 * @src: source particle.
 * @dst: destination particle.
 *
 * Function provided to compute the particle/particle direct interaction.
 */

/**
 * AranParticle2MultipoleFunc3d:
 * @src: source particle.
 * @dst_node: @dst tree node info.
 * @dst: destination development.
 *
 * Function provided to accumulate @src contribution into @dst multipole
 * expansion.
 */

/**
 * AranMultipole2MultipoleFunc3d:
 * @src_node: @src tree node info.
 * @src: source development.
 * @dst_node: @dst tree node info.
 * @dst: destination development.
 *
 * Function provided to translate a multipole expansion from @src to @dst.
 */

/**
 * AranMultipole2LocalFunc3d:
 * @src_node: @src tree node info.
 * @src: source development.
 * @dst_node: @dst tree node info.
 * @dst: destination development.
 *
 * Function provided to translate a multipole expansion from @src to @dst local
 * expansion.
 *
 * If the function returns %FALSE, LibAran will consider that the translation
 * is to be avoided. In this case, the near interaction p2p will be called for
 * each pair of src and dst particles.
 *
 * WARNING: The function _must_ be symmetric. Otherwise, the behaviour of
 * LibAran is undefined.
 *
 * Returns: %TRUE if the translation took place.
 */

/**
 * AranLocal2LocalFunc3d:
 * @src_node: @src tree node info.
 * @src: source development.
 * @dst_node: @dst tree node info.
 * @dst: destination development.
 *
 * Function provided to translate a local expansion from @src to @dst.
 */

/**
 * AranLocal2ParticleFunc3d:
 * @src_node: @src tree node info.
 * @src: source development.
 * @dst: destination particle.
 *
 * Function provided to evaluate @src contribuition and accumulate it into
 * particle @dst.
 */

static GMemChunk *aran_solver3d_mem_chunk = 0;
static guint aran_solver3d_instances_count = 0;

/* static functions: */
static void aran_solver3d_finalize ();
static AranSolver3d *_solver3d_alloc ();

static void aran_solver3d_finalize ()
{
  if (aran_solver3d_mem_chunk)
    {
      g_mem_chunk_destroy (aran_solver3d_mem_chunk);
      aran_solver3d_mem_chunk = 0;
    }
}

static AranSolver3d *_solver3d_alloc ()
{
  AranSolver3d *solver;

  if (!aran_solver3d_mem_chunk)
    {
      aran_solver3d_mem_chunk = g_mem_chunk_create (AranSolver3d,
                                                    ARAN_SOLVER3D_PREALLOC,
                                                    G_ALLOC_ONLY);
    }

  aran_solver3d_instances_count ++;

  solver = g_chunk_new (AranSolver3d, aran_solver3d_mem_chunk);

  solver->devel = NULL;
  solver->zero = NULL;

  solver->p2p = NULL;

  solver->p2m = NULL;
  solver->m2m = NULL;
  solver->m2l = NULL;
  solver->l2l = NULL;
  solver->l2p = NULL;

  aran_solver3d_reinit_stats (solver);

  return solver;
}

static void _solver3d_dealloc (AranSolver3d *solver)
{
  g_chunk_free (solver, aran_solver3d_mem_chunk);

  aran_solver3d_instances_count --;

  if (aran_solver3d_instances_count == 0)
    aran_solver3d_finalize ();
}

/* private function */
void aran_solver3d_init ()
{
  static gboolean wasinit = FALSE;

  if (! wasinit)
    {
      wasinit = TRUE;
      g_atexit (aran_solver3d_finalize);
    }
}


/*----------------------------------------------------*/

static void near_func (const VsgPRTree3dNodeInfo *one_info,
                       const VsgPRTree3dNodeInfo *other_info,
                       AranSolver3d *solver)
{
  GSList *one_list = one_info->point_list;

  if (solver->p2p == NULL) return;

  while (one_list)
    {
      VsgPoint3 one_point = (VsgPoint3) one_list->data;
      GSList *other_list =
        (one_info != other_info)? other_info->point_list : one_list;

      while (other_list)
        {
          VsgPoint3 other_point = (VsgPoint3) other_list->data;

          /* Particle to Particle interaction */
          solver->p2p (one_point, other_point);
          solver->p2p_counter ++;

          other_list = other_list->next;
        }

      one_list = one_list->next;
    }
}

static gboolean far_func (const VsgPRTree3dNodeInfo *one_info,
                          const VsgPRTree3dNodeInfo *other_info,
                          AranSolver3d *solver)
{
  gboolean done;

  /* Multipole to Local transformation */
  if (solver->m2l != NULL)
    {
      gpointer one_dev = one_info->user_data;
      gpointer other_dev = other_info->user_data;

      /* both ways in order to get symmetric exchange */
      done = solver->m2l (one_info, one_dev,
                          other_info, other_dev);

      if (!done) return FALSE;

      solver->m2l_counter ++;

      done = solver->m2l (other_info, other_dev, one_info, one_dev);

      if (!done)
        {
          /* should _NOT_ happen : user error */
          g_error ("m2l function (0x%p) return status not symmetric.",
                   solver->m2l);
          return FALSE;
        }

      solver->m2l_counter ++;
    }

  return TRUE;
}



static void clear_func (const VsgPRTree3dNodeInfo *node_info,
                        AranSolver3d *solver)
{
  gpointer node_dev = node_info->user_data;

#ifdef VSG_HAVE_MPI
  if (VSG_PRTREE3D_NODE_INFO_IS_REMOTE (node_info)) return;
#endif

  solver->zero (node_dev);
  solver->zero_counter ++;
}

static void up_func (const VsgPRTree3dNodeInfo *node_info,
                     AranSolver3d *solver)
{
  gpointer node_dev = node_info->user_data;

#ifdef VSG_HAVE_MPI
  if (VSG_PRTREE3D_NODE_INFO_IS_REMOTE (node_info)) return;
#endif

  if (node_info->isleaf)
    {
      if (solver->p2m != NULL)
        {
          GSList *node_list = node_info->point_list;

          while (node_list)
            {
              VsgPoint3 node_point = (VsgPoint3) node_list->data;

              /* Particle to Multipole gathering */
              solver->p2m (node_point, node_info, node_dev);
              solver->p2m_counter ++;

              node_list = node_list->next;
            }
        }
    }

  if (solver->m2m != NULL && node_info->point_count != 0 &&
      node_info->father_info)
    {
      /* Multipole to Multipole translation */
      solver->m2m (node_info,
                   node_dev,
                   node_info->father_info,
                   node_info->father_info->user_data);
      solver->m2m_counter ++;
    }
}

static void down_func (const VsgPRTree3dNodeInfo *node_info,
                       AranSolver3d *solver)
{
  gpointer node_dev = node_info->user_data;

#ifdef VSG_HAVE_MPI
  if (VSG_PRTREE3D_NODE_INFO_IS_REMOTE (node_info)) return;
#endif

  if (solver->l2l != NULL && node_info->point_count != 0 &&
      node_info->father_info)
    {
      /* Local to Local translation */
      solver->l2l (node_info->father_info,
                   node_info->father_info->user_data,
                   node_info,
                   node_dev);
      solver->l2l_counter ++;
    }

  if ((node_info->isleaf))
    {
      if (solver->l2p != NULL)
        {
          GSList *node_list = node_info->point_list;

          while (node_list)
            {
              VsgPoint3 node_point = (VsgPoint3) node_list->data;

              /* Local to Particle distribution */
              solver->l2p (node_info, node_dev, node_point);
              solver->l2p_counter ++;

              node_list = node_list->next;
            }
        }
    }
}



/* public functions */

/**
 * aran_solver3d_new:
 * @prtree: a #VsgPRTree3d.
 * @devel_type: #GType of development to be used by the solver.
 * @devel: instance of type @devel_type for cloning.
 * @zero: @devel_type zeroing function.
 *
 * Creates a new #AranSolver3d from @prtree and a development type. @prtree
 * is considered as owned by the new #AranSolver3d (ie. #VsgPRTree3d cannot be
 * shared between different #AranSolver3d nor being reused elsewhere).
 *
 * Returns: newly allocated structure.
 */
AranSolver3d *aran_solver3d_new (VsgPRTree3d *prtree,
                                 GType devel_type,
                                 gpointer devel,
                                 AranZeroFunc zero)
{
  VsgVector3d lbound = {-1., -1., -1.};
  VsgVector3d ubound = {1., 1., 1.};
  AranSolver3d *solver;

  solver = _solver3d_alloc ();

  if (prtree == NULL)
    solver->prtree = vsg_prtree3d_new (&lbound, &ubound, NULL, 0);
  else
    solver->prtree = prtree;

  aran_solver3d_set_development (solver, devel_type, devel, zero);

  return solver;
}

/**
 * aran_solver3d_free:
 * @solver: an #AranSolver3d.
 *
 * Deallocates all memory associated with @solver (Even the #VsgPRTree3d is
 * freed).
 */
void aran_solver3d_free (AranSolver3d *solver)
{
  if (solver == NULL) return;

  vsg_prtree3d_free (solver->prtree);

  if (solver->devel != NULL)
    g_boxed_free (solver->devel_type, solver->devel);

  _solver3d_dealloc (solver);
}

/**
 * aran_solver3d_set_development:
 * @solver: an #AranSolver3d.
 * @devel_type: #GType of development to be used by the solver.
 * @devel: instance of type @devel_type for cloning.
 * @zero: @devel_type zeroing function.
 *
 * Associates @solver with a new development definition.
 */
void aran_solver3d_set_development (AranSolver3d *solver,
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
      vsg_prtree3d_set_node_data (solver->prtree, devel_type, devel);
    }
  else
    {
      vsg_prtree3d_set_node_data (solver->prtree,
                                  G_TYPE_NONE,
                                  NULL);
    }
}

/**
 * aran_solver3d_set_functions:
 * @solver: an #AranSolver3d.
 * @p2p: particle 2 particle function.
 * @p2m: particle 2 multipole function.
 * @m2m: multipole 2 multipole function.
 * @m2l: multipole 2 local function.
 * @l2l: local 2 local function.
 * @l2p: local 2 particle function.
 *
 * Associates @solver with a set of FMM functions.
 */
void aran_solver3d_set_functions (AranSolver3d *solver,
                                  AranParticle2ParticleFunc3d p2p,
                                  AranParticle2MultipoleFunc3d p2m,
                                  AranMultipole2MultipoleFunc3d m2m,
                                  AranMultipole2LocalFunc3d m2l,
                                  AranLocal2LocalFunc3d l2l,
                                  AranLocal2ParticleFunc3d l2p)
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
 * aran_solver3d_reinit_stats:
 * @solver: an #AranSolver3d.
 *
 * Resets @solver Near/Far counters to zero.
 */
void aran_solver3d_reinit_stats (AranSolver3d *solver)
{
  g_return_if_fail (solver != NULL);

  solver->zero_counter = 0;
  solver->p2p_counter = 0;

  solver->p2m_counter = 0;
  solver->m2m_counter = 0;
  solver->m2l_counter = 0;
  solver->l2l_counter = 0;
  solver->l2p_counter = 0;
}

/**
 * aran_solver3d_get_stats:
 * @solver: an #AranSolver3d.
 * @zero: Clear function count result.
 * @p2p: particle 2 particle function count result.
 * @p2m: particle 2 multipole function count result.
 * @m2m: multipole 2 multipole function count result.
 * @m2l: multipole 2 local function count result.
 * @l2l: local 2 local function count result.
 * @l2p: local 2 particle function count result.
 *
 * Retrieves count values for @solver differents functions. Each value
 * represents the number of calls of the corresponding functino since
 * last call to aran_solver3d_reinit_stats().
 */
void aran_solver3d_get_stats (AranSolver3d *solver, glong *zero_count,
                              glong *p2p_count, glong *p2m_count,
                              glong *m2m_count, glong *m2l_count,
                              glong *l2l_count, glong *l2p_count)
{
  g_return_if_fail (solver != NULL);

  *zero_count = solver->zero_counter;
  *p2p_count = solver->p2p_counter;
  *p2m_count = solver->p2m_counter;
  *m2m_count = solver->m2m_counter;
  *m2l_count = solver->m2l_counter;
  *l2l_count = solver->l2l_counter;
  *l2p_count = solver->l2p_counter;
}

/**
 * aran_solver3d_get_tolerance:
 * @solver: an #AranSolver3d.
 *
 * Inquiry of the spatial tolerance used in the associated #VsgPRTree3d.
 *
 * Returns: spatial tolerance.
 */
gdouble aran_solver3d_get_tolerance (AranSolver3d *solver)
{
  g_return_val_if_fail (solver != NULL, -1.);

  return vsg_prtree3d_get_tolerance (solver->prtree);
}

/**
 * aran_solver3d_set_tolerance:
 * @solver: an #AranSolver3d.
 * @tolerance: a positive #gdouble.
 *
 * Sets the spatial tolerance used in the associated #VsgPRTree3d.
 */
void aran_solver3d_set_tolerance (AranSolver3d *solver,
                                  gdouble tolerance)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_set_tolerance (solver->prtree, tolerance);
}

/**
 * aran_solver3d_get_bounds:
 * @solver: an #AranSolver3d.
 * @lbound: bbox lower bound.
 * @ubound: bbox upper bound.
 *
 * Gets the bounding box for the associated #VsgPRTree3d.
 */
void aran_solver3d_get_bounds (AranSolver3d *solver,
                               VsgVector3d *lbound,
                               VsgVector3d *ubound)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_get_bounds (solver->prtree, lbound, ubound);
}

/**
 * aran_solver3d_depth:
 * @solver: an #AranSolver3d.
 *
 * Computes the depth of the associated #VsgPRTree3d.
 *
 * Returns: tree depth.
 */
guint aran_solver3d_depth (const AranSolver3d *solver)
{
  g_return_val_if_fail (solver != NULL, 0);

  return vsg_prtree3d_depth (solver->prtree);
}

/**
 * aran_solver3d_point_count:
 * @solver: an #AranSolver3d.
 *
 * Counts the number of points (particles) in the associated #VsgPRTree3d.
 *
 * Returns: point count.
 */
guint aran_solver3d_point_count (const AranSolver3d *solver)
{
  g_return_val_if_fail (solver != NULL, 0);

  return vsg_prtree3d_point_count (solver->prtree);
}

/**
 * aran_solver3d_insert_point:
 * @solver: an #AranSolver3d.
 * @point: a particle.
 *
 * Inserts a particle @point into the associated #VsgPRTree3d.
 */
void aran_solver3d_insert_point (AranSolver3d *solver,
                                 VsgPoint3 point)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_insert_point (solver->prtree, point);
}

/**
 * aran_solver3d_insert_point_local:
 * @solver: an #AranSolver3d.
 * @point: a particle.
 *
 * Inserts a particle @point into the associated #VsgPRTree3d only if
 * @point falls in a local region of @solver's tree.
 *
 * Returns: #TRUE upon succesfull insertion.
 */
gboolean aran_solver3d_insert_point_local (AranSolver3d *solver,
                                           VsgPoint2 point)
{
  g_return_val_if_fail (solver != NULL, FALSE);

  return vsg_prtree3d_insert_point_local (solver->prtree, point);
}

/**
 * aran_solver3d_remove_point:
 * @solver: an #AranSolver3d.
 * @point: a particle.
 *
 * Removes a particle from the associated #VsgPRTree3d.
 *
 * Returns: %TRUE on success.
 */
gboolean aran_solver3d_remove_point (AranSolver3d *solver,
                                     VsgPoint3 point)
{
  g_return_val_if_fail (solver != NULL, FALSE);

  return vsg_prtree3d_remove_point (solver->prtree, point);
}

/**
 * aran_solver3d_find_point:
 * @solver: an #AranSolver3d.
 * @selector: a particle.
 *
 * See vsg_prtree3d_find_point() in Vsg API docs.
 *
 * Returns: found #VsgPoint3 or NULL.
 */
VsgPoint3 aran_solver3d_find_point (AranSolver3d *solver,
                                    VsgPoint3 selector)
{
  g_return_val_if_fail (solver != NULL, NULL);

  return vsg_prtree3d_find_point (solver->prtree, selector);
}


/**
 * aran_solver3d_foreach_point:
 * @solver: an #AranSolver3d.
 * @func: foreach function.
 * @user_data: pointer to pass to @func.
 *
 * See vsg_prtree3d_foreach_point() in Vsg API docs.
 */
void aran_solver3d_foreach_point (AranSolver3d *solver,
                                  GFunc func,
                                  gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_foreach_point (solver->prtree, func, user_data);
}

/**
 * aran_solver3d_traverse:
 * @solver: an #AranSolver3d.
 * @func: traverse function.
 * @user_data: pointer to pass to @func.
 *
 * See vsg_prtree3d_traverse() in Vsg API docs.
 */
void aran_solver3d_traverse (AranSolver3d *solver,
                             GTraverseType order,
                             VsgPRTree3dFunc func,
                             gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_traverse (solver->prtree, order, func, user_data);
}

/**
 * aran_solver3d_foreach_point_custom:
 * @solver: an #AranSolver3d.
 * @selector: a #VsgRegion3.
 * @locfunc: localization function for @selector.
 * @func: foreach function.
 * @user_data: pointer to pass to @func.
 *
 * See vsg_prtree3d_foreach_point_custom() in Vsg API docs.
 */
void aran_solver3d_foreach_point_custom (AranSolver3d *solver,
                                         VsgRegion3 selector,
                                         VsgRegion3Point3LocFunc locfunc,
                                         GFunc func,
                                         gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_foreach_point_custom (solver->prtree,
                                     selector, locfunc,
                                     func, user_data);
}

/**
 * aran_solver3d_solve:
 * @solver: an #AranSolver3d.
 *
 * Solves the FMM problem for @solver.
 */
void aran_solver3d_solve (AranSolver3d *solver)
{
  g_return_if_fail (solver != NULL);

#ifdef VSG_HAVE_MPI
  {
    VsgPRTreeParallelConfig pconfig;

    vsg_prtree3d_get_parallel (solver->prtree, &pconfig);
  }
#endif

  /* clear multipole and local developments before the big work */
  vsg_prtree3d_traverse (solver->prtree, G_POST_ORDER,
                         (VsgPRTree3dFunc) clear_func,
                         solver);

  /* gather information in Multipole development */
  vsg_prtree3d_traverse (solver->prtree, G_POST_ORDER,
                         (VsgPRTree3dFunc) up_func,
                         solver);

#ifdef VSG_HAVE_MPI
  /* gather shared in_counts */
  {
    VsgPRTreeParallelConfig pc;

    vsg_prtree3d_get_parallel (solver->prtree, &pc);
    vsg_prtree3d_shared_nodes_allreduce (solver->prtree,
                                         &pc.node_data.visit_forward);
  }
#endif /* VSG_HAVE_MPI */

  /* transmit info from Multipole to Local developments */
  vsg_prtree3d_near_far_traversal (solver->prtree,
                                   (VsgPRTree3dFarInteractionFunc) far_func,
                                   (VsgPRTree3dInteractionFunc) near_func,
                                   solver);

  /* distribute information through Local developments towards particles */
  vsg_prtree3d_traverse (solver->prtree, G_PRE_ORDER,
                         (VsgPRTree3dFunc) down_func,
                         solver);
}

void aran_solver3d_set_children_order_hilbert (AranSolver3d *solver)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_set_children_order_hilbert (solver->prtree);

}

void aran_solver3d_set_children_order_default (AranSolver3d *solver)
{
  g_return_if_fail (solver != NULL);

 vsg_prtree3d_set_children_order_default (solver->prtree);
}

#ifdef VSG_HAVE_MPI

void aran_solver3d_set_parallel (AranSolver3d *solver,
                                 VsgPRTreeParallelConfig *pconfig)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_set_parallel (solver->prtree, pconfig);
}

void aran_solver3d_get_parallel (AranSolver3d *solver,
                                 VsgPRTreeParallelConfig *pconfig)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_get_parallel (solver->prtree, pconfig);
}

void aran_solver3d_migrate_flush (AranSolver3d *solver)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_migrate_flush (solver->prtree);
}

void aran_solver3d_distribute_nodes (AranSolver3d *solver,
                                     VsgPRTree3dDistributionFunc func,
                                     gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_distribute_nodes (solver->prtree, func, user_data);
}

void aran_solver3d_distribute_contiguous_leaves (AranSolver3d *solver)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_distribute_contiguous_leaves (solver->prtree);
}

#endif

