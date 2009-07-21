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

#include "aransolver2d.h"

#include <glib-object.h>
/* #include <vsg/vsgd-inline.h> */

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
 * @dst_node: @dst tree node info.
 * @dst: destination development.
 *
 * Function provided to accumulate @src contribution into @dst multipole
 * expansion.
 */

/**
 * AranMultipole2MultipoleFunc2d:
 * @src_node: @src tree node info.
 * @src: source development.
 * @dst_node: @dst tree node info.
 * @dst: destination development.
 *
 * Function provided to translate a multipole expansion from @src to @dst.
 */

/**
 * AranMultipole2LocalFunc2d:
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
 * AranLocal2LocalFunc2d:
 * @src_node: @src tree node info.
 * @src: source development.
 * @dst_node: @dst tree node info.
 * @dst: destination development.
 *
 * 
 * Function provided to translate a local expansion from @src to @dst.
 */

/**
 * AranLocal2ParticleFunc2d:
 * @src_node: @src tree node info.
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
static glong _p2p_count = 0;
static glong _m2l_count = 0;
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
          _p2p_count ++;
	  other_list = other_list->next;
	}

      one_list = one_list->next;
    }
}

static gboolean far_func (const VsgPRTree2dNodeInfo *one_info,
			  const VsgPRTree2dNodeInfo *other_info,
			  AranSolver2d *solver)
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

      done = solver->m2l (other_info, other_dev,
			  one_info, one_dev);
      _m2l_count ++;

      if (!done)
	{
	  /* should _NOT_ happen : user error */
	  g_error ("m2l function (0x%p) return status not symmetric.",
		   solver->m2l);
	  return FALSE;
	}
    }

  return TRUE;
}


static void clear_func (const VsgPRTree2dNodeInfo *node_info,
                        AranSolver2d *solver)
{
  gpointer node_dev = node_info->user_data;

#ifdef VSG_HAVE_MPI
  if (VSG_PRTREE2D_NODE_INFO_IS_REMOTE (node_info)) return;
#endif

  aran_solver2d_clear_development (solver, node_dev);
}

static void up_func (const VsgPRTree2dNodeInfo *node_info,
		     AranSolver2d *solver)
{
  gpointer node_dev = node_info->user_data;

#ifdef VSG_HAVE_MPI
  if (VSG_PRTREE2D_NODE_INFO_IS_REMOTE (node_info)) return;
#endif

  if (node_info->isleaf)
    {
      if (solver->p2m != NULL)
	{
	  GSList *node_list = node_info->point_list;

	  while (node_list)
	    {
	      VsgPoint2 node_point = (VsgPoint2) node_list->data;

	      /* Particle to Multipole gathering */
	      solver->p2m (node_point, node_info, node_dev);

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
    }
}

static void down_func (const VsgPRTree2dNodeInfo *node_info,
		       AranSolver2d *solver)
{
  gpointer node_dev = node_info->user_data;

#ifdef VSG_HAVE_MPI
  if (VSG_PRTREE2D_NODE_INFO_IS_REMOTE (node_info)) return;
#endif

  if (solver->l2l != NULL && node_info->point_count != 0 &&
      node_info->father_info)
    {
      /* Local to Local translation */
      solver->l2l (node_info->father_info,
                   node_info->father_info->user_data,
                   node_info,
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
	      solver->l2p (node_info, node_dev, node_point);

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
    solver->prtree = vsg_prtree2d_new (&lbound, &ubound, NULL, 0);
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

  if (solver->devel != NULL)
    g_boxed_free (solver->devel_type, solver->devel);

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
 * aran_solver2d_insert_point:
 * @solver: an #AranSolver2d.
 * @point: a particle.
 *
 * Inserts a particle @point into the associated #VsgPRTree2d only if
 * @point falls in a local region of @solver's tree.
 *
 * Returns: #TRUE upon succesfull insertion.
 */
gboolean aran_solver2d_insert_point_local (AranSolver2d *solver,
                                           VsgPoint2 point)
{
  g_return_val_if_fail (solver != NULL, FALSE);

  return vsg_prtree2d_insert_point_local (solver->prtree, point);
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

void _traverse_print_dev (VsgPRTree2dNodeInfo *node_info, FILE *f)
{
  if (! VSG_PRTREE2D_NODE_INFO_IS_REMOTE (node_info))
    {
      AranDevelopment2d *dev = (AranDevelopment2d *) node_info->user_data;

      if (VSG_PRTREE2D_NODE_INFO_IS_LOCAL (node_info))
        fprintf(f, "* ");
      vsg_prtree_key2d_write (&node_info->id, f);
      fprintf(f, " ");
      aran_laurent_seriesd_write (dev->multipole, f);
      fprintf (f, "\n");
    }
}

/**
 * aran_solver2d_solve:
 * @solver: an #AranSolver2d.
 *
 * Solves the FMM problem for @solver.
 */
void aran_solver2d_solve (AranSolver2d *solver)
{
  GTimer *timer = NULL;
  gdouble t1, t2;
  gint rk = 0;

  g_return_if_fail (solver != NULL);

  _p2p_count = 0;
  _m2l_count = 0;
  timer = g_timer_new ();
#ifdef VSG_HAVE_MPI
  MPI_Comm_rank (MPI_COMM_WORLD, &rk);
#endif

  /* clear multipole and local developments before the big work */
  vsg_prtree2d_traverse (solver->prtree, G_POST_ORDER,
                         (VsgPRTree2dFunc) clear_func,
                         solver);

  t1 = g_timer_elapsed (timer, NULL);
  g_printerr ("%d : zero elapsed=%f seconds\n", rk, t1);

  /* gather information in Multipole development */
  vsg_prtree2d_traverse (solver->prtree, G_POST_ORDER,
                         (VsgPRTree2dFunc) up_func,
                         solver);

#ifdef VSG_HAVE_MPI
  /* gather shared in_counts */
  {
    VsgPRTreeParallelConfig pc;

    vsg_prtree2d_get_parallel (solver->prtree, &pc);
    vsg_prtree2d_shared_nodes_allreduce (solver->prtree,
                                         &pc.node_data.visit_forward);
  }
#endif /* VSG_HAVE_MPI */

  t2 = g_timer_elapsed (timer, NULL);
  g_printerr ("%d : up elapsed=%f seconds\n", rk, t2-t1);
  t1 = t2;

  /* transmit info from Multipole to Local developments */
  vsg_prtree2d_near_far_traversal (solver->prtree,
                                   (VsgPRTree2dFarInteractionFunc) far_func,
                                   (VsgPRTree2dInteractionFunc) near_func,
                                   solver);

  t2 = g_timer_elapsed (timer, NULL);
  g_printerr ("%d : nf elapsed=%f seconds\n", rk, t2-t1);
  t1 = t2;

  /* distribute information through Local developments towards particles */
  vsg_prtree2d_traverse (solver->prtree, G_PRE_ORDER,
                         (VsgPRTree2dFunc) down_func,
                         solver);

  t2 = g_timer_elapsed (timer, NULL);
  g_printerr ("%d : down elapsed=%f seconds\n", rk, t2-t1);

  g_printerr ("%d : p2p count=%ld\n", rk, _p2p_count);
  g_printerr ("%d : m2l count=%ld\n", rk, _m2l_count);

  g_timer_destroy (timer);
}

void aran_solver2d_set_children_order_hilbert (AranSolver2d *solver)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_set_children_order_hilbert (solver->prtree);

}

void aran_solver2d_set_children_order_default (AranSolver2d *solver)
{
  g_return_if_fail (solver != NULL);

 vsg_prtree2d_set_children_order_default (solver->prtree);
}

#ifdef VSG_HAVE_MPI

void aran_solver2d_set_parallel (AranSolver2d *solver,
                                 VsgPRTreeParallelConfig *pconfig)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_set_parallel (solver->prtree, pconfig);
}

void aran_solver2d_migrate_flush (AranSolver2d *solver)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_migrate_flush (solver->prtree);
}

void aran_solver2d_distribute_nodes (AranSolver2d *solver,
                                     VsgPRTree2dDistributionFunc func,
                                     gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_distribute_nodes (solver->prtree, func, user_data);
}

void aran_solver2d_distribute_contiguous_leaves (AranSolver2d *solver)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_distribute_contiguous_leaves (solver->prtree);
}

#endif

