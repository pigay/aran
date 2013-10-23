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

#include "aran-config.h"

#include <string.h>

#include <glib-object.h>
/* #include <vsg/vsgd-inline.h> */

#include <vsg/vsgtiming.h>

#include "aransolver2d.h"
#include "aranprofile.h"
#include "aranprofiledb.h"

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
  AranParticle2LocalFunc2d p2l;
  AranMultipole2ParticleFunc2d m2p;

  guint semifar_threshold;

  glong zero_counter;
  glong p2p_counter, p2p_remote_counter;
  glong p2m_counter;
  glong m2m_counter;
  glong m2l_counter, m2l_remote_counter;
  glong l2l_counter;
  glong l2p_counter;
  glong p2l_counter;
  glong m2p_counter;

  gdouble p2p_time;
  gdouble p2m_time;
  gdouble m2m_time;
  gdouble m2l_time;
  gdouble l2l_time;
  gdouble l2p_time;
  gdouble p2l_time;
  gdouble m2p_time;
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

#define _USE_G_SLICES GLIB_CHECK_VERSION (2, 10, 0)

#if ! _USE_G_SLICES

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
#endif /* ! _USE_G_SLICES */

static AranSolver2d *_solver2d_alloc ()
{
  AranSolver2d *solver;

#if _USE_G_SLICES
  solver = g_slice_new(AranSolver2d);
#else
  if (!aran_solver2d_mem_chunk)
    {
      aran_solver2d_mem_chunk = g_mem_chunk_create (AranSolver2d,
                                                    ARAN_SOLVER2D_PREALLOC,
                                                    G_ALLOC_ONLY);
    }

  aran_solver2d_instances_count ++;

  solver = g_chunk_new (AranSolver2d, aran_solver2d_mem_chunk);
#endif /* _USE_G_SLICES */

  solver->devel = NULL;
  solver->zero = NULL;

  solver->p2p = NULL;

  solver->p2m = NULL;
  solver->m2m = NULL;
  solver->m2l = NULL;
  solver->l2l = NULL;
  solver->l2p = NULL;
  solver->p2l = NULL;
  solver->m2p = NULL;

  solver->semifar_threshold = 0;

  aran_solver2d_reinit_stats (solver);

  solver->p2p_time = -1.;
  solver->p2m_time = -1.;
  solver->m2m_time = -1.;
  solver->m2l_time = -1.;
  solver->l2l_time = -1.;
  solver->l2p_time = -1.;
  solver->p2l_time = -1.;
  solver->m2p_time = -1.;

  return solver;
}

static void _solver2d_dealloc (AranSolver2d *solver)
{
#if _USE_G_SLICES
  g_slice_free (AranSolver2d, solver);
#else
  g_chunk_free (solver, aran_solver2d_mem_chunk);

  aran_solver2d_instances_count --;

  if (aran_solver2d_instances_count == 0)
    aran_solver2d_finalize ();
#endif /* _USE_G_SLICES */
}

/* private function */
void aran_solver2d_init ()
{
#if ! _USE_G_SLICES
  static gboolean wasinit = FALSE;

  if (! wasinit)
    {
      wasinit = TRUE;
      g_atexit (aran_solver2d_finalize);
    }
#endif /* ! _USE_G_SLICES */
}


/*----------------------------------------------------*/
static void nop_near_func (const VsgPRTree2dNodeInfo *one_info,
                       const VsgPRTree2dNodeInfo *other_info,
                       AranSolver2d *solver)
{
}

/* general case near_func algorithm */
static void near_func_default (const VsgPRTree2dNodeInfo *one_info,
                               const VsgPRTree2dNodeInfo *other_info,
                               AranSolver2d *solver)
{
  GSList *one_list = one_info->point_list;

  while (one_list)
    {
      VsgPoint2 one_point = (VsgPoint2) one_list->data;
      GSList *other_list = other_info->point_list;

      while (other_list)
        {
          VsgPoint2 other_point = (VsgPoint2) other_list->data;

          /* Particle to Particle interaction */
          solver->p2p (one_point, other_point);

          other_list = other_list->next;
        }

      one_list = one_list->next;
    }

  solver->p2p_counter +=one_info->point_count * other_info->point_count;
  if (VSG_PRTREE2D_NODE_INFO_IS_REMOTE (one_info) ||
      VSG_PRTREE2D_NODE_INFO_IS_REMOTE (other_info))
    solver->p2p_remote_counter +=one_info->point_count * other_info->point_count;
}

/* near_func algorithm for reflexive interaction (one_info == other_info */
static void near_func_reflexive (const VsgPRTree2dNodeInfo *one_info,
                                 const VsgPRTree2dNodeInfo *other_info,
                                 AranSolver2d *solver)
{
  GSList *one_list = one_info->point_list;

  while (one_list)
    {
      VsgPoint2 one_point = (VsgPoint2) one_list->data;
      GSList *other_list = one_list;

      while (other_list)
        {
          VsgPoint2 other_point = (VsgPoint2) other_list->data;

          /* Particle to Particle interaction */
          solver->p2p (one_point, other_point);

          other_list = other_list->next;
        }

      one_list = one_list->next;
    }

  solver->p2p_counter +=
    (one_info->point_count * (one_info->point_count+1)) / 2;

  if (VSG_PRTREE2D_NODE_INFO_IS_REMOTE (one_info))
    solver->p2p_remote_counter +=one_info->point_count * other_info->point_count;
}

static void near_func (const VsgPRTree2dNodeInfo *one_info,
                       const VsgPRTree2dNodeInfo *other_info,
                       AranSolver2d *solver)
{
  if (one_info == other_info)
    near_func_reflexive (one_info, other_info, solver);
  else
    near_func_default (one_info, other_info, solver);
}


static void nop_far_func (const VsgPRTree2dNodeInfo *one_info,
                          const VsgPRTree2dNodeInfo *other_info,
                          AranSolver2d *solver)
{
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
      solver->m2l (one_info, one_dev, other_info, other_dev);

      solver->m2l_counter ++;

      solver->m2l (other_info, other_dev, one_info, one_dev);

      solver->m2l_counter ++;

      if (VSG_PRTREE3D_NODE_INFO_IS_REMOTE (one_info) ||
          VSG_PRTREE3D_NODE_INFO_IS_REMOTE (other_info))
        solver->m2l_remote_counter += 2;
    }
}

static void semifar_func (const VsgPRTree2dNodeInfo *one_info,
                          const VsgPRTree2dNodeInfo *other_info,
                          AranSolver2d *solver)
{
  const VsgPRTree2dNodeInfo *info;
  GSList *list;
  gpointer dev = one_info->user_data;
  guint point_count = 0;

  if (one_info->depth > other_info->depth)
    {
      list = other_info->point_list;
      info = one_info;
      dev = one_info->user_data;
    }
  else
    {
      list = one_info->point_list;
      info = other_info;
      dev = other_info->user_data;
    }

  /* g_printerr ("semifar one[%#lx %#lx %d] other[%#lx %#lx %d]\n", */
  /*             one_info->id.x, one_info->id.y, one_info->depth, */
  /*             other_info->id.x, other_info->id.y, other_info->depth); */

  while (list)
    {
      VsgPoint3 point = (VsgPoint3) list->data;

      solver->p2l (point, info, dev);
      solver->m2p (info, dev, point);

      point_count ++;
      list = g_slist_next (list);
    }

  solver->p2l_counter += point_count;
  solver->m2p_counter += point_count;
}


static void clear_func (const VsgPRTree2dNodeInfo *node_info,
                        AranSolver2d *solver)
{
  gpointer node_dev = node_info->user_data;

#ifdef VSG_HAVE_MPI
  if (VSG_PRTREE2D_NODE_INFO_IS_REMOTE (node_info)) return;
#endif

  solver->zero (node_dev);
  solver->zero_counter ++;
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
      solver->l2l_counter ++;
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
              solver->l2p_counter ++;

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
 * Associates @solver with a minimal set of FMM functions.
 */
void aran_solver2d_set_functions (AranSolver2d *solver,
                                  AranParticle2ParticleFunc2d p2p,
                                  AranParticle2MultipoleFunc2d p2m,
                                  AranMultipole2MultipoleFunc2d m2m,
                                  AranMultipole2LocalFunc2d m2l,
                                  AranLocal2LocalFunc2d l2l,
                                  AranLocal2ParticleFunc2d l2p)
{
  aran_solver2d_set_functions_full (solver, p2p, p2m, m2m, m2l, l2l, l2p,
                                    NULL, NULL, 0);
}

/**
 * aran_solver2d_set_functions_full:
 * @solver: an #AranSolver2d.
 * @p2p: particle 2 particle function.
 * @p2m: particle 2 multipole function.
 * @m2m: multipole 2 multipole function.
 * @m2l: multipole 2 local function.
 * @l2l: local 2 local function.
 * @l2p: local 2 particle function.
 * @p2l: particle 2 local function.
 * @m2p: multipole 2 particle function.
 * @semifar_threshold: beyond this number of particles, @p2l and @m2p
 * can be called instead of @p2p. If @semifar_threshold is #G_MAXUINT,
 * @p2p will always be called.
 *
 * Associates @solver with a complete set of FMM functions.
 */
void aran_solver2d_set_functions_full (AranSolver2d *solver,
                                       AranParticle2ParticleFunc2d p2p,
                                       AranParticle2MultipoleFunc2d p2m,
                                       AranMultipole2MultipoleFunc2d m2m,
                                       AranMultipole2LocalFunc2d m2l,
                                       AranLocal2LocalFunc2d l2l,
                                       AranLocal2ParticleFunc2d l2p,
                                       AranParticle2LocalFunc2d p2l,
                                       AranMultipole2ParticleFunc2d m2p,
                                       guint semifar_threshold)
{
  g_return_if_fail (solver != NULL);

  if (solver->p2p != p2p) solver->p2p_time = -1.;
  if (solver->p2m != p2m) solver->p2m_time = -1.;
  if (solver->m2m != m2m) solver->m2m_time = -1.;
  if (solver->m2l != m2l) solver->m2l_time = -1.;
  if (solver->l2l != l2l) solver->l2l_time = -1.;
  if (solver->l2p != l2p) solver->l2p_time = -1.;
  if (solver->p2l != p2l) solver->p2l_time = -1.;
  if (solver->m2p != m2p) solver->m2p_time = -1.;

  solver->p2p = p2p;

  solver->p2m = p2m;
  solver->m2m = m2m;
  solver->m2l = m2l;
  solver->l2l = l2l;
  solver->l2p = l2p;
  solver->p2l = p2l;
  solver->m2p = m2p;

  solver->semifar_threshold = semifar_threshold;
}

void aran_solver2d_get_functions_full (AranSolver2d *solver,
                                       AranParticle2ParticleFunc2d *p2p,
                                       AranParticle2MultipoleFunc2d *p2m,
                                       AranMultipole2MultipoleFunc2d *m2m,
                                       AranMultipole2LocalFunc2d *m2l,
                                       AranLocal2LocalFunc2d *l2l,
                                       AranLocal2ParticleFunc2d *l2p,
                                       AranParticle2LocalFunc2d *p2l,
                                       AranMultipole2ParticleFunc2d *m2p,
                                       guint *semifar_threshold)
{
  g_return_if_fail (solver != NULL);

  if (p2m) *p2m = solver->p2m;
  if (m2m) *m2m = solver->m2m;
  if (m2l) *m2l = solver->m2l;
  if (l2l) *l2l = solver->l2l;
  if (l2p) *l2p = solver->l2p;
  if (p2l) *p2l = solver->p2l;
  if (m2p) *m2p = solver->m2p;

  if (semifar_threshold) *semifar_threshold = solver->semifar_threshold;
}

void aran_solver2d_db_profile_operators (AranSolver2d *solver, gdouble order)
{
  g_return_if_fail (solver != NULL);

  if (solver->p2p)
    {
      solver->p2p_time = aran_profile_db_address_eval (solver->p2p, 0.);
    }

  if (solver->p2m)
    {
      solver->p2m_time = aran_profile_db_address_eval (solver->p2m, order);
    }

  if (solver->p2l)
    {
      solver->p2l_time = aran_profile_db_address_eval (solver->p2l, order);
    }

  if (solver->l2p)
    {
      solver->l2p_time = aran_profile_db_address_eval (solver->p2p, order);
    }

  if (solver->m2p)
    {
      solver->m2p_time = aran_profile_db_address_eval (solver->m2p, order);
    }

  if (solver->m2m)
    {
      solver->m2m_time = aran_profile_db_address_eval (solver->m2m, order);
    }

  if (solver->m2l)
    {
      solver->m2l_time = aran_profile_db_address_eval (solver->m2l, order);
    }

  if (solver->l2l)
    {
      solver->l2l_time = aran_profile_db_address_eval (solver->l2l, order);
    }

  g_printerr ("db_p2p_time %p = %g\n", solver->p2p, solver->p2p_time);
  g_printerr ("db_p2m_time %p = %g\n", solver->p2m, solver->p2m_time);
  g_printerr ("db_p2l_time %p = %g\n", solver->p2l, solver->p2l_time);
  g_printerr ("db_m2m_time %p = %g\n", solver->m2m, solver->m2m_time);
  g_printerr ("db_m2l_time %p = %g\n", solver->m2l, solver->m2l_time);
  g_printerr ("db_l2l_time %p = %g\n", solver->l2l, solver->l2l_time);
  g_printerr ("db_l2p_time %p = %g\n", solver->l2p, solver->l2p_time);
  g_printerr ("db_m2p_time %p = %g\n", solver->m2p, solver->m2p_time);

}
void aran_solver2d_profile_operators (AranSolver2d *solver,
                                      AranParticleInitFunc2d point_init,
                                      gpointer p1, gpointer p2)
{
  gdouble t;
  guint maxbox;
  VsgPRTree2dNodeInfo nodeinfo1 = {{0., 0.},};
  VsgPRTree2dNodeInfo nodeinfo2 = {{1., 1.},};
  gpointer dev1, dev2;

  g_return_if_fail (solver != NULL);

  maxbox = vsg_prtree2d_get_max_point (solver->prtree);

  if (solver->p2p)
    {
      t = aran_profile_nearfunc_2d (solver->p2p, point_init, p1, p2, maxbox);
      solver->p2p_time = t / (maxbox * maxbox);;
    }

  dev1 = g_boxed_copy (solver->devel_type, solver->devel);

  if (solver->p2m)
    {
      t = aran_profile_p2m_2d (solver->p2m, solver->zero, p1, &nodeinfo1, dev1);
      solver->p2m_time = t;
    }

  if (solver->l2p)
    {
      t = aran_profile_l2p_2d (solver->l2p, point_init, &nodeinfo1, dev1, p1);
      solver->l2p_time = t;
    }
  if (solver->p2l)
    {
      t = aran_profile_p2l_2d (solver->p2l, solver->zero, p1, &nodeinfo1, dev1);
      solver->p2l_time = t;
    }

  if (solver->m2p)
    {
      t = aran_profile_m2p_2d (solver->m2p, point_init, &nodeinfo1, dev1, p1);
      solver->m2p_time = t;
    }

  dev2 = g_boxed_copy (solver->devel_type, solver->devel);

  if (solver->m2m)
    {
      t = aran_profile_m2m_2d (solver->m2m, solver->zero, &nodeinfo1, dev1,
                               &nodeinfo2, dev2);
      solver->m2m_time = t;
    }

  if (solver->m2l)
    {
      t = aran_profile_m2l_2d (solver->m2l, solver->zero, &nodeinfo1, dev1,
                               &nodeinfo2, dev2);
      solver->m2l_time = t;
    }

  if (solver->l2l)
    {
      t = aran_profile_l2l_2d (solver->l2l, solver->zero, &nodeinfo1, dev1,
                               &nodeinfo2, dev2);
      solver->l2l_time = t;
    }

  g_boxed_free (solver->devel_type, dev1);
  g_boxed_free (solver->devel_type, dev2);

  /* g_printerr ("p2p_time = %g\n", solver->p2p_time); */
  /* g_printerr ("p2m_time = %g\n", solver->p2m_time); */
  /* g_printerr ("p2l_time = %g\n", solver->p2l_time); */
  /* g_printerr ("m2m_time = %g\n", solver->m2m_time); */
  /* g_printerr ("m2l_time = %g\n", solver->m2l_time); */
  /* g_printerr ("l2l_time = %g\n", solver->l2l_time); */
  /* g_printerr ("l2p_time = %g\n", solver->l2p_time); */
  /* g_printerr ("m2p_time = %g\n", solver->m2p_time); */

}

/**
 * aran_solver2d_reinit_stats:
 * @solver: an #AranSolver2d.
 *
 * Resets @solver Near/Far counters to zero.
 */
void aran_solver2d_reinit_stats (AranSolver2d *solver)
{
  g_return_if_fail (solver != NULL);

  solver->zero_counter = 0;
  solver->p2p_counter = 0;
  solver->p2p_remote_counter = 0;

  solver->p2m_counter = 0;
  solver->m2m_counter = 0;
  solver->m2l_counter = 0;
  solver->m2l_remote_counter = 0;
  solver->l2l_counter = 0;
  solver->l2p_counter = 0;
  solver->p2l_counter = 0;
  solver->m2p_counter = 0;
}

/**
 * aran_solver2d_get_stats:
 * @solver: an #AranSolver2d.
 * @zero: Clear function count result.
 * @p2p: particle 2 particle function count result.
 * @p2m: particle 2 multipole function count result.
 * @m2m: multipole 2 multipole function count result.
 * @m2l: multipole 2 local function count result.
 * @l2l: local 2 local function count result.
 * @l2p: local 2 particle function count result.
 * @p2l: particle 2 localfunction count result.
 * @m2p: multipole 2 particle function count result.
 * @p2p_remote: number of p2p calls with remote nodes.
 * @m2l_remote: number of m2l calls with remote nodes.
 *
 * Retrieves count values for @solver differents functions. Each value
 * represents the number of calls of the corresponding functino since
 * last call to aran_solver2d_reinit_stats().
 */
void aran_solver2d_get_stats (AranSolver2d *solver, glong *zero_count,
                              glong *p2p_count, glong *p2m_count,
                              glong *m2m_count, glong *m2l_count,
                              glong *l2l_count, glong *l2p_count,
                              glong *p2l_count, glong *m2p_count,
                              glong *p2p_remote_count, glong *m2l_remote_count)
{
  g_return_if_fail (solver != NULL);

  *zero_count = solver->zero_counter;
  *p2p_count = solver->p2p_counter;
  *p2p_remote_count = solver->p2p_remote_counter;
  *p2m_count = solver->p2m_counter;
  *m2m_count = solver->m2m_counter;
  *m2l_count = solver->m2l_counter;
  *m2l_remote_count = solver->m2l_remote_counter;
  *l2l_count = solver->l2l_counter;
  *l2p_count = solver->l2p_counter;
  *p2l_count = solver->p2l_counter;
  *m2p_count = solver->m2p_counter;
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
 * aran_solver2d_insert_point_local:
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
 * aran_solver2d_traverse:
 * @solver: an #AranSolver2d.
 * @func: traverse function.
 * @user_data: pointer to pass to @func.
 *
 * See vsg_prtree2d_traverse() in Vsg API docs.
 */
void aran_solver2d_traverse (AranSolver2d *solver,
                             GTraverseType order,
                             VsgPRTree2dFunc func,
                             gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_traverse (solver->prtree, order, func, user_data);
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

guint aran_solver2d_optimal_semifar_threshold (AranSolver2d *solver)
{
  guint semifar_threshold;

  g_return_val_if_fail (solver != NULL, G_MAXUINT);

  if (solver->p2l_time < 0. || solver->m2p_time < 0. || solver->p2p_time <= 0.)
    {
      semifar_threshold = G_MAXUINT;
    }
  else
    {
      semifar_threshold = (solver->p2l_time + solver->m2p_time) / solver->p2p_time;
    }

#ifdef VSG_HAVE_MPI
  /* semifar threshold values are to be consistent across all processors */
  {
    MPI_Comm communicator = vsg_prtree2d_get_communicator (solver->prtree);

    if (communicator != MPI_COMM_NULL)
      {
        guint tmp;
        gint sz;

        MPI_Comm_size (communicator, &sz);

        if (sz > 1)
          {
            MPI_Allreduce (&semifar_threshold, &tmp, 1, MPI_UNSIGNED, MPI_MAX, communicator);
            semifar_threshold = tmp;
          }
      }
  }
#endif

  return semifar_threshold;
}

/**
 * aran_solver2d_solve:
 * @solver: an #AranSolver2d.
 *
 * Solves the FMM problem for @solver.
 */
void aran_solver2d_solve (AranSolver2d *solver)
{
  VsgPRTree2dFarInteractionFunc far;
  VsgPRTree2dInteractionFunc near;
  VsgPRTree2dSemifarInteractionFunc semifar;

  g_return_if_fail (solver != NULL);

  VSG_TIMING_START (solve, vsg_prtree2d_get_communicator (solver->prtree));

  /*set interaction functions from solevr configuration */
  far = (VsgPRTree2dFarInteractionFunc)
    ((solver->m2l != NULL) ? far_func : nop_far_func);

  near = (VsgPRTree2dInteractionFunc)
    ((solver->p2p != NULL) ? near_func : nop_near_func);

  semifar = (VsgPRTree2dSemifarInteractionFunc)
    ((solver->p2l != NULL) && (solver->m2p != NULL) ? semifar_func : NULL);

  if (semifar != NULL && solver->semifar_threshold == 0)
    {
      /* find optimal semifar_threshold */
      solver->semifar_threshold = aran_solver2d_optimal_semifar_threshold (solver);

      g_printerr ("semifar threshold: %u\n", solver->semifar_threshold);
    }

  /* clear multipole and local developments before the big work */
  vsg_prtree2d_traverse (solver->prtree, G_POST_ORDER,
                         (VsgPRTree2dFunc) clear_func,
                         solver);

  VSG_TIMING_START (up, vsg_prtree2d_get_communicator (solver->prtree));

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

  VSG_TIMING_END (up, stderr);

  /* transmit info from Multipole to Local developments */
  vsg_prtree2d_near_far_traversal_full (solver->prtree, far, near, semifar,
                                        solver->semifar_threshold, solver);

  VSG_TIMING_START (down, vsg_prtree2d_get_communicator (solver->prtree));

  /* distribute information through Local developments towards particles */
  vsg_prtree2d_traverse (solver->prtree, G_PRE_ORDER,
                         (VsgPRTree2dFunc) down_func,
                         solver);

  VSG_TIMING_END (down, stderr);

  VSG_TIMING_END (solve, stderr);
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

MPI_Comm
aran_solver2d_get_communicator (AranSolver2d *solver)
{
  g_return_val_if_fail (solver != NULL, NULL);

  return vsg_prtree2d_get_communicator (solver->prtree);
}


void aran_solver2d_set_parallel (AranSolver2d *solver,
                                 VsgPRTreeParallelConfig *pconfig)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_set_parallel (solver->prtree, pconfig);
}

void aran_solver2d_get_parallel (AranSolver2d *solver,
                                 VsgPRTreeParallelConfig *pconfig)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_get_parallel (solver->prtree, pconfig);
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

void aran_solver2d_set_nf_isleaf (AranSolver2d *solver,
                                  VsgPRTree2dNFIsleafFunc isleaf,
                                  gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree2d_set_nf_isleaf (solver->prtree, isleaf, user_data);
}

