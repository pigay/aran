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

#include "aransolver3d.h"
#include "aranprofile.h"
#include "aranprofiledb.h"

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
  AranParticle2LocalFunc3d p2l;
  AranMultipole2ParticleFunc3d m2p;

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

#define _USE_G_SLICES GLIB_CHECK_VERSION (2, 10, 0)

#if ! _USE_G_SLICES

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
#endif /* ! _USE_G_SLICES */

static AranSolver3d *_solver3d_alloc ()
{
  AranSolver3d *solver;

#if _USE_G_SLICES
  solver = g_slice_new(AranSolver3d);
#else
  if (!aran_solver3d_mem_chunk)
    {
      aran_solver3d_mem_chunk = g_mem_chunk_create (AranSolver3d,
                                                    ARAN_SOLVER3D_PREALLOC,
                                                    G_ALLOC_ONLY);
    }

  aran_solver3d_instances_count ++;

  solver = g_chunk_new (AranSolver3d, aran_solver3d_mem_chunk);
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

  aran_solver3d_reinit_stats (solver);

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

static void _solver3d_dealloc (AranSolver3d *solver)
{
#if _USE_G_SLICES
  g_slice_free (AranSolver3d, solver);
#else
  g_chunk_free (solver, aran_solver3d_mem_chunk);

  aran_solver3d_instances_count --;

  if (aran_solver3d_instances_count == 0)
    aran_solver3d_finalize ();
#endif /* _USE_G_SLICES */
}

/* private function */
void aran_solver3d_init ()
{
#if ! _USE_G_SLICES
  static gboolean wasinit = FALSE;

  if (! wasinit)
    {
      wasinit = TRUE;
      g_atexit (aran_solver3d_finalize);
    }
#endif /* ! _USE_G_SLICES */
}


/*----------------------------------------------------*/
static void nop_near_func (const VsgPRTree3dNodeInfo *one_info,
		       const VsgPRTree3dNodeInfo *other_info,
		       AranSolver3d *solver)
{
}

/* general case near_func algorithm */
static void near_func_default (const VsgPRTree3dNodeInfo *one_info,
                               const VsgPRTree3dNodeInfo *other_info,
                               AranSolver3d *solver)
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
  if (VSG_PRTREE3D_NODE_INFO_IS_PRIVATE_REMOTE (one_info) ||
      VSG_PRTREE3D_NODE_INFO_IS_PRIVATE_REMOTE (other_info))
    solver->p2p_remote_counter +=one_info->point_count * other_info->point_count;
}

/* near_func algorithm for reflexive interaction (one_info == other_info) */
static void near_func_reflexive (const VsgPRTree3dNodeInfo *one_info,
                                 const VsgPRTree3dNodeInfo *other_info,
                                 AranSolver3d *solver)
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

  if (VSG_PRTREE3D_NODE_INFO_IS_PRIVATE_REMOTE (one_info))
    solver->p2p_remote_counter +=one_info->point_count * other_info->point_count;
}

static void near_func (const VsgPRTree3dNodeInfo *one_info,
                       const VsgPRTree3dNodeInfo *other_info,
                       AranSolver3d *solver)
{
  if (vsg_prtree_key3d_equals (&one_info->id, &other_info->id))
    near_func_reflexive (one_info, other_info, solver);
  else
    near_func_default (one_info, other_info, solver);
}


static void nop_far_func (const VsgPRTree3dNodeInfo *one_info,
                          const VsgPRTree3dNodeInfo *other_info,
                          AranSolver3d *solver)
{
}

static void far_func (const VsgPRTree3dNodeInfo *one_info,
                      const VsgPRTree3dNodeInfo *other_info,
                      AranSolver3d *solver)
{
  /* Multipole to Local transformation */
  if (solver->m2l != NULL)
    {
      gpointer one_dev = one_info->user_data;
      gpointer other_dev = other_info->user_data;

      /* both ways in order to get symmetric exchange */
      solver->m2l (one_info, one_dev, other_info, other_dev);

      solver->m2l (other_info, other_dev, one_info, one_dev);

      solver->m2l_counter += 2;

      if (VSG_PRTREE3D_NODE_INFO_IS_PRIVATE_REMOTE (one_info) ||
          VSG_PRTREE3D_NODE_INFO_IS_PRIVATE_REMOTE (other_info))
        solver->m2l_remote_counter += 2;
    }
}


static void semifar_func (const VsgPRTree3dNodeInfo *one_info,
                          const VsgPRTree3dNodeInfo *other_info,
                          AranSolver3d *solver)
{
  const VsgPRTree3dNodeInfo *info;
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

  /* g_printerr ("semifar one[%#lx %#lx %#lx %d] other[%#lx %#lx %#lx %d]\n", */
  /*             one_info->id.x, one_info->id.y, one_info->id.z, one_info->depth, */
  /*             other_info->id.x, other_info->id.y, other_info->id.z, other_info->depth); */

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


static void clear_func (const VsgPRTree3dNodeInfo *node_info,
                        AranSolver3d *solver)
{
  gpointer node_dev = node_info->user_data;

#ifdef VSG_HAVE_MPI
  if (VSG_PRTREE3D_NODE_INFO_IS_PRIVATE_REMOTE (node_info)) return;
#endif

  solver->zero (node_dev);
  solver->zero_counter ++;
}

static void up_func (const VsgPRTree3dNodeInfo *node_info,
                     AranSolver3d *solver)
{
  gpointer node_dev = node_info->user_data;

#ifdef VSG_HAVE_MPI
  if (VSG_PRTREE3D_NODE_INFO_IS_PRIVATE_REMOTE (node_info)) return;
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
  if (VSG_PRTREE3D_NODE_INFO_IS_PRIVATE_REMOTE (node_info)) return;
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
 * Associates @solver with a minimal set of FMM functions.
 */
void aran_solver3d_set_functions (AranSolver3d *solver,
                                  AranParticle2ParticleFunc3d p2p,
                                  AranParticle2MultipoleFunc3d p2m,
                                  AranMultipole2MultipoleFunc3d m2m,
                                  AranMultipole2LocalFunc3d m2l,
                                  AranLocal2LocalFunc3d l2l,
                                  AranLocal2ParticleFunc3d l2p)
{
  aran_solver3d_set_functions_full (solver, p2p, p2m, m2m, m2l, l2l, l2p,
                                    NULL, NULL, 0);
}

/**
 * aran_solver3d_set_functions_full:
 * @solver: an #AranSolver3d.
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
void aran_solver3d_set_functions_full (AranSolver3d *solver,
                                       AranParticle2ParticleFunc3d p2p,
                                       AranParticle2MultipoleFunc3d p2m,
                                       AranMultipole2MultipoleFunc3d m2m,
                                       AranMultipole2LocalFunc3d m2l,
                                       AranLocal2LocalFunc3d l2l,
                                       AranLocal2ParticleFunc3d l2p,
                                       AranParticle2LocalFunc3d p2l,
                                       AranMultipole2ParticleFunc3d m2p,
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

void aran_solver3d_get_functions_full (AranSolver3d *solver,
                                       AranParticle2ParticleFunc3d *p2p,
                                       AranParticle2MultipoleFunc3d *p2m,
                                       AranMultipole2MultipoleFunc3d *m2m,
                                       AranMultipole2LocalFunc3d *m2l,
                                       AranLocal2LocalFunc3d *l2l,
                                       AranLocal2ParticleFunc3d *l2p,
                                       AranParticle2LocalFunc3d *p2l,
                                       AranMultipole2ParticleFunc3d *m2p,
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

void aran_solver3d_db_profile_operators (AranSolver3d *solver, gdouble order)
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
void aran_solver3d_profile_operators (AranSolver3d *solver,
                                      AranParticleInitFunc3d point_init,
                                      gpointer p1, gpointer p2)
{
  gdouble t;
  guint maxbox;
  VsgPRTree3dNodeInfo nodeinfo1 = {{0., 0., 0.},};
  VsgPRTree3dNodeInfo nodeinfo2 = {{1., 1., 1.},};
  gpointer dev1, dev2;

  g_return_if_fail (solver != NULL);

  maxbox = vsg_prtree3d_get_max_point (solver->prtree);

  if (solver->p2p)
    {
      t = aran_profile_nearfunc_3d (solver->p2p, point_init, p1, p2, maxbox);
      solver->p2p_time = t / (maxbox * maxbox);;
    }

  dev1 = g_boxed_copy (solver->devel_type, solver->devel);

  if (solver->p2m)
    {
      t = aran_profile_p2m_3d (solver->p2m, solver->zero, p1, &nodeinfo1, dev1);
      solver->p2m_time = t;
    }

  if (solver->l2p)
    {
      t = aran_profile_l2p_3d (solver->l2p, point_init, &nodeinfo1, dev1, p1);
      solver->l2p_time = t;
    }
  if (solver->p2l)
    {
      t = aran_profile_p2l_3d (solver->p2l, solver->zero, p1, &nodeinfo1, dev1);
      solver->p2l_time = t;
    }

  if (solver->m2p)
    {
      t = aran_profile_m2p_3d (solver->m2p, point_init, &nodeinfo1, dev1, p1);
      solver->m2p_time = t;
    }

  dev2 = g_boxed_copy (solver->devel_type, solver->devel);

  if (solver->m2m)
    {
      t = aran_profile_m2m_3d (solver->m2m, solver->zero, &nodeinfo1, dev1,
                               &nodeinfo2, dev2);
      solver->m2m_time = t;
    }

  if (solver->m2l)
    {
      t = aran_profile_m2l_3d (solver->m2l, solver->zero, &nodeinfo1, dev1,
                               &nodeinfo2, dev2);
      solver->m2l_time = t;
    }

  if (solver->l2l)
    {
      t = aran_profile_l2l_3d (solver->l2l, solver->zero, &nodeinfo1, dev1,
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
 * aran_solver3d_get_stats:
 * @solver: an #AranSolver3d.
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
 * last call to aran_solver3d_reinit_stats().
 */
void aran_solver3d_get_stats (AranSolver3d *solver, glong *zero_count,
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

guint aran_solver3d_optimal_semifar_threshold (AranSolver3d *solver)
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
    MPI_Comm communicator = vsg_prtree3d_get_communicator (solver->prtree);

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
 * aran_solver3d_solve:
 * @solver: an #AranSolver3d.
 *
 * Solves the FMM problem for @solver.
 */
void aran_solver3d_solve (AranSolver3d *solver)
{
  VsgPRTree3dFarInteractionFunc far;
  VsgPRTree3dInteractionFunc near;
  VsgPRTree3dSemifarInteractionFunc semifar;

  g_return_if_fail (solver != NULL);

  VSG_TIMING_START (solve, vsg_prtree3d_get_communicator (solver->prtree));

  /*set interaction functions from solevr configuration */
  far = (VsgPRTree3dFarInteractionFunc)
    ((solver->m2l != NULL) ? far_func : nop_far_func);

  near = (VsgPRTree3dInteractionFunc)
    ((solver->p2p != NULL) ? near_func : nop_near_func);

  semifar = (VsgPRTree3dSemifarInteractionFunc)
    ((solver->p2l != NULL) && (solver->m2p != NULL) ? semifar_func : NULL);

  if (semifar != NULL && solver->semifar_threshold == 0)
    {
      /* find optimal semifar_threshold */
      solver->semifar_threshold = aran_solver3d_optimal_semifar_threshold (solver);

      g_printerr ("semifar threshold: %u\n", solver->semifar_threshold);
    }

  /* clear multipole and local developments before the big work */
  vsg_prtree3d_traverse (solver->prtree, G_POST_ORDER,
                         (VsgPRTree3dFunc) clear_func,
                         solver);

  VSG_TIMING_START (up, vsg_prtree3d_get_communicator (solver->prtree));

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

  VSG_TIMING_END (up, stderr);

  /* transmit info from Multipole to Local developments */
  vsg_prtree3d_near_far_traversal_full (solver->prtree, far, near, semifar,
                                        solver->semifar_threshold, solver);

  VSG_TIMING_START (down, vsg_prtree3d_get_communicator (solver->prtree));

  /* distribute information through Local developments towards particles */
  vsg_prtree3d_traverse (solver->prtree, G_PRE_ORDER,
                         (VsgPRTree3dFunc) down_func,
                         solver);

  VSG_TIMING_END (down, stderr);

  VSG_TIMING_END (solve, stderr);
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

MPI_Comm
aran_solver3d_get_communicator (AranSolver3d *solver)
{
  g_return_val_if_fail (solver != NULL, NULL);

  return vsg_prtree3d_get_communicator (solver->prtree);
}


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


void aran_solver3d_set_nf_isleaf (AranSolver3d *solver,
                                  VsgPRTree3dNFIsleafFunc isleaf,
                                  gpointer user_data)
{
  g_return_if_fail (solver != NULL);

  vsg_prtree3d_set_nf_isleaf (solver->prtree, isleaf, user_data);
}
