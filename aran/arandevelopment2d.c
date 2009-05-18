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

#include "arandevelopment2d.h"

/**
 * ARAN_TYPE_DEVELOPMENT2D:
 *
 * #AranDevelopment2d #GBoxed #GType.
 */

/**
 * AranDevelopment2d:
 * @multipole: Multipole expansion.
 * @local: Local expansion.
 *
 * A structure used as #VsgPRTree2d node_data within an #AranSolver2d.
 */

/* functions */

GType aran_development2d_get_type ()
{
  static GType devtype = G_TYPE_NONE;

  if (G_UNLIKELY (devtype == G_TYPE_NONE))
    {
      devtype =
	g_boxed_type_register_static ("AranDevelopment2d",
				      (GBoxedCopyFunc) aran_development2d_clone,
				      (GBoxedFreeFunc) aran_development2d_free);
    }

  return devtype;
}

/**
 * aran_development2d_new:
 * @posdeg: a #guint8.
 * @negdeg: a #guint8.
 *
 * Allocates a new #AranDeveloment2d structure.
 *
 * Returns: newly allocated structure.
 */
AranDevelopment2d *aran_development2d_new (guint8 posdeg, guint8 negdeg)
{
  AranDevelopment2d *result =
    (AranDevelopment2d *) g_malloc (sizeof (AranDevelopment2d));

  result->multipole = aran_laurent_seriesd_new (posdeg, negdeg);
  result->local = aran_laurent_seriesd_new (MAX (posdeg, negdeg), 0);

  return result;
}

/**
 * aran_development2d_free:
 * @ad: an #AranDevelopment2d.
 *
 * Deallocates @ad and all associated memory.
 */
void aran_development2d_free (AranDevelopment2d *ad)
{
  aran_laurent_seriesd_free (ad->multipole);
  aran_laurent_seriesd_free (ad->local);
  g_free (ad);
}

/**
 * aran_development2d_copy:
 * @src: an #AranDevelopment2d.
 * @dst: an #AranDevelopment2d.
 *
 * Copies @src into @dst. Development degrees must coincide.
 */
void aran_development2d_copy (const AranDevelopment2d *src,
                              AranDevelopment2d *dst)
{
  g_return_if_fail (src != NULL);
  g_return_if_fail (dst != NULL);

  aran_laurent_seriesd_copy (src->multipole, dst->multipole);
  aran_laurent_seriesd_copy (src->local, dst->local);
}

/**
 * aran_development2d_clone:
 * @src: an #AranDevelopment2d.
 *
 * Duplicates @src.
 *
 * Returns: newly allocated copy of @src.
 */
AranDevelopment2d *aran_development2d_clone (AranDevelopment2d *src)
{
  AranDevelopment2d *dst;

  g_return_val_if_fail (src != NULL, NULL);

  dst =
    aran_development2d_new (aran_laurent_seriesd_get_posdeg (src->multipole),
			    aran_laurent_seriesd_get_negdeg (src->multipole));

  aran_development2d_copy (src, dst);

  return dst;
}

/**
 * aran_development2d_set_zero:
 * @ad: an #AranDevelopment2d.
 *
 * Sets all @ad coefficients to zero.
 */
void aran_development2d_set_zero (AranDevelopment2d *ad)
{
  g_return_if_fail (ad != NULL);

  aran_laurent_seriesd_set_zero (ad->multipole);
  aran_laurent_seriesd_set_zero (ad->local);
}

/**
 * aran_development2d_write:
 * @ad: an #AranDevelopment2d.
 * @file: output file.
 *
 * Writes @ad to @file.
 */
void aran_development2d_write (AranDevelopment2d *ad, FILE *file)
{
  g_return_if_fail (ad != NULL);

  fprintf (file, "{multipole= ");
  aran_laurent_seriesd_write (ad->multipole, file);
  fprintf (file, ", local= ");
  aran_laurent_seriesd_write (ad->local, file);
  fprintf (file, "}");
}

/**
 * aran_development2d_m2m:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment2d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment2d.
 *
 * Performs multipole 2 multipole translation between @src and @dst.
 */
void aran_development2d_m2m (const VsgPRTree2dNodeInfo *src_node,
			     AranDevelopment2d *src,
			     const VsgPRTree2dNodeInfo *dst_node,
			     AranDevelopment2d *dst)
{
  gcomplex128 zsrc = src_node->center.x + G_I*src_node->center.y;
  gcomplex128 zdst = dst_node->center.x + G_I*dst_node->center.y;

  aran_laurent_seriesd_translate (src->multipole, zsrc, dst->multipole, zdst);
}

/**
 * aran_development2d_m2l:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment2d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment2d.
 *
 * Performs multipole 2 local translation between @src and @dst.
 *
 * Returns: %TRUE.
 */
gboolean aran_development2d_m2l (const VsgPRTree2dNodeInfo *src_node,
                                 AranDevelopment2d *src,
                                 const VsgPRTree2dNodeInfo *dst_node,
                                 AranDevelopment2d *dst)
{
  gcomplex128 zsrc = src_node->center.x + G_I*src_node->center.y;
  gcomplex128 zdst = dst_node->center.x + G_I*dst_node->center.y;

  aran_laurent_seriesd_to_taylor (src->multipole, zsrc, dst->local, zdst);

  return TRUE;
}

/**
 * aran_development2d_l2l:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment2d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment2d.
 *
 * Performs local 2 local translation between @src and @dst.
 */
void aran_development2d_l2l (const VsgPRTree2dNodeInfo *src_node,
			     AranDevelopment2d *src,
			     const VsgPRTree2dNodeInfo *dst_node,
			     AranDevelopment2d *dst)
{
  gcomplex128 zsrc = src_node->center.x + G_I*src_node->center.y;
  gcomplex128 zdst = dst_node->center.x + G_I*dst_node->center.y;

  aran_laurent_seriesd_translate (src->local, zsrc, dst->local, zdst);
}

/**
 * aran_development2d_multipole_evaluate:
 * @devel_node: tree node info of @devel.
 * @devel: an #AranDevelopment2d.
 * @pos: evaluation position.
 *
 * Evaluates the multipole part of @devel at @pos.
 *
 * Returns: value of multipole part of @devel(@pos).
 */
gcomplex128
aran_development2d_multipole_evaluate (const VsgPRTree2dNodeInfo *devel_node,
                                       AranDevelopment2d *devel,
                                       const VsgVector2d *pos)
{
  gcomplex128 zl = devel_node->center.x + G_I*devel_node->center.y;
  gcomplex128 zp = pos->x + G_I*pos->y;

  return aran_laurent_seriesd_evaluate (devel->multipole, zp-zl);
}

/**
 * aran_development2d_local_evaluate:
 * @devel_node: tree node info of @devel.
 * @devel: an #AranDevelopment2d.
 * @pos: evaluation position.
 *
 * Evaluates the local part of @devel at @pos.
 *
 * Returns: value of local part of @devel(@pos).
 */
gcomplex128
aran_development2d_local_evaluate (const VsgPRTree2dNodeInfo *devel_node,
                                   AranDevelopment2d *devel,
                                   const VsgVector2d *pos)
{
  gcomplex128 zl = devel_node->center.x + G_I*devel_node->center.y;
  gcomplex128 zp = pos->x + G_I*pos->y;

  return aran_laurent_seriesd_evaluate (devel->local, zp-zl);
}

#ifdef VSG_HAVE_MPI

void aran_development2d_vtable_init (VsgParallelVTable *vtable, guint8 posdeg,
                                     guint8 negdeg)
{
  vtable->alloc = (VsgMigrableAllocDataFunc) aran_development2d_alloc;
  vtable->alloc_data = aran_development2d_new ((posdeg), (negdeg));

  vtable->destroy = aran_development2d_destroy;
  vtable->destroy_data = NULL;

  vtable->migrate.pack =
    (VsgMigrablePackDataFunc) aran_development2d_migrate_pack;
  vtable->migrate.pack_data = NULL;

  vtable->migrate.unpack =
    (VsgMigrablePackDataFunc) aran_development2d_migrate_unpack;
  vtable->migrate.unpack_data = NULL;

  vtable->visit_forward.pack =
    (VsgMigrablePackDataFunc) aran_development2d_visit_fw_pack;
  vtable->visit_forward.pack_data = NULL;

  vtable->visit_forward.unpack =
    (VsgMigrablePackDataFunc) aran_development2d_visit_fw_unpack;
  vtable->visit_forward.unpack_data = NULL;

  vtable->visit_forward.reduce =
    (VsgMigrableReductionDataFunc) aran_development2d_visit_fw_reduce;
  vtable->visit_forward.reduce_data = NULL;

  vtable->visit_backward.pack =
    (VsgMigrablePackDataFunc) aran_development2d_visit_bw_pack;
  vtable->visit_backward.pack_data = NULL;

  vtable->visit_backward.unpack =
    (VsgMigrablePackDataFunc) aran_development2d_visit_bw_unpack;
  vtable->visit_backward.unpack_data = NULL;

  vtable->visit_backward.reduce =
    (VsgMigrableReductionDataFunc) aran_development2d_visit_bw_reduce;
  vtable->visit_backward.reduce_data = NULL;
}

void aran_development2d_vtable_clear (VsgParallelVTable *vtable)
{
  g_return_if_fail (vtable != NULL);
  g_return_if_fail (vtable->alloc_data != NULL);

  aran_development2d_free (vtable->alloc_data);
}

#define ARAN_DEVELOPMENT2D_VTABLE_DESTROY(vtable) \
  {aran_development2d_free ((vtable).alloc_data);}

/**
 * aran_development2d_alloc:
 * @resident: unused.
 * @src: an example #AranDevelopment2d to copy from.
 *
 * Allocates a new #AranDevelopment2d by clonig @src.
 *
 * Returns: a copy of @src.
 */
gpointer aran_development2d_alloc (gboolean resident, AranDevelopment2d *src)
{
  gpointer ret;

  ret = g_boxed_copy (ARAN_TYPE_DEVELOPMENT2D, src);

  return ret;
}

/**
 * aran_development2d_destroy:
 * @data: A #AranDevelopment2d.
 * @resident: unused.
 * @user_data: unused.
 *
 * Deletes @data from memory.
 */
void aran_development2d_destroy (gpointer data, gboolean resident,
                                gpointer user_data)
{
  g_assert (data != NULL);
  g_boxed_free (ARAN_TYPE_DEVELOPMENT2D, data);
}

/**
 * aran_development2d_local_evaluate:
 * @devel: an #AranDevelopment2d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Performs complete packing of @devel into @pm for a migration
 * between processors.
 */
void aran_development2d_migrate_pack (AranDevelopment2d *devel,
                                      VsgPackedMsg *pm,
                                      gpointer user_data)
{
  aran_laurent_seriesd_pack (devel->multipole, pm);
  aran_laurent_seriesd_pack (devel->local, pm);
}

/**
 * aran_development2d_migrate_unpack:
 * @devel: an #AranDevelopment2d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Unpacks information from @pm in a migration between processors and
 * stores it in @devel.
 */
void aran_development2d_migrate_unpack (AranDevelopment2d *devel,
                                        VsgPackedMsg *pm,
                                        gpointer user_data)
{
  aran_laurent_seriesd_unpack (devel->multipole, pm);
  aran_laurent_seriesd_unpack (devel->local, pm);
}

/**
 * aran_development2d_visit_fw_pack:
 * @devel: an #AranDevelopment2d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Performs packing of @devel into @pm for a near/far forward visit
 * of a local node to another processor.
 */
void aran_development2d_visit_fw_pack (AranDevelopment2d *devel,
                                       VsgPackedMsg *pm,
                                       gpointer user_data)
{
  aran_laurent_seriesd_pack (devel->multipole, pm);
}

/**
 * aran_development2d_visit_fw_unpack:
 * @devel: an #AranDevelopment2d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Unpacks information from @pm for a near/far forward visit of a
 * remote node.
 */
void aran_development2d_visit_fw_unpack (AranDevelopment2d *devel,
                                         VsgPackedMsg *pm,
                                         gpointer user_data)
{
  aran_laurent_seriesd_unpack (devel->multipole, pm);
}

/**
 * aran_development2d_visit_fw_reduce:
 * @a: source #AranDevelopment2d.
 * @b: destination #AranDevelopment2d.
 * @user_data: unused.
 *
 * Forward visit reduction operator for #AranDevelopment2d.
 */
void aran_development2d_visit_fw_reduce (AranDevelopment2d *a,
                                         AranDevelopment2d *b,
                                         gpointer user_data)
{
  aran_laurent_seriesd_add (a->multipole, b->multipole, b->multipole);
}

/**
 * aran_development2d_visit_bw_pack:
 * @devel: an #AranDevelopment2d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Performs packing of @devel into @pm for a near/far backward visit
 * of a remote node to its original processor.
 */
void aran_development2d_visit_bw_pack (AranDevelopment2d *devel,
                                       VsgPackedMsg *pm,
                                       gpointer user_data)

{
  aran_laurent_seriesd_pack (devel->local, pm);
}

/**
 * aran_development2d_visit_bw_unpack:
 * @devel: an #AranDevelopment2d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Unpacks information from @pm for a near/far backward visit of a
 * remote node.
 */
void aran_development2d_visit_bw_unpack (AranDevelopment2d *devel,
                                         VsgPackedMsg *pm,
                                         gpointer user_data)
{
  aran_laurent_seriesd_unpack (devel->local, pm);
}

/**
 * aran_development2d_visit_bw_reduce:
 * @a: source #AranDevelopment2d.
 * @b: destination #AranDevelopment2d.
 * @user_data: unused.
 *
 * Backrward visit reduction operator for #AranDevelopment2d.
 */
void aran_development2d_visit_bw_reduce (AranDevelopment2d *a,
                                         AranDevelopment2d *b,
                                         gpointer user_data)
{
  aran_laurent_seriesd_add (a->local, b->local, b->local);

}

#endif /* VSG_HAVE_MPI */
