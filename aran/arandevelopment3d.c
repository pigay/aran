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

#include "arandevelopment3d.h"

/* #include <vsg/vsgd-inline.h> */

/**
 * ARAN_TYPE_DEVELOPMENT3D:
 *
 * #AranDevelopment3d #GBoxed #GType.
 */

/**
 * AranDevelopment3d:
 * @multipole: Multipole expansion.
 * @local: Local expansion.
 *
 * A structure used as #VsgPRTree3d node_data within an #AranSolver3d.
 */

/* functions */

GType aran_development3d_get_type ()
{
  static GType devtype = G_TYPE_NONE;

  if (G_UNLIKELY (devtype == G_TYPE_NONE))
    {
      devtype =
	g_boxed_type_register_static ("AranDevelopment3d",
				      (GBoxedCopyFunc) aran_development3d_clone,
				      (GBoxedFreeFunc) aran_development3d_free);
    }

  return devtype;
}

/**
 * aran_development3d_new:
 * @posdeg: a #guint8.
 * @negdeg: a #guint8.
 *
 * Allocates a new #AranDeveloment3d structure.
 *
 * Returns: newly allocated structure.
 */
AranDevelopment3d *aran_development3d_new (guint8 posdeg, guint8 negdeg)
{
  AranDevelopment3d *result =
    (AranDevelopment3d *) g_malloc (sizeof (AranDevelopment3d));

  result->multipole = aran_spherical_seriesd_new (posdeg, negdeg);
  result->local = aran_spherical_seriesd_new (MAX (posdeg, negdeg), 0);

  return result;
}

/**
 * aran_development3d_free:
 * @ad: an #AranDevelopment3d.
 *
 * Deallocates @ad and all associated memory.
 */
void aran_development3d_free (AranDevelopment3d *ad)
{
  aran_spherical_seriesd_free (ad->multipole);
  aran_spherical_seriesd_free (ad->local);
  g_free (ad);
}

/**
 * aran_development3d_copy:
 * @src: an #AranDevelopment3d.
 * @dst: an #AranDevelopment3d.
 *
 * Copies @src into @dst. Development degrees must coincide.
 */
void aran_development3d_copy (const AranDevelopment3d *src,
                              AranDevelopment3d *dst)
{
  g_return_if_fail (src != NULL);
  g_return_if_fail (dst != NULL);

  aran_spherical_seriesd_copy (src->multipole, dst->multipole);
  aran_spherical_seriesd_copy (src->local, dst->local);
}

/**
 * aran_development3d_clone:
 * @src: an #AranDevelopment3d.
 *
 * Duplicates @src.
 *
 * Returns: newly allocated copy of @src.
 */
AranDevelopment3d *aran_development3d_clone (AranDevelopment3d *src)
{
  AranDevelopment3d *dst;

  g_return_val_if_fail (src != NULL, NULL);

  dst =
    aran_development3d_new (aran_spherical_seriesd_get_posdeg (src->multipole),
			    aran_spherical_seriesd_get_negdeg (src->multipole));

  aran_development3d_copy (src, dst);

  return dst;
}

/**
 * aran_development3d_set_zero:
 * @ad: an #AranDevelopment3d.
 *
 * Sets all @ad coefficients to zero.
 */
void aran_development3d_set_zero (AranDevelopment3d *ad)
{
  g_return_if_fail (ad != NULL);

  aran_spherical_seriesd_set_zero (ad->multipole);
  aran_spherical_seriesd_set_zero (ad->local);
}

/**
 * aran_development3d_write:
 * @ad: an #AranDevelopment3d.
 * @file: output file.
 *
 * Writes @ad to @file.
 */
void aran_development3d_write (AranDevelopment3d *ad, FILE *file)
{
  g_return_if_fail (ad != NULL);

  fprintf (file, "{multipole= ");
  aran_spherical_seriesd_write (ad->multipole, file);
  fprintf (file, ", local= ");
  aran_spherical_seriesd_write (ad->local, file);
  fprintf (file, "}");
}

/**
 * aran_development3d_m2m:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment3d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 multipole translation between @src and @dst.
 */
void aran_development3d_m2m (const VsgPRTree3dNodeInfo *src_node,
			     AranDevelopment3d *src,
			     const VsgPRTree3dNodeInfo *dst_node,
			     AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate (src->multipole, &src_node->center,
				    dst->multipole, &dst_node->center);
}

/**
 * aran_development3d_m2l:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment3d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 local translation between @src and @dst.
 *
 * Returns: %TRUE.
 */
gboolean aran_development3d_m2l (const VsgPRTree3dNodeInfo *src_node,
				 AranDevelopment3d *src,
				 const VsgPRTree3dNodeInfo *dst_node,
				 AranDevelopment3d *dst)
{
  aran_spherical_seriesd_to_local (src->multipole, &src_node->center,
				   dst->local, &dst_node->center);

  return TRUE;
}

/**
 * aran_development3d_l2l:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment3d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment3d.
 *
 * Performs local 2 local translation between @src and @dst.
 */
void aran_development3d_l2l (const VsgPRTree3dNodeInfo *src_node,
			     AranDevelopment3d *src,
			     const VsgPRTree3dNodeInfo *dst_node,
			     AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate (src->local, &src_node->center,
				    dst->local, &dst_node->center);
}

/**
 * aran_development3d_m2m_kkylin:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment3d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 multipole translation between @src and @dst with the
 * K. Kylin formulas.
 */
void aran_development3d_m2m_kkylin (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate_kkylin (src->multipole, &src_node->center,
                                           dst->multipole, &dst_node->center);
}

/**
 * aran_development3d_m2l_kkylin:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment3d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 local translation between @src and @dst with the
 * K. Kylin formulas.
 *
 * Returns: %TRUE.
 */
gboolean aran_development3d_m2l_kkylin (const VsgPRTree3dNodeInfo *src_node,
					AranDevelopment3d *src,
					const VsgPRTree3dNodeInfo *dst_node,
					AranDevelopment3d *dst)
{
  aran_spherical_seriesd_to_local_kkylin (src->multipole, &src_node->center,
                                          dst->local, &dst_node->center);

  return TRUE;
}

/**
 * aran_development3d_l2l_kkylin:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment3d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment3d.
 *
 * Performs local 2 local translation between @src and @dst with the
 * K. Kylin formulas.
 */
void aran_development3d_l2l_kkylin (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate_kkylin (src->local, &src_node->center,
                                           dst->local, &dst_node->center);
}

/**
 * aran_development3d_m2m_rotate:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment3d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 multipole translation between @src and @dst with the
 * "Point and Shoot" technique: rotate @src to @dst_node->center direction and
 * then translate the rotated series along the Z axis. Rotate it back and
 * accumulate it into dst. These operations lead to a O(p^3) scheme.
 */
void aran_development3d_m2m_rotate (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate_rotate (src->multipole, &src_node->center,
                                           dst->multipole, &dst_node->center);
}

/**
 * aran_development3d_m2l_rotate:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment3d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 local translation between @src and @dst with the
 * "Point and Shoot" technique: rotate @src to @dst_node->center direction and
 * then translate the rotated series along the Z axis. Rotate it back and
 * accumulate it into dst. These operations lead to a O(p^3) scheme.
 *
 * Returns: %TRUE.
 */
gboolean aran_development3d_m2l_rotate (const VsgPRTree3dNodeInfo *src_node,
					AranDevelopment3d *src,
					const VsgPRTree3dNodeInfo *dst_node,
					AranDevelopment3d *dst)
{
  aran_spherical_seriesd_to_local_rotate (src->multipole, &src_node->center,
                                          dst->local, &dst_node->center);

  return TRUE;
}

/**
 * aran_development3d_l2l_rotate:
 * @src_node: @src tree node info.
 * @src: an #AranDevelopment3d.
 * @dst_node: @dst tree node info.
 * @dst: an #AranDevelopment3d.
 *
 * Performs local 2 local translation between @src and @dst with the
 * "Point and Shoot" technique: rotate @src to @dst_node->center direction and
 * then translate the rotated series along the Z axis. Rotate it back and
 * accumulate it into dst. These operations lead to a O(p^3) scheme.
 */
void aran_development3d_l2l_rotate (const VsgPRTree3dNodeInfo *src_node,
                                    AranDevelopment3d *src,
                                    const VsgPRTree3dNodeInfo *dst_node,
                                    AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate_rotate (src->local, &src_node->center,
                                           dst->local, &dst_node->center);
}

/**
 * aran_development3d_multipole_evaluate:
 * @devel_node: tree node info of @devel.
 * @devel: an #AranDevelopment3d.
 * @pos: evaluation position.
 *
 * Evaluates the multipole part of @devel at @pos.
 *
 * Returns: value of multipole part of @devel(@pos).
 */
gcomplex128
aran_development3d_multipole_evaluate (const VsgPRTree3dNodeInfo *devel_node,
                                       AranDevelopment3d *devel,
                                       const VsgVector3d *pos)
{
  VsgVector3d tmp;

  vsg_vector3d_sub (pos, &devel_node->center, &tmp);

  return aran_spherical_seriesd_evaluate (devel->multipole, &tmp);
}

/**
 * aran_development3d_local_evaluate:
 * @devel_node: tree node info of @devel.
 * @devel: an #AranDevelopment3d.
 * @pos: evaluation position.
 *
 * Evaluates the local part of @devel at @pos.
 *
 * Returns: value of local part of @devel(@pos).
 */
gcomplex128 aran_development3d_local_evaluate (const VsgPRTree3dNodeInfo *devel_node,
					       AranDevelopment3d *devel,
					       const VsgVector3d *pos)
{
  VsgVector3d tmp;

  vsg_vector3d_sub (pos, &devel_node->center, &tmp);

  return aran_spherical_seriesd_evaluate (devel->local, &tmp);
}

/**
 * aran_development3d_local_gradient_evaluate:
 * @devel_node: tree node info of @devel.
 * @devel: an #AranDevelopment3d.
 * @pos: evaluation position.
 * @grad: result gradient.
 *
 * Evaluates the gradient of the local part of @devel at @pos.
 */
void
aran_development3d_local_gradient_evaluate (const VsgPRTree3dNodeInfo *devel_node,
                                            AranDevelopment3d *devel,
                                            const VsgVector3d *pos,
                                            VsgVector3d *grad)
{
  VsgVector3d tmp;

  vsg_vector3d_sub (pos, &devel_node->center, &tmp);

  aran_spherical_seriesd_local_gradient_evaluate (devel->local, &tmp, grad);
}

#ifdef VSG_HAVE_MPI

#include <aran/aransphericalseriesd-private.h>

void aran_development3d_vtable_init (VsgParallelVTable *vtable, guint8 posdeg,
                                     guint8 negdeg)
{
  vtable->alloc = (VsgMigrableAllocDataFunc) aran_development3d_alloc;
  vtable->alloc_data = aran_development3d_new ((posdeg), (negdeg));

  vtable->destroy = aran_development3d_destroy;
  vtable->destroy_data = NULL;

  vtable->migrate.pack =
    (VsgMigrablePackDataFunc) aran_development3d_migrate_pack;
  vtable->migrate.pack_data = NULL;

  vtable->migrate.unpack =
    (VsgMigrablePackDataFunc) aran_development3d_migrate_unpack;
  vtable->migrate.unpack_data = NULL;

  vtable->visit_forward.pack =
    (VsgMigrablePackDataFunc) aran_development3d_visit_fw_pack;
  vtable->visit_forward.pack_data = NULL;

  vtable->visit_forward.unpack =
    (VsgMigrablePackDataFunc) aran_development3d_visit_fw_unpack;
  vtable->visit_forward.unpack_data = NULL;

  vtable->visit_forward.reduce =
    (VsgMigrableReductionDataFunc) aran_development3d_visit_fw_reduce;
  vtable->visit_forward.reduce_data = NULL;

  vtable->visit_backward.pack =
    (VsgMigrablePackDataFunc) aran_development3d_visit_bw_pack;
  vtable->visit_backward.pack_data = NULL;

  vtable->visit_backward.unpack =
    (VsgMigrablePackDataFunc) aran_development3d_visit_bw_unpack;
  vtable->visit_backward.unpack_data = NULL;

  vtable->visit_backward.reduce =
    (VsgMigrableReductionDataFunc) aran_development3d_visit_bw_reduce;
  vtable->visit_backward.reduce_data = NULL;

/*   { */
/*     AranSphericalSeriesd *multipole = */
/*       ((AranDevelopment3d *) vtable->alloc_data)->multipole; */
/*     AranSphericalSeriesd *local = */
/*       ((AranDevelopment3d *) vtable->alloc_data)->local; */
/*     g_printerr ("dev3d msg size M=%lu L=%lu\n",  */
/* 		_spherical_seriesd_size (multipole->posdeg, multipole->negdeg), */
/* 		_spherical_seriesd_size (local->posdeg, local->negdeg)); */
/*   } */
}

void aran_development3d_vtable_clear (VsgParallelVTable *vtable)
{
  g_return_if_fail (vtable != NULL);
  g_return_if_fail (vtable->alloc_data != NULL);

  aran_development3d_free (vtable->alloc_data);
}

#define ARAN_DEVELOPMENT3D_VTABLE_DESTROY(vtable) \
  {aran_development3d_free ((vtable).alloc_data);}

/**
 * aran_development3d_alloc:
 * @resident: unused.
 * @src: an example #AranDevelopment3d to copy from.
 *
 * Allocates a new #AranDevelopment3d by clonig @src.
 *
 * Returns: a copy of @src.
 */
gpointer aran_development3d_alloc (gboolean resident, AranDevelopment3d *src)
{
  gpointer ret;

  ret = g_boxed_copy (ARAN_TYPE_DEVELOPMENT3D, src);

  return ret;
}

/**
 * aran_development3d_destroy:
 * @data: A #AranDevelopment3d.
 * @resident: unused.
 * @user_data: unused.
 *
 * Deletes @data from memory.
 */
void aran_development3d_destroy (gpointer data, gboolean resident,
                                gpointer user_data)
{
  g_assert (data != NULL);
  g_boxed_free (ARAN_TYPE_DEVELOPMENT3D, data);
}

/**
 * aran_development3d_local_evaluate:
 * @devel: an #AranDevelopment3d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Performs complete packing of @devel into @pm for a migration
 * between processors.
 */
void aran_development3d_migrate_pack (AranDevelopment3d *devel,
                                      VsgPackedMsg *pm,
                                      gpointer user_data)
{
  aran_spherical_seriesd_pack (devel->multipole, pm);
  aran_spherical_seriesd_pack (devel->local, pm);
}

/**
 * aran_development3d_migrate_unpack:
 * @devel: an #AranDevelopment3d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Unpacks information from @pm in a migration between processors and
 * stores it in @devel.
 */
void aran_development3d_migrate_unpack (AranDevelopment3d *devel,
                                        VsgPackedMsg *pm,
                                        gpointer user_data)
{
  aran_spherical_seriesd_unpack (devel->multipole, pm);
  aran_spherical_seriesd_unpack (devel->local, pm);
}

/**
 * aran_development3d_visit_fw_pack:
 * @devel: an #AranDevelopment3d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Performs packing of @devel into @pm for a near/far forward visit
 * of a local node to another processor.
 */
void aran_development3d_visit_fw_pack (AranDevelopment3d *devel,
                                       VsgPackedMsg *pm,
                                       gpointer user_data)
{
  aran_spherical_seriesd_pack (devel->multipole, pm);
}

/**
 * aran_development3d_visit_fw_unpack:
 * @devel: an #AranDevelopment3d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Unpacks information from @pm for a near/far forward visit of a
 * remote node.
 */
void aran_development3d_visit_fw_unpack (AranDevelopment3d *devel,
                                         VsgPackedMsg *pm,
                                         gpointer user_data)
{
  aran_spherical_seriesd_unpack (devel->multipole, pm);
}

/**
 * aran_development3d_visit_fw_reduce:
 * @a: source #AranDevelopment3d.
 * @b: destination #AranDevelopment3d.
 * @user_data: unused.
 *
 * Forward visit reduction operator for #AranDevelopment3d.
 */
void aran_development3d_visit_fw_reduce (AranDevelopment3d *a,
                                         AranDevelopment3d *b,
                                         gpointer user_data)
{
  aran_spherical_seriesd_add (a->multipole, b->multipole, b->multipole);
}

/**
 * aran_development3d_visit_bw_pack:
 * @devel: an #AranDevelopment3d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Performs packing of @devel into @pm for a near/far backward visit
 * of a remote node to its original processor.
 */
void aran_development3d_visit_bw_pack (AranDevelopment3d *devel,
                                       VsgPackedMsg *pm,
                                       gpointer user_data)

{
  aran_spherical_seriesd_pack (devel->local, pm);
}

/**
 * aran_development3d_visit_bw_unpack:
 * @devel: an #AranDevelopment3d.
 * @pm: a #VsgPackedMsg.
 * @user_data: unused.
 *
 * Unpacks information from @pm for a near/far backward visit of a
 * remote node.
 */
void aran_development3d_visit_bw_unpack (AranDevelopment3d *devel,
                                         VsgPackedMsg *pm,
                                         gpointer user_data)
{
  aran_spherical_seriesd_unpack (devel->local, pm);
}

/**
 * aran_development3d_visit_bw_reduce:
 * @a: source #AranDevelopment3d.
 * @b: destination #AranDevelopment3d.
 * @user_data: unused.
 *
 * Backrward visit reduction operator for #AranDevelopment3d.
 */
void aran_development3d_visit_bw_reduce (AranDevelopment3d *a,
                                         AranDevelopment3d *b,
                                         gpointer user_data)
{
  aran_spherical_seriesd_add (a->local, b->local, b->local);

}

#endif /* VSG_HAVE_MPI */
