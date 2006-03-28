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

#include "arandevelopment3d.h"


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
 * @src_center: @src center.
 * @src: an #AranDevelopment3d.
 * @dst_center: @dst center.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 multipole translation between @src and @dst.
 */
void aran_development3d_m2m (const VsgVector3d *src_center,
			     AranDevelopment3d *src,
			     const VsgVector3d *dst_center,
			     AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate (src->multipole, src_center,
				    dst->multipole, dst_center);
}

/**
 * aran_development3d_m2l:
 * @src_center: @src center.
 * @src: an #AranDevelopment3d.
 * @dst_center: @dst center.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 local translation between @src and @dst.
 */
void aran_development3d_m2l (const VsgVector3d *src_center,
			     AranDevelopment3d *src,
			     const VsgVector3d *dst_center,
			     AranDevelopment3d *dst)
{
  aran_spherical_seriesd_to_local (src->multipole, src_center,
				   dst->local, dst_center);
}

/**
 * aran_development3d_l2l:
 * @src_center: @src center.
 * @src: an #AranDevelopment3d.
 * @dst_center: @dst center.
 * @dst: an #AranDevelopment3d.
 *
 * Performs local 2 local translation between @src and @dst.
 */
void aran_development3d_l2l (const VsgVector3d *src_center,
			     AranDevelopment3d *src,
			     const VsgVector3d *dst_center,
			     AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate (src->local, src_center,
				    dst->local, dst_center);
}

/**
 * aran_development3d_m2m_kkylin:
 * @src_center: @src center.
 * @src: an #AranDevelopment3d.
 * @dst_center: @dst center.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 multipole translation between @src and @dst with the
 * K. Kylin formulas.
 */
void aran_development3d_m2m_kkylin (const VsgVector3d *src_center,
                                    AranDevelopment3d *src,
                                    const VsgVector3d *dst_center,
                                    AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate_kkylin (src->multipole, src_center,
                                           dst->multipole, dst_center);
}

/**
 * aran_development3d_m2l_kkylin:
 * @src_center: @src center.
 * @src: an #AranDevelopment3d.
 * @dst_center: @dst center.
 * @dst: an #AranDevelopment3d.
 *
 * Performs multipole 2 local translation between @src and @dst with the
 * K. Kylin formulas.
 */
void aran_development3d_m2l_kkylin (const VsgVector3d *src_center,
                                    AranDevelopment3d *src,
                                    const VsgVector3d *dst_center,
                                    AranDevelopment3d *dst)
{
  aran_spherical_seriesd_to_local_kkylin (src->multipole, src_center,
                                          dst->local, dst_center);
}

/**
 * aran_development3d_l2l_kkylin:
 * @src_center: @src center.
 * @src: an #AranDevelopment3d.
 * @dst_center: @dst center.
 * @dst: an #AranDevelopment3d.
 *
 * Performs local 2 local translation between @src and @dst with the
 * K. Kylin formulas.
 */
void aran_development3d_l2l_kkylin (const VsgVector3d *src_center,
                                    AranDevelopment3d *src,
                                    const VsgVector3d *dst_center,
                                    AranDevelopment3d *dst)
{
  aran_spherical_seriesd_translate_kkylin (src->local, src_center,
                                           dst->local, dst_center);
}

/**
 * aran_development3d_multipole_evaluate:
 * @devel_center: center of @devel.
 * @devel: an #AranDevelopment3d.
 * @pos: evaluation position.
 *
 * Evaluates the multipole part of @devel at @pos.
 *
 * Returns: value of multipole part of @devel(@pos).
 */
gcomplex128
aran_development3d_multipole_evaluate (const VsgVector3d *devel_center,
                                       AranDevelopment3d *devel,
                                       const VsgVector3d *pos)
{
  VsgVector3d tmp;

  vsg_vector3d_sub (pos, devel_center, &tmp);

  return aran_spherical_seriesd_evaluate (devel->multipole, &tmp);
}

/**
 * aran_development3d_local_evaluate:
 * @devel_center: center of @devel.
 * @devel: an #AranDevelopment3d.
 * @pos: evaluation position.
 *
 * Evaluates the local part of @devel at @pos.
 *
 * Returns: value of local part of @devel(@pos).
 */
gcomplex128 aran_development3d_local_evaluate (const VsgVector3d *devel_center,
					       AranDevelopment3d *devel,
					       const VsgVector3d *pos)
{
  VsgVector3d tmp;

  vsg_vector3d_sub (pos, devel_center, &tmp);

  return aran_spherical_seriesd_evaluate (devel->local, &tmp);
}

/**
 * aran_development3d_local_gradient_evaluate:
 * @devel_center: center of @devel.
 * @devel: an #AranDevelopment3d.
 * @pos: evaluation position.
 * @grad: result gradient.
 *
 * Evaluates the gradient of the local part of @devel at @pos.
 */
void
aran_development3d_local_gradient_evaluate (const VsgVector3d *devel_center,
                                            AranDevelopment3d *devel,
                                            const VsgVector3d *pos,
                                            VsgVector3d *grad)
{
  VsgVector3d tmp;

  vsg_vector3d_sub (pos, devel_center, &tmp);

  aran_spherical_seriesd_local_gradient_evaluate (devel->local, &tmp, grad);
}

