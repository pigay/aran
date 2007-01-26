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
 * @src_center: @src center.
 * @src: an #AranDevelopment2d.
 * @dst_center: @dst center.
 * @dst: an #AranDevelopment2d.
 *
 * Performs multipole 2 multipole translation between @src and @dst.
 */
void aran_development2d_m2m (const VsgVector2d *src_center,
			     AranDevelopment2d *src,
			     const VsgVector2d *dst_center,
			     AranDevelopment2d *dst)
{
  gcomplex128 zsrc = src_center->x + G_I*src_center->y;
  gcomplex128 zdst = dst_center->x + G_I*dst_center->y;

  aran_laurent_seriesd_translate (src->multipole, zsrc, dst->multipole, zdst);
}

/**
 * aran_development2d_m2l:
 * @src_center: @src center.
 * @src: an #AranDevelopment2d.
 * @dst_center: @dst center.
 * @dst: an #AranDevelopment2d.
 *
 * Performs multipole 2 local translation between @src and @dst.
 */
void aran_development2d_m2l (const VsgVector2d *src_center,
			     AranDevelopment2d *src,
			     const VsgVector2d *dst_center,
			     AranDevelopment2d *dst)
{
  gcomplex128 zsrc = src_center->x + G_I*src_center->y;
  gcomplex128 zdst = dst_center->x + G_I*dst_center->y;

  aran_laurent_seriesd_to_taylor (src->multipole, zsrc, dst->local, zdst);
}

/**
 * aran_development2d_l2l:
 * @src_center: @src center.
 * @src: an #AranDevelopment2d.
 * @dst_center: @dst center.
 * @dst: an #AranDevelopment2d.
 *
 * Performs local 2 local translation between @src and @dst.
 */
void aran_development2d_l2l (const VsgVector2d *src_center,
			     AranDevelopment2d *src,
			     const VsgVector2d *dst_center,
			     AranDevelopment2d *dst)
{
  gcomplex128 zsrc = src_center->x + G_I*src_center->y;
  gcomplex128 zdst = dst_center->x + G_I*dst_center->y;

  aran_laurent_seriesd_translate (src->local, zsrc, dst->local, zdst);
}

/**
 * aran_development2d_multipole_evaluate:
 * @devel_center: center of @devel.
 * @devel: an #AranDevelopment2d.
 * @pos: evaluation position.
 *
 * Evaluates the multipole part of @devel at @pos.
 *
 * Returns: value of multipole part of @devel(@pos).
 */
gcomplex128
aran_development2d_multipole_evaluate (const VsgVector2d *devel_center,
                                       AranDevelopment2d *devel,
                                       const VsgVector2d *pos)
{
  gcomplex128 zl = devel_center->x + G_I*devel_center->y;
  gcomplex128 zp = pos->x + G_I*pos->y;

  return aran_laurent_seriesd_evaluate (devel->multipole, zp-zl);
}

/**
 * aran_development2d_local_evaluate:
 * @devel_center: center of @devel.
 * @devel: an #AranDevelopment2d.
 * @pos: evaluation position.
 *
 * Evaluates the local part of @devel at @pos.
 *
 * Returns: value of local part of @devel(@pos).
 */
gcomplex128 aran_development2d_local_evaluate (const VsgVector2d *devel_center,
					       AranDevelopment2d *devel,
					       const VsgVector2d *pos)
{
  gcomplex128 zl = devel_center->x + G_I*devel_center->y;
  gcomplex128 zp = pos->x + G_I*pos->y;

  return aran_laurent_seriesd_evaluate (devel->local, zp-zl);
}
