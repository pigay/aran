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

#ifndef __ARAN_DEVELOPMENT2D_H__
#define __ARAN_DEVELOPMENT2D_H__

#include <glib-object.h>

#include <vsg/vsgd.h>
#include <aran/aranlaurentseriesd.h>

G_BEGIN_DECLS;

/* macros */
#define ARAN_TYPE_DEVELOPMENT2D (aran_development2d_get_type ())

/* typedefs */

typedef struct _AranDevelopment2d AranDevelopment2d;

struct _AranDevelopment2d
{
  AranLaurentSeriesd *multipole;
  AranLaurentSeriesd *local;
};

/* functions */
GType aran_development2d_get_type ();

AranDevelopment2d *aran_development2d_new (guint8 posdeg, guint8 negdeg);

void aran_development2d_free (AranDevelopment2d *ad);

void aran_development2d_copy (const AranDevelopment2d *src,
                              AranDevelopment2d *dst);

AranDevelopment2d *aran_development2d_clone (AranDevelopment2d *src);

void aran_development2d_set_zero (AranDevelopment2d *ad);

void aran_development2d_write (AranDevelopment2d *ad, FILE *file);

void aran_development2d_m2m (const VsgVector2d *src_center,
			     AranDevelopment2d *src,
			     const VsgVector2d *dst_center,
			     AranDevelopment2d *dst);

void aran_development2d_m2l (const VsgVector2d *src_center,
			     AranDevelopment2d *src,
			     const VsgVector2d *dst_center,
			     AranDevelopment2d *dst);

void aran_development2d_l2l (const VsgVector2d *src_center,
			     AranDevelopment2d *src,
			     const VsgVector2d *dst_center,
			     AranDevelopment2d *dst);

gcomplex128
aran_development2d_multipole_evaluate (const VsgVector2d *devel_center,
                                       AranDevelopment2d *devel,
                                       const VsgVector2d *pos);

gcomplex128 aran_development2d_local_evaluate (const VsgVector2d *devel_center,
					       AranDevelopment2d *devel,
					       const VsgVector2d *pos);

G_END_DECLS;

#endif /* __ARAN_DEVELOPMENT2D_H__ */
