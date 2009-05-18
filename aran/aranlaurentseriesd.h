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

#ifndef __ARAN_LAURENT_SERIESD_H__
#define __ARAN_LAURENT_SERIESD_H__

#include <glib.h>

#include <vsg/vsgd.h>
#ifdef VSG_HAVE_MPI
#include <vsg/vsgpackedmsg.h>
#endif


#include <aran/arancomplex.h>

G_BEGIN_DECLS;

/* typedefs */

typedef struct _AranLaurentSeriesd AranLaurentSeriesd;

/* functions */

AranLaurentSeriesd *aran_laurent_seriesd_new (guint8 posdeg, guint8 negdeg);

void aran_laurent_seriesd_free (AranLaurentSeriesd *als);

void aran_laurent_seriesd_copy (const AranLaurentSeriesd *src,
                                AranLaurentSeriesd *dst);

AranLaurentSeriesd *aran_laurent_seriesd_clone (AranLaurentSeriesd *src);

gcomplex128 *aran_laurent_seriesd_get_term (AranLaurentSeriesd *als, gint i);

guint8 aran_laurent_seriesd_get_posdeg (AranLaurentSeriesd *als);

guint8 aran_laurent_seriesd_get_negdeg (AranLaurentSeriesd *als);

void aran_laurent_seriesd_set_zero (AranLaurentSeriesd *als);

void aran_laurent_seriesd_add (AranLaurentSeriesd *one,
                               AranLaurentSeriesd *other,
                               AranLaurentSeriesd *result);

void aran_laurent_seriesd_write (AranLaurentSeriesd *als, FILE *file);

#ifdef VSG_HAVE_MPI
void aran_laurent_seriesd_pack (AranLaurentSeriesd *als, VsgPackedMsg *pm);
void aran_laurent_seriesd_unpack (AranLaurentSeriesd *als, VsgPackedMsg *pm);
#endif

gcomplex128 aran_laurent_seriesd_evaluate (AranLaurentSeriesd *als,
					   gcomplex128 z);

void aran_laurent_seriesd_translate (AranLaurentSeriesd *src,
				     gcomplex128 zsrc,
				     AranLaurentSeriesd *dst,
				     gcomplex128 zdst);

void aran_laurent_seriesd_to_taylor (AranLaurentSeriesd *src,
				     gcomplex128 zsrc,
				     AranLaurentSeriesd *dst,
				     gcomplex128 zdst);

G_END_DECLS;

#endif /* __ARAN_LAURENT_SERIESD_H__ */
