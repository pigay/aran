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

#include "aranlaurentseriesd.h"

#include "aranbinomial.h"

#include <string.h>

/**
 * AranLaurentSeriesd:
 *
 * Opaque structure. Possesses only private data.
 */
struct _AranLaurentSeriesd
{
  guint8 posdeg;
  guint8 negdeg;
};

/* macros */
#define ARAN_LAURENT_SERIESD_TERM(als, i) ( \
(((als) == NULL) || ((i) < -(als)->negdeg) || \
 (-(i) < -(als)->posdeg)) ? NULL : \
((gcomplex128 *) (((AranLaurentSeriesd *) (als)) + 1) + \
((als)->posdeg - (i))) \
)

/* static functions */
static inline
void aran_taylor_translate (AranLaurentSeriesd *src,
			    AranLaurentSeriesd *dst,
			    gcomplex128 zd_m_zs)
{
  gint16 i, j;
  gcomplex128 *srcterm, *dstterm;
  gcomplex128 sum;

  if (src->posdeg > dst->posdeg)
    g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__);

  dstterm = ARAN_LAURENT_SERIESD_TERM (dst, dst->posdeg);

  for (i=dst->posdeg; i>=0; i --)
    {
      srcterm = ARAN_LAURENT_SERIESD_TERM (src, src->posdeg);
      sum = *srcterm;

      for (j=src->posdeg-1; j>=i; j --)
	{
	  srcterm ++;
	  sum = sum*zd_m_zs + aran_fast_binomial (j, i)*(*srcterm);
	}
      *dstterm += sum;
      dstterm ++;
    }
}

/* functions */

/**
 * aran_laurent_seriesd_new:
 * @posdeg: positive degree of the series.
 * @negdeg: negative degree of the series.
 *
 * Allocates an #AranLaurentSeriesd with the specified degrees. Polynomial
 * terms will go from z^posdeg to 1/z^(negdeg+1).
 *
 * Returns: newly allocated structure.
 */
AranLaurentSeriesd *aran_laurent_seriesd_new (guint8 posdeg, guint8 negdeg)
{
  AranLaurentSeriesd *als =
    (AranLaurentSeriesd *) g_malloc0 (sizeof (AranLaurentSeriesd) +
				      (posdeg + negdeg + 1) *
				      sizeof (gcomplex128));

  als->posdeg = posdeg;
  als->negdeg = negdeg;

  return als;
}

/**
 * aran_laurent_seriesd_free:
 * @als: an #AranLaurentSeriesd.
 *
 * Deallocates all memory associated to @als.
 */
void aran_laurent_seriesd_free (AranLaurentSeriesd *als)
{
  g_free (als);
}

/**
 * aran_laurent_seriesd_copy:
 * @src: an #AranLaurentSeriesd.
 * @dst: an #AranLaurentSeriesd.
 *
 * Copies @src into @dst. Conditions @src->posdeg > @dst->posdeg and
 * @src->negdeg > @dst->negdeg must hold.
 */
void aran_laurent_seriesd_copy (const AranLaurentSeriesd *src,
                                AranLaurentSeriesd *dst)
{
  guint8 posdeg = MIN (src->posdeg, dst->posdeg);
  guint8 size = posdeg + MIN (src->negdeg, dst->negdeg) + 1;
  void *dstp = ARAN_LAURENT_SERIESD_TERM (dst, posdeg);
  void *srcp = ARAN_LAURENT_SERIESD_TERM (src, posdeg);

  if (src->posdeg > dst->posdeg || src->negdeg > dst->negdeg)
    g_warning ("Loosing precision in \"%s\"", __PRETTY_FUNCTION__);

  aran_laurent_seriesd_set_zero (dst);

  memcpy (dstp, srcp, size * sizeof (gcomplex128));
}

/**
 * aran_laurent_seriesd_get_term:
 * @als: an #AranLaurentSeriesd.
 * @i: a #guint.
 *
 * Computes the localization of term @i in @als.
 *
 * Returns: address of term @i.
 */
gcomplex128 *aran_laurent_seriesd_get_term (AranLaurentSeriesd *als, gint i)
{
  return ARAN_LAURENT_SERIESD_TERM (als, i);
}

/**
 * aran_laurent_seriesd_get_posdeg:
 * @als: an #AranLaurentSeriesd.
 *
 * Returns: @als positive degree.
 */
guint8 aran_laurent_seriesd_get_posdeg (AranLaurentSeriesd *als)
{
  g_return_val_if_fail (als != NULL, 0);
  return als->posdeg;
}

/**
 * aran_laurent_seriesd_get_negdeg:
 * @als: an #AranLaurentSeriesd.
 *
 * Returns: @als negative degree.
 */
guint8 aran_laurent_seriesd_get_negdeg (AranLaurentSeriesd *als)
{
  g_return_val_if_fail (als != NULL, 0);
  return als->negdeg;
}

/**
 * aran_laurent_seriesd_clone:
 * @src: an #AranLaurentSeriesd.
 *
 * Duplicates @src.
 *
 * Returns: a newly allocated #AranLaurentSeriesd identical to @src.
 */
AranLaurentSeriesd *aran_laurent_seriesd_clone (AranLaurentSeriesd *src)
{
  AranLaurentSeriesd *dst;

  g_return_val_if_fail (src != NULL, NULL);

  dst = aran_laurent_seriesd_new (src->posdeg, src->negdeg);

  aran_laurent_seriesd_copy (src, dst);
/*   dst = g_memdup (src, sizeof (AranLaurentSeriesd) + */
/* 		  (src->posdeg + src->negdeg + 1) * */
/* 		  sizeof (gcomplex128)); */

  return dst;
}

/**
 * aran_laurent_seriesd_set_zero:
 * @als: an #AranLaurentSeriesd.
 *
 * Nullifies all polynomial coefficients in @als.
 */
void aran_laurent_seriesd_set_zero (AranLaurentSeriesd *als)
{
  void *ptr = ARAN_LAURENT_SERIESD_TERM (als, als->posdeg);
  guint16 size;

  g_return_if_fail (ptr != NULL);

  size = (als->posdeg+als->negdeg+1) * sizeof (gcomplex128);

  memset (ptr, 0, size);
}

/**
 * aran_laurent_seriesd_write:
 * @als: an #AranLaurentSeriesd.
 * @file: output file.
 *
 * Writes @als to @file.
 */
void aran_laurent_seriesd_write (AranLaurentSeriesd *als, FILE *file)
{
  gint16 i;
  gcomplex128 *term;

  g_return_if_fail (als != NULL);

  fprintf (file, "[");

  for (i=-als->negdeg; i<als->posdeg; i++)
    {
      term = ARAN_LAURENT_SERIESD_TERM (als, i);
      fprintf (file, "%d:(%e,%e), ", i, creal (*term), cimag (*term));
    }

  term = ARAN_LAURENT_SERIESD_TERM (als, als->posdeg);
  fprintf (file, "%d:(%e,%e)] ", als->posdeg, creal (*term), cimag (*term));
}

/**
 * aran_laurent_seriesd_evaluate:
 * @als: an #AranLaurentSeriesd.
 * @z: evaluation point.
 *
 * Evaluates @als at @z location.
 *
 * Returns: value of @als(@z).
 */
gcomplex128 aran_laurent_seriesd_evaluate (AranLaurentSeriesd *als,
					   gcomplex128 z)
{
  guint8 i;
  gcomplex128 *term;
  gcomplex128 pos, neg, invz;

  g_return_val_if_fail (als != NULL, 0.);

  term = ARAN_LAURENT_SERIESD_TERM (als, als->posdeg);
  pos = *term;

  for (i=1; i<=als->posdeg; i++)
    {
      term ++;
      pos = pos*z + *term;
    }

  neg = 0.;

  if (als->negdeg != 0)
    {
      invz = 1. / z;

      term = ARAN_LAURENT_SERIESD_TERM (als, - als->negdeg);
      neg = *term;

      for (i=als->negdeg-1; i>0; i --)
	{
	  term --;
	  neg = neg*invz + *term;
	}
      neg *= invz; /* negative terms begin with degree -1 */
    }

  return pos + neg;
}

/**
 * aran_laurent_seriesd_translate:
 * @src: an #AranLaurentSeriesd.
 * @zsrc: @src center.
 * @dst: an #AranLaurentSeriesd.
 * @zdst: @dst center.
 *
 * Computes the translation of @src to @zdst and accumulates it into @dst.
 */
void aran_laurent_seriesd_translate (AranLaurentSeriesd *src,
				     gcomplex128 zsrc,
				     AranLaurentSeriesd *dst,
				     gcomplex128 zdst)
{
  gint16 i, j;
  gcomplex128 zd_m_zs = zdst-zsrc;

  g_return_if_fail (src != NULL);
  g_return_if_fail (dst != NULL);

  aran_taylor_translate (src, dst, zd_m_zs);

  if (src->negdeg > 0)
    {
      gcomplex128 sum;
      gcomplex128 zs_m_zd = -zd_m_zs;
      gcomplex128 *srcterm, *dstterm;

      if (src->negdeg > dst->negdeg)
	g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__);

      dstterm = ARAN_LAURENT_SERIESD_TERM (dst, -1);

      for (i=1; i<=dst->negdeg; i ++)
	{
	  srcterm = ARAN_LAURENT_SERIESD_TERM (src, -1);
	  sum = 0.;

	  for (j=1; j<MIN (i, src->negdeg); j ++)
	    {
	      sum = (sum + aran_fast_binomial (i-1, j-1)*(*srcterm)) * zs_m_zd;
	      srcterm ++;
	    }
	  /* negative terms begin with degree -1 */
	  *dstterm += sum + *srcterm; 
	  dstterm ++;
	}
    }
}

/**
 * aran_laurent_seriesd_to_taylor:
 * @src: an #AranLaurentSeriesd.
 * @zsrc: @src center.
 * @dst: an #AranLaurentSeriesd.
 * @zdst: @dst center.
 *
 * Computes the transformation of @src to a Taylor series at @zdst and
 * accumulates it into @dst.
 */
void aran_laurent_seriesd_to_taylor (AranLaurentSeriesd *src,
				     gcomplex128 zsrc,
				     AranLaurentSeriesd *dst,
				     gcomplex128 zdst)
{
  gint16 i, j;
  gcomplex128 zd_m_zs = zdst-zsrc;

  g_return_if_fail (src != NULL);
  g_return_if_fail (dst != NULL);

  aran_taylor_translate (src, dst, zd_m_zs);

  if (src->negdeg > 0)
    {
      gcomplex128 sum, pow;
      gcomplex128 inv_zd_m_zs = 1. / zd_m_zs;
      gcomplex128 *srcterm, *dstterm;

      if (src->negdeg > dst->posdeg)
	g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__);

      dstterm = ARAN_LAURENT_SERIESD_TERM (dst, 0);

      pow = 1.;

      for (i=0; i<=dst->posdeg; i ++)
	{
	  srcterm = ARAN_LAURENT_SERIESD_TERM (src, -src->negdeg+1);
	  sum = 0;

	  for (j=src->negdeg-1; j>=1; j --)
	    {
	      sum = (sum + aran_fast_binomial (i+j-1, j-1) * (*srcterm)) *
		inv_zd_m_zs;
	      srcterm --;
	    }

	  *dstterm += sum * pow; 
	  dstterm --;
	  pow *= -inv_zd_m_zs;
	}
    }
}
 
