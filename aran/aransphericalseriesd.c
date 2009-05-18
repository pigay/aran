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

#include "aransphericalseriesd.h"
#include "aransphericalseriesd-private.h"

#include "aransphericalharmonic.h"
#include "arancoefficientbufferd.h"
#include "aranbinomialbufferd.h"

/* #include <vsg/vsgd-inline.h> */

#include <string.h>
#include <math.h>

/**
 * AranSphericalSeriesd:
 *
 * Opaque structure. Possesses only private data.
 */
static AranBinomialBufferd *_alpha_buffer = NULL;
static AranCoefficientBufferd *_beta_buffer = NULL;

static gdouble _beta_generator (guint l, AranCoefficientBufferd * buf)
{
  return sqrt ((4. * G_PI) / (l + l + 1.));
}

static gdouble _alpha_generator (guint l, guint m, AranBinomialBufferd * buf)
{
  if (l == 0)
    return 1.;

  if (m == 0)
    return *aran_binomial_bufferd_get (buf, l - 1, m) / l;

  if (l == m)
    return *aran_binomial_bufferd_get (buf, l - 1, l - 1) /
      sqrt ((l + l - 1.) * (l + l));

  return sqrt ((l - m + 1.) / (l + m)) * *aran_binomial_bufferd_get (buf, l,
                                                                     m - 1);

}

gdouble aran_spherical_seriesd_beta (gint l)
{
/*   return sqrt ((4.*G_PI)/(l+l+1.)); */
  return *aran_coefficient_bufferd_get_unsafe (_beta_buffer, l);
}

gdouble aran_spherical_seriesd_alpha (guint n, guint p)
{
  return *aran_binomial_bufferd_get_unsafe (_alpha_buffer, n, p);
}

static void _atexit ()
{
  aran_coefficient_bufferd_free (_beta_buffer);
  _beta_buffer = NULL;

  aran_binomial_bufferd_free (_alpha_buffer);
  _alpha_buffer = NULL;
}

void aran_spherical_seriesd_beta_require (guint deg)
{
  aran_coefficient_bufferd_require (_beta_buffer, deg);
}

void aran_spherical_seriesd_alpha_require (guint deg)
{
  aran_binomial_bufferd_require (_alpha_buffer, deg);
}

/* functions */

/**
 * aran_spherical_seriesd_new:
 * @posdeg: the positive expansion degree.
 * @negdeg: the negative expansion degree.
 *
 * Creates a new #AranSphericalSeriesd with required positive (local) and
 * negative (multipole) degrees.
 *
 * Returns: Newly allocated #AranSphericalSeriesd structure.
 */
AranSphericalSeriesd *aran_spherical_seriesd_new (guint8 posdeg, guint8 negdeg)
{
  AranSphericalSeriesd *ass;

  if (_beta_buffer == NULL)
    {
      _beta_buffer = aran_coefficient_bufferd_new (_beta_generator,
                                                   posdeg + negdeg);
      _alpha_buffer = aran_binomial_bufferd_new (_alpha_generator,
                                                 posdeg + negdeg);
      g_atexit (_atexit);
    }
  else
    {
      aran_spherical_seriesd_beta_require (posdeg + negdeg);
      aran_spherical_seriesd_alpha_require (posdeg + negdeg);
    }

  ass = (AranSphericalSeriesd *)
    g_malloc0 (ARAN_SPHERICAL_SERIESD_SIZE (posdeg, negdeg));

  ass->posdeg = posdeg;
  ass->negdeg = negdeg;

  return ass;
}

/**
 * aran_spherical_seriesd_free:
 * @ass: an #AranSphericalSeriesd.
 *
 * Frees memory allocated for @ass.
 */
void aran_spherical_seriesd_free (AranSphericalSeriesd * ass)
{
  g_free (ass);
}

/**
 * aran_spherical_seriesd_copy:
 * @src: an #AranSphericalSeriesd.
 * @dst: an #AranSphericalSeriesd.
 *
 * Copies @src to @dst up to @dst degrees: can loose precision when
 * @src->posdeg > @dst->posdeg for example.
 */
void aran_spherical_seriesd_copy (const AranSphericalSeriesd * src,
                                  AranSphericalSeriesd * dst)
{
  guint8 nd = MIN (dst->negdeg, src->negdeg);
  guint8 pd = MIN (dst->posdeg, src->posdeg);

  if (src->posdeg > dst->posdeg || src->negdeg > dst->negdeg)
    g_critical ("Could loose precision in \"%s\"", __PRETTY_FUNCTION__);

  aran_spherical_seriesd_set_zero (dst);

  memcpy (_spherical_seriesd_get_pos_term (dst, 0, 0),
          _spherical_seriesd_get_pos_term (src, 0, 0),
          _spherical_seriesd_size (pd, 0) * sizeof (gcomplex128));

  if (nd != 0)
    memcpy (_spherical_seriesd_get_neg_term (dst, 0, 0),
            _spherical_seriesd_get_neg_term (src, 0, 0),
            (_spherical_seriesd_size (0, nd) -
             _spherical_seriesd_size (0, 0)) * sizeof (gcomplex128));
}

/**
 * aran_spherical_seriesd_clone:
 * @src: an #AranSphericalSeriesd.
 *
 * Duplicates @src.
 *
 * Returns: a newly allocated #AranSphericalSeriesd identical to @src.
 */
AranSphericalSeriesd *
aran_spherical_seriesd_clone (const AranSphericalSeriesd * src)
{
  g_return_val_if_fail (src != NULL, NULL);

  return g_memdup (src,
                   ARAN_SPHERICAL_SERIESD_SIZE (src->posdeg, src->negdeg));
}

/**
 * aran_spherical_seriesd_get_term:
 * @ass: an #AranSphericalSeriesd.
 * @i: a #gint. Condition -@ass->negdeg <= @i <= @ass->posdeg must hold.
 * @j: a #gint. Condition -@i <= @j <= @i must hold.
 *
 * Provides access to a specified term in @ass.
 * Returns: address of the coefficient.
 */
gcomplex128 *aran_spherical_seriesd_get_term (AranSphericalSeriesd * ass,
                                              gint i, gint j)
{
  return _spherical_seriesd_get_term (ass, i, j);
}

/**
 * aran_spherical_seriesd_get_posdeg:
 * @ass: an #AranSphericalSeriesd.
 *
 * Returns: positive degree of @ass.
 */
guint8 aran_spherical_seriesd_get_posdeg (const AranSphericalSeriesd * ass)
{
  g_return_val_if_fail (ass != NULL, 0);
  return ass->posdeg;
}

/**
 * aran_spherical_seriesd_get_negdeg:
 * @ass: an #AranSphericalSeriesd.
 *
 * Returns: negative degree of @ass.
 */
guint8 aran_spherical_seriesd_get_negdeg (const AranSphericalSeriesd * ass)
{
  g_return_val_if_fail (ass != NULL, 0);
  return ass->negdeg;
}

/**
 * aran_spherical_seriesd_set_zero:
 * @ass: an #AranSphericalSeriesd.
 *
 * Sets all @ass coefficients to zero. @ass degrees are unchanged.
 */
void aran_spherical_seriesd_set_zero (AranSphericalSeriesd * ass)
{
  void *ptr = _spherical_seriesd_get_pos_term (ass, 0, 0);
  guint16 size;

  g_return_if_fail (ptr != NULL);

  size = _spherical_seriesd_size (ass->posdeg, ass->negdeg) *
    sizeof (gcomplex128);

  memset (ptr, 0, size);
}

/**
 * aran_spherical_seriesd_write:
 * @ass: an #AranSphericalSeriesd.
 * @file: output file.
 *
 * Writes @ass to @file.
 */
void aran_spherical_seriesd_write (const AranSphericalSeriesd * ass,
                                   FILE * file)
{
  gint16 i, j;
  gcomplex128 *term;

  g_return_if_fail (ass != NULL);

  fprintf (file, "[");

  for (i = -ass->negdeg; i <= ass->posdeg; i++)
    {
      gint16 maxj = ABS (i);

      if (i < 0)
        maxj--;

      for (j = 0; j <= maxj; j++)
        {
          term = _spherical_seriesd_get_term (ass, i, j);
          fprintf (file, "[%d,%d]:(%e,%e), ",
                   i, j, creal (*term), cimag (*term));
        }
    }

  fprintf (file, "]");
}

/**
 * aran_spherical_seriesd_evaluate_internal:
 * @ass: an #AranSphericalSeriesd.
 * @r: radius.
 * @cost: cos (theta).
 * @sint: sin (theta).
 * @cosp: cos (phi).
 * @sinp: sin (phi).
 *
 * Evaluates @ass at a point defined by its spherical coordinates.
 * Returns: value of @ass at the specified location.
 */
gcomplex128
aran_spherical_seriesd_evaluate_internal (const AranSphericalSeriesd * ass,
                                          gdouble r,
                                          gdouble cost, gdouble sint,
                                          gdouble cosp, gdouble sinp)
{
  gint l, m;
  gint n = MAX (ass->posdeg, ((gint) ass->negdeg) - 1);
  gcomplex128 harmonics[((n + 1) * (n + 2)) / 2];
  gcomplex128 *coefficient, *hterm;
  gcomplex128 expp = cosp + G_I * sinp;
  gcomplex128 res = 0.;

  aran_spherical_harmonic_evaluate_multiple_internal (n, cost, sint, expp,
                                                      harmonics);

  for (l = ass->posdeg; l >= 0; l--)
    {
      gcomplex128 sum = 0.;

      coefficient = _spherical_seriesd_get_pos_term (ass, l, 0);
      hterm = aran_spherical_harmonic_multiple_get_term (l, 0, harmonics);

      for (m = 1; m <= l; m++)
        {
          sum += creal (coefficient[m] * hterm[m]);
        }

      sum = 2. * sum + coefficient[0] * hterm[0];

      res = res * r + sum;
    }

  if (ass->negdeg != 0)
    {
      gcomplex128 negres = 0.;
      gdouble invr = 1. / r;

      for (l = ass->negdeg - 1; l >= 0; l--)
        {
          gcomplex128 sum = 0.;

          coefficient = _spherical_seriesd_get_neg_term (ass, l, 0);
          hterm = aran_spherical_harmonic_multiple_get_term (l, 0, harmonics);

          for (m = 1; m <= l; m++)
            {
              sum += creal (coefficient[m] * hterm[m]);
            }

          sum = 2. * sum + coefficient[0] * hterm[0];

          negres = (negres + sum) * invr;
        }

      res += negres;
    }

  return res;
}

/**
 * aran_spherical_seriesd_evaluate:
 * @ass: an #AranSphericalSeriesd.
 * @x: evaluation location.
 *
 * Evaluates @ass at @x. This is a convenience wrapper around
 * aran_spherical_seriesd_evaluate_internal().
 *
 * Returns: value of @ass at @x.
 */
gcomplex128 aran_spherical_seriesd_evaluate (const AranSphericalSeriesd * ass,
                                             const VsgVector3d * x)
{
  gdouble r, cost, sint, cosp, sinp;

  vsg_vector3d_to_spherical_internal (x, &r, &cost, &sint, &cosp, &sinp);

  return aran_spherical_seriesd_evaluate_internal (ass, r, cost, sint, cosp,
                                                   sinp);
}

/**
 * aran_spherical_seriesd_local_gradient_evaluate_internal:
 * @ass: an #AranSphericalSeriesd.
 * @r: radius
 * @cost: cos (theta).
 * @sint: sin (theta).
 * @cosp: cos (phi).
 * @sinp: sin (phi).
 * @dr: gradient radius coordinate.
 * @dt: gradient theta coordinate.
 * @dp: gradient phi coordinate.
 *
 * Evaluates gradient of local part of @ass at a point defined by its spherical
 * coordinates.
 */
void aran_spherical_seriesd_local_gradient_evaluate_internal
  (const AranSphericalSeriesd * ass,
   gdouble r,
   gdouble cost, gdouble sint,
   gdouble cosp, gdouble sinp,
   gcomplex128 * dr, gcomplex128 * dt, gcomplex128 * dp)
{
  gint l, m;
  gcomplex128 expp = cosp + G_I * sinp;
  gcomplex128 harmonics[((ass->posdeg + 1) * (ass->posdeg + 2)) / 2];
  gcomplex128 special_harmonics[((ass->posdeg + 1) * (ass->posdeg + 2)) / 2];
  gcomplex128 *hterm, *srcterm, *shterm;

  aran_spherical_harmonic_pre_gradient_multiple_internal (ass->posdeg,
                                                          cost, sint, expp,
                                                          harmonics,
                                                          special_harmonics);

  *dr = 0.;
  *dt = 0.;
  *dp = 0.;

  for (l = ass->posdeg; l >= 1; l--)
    {
      gcomplex128 sumr = 0., sumt = 0., sump = 0.;

      hterm = aran_spherical_harmonic_multiple_get_term (l, 0, harmonics);
      shterm =
        aran_spherical_harmonic_multiple_get_term (l, 0, special_harmonics);
      srcterm = _spherical_seriesd_get_pos_term (ass, l, 0);

      for (m = 1; m <= l; m++)
        {
          sumr += creal (hterm[m] * srcterm[m]);
          sumt += creal (shterm[m] * srcterm[m]) * (cost * m);
          sump += cimag (shterm[m] * srcterm[m]) * m;
        }
      for (m = 1; m <= l - 1; m++)
        {
          sumt += creal (hterm[m + 1] * srcterm[m] * conj (expp)) *
            sqrt ((l - m) * (l + m + 1));
        }
      sumr = 2. * sumr + creal (hterm[0] * srcterm[0]);
      sumt = 2. * sumt + creal (hterm[1] * srcterm[0] * conj (expp)) *
        sqrt (l * (l + 1));
      sump *= -2.;

      *dr = (*dr) * r + sumr * l;
      *dt = (*dt) * r + sumt;
      *dp = (*dp) * r + sump;
    }
}

static void local_to_cartesian (gdouble r,
                                gdouble cost, gdouble sint,
                                gdouble cosp, gdouble sinp,
                                gdouble dr, gdouble dt, gdouble dp,
                                VsgVector3d * grad)
{
  grad->x = sint * cosp * dr + cost * cosp * dt - sinp * dp;
  grad->y = sint * sinp * dr + cost * sinp * dt + cosp * dp;
  grad->z = cost * dr - sint * dt;
}

/**
 * aran_spherical_seriesd_local_gradient_evaluate:
 * @ass: an #AranSphericalSeriesd.
 * @x: location.
 * @grad: result
 *
 * Evaluates the gradient local part of @ass at @x. This is a convenience
 * wrapper around aran_spherical_seriesd_local_gradient_evaluate_internal().
 */
void aran_spherical_seriesd_local_gradient_evaluate
  (const AranSphericalSeriesd * ass, const VsgVector3d * x, VsgVector3d * grad)
{
  gdouble r, cost, sint, cosp, sinp;
  gcomplex128 dr, dt, dp;

  vsg_vector3d_to_spherical_internal (x, &r, &cost, &sint, &cosp, &sinp);

  aran_spherical_seriesd_local_gradient_evaluate_internal (ass,
                                                           r,
                                                           cost, sint,
                                                           cosp, sinp,
                                                           &dr, &dt, &dp);

  local_to_cartesian (r, cost, sint, cosp, sinp,
                      creal (dr), creal (dt), creal (dp), grad);
}
