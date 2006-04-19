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

#include "aransphericalseriesd.h"
#include "aransphericalseriesd-private.h"

#include "aransphericalharmonic.h"

#include <string.h>
#include <math.h>

/* functions */


static void aran_local_translate_vertical (const AranSphericalSeriesd * src,
                                           AranSphericalSeriesd * dst,
                                           gdouble r,
                                           gdouble cost,
                                           gdouble cosp, gdouble sinp)
{
  gint l, m;
  gint n;
  gdouble rpow[src->posdeg + 1];
  gdouble pow;
  gcomplex128 *srcterm, *dstterm, *hterm;
  gint d = MAX (src->posdeg, dst->posdeg);
  gcomplex128 harmonics[((d + 1) * (d + 2)) / 2];
  gcomplex128 expp = cosp + G_I * sinp;

  if (src->posdeg > dst->posdeg)
    g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__);

  aran_spherical_seriesd_beta_require (d);
  aran_spherical_seriesd_alpha_require (d);

  aran_spherical_harmonic_evaluate_multiple_internal (src->posdeg, cost, 0.,
                                                      expp, harmonics);

  pow = 1.;
  for (l = 0; l <= src->posdeg; l++)
    {
      rpow[l] = pow;
      pow *= r;
    }

  for (l = 0; l <= dst->posdeg; l++)
    {
      for (m = 0; m <= l; m++)
        {
          dstterm = _spherical_seriesd_get_pos_term (dst, l, m);

          for (n = l; n <= src->posdeg; n++)
            {
              gdouble normaliz = aran_spherical_seriesd_beta (n - l) *
                aran_spherical_seriesd_beta (l) /
                aran_spherical_seriesd_beta (n);
              gdouble factor;

              srcterm = _spherical_seriesd_get_pos_term (src, n, 0);
              hterm = aran_spherical_harmonic_multiple_get_term (n - l, 0,
                                                                 harmonics);

              /* o-m = 0 */
              factor =
                aran_spherical_seriesd_alpha (n - l, 0) *
                aran_spherical_seriesd_alpha (l, m) /
                aran_spherical_seriesd_alpha (n, m);

              *dstterm += hterm[0] * srcterm[m] * factor * normaliz *
                rpow[n - l];
            }
        }
    }
}

static void aran_local_translate (const AranSphericalSeriesd * src,
                                  AranSphericalSeriesd * dst,
                                  gdouble r,
                                  gdouble cost, gdouble sint,
                                  gdouble cosp, gdouble sinp)
{
  gint l, m;
  gint n, o;
  gdouble rpow[src->posdeg + 1];
  gdouble pow;
  gcomplex128 *srcterm, *dstterm, *hterm;
  gint d = MAX (src->posdeg, dst->posdeg);
  gcomplex128 harmonics[((d + 1) * (d + 2)) / 2];
  gcomplex128 expp = cosp + G_I * sinp;

  if (src->posdeg > dst->posdeg)
    g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__);

  aran_spherical_seriesd_beta_require (d);
  aran_spherical_seriesd_alpha_require (d);

  aran_spherical_harmonic_evaluate_multiple_internal (src->posdeg, cost, sint,
                                                      expp, harmonics);

  pow = 1.;
  for (l = 0; l <= src->posdeg; l++)
    {
      rpow[l] = pow;
      pow *= r;
    }

  for (l = 0; l <= dst->posdeg; l++)
    {
      for (m = 0; m <= l; m++)
        {
          dstterm = _spherical_seriesd_get_pos_term (dst, l, m);

          for (n = l; n <= src->posdeg; n++)
            {
              gdouble normaliz = aran_spherical_seriesd_beta (n - l) *
                aran_spherical_seriesd_beta (l) /
                aran_spherical_seriesd_beta (n);
              gcomplex128 sum = 0.;

              srcterm = _spherical_seriesd_get_pos_term (src, n, 0);
              hterm = aran_spherical_harmonic_multiple_get_term (n - l, 0,
                                                                 harmonics);
              for (o = l + m - n; o <= m + n - l; o++)
                {
                  gint o_m_m = o - m;
                  guint abs_o_m_m = ABS (o_m_m);

                  gdouble factor =
                    aran_spherical_seriesd_alpha (n - l, abs_o_m_m) *
                    aran_spherical_seriesd_alpha (l, ABS (m)) /
                    aran_spherical_seriesd_alpha (n, ABS (o));

                  gcomplex128 h = hterm[abs_o_m_m];

                  /* h = Y_(n-l)^(o-m) */
                  if (o_m_m < 0)
                    h = _sph_sym (h, abs_o_m_m);

                  if (o >= 0)
                    h *= srcterm[o];
                  else
                    h *= _sph_sym (srcterm[-o], -o);

                  sum += h * factor;
                }

              *dstterm += sum * (normaliz * rpow[n - l]);
            }
        }
    }
}

static void
aran_multipole_translate_vertical (const AranSphericalSeriesd * src,
                                   AranSphericalSeriesd * dst,
                                   gdouble r,
                                   gdouble cost,
                                   gdouble cosp, gdouble sinp)
{
  gint l, m;
  gint n;
  gdouble rpow[dst->negdeg];
  gdouble pow;
  gcomplex128 *srcterm, *dstterm, *hterm;
  gint d = MAX (src->negdeg, dst->negdeg) - 1;
  gcomplex128 harmonics[((d + 1) * (d + 2)) / 2];
  gcomplex128 expp = cosp + G_I * sinp;

  if (src->negdeg > dst->negdeg)
    g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__);

  aran_spherical_seriesd_alpha_require (d);
  aran_spherical_seriesd_beta_require (d);

  aran_spherical_harmonic_evaluate_multiple_internal (dst->negdeg - 1, cost,
                                                      0., expp, harmonics);

  pow = 1.;
  for (l = 0; l < dst->negdeg; l++)
    {
      rpow[l] = pow;
      pow *= r;
      for (m = 0; m <= l; m++)
        {
          dstterm = _spherical_seriesd_get_neg_term (dst, l, m);

          for (n = m; n <= MIN (l, src->negdeg - 1); n++)
            {
              gdouble normaliz = aran_spherical_seriesd_beta (l - n) *
                aran_spherical_seriesd_beta (l) /
                aran_spherical_seriesd_beta (n);
              gdouble factor;

              srcterm = _spherical_seriesd_get_neg_term (src, n, 0);
              hterm = aran_spherical_harmonic_multiple_get_term (l - n, 0,
                                                                 harmonics);

              /* m-o = 0 */
              factor =
                aran_spherical_seriesd_alpha (l - n, 0) *
                aran_spherical_seriesd_alpha (n, m) /
                aran_spherical_seriesd_alpha (l, m);

              *dstterm += conj (hterm[0]) * srcterm[m] * factor * normaliz *
                rpow[l - n];
            }
        }
    }
}

static void aran_multipole_translate (const AranSphericalSeriesd * src,
                                      AranSphericalSeriesd * dst,
                                      gdouble r,
                                      gdouble cost, gdouble sint,
                                      gdouble cosp, gdouble sinp)
{
  gint l, m;
  gint n, o;
  gdouble rpow[dst->negdeg];
  gdouble pow;
  gcomplex128 *srcterm, *dstterm, *hterm;
  gint d = MAX (src->negdeg, dst->negdeg) - 1;
  gcomplex128 harmonics[((d + 1) * (d + 2)) / 2];
  gcomplex128 expp = cosp + G_I * sinp;

  if (src->negdeg > dst->negdeg)
    g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__);

  aran_spherical_seriesd_alpha_require (d);
  aran_spherical_seriesd_beta_require (d);

  aran_spherical_harmonic_evaluate_multiple_internal (dst->negdeg - 1, cost,
                                                      sint, expp, harmonics);

  pow = 1.;
  for (l = 0; l < dst->negdeg; l++)
    {
      rpow[l] = pow;
      pow *= r;
      for (m = 0; m <= l; m++)
        {
          dstterm = _spherical_seriesd_get_neg_term (dst, l, m);

          for (n = 0; n <= MIN (l, src->negdeg - 1); n++)
            {
              gdouble normaliz = aran_spherical_seriesd_beta (l - n) *
                aran_spherical_seriesd_beta (l) /
                aran_spherical_seriesd_beta (n);
              gcomplex128 sum = 0.;

              srcterm = _spherical_seriesd_get_neg_term (src, n, 0);
              hterm = aran_spherical_harmonic_multiple_get_term (l - n, 0,
                                                                 harmonics);
              for (o = MAX (-n, m + n - l); o <= MIN (n, l + m - n); o++)
                {
                  guint abs_m_m_o = ABS (m - o);
                  gdouble factor =
                    aran_spherical_seriesd_alpha (l - n, abs_m_m_o) *
                    aran_spherical_seriesd_alpha (n, ABS (o)) /
                    aran_spherical_seriesd_alpha (l, ABS (m));

                  gcomplex128 h = conj (hterm[abs_m_m_o]);

                  /* h=conj ( Y_(l-n)^(m-o) ) */
                  if ((m - o) < 0)
                    h = _sph_sym (h, abs_m_m_o);

                  if (o >= 0)
                    h *= srcterm[o];
                  else
                    h *= _sph_sym (srcterm[-o], -o);

                  sum += h * factor;
                }

              *dstterm += sum * normaliz * rpow[l - n];
            }
        }
    }
}

void
aran_spherical_seriesd_translate_vertical (const AranSphericalSeriesd * src,
                                           AranSphericalSeriesd * dst,
                                           gdouble r,
                                           gdouble cost,
                                           gdouble cosp, gdouble sinp)
{
  aran_local_translate_vertical (src, dst, r, cost, cosp, sinp);

  if (src->negdeg > 0)
    {
      aran_multipole_translate_vertical (src, dst, r, -cost, -cosp, -sinp);
    }
}

void
aran_spherical_seriesd_multipole_to_local_vertical
(const AranSphericalSeriesd * src,
 AranSphericalSeriesd * dst,
 gdouble r,
 gdouble cost,
 gdouble cosp, gdouble sinp)
{
  gint l, m;
  gint n;
  gint d = dst->posdeg + src->negdeg;
  gdouble rpow[d + 1];
  gdouble pow, inv_r;
  gcomplex128 *srcterm, *dstterm, *hterm;
  gcomplex128 harmonics[((d + 1) * (d + 2)) / 2];
  gcomplex128 expp = cosp + G_I * sinp;
  gdouble sign;

  aran_spherical_seriesd_alpha_require (d);
  aran_spherical_seriesd_beta_require (d);

  aran_spherical_harmonic_evaluate_multiple_internal (d, cost, 0., expp,
                                                      harmonics);

  inv_r = 1. / r;
  pow = 1.;
  for (l = 0; l <= d; l++)
    {
      rpow[l] = pow;
      pow *= inv_r;
    }

  sign = 1.;
  for (l = 0; l <= dst->posdeg; l++)
    {
      for (m = 0; m <= l; m++)
        {
          dstterm = _spherical_seriesd_get_pos_term (dst, l, m);

          for (n = m; n < src->negdeg; n++)
            {
              gdouble normaliz = aran_spherical_seriesd_beta (l + n) *
                aran_spherical_seriesd_beta (l) /
                aran_spherical_seriesd_beta (n);
              gcomplex128 sum;
              gdouble factor;
              gcomplex128 h;

              srcterm = _spherical_seriesd_get_neg_term (src, n, 0);
              hterm = aran_spherical_harmonic_multiple_get_term (l + n, 0,
                                                                 harmonics);

              /* o == -m */
              factor = aran_spherical_seriesd_alpha (l, m) *
                aran_spherical_seriesd_alpha (n, m);

              /* h= Y_(l+n)^(m+o) */
              h = hterm[0];

              sum = h * _sph_sym (srcterm[m], m) * factor /
                aran_spherical_seriesd_alpha (l + n, 0);

              *dstterm += conj (sum) * sign * normaliz * rpow[l + n + 1];
            }
        }
      sign = -sign;
    }
}

static void aran_multipole_to_local (const AranSphericalSeriesd * src,
                                     AranSphericalSeriesd * dst,
                                     gdouble r,
                                     gdouble cost, gdouble sint,
                                     gdouble cosp, gdouble sinp)
{
  gint l, m;
  gint n, o;
  gint d = dst->posdeg + src->negdeg;
  gdouble rpow[d + 1];
  gdouble pow, inv_r;
  gcomplex128 *srcterm, *dstterm, *hterm;
  gcomplex128 harmonics[((d + 1) * (d + 2)) / 2];
  gcomplex128 expp = cosp + G_I * sinp;
  gdouble sign;

/*   if ((src->negdeg-1) > dst->posdeg) */
/*     g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__); */

  /* vertical translation */
  if (ABS (sint) < 1.e-5) /* sint in m2l translation is zero or high values */
    {
      aran_spherical_seriesd_multipole_to_local_vertical (src, dst, r,
                                                          cost, cosp, sinp);
      return;
    }

  aran_spherical_seriesd_alpha_require (d);
  aran_spherical_seriesd_beta_require (d);

  aran_spherical_harmonic_evaluate_multiple_internal (d, cost, sint, expp,
                                                      harmonics);

  inv_r = 1. / r;
  pow = 1.;
  for (l = 0; l <= d; l++)
    {
      rpow[l] = pow;
      pow *= inv_r;
    }

  sign = 1.;
  for (l = 0; l <= dst->posdeg; l++)
    {
      for (m = 0; m <= l; m++)
        {
          dstterm = _spherical_seriesd_get_pos_term (dst, l, m);

          for (n = 0; n < src->negdeg; n++)
            {
              gdouble normaliz = aran_spherical_seriesd_beta (l + n) *
                aran_spherical_seriesd_beta (l) /
                aran_spherical_seriesd_beta (n);
              gcomplex128 sum = 0.;

              srcterm = _spherical_seriesd_get_neg_term (src, n, 0);
              hterm = aran_spherical_harmonic_multiple_get_term (l + n, 0,
                                                                 harmonics);

              sum = aran_spherical_seriesd_alpha (l, m) *
                aran_spherical_seriesd_alpha (n, 0) /
                aran_spherical_seriesd_alpha (l + n, m) * hterm[m] *
                srcterm[0];

              for (o = 1; o <= n; o++)
                {
                  guint m_p_o = m + o;

                  gdouble factor = aran_spherical_seriesd_alpha (l, m) *
                    aran_spherical_seriesd_alpha (n, o);

                  /* h= Y_(l+n)^(m+o) */
                  gcomplex128 h = hterm[m_p_o];

                  sum += h * srcterm[o] * factor /
                    aran_spherical_seriesd_alpha (l + n, m_p_o);

                  m_p_o = ABS (m - o);

                  /* h= Y_(l+n)^(m-o) */
                  h = hterm[m_p_o];
                  if ((m - o) < 0)
                    h = _sph_sym (h, m_p_o);

                  sum += h * _sph_sym (srcterm[o], o) * factor /
                    aran_spherical_seriesd_alpha (l + n, m_p_o);
                }

              *dstterm += conj (sum) * sign * normaliz * rpow[l + n + 1];
            }
        }
      sign = -sign;
    }
}

/**
 * aran_spherical_seriesd_translate:
 * @src: source expansion series.
 * @xsrc: @src center.
 * @dst: destination expansion series.
 * @xdst: @dst center.
 *
 * Computes the translation of @src expansion series and accumulates it into
 * @dst. This is the #AranSphericalSeriesd equivalent of @dst += T*@src. Where
 * T is  the translation operator.
 */
void aran_spherical_seriesd_translate (const AranSphericalSeriesd * src,
                                       const VsgVector3d * xsrc,
                                       AranSphericalSeriesd * dst,
                                       const VsgVector3d * xdst)
{
  VsgVector3d tmp;
  gdouble r, cost, sint, cosp, sinp;

  vsg_vector3d_sub (xdst, xsrc, &tmp);

  vsg_vector3d_to_spherical_internal (&tmp, &r, &cost, &sint, &cosp, &sinp);

  aran_local_translate (src, dst, r, cost, sint, cosp, sinp);

  if (src->negdeg > 0)
    {
      aran_multipole_translate (src, dst, r, -cost, sint, -cosp, -sinp);
    }
}

/**
 * aran_spherical_seriesd_to_local:
 * @src: source expansion series.
 * @xsrc: @src center.
 * @dst: destination expansion series.
 * @xdst: @dst center.
 *
 * Like aran_spherical_seriesd_translate() except the multipole part of
 * @src is transformed into a local expansion in @dst.
 */
void aran_spherical_seriesd_to_local (const AranSphericalSeriesd * src,
                                      const VsgVector3d * xsrc,
                                      AranSphericalSeriesd * dst,
                                      const VsgVector3d * xdst)
{
  VsgVector3d tmp;
  gdouble r, cost, sint, cosp, sinp;

  vsg_vector3d_sub (xdst, xsrc, &tmp);

  vsg_vector3d_to_spherical_internal (&tmp, &r, &cost, &sint, &cosp, &sinp);

  aran_local_translate (src, dst, r, cost, sint, cosp, sinp);

  if (src->negdeg > 0)
    {
      aran_multipole_to_local (src, dst, r, cost, sint, cosp, sinp);
    }
}
