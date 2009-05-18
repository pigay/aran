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

/* #include <vsg/vsgd-inline.h> */

#include <string.h>
#include <math.h>

/* functions */

static void aran_multipole_translate_kkylin (const AranSphericalSeriesd * src,
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
  gcomplex128 buf1[(d + 1) * (d + 1)], buf2[(d + 1) * (d + 1)];
  gcomplex128 *row1[d + 1], *row2[d + 1];
  gcomplex128 **partial, **partial_m_1, **flip;
  gcomplex128 f_sint_expp;

/*   if (src->negdeg > dst->negdeg) */
/*     g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__); */

  g_return_if_fail (src->negdeg <= dst->negdeg);

  /* partial sums buffers */
  partial = row1;
  partial_m_1 = row2;

  for (m = 0; m <= d; m++)
    {
      row1[m] = buf1 + m * (d + 1);
      row2[m] = buf2 + m * (d + 1);
    }

  aran_spherical_seriesd_alpha_require (d);
  aran_spherical_seriesd_beta_require (d);

  aran_spherical_harmonic_evaluate_multiple_internal (dst->negdeg - 1, cost,
                                                      sint, expp, harmonics);

  pow = 1.;
  rpow[0] = pow;
  pow *= r;

  srcterm = _spherical_seriesd_get_neg_term (src, 0, 0);
  hterm = aran_spherical_harmonic_multiple_get_term (0, 0, harmonics);

  partial[0][0] = conj (hterm[0]) * srcterm[0];

  dstterm = _spherical_seriesd_get_neg_term (dst, 0, 0);
  dstterm[0] += partial[0][0] * aran_spherical_seriesd_beta (0);

  for (l = 1; l < dst->negdeg; l++)
    {
      gdouble fact;

      rpow[l] = pow;
      pow *= r;

      flip = partial;
      partial = partial_m_1;
      partial_m_1 = flip;

      dstterm = _spherical_seriesd_get_neg_term (dst, l, 0);

      for (n = 0; n <= l - 1; n++)
        {
          hterm = aran_spherical_harmonic_multiple_get_term (l - n, 0,
                                                             harmonics);
          srcterm = _spherical_seriesd_get_neg_term (src, n, 0);

          fact = aran_spherical_seriesd_beta (l - n) /
            aran_spherical_seriesd_beta (n);

          /* partial sum: l, 0, n */
          partial[0][n] = 0.;

          for (o = MAX (-n, n - l); o <= MIN (n, l - n); o++)
            {
              guint abs_o = ABS (o);
              gcomplex128 h = hterm[abs_o];     /* h = Y_(l-n)^(-o) */
              gcomplex128 s = srcterm[abs_o];

              if (o > 0)
                h = _sph_sym (h, abs_o);
              else if (o < 0)
                s = _sph_sym (s, abs_o);

              partial[0][n] +=
                (rpow[l - n] * aran_spherical_seriesd_alpha (n, abs_o) *
                 aran_spherical_seriesd_alpha (l - n, abs_o)) * conj (h) * s;
            }

          partial[0][n] *= fact;

          if (l > 1)            /* don't overwrite partial[0][n] */
            {
              /* partial sum: l, l-1, n */
              partial[l - 1][n] = 0.;

              for (o = MAX (-n, n - 1); o <= MIN (n, l + l - 1 - n); o++)
                {
                  guint abs_o = ABS (o);
                  gcomplex128 h = hterm[ABS (l - o - 1)]; /* h = Y_(l-n)^(l-o-1) */
                  gcomplex128 s = srcterm[abs_o];

                  if ((l - o - 1) < 0)
                    h = _sph_sym (h, ABS (l - o - 1));

                  if (o < 0)
                    s = _sph_sym (s, abs_o);

                  partial[l - 1][n] +=
                    (rpow[l - n] * aran_spherical_seriesd_alpha (n, abs_o) *
                     aran_spherical_seriesd_alpha (l - n, ABS (l - o - 1))) *
                    conj (h) * s;
                }
              partial[l - 1][n] *= fact;
            }

          /* partial sum: l, l, n */
          partial[l][n] = 0.;

          for (o = MAX (-n, n); o <= MIN (n, l + l - n); o++)
            {
              guint abs_o = ABS (o);
              gcomplex128 h = hterm[ABS (l - o)];       /* h = Y_(l-n)^(l-o) */
              gcomplex128 s = srcterm[abs_o];

              if ((l - o) < 0)
                h = _sph_sym (h, ABS (l - o));

              if (o < 0)
                s = _sph_sym (s, abs_o);

              partial[l][n] +=
                (rpow[l - n] * aran_spherical_seriesd_alpha (n, abs_o) *
                 aran_spherical_seriesd_alpha (l - n, ABS (l - o))) *
                conj (h) * s;
            }

          partial[l][n] *= fact;

          /* recurrence formula on partial sums */
          fact = 0.5 * r / (l - n);
          f_sint_expp = fact * sint * expp;
          fact *= 2. * cost;
          for (m = 1; m < l - 1; m++)
            {
              partial[m][n] = f_sint_expp * partial_m_1[m + 1][n] -
                conj (f_sint_expp) * partial_m_1[m - 1][n] +
                fact * partial_m_1[m][n];
            }
        }

      /* partial[*][l] direct formulas */
      fact = aran_spherical_seriesd_beta (0) / aran_spherical_seriesd_beta (l);
      hterm = aran_spherical_harmonic_multiple_get_term (0, 0, harmonics);
      for (m = 0; m <= l; m++)
        {
          partial[m][l] = (aran_spherical_seriesd_alpha (l, m) * fact) *
            conj (hterm[0]) *
            *_spherical_seriesd_get_neg_term (src, l, m);
        }

      /* set multipole terms from partial sums */
      for (m = 0; m <= l; m++)
        {
          gcomplex128 sum = 0.;

          for (n = 0; n <= l; n++)
            sum += partial[m][n];

          dstterm[m] += sum * (aran_spherical_seriesd_beta (l) /
                               aran_spherical_seriesd_alpha (l, m));
        }
    }
}

static void aran_local_translate_kkylin (const AranSphericalSeriesd * src,
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
  gcomplex128 buf1[(d + 1) * (d + 1)], buf2[(d + 1) * (d + 1)];
  gcomplex128 *row1[d + 1], *row2[d + 1];
  gcomplex128 **partial, **partial_p_1, **flip;
  gcomplex128 f_sint_expp, y00;

  g_return_if_fail (src->posdeg <= dst->posdeg);

  /* partial sums buffers */
  partial = row1;
  partial_p_1 = row2;

  for (m = 0; m <= d; m++)
    {
      row1[m] = buf1 + m * (d + 1);
      row2[m] = buf2 + m * (d + 1);
    }

  aran_spherical_seriesd_alpha_require (d);
  aran_spherical_seriesd_beta_require (d);

  aran_spherical_harmonic_evaluate_multiple_internal (src->posdeg, cost,
                                                      sint, expp, harmonics);

  y00 = *aran_spherical_harmonic_multiple_get_term (0, 0, harmonics);

  pow = 1.;
  for (l = 0; l <= src->posdeg; l++)
    {
      rpow[l] = pow;
      pow *= r;
    }

  for (l = dst->posdeg; l >= 0; l--)
    {
      gdouble fact;

      dstterm = _spherical_seriesd_get_pos_term (dst, l, 0);

      if (l <= src->posdeg)
        {
          /* partial[*][l] direct formulas */
          srcterm = _spherical_seriesd_get_pos_term (src, l, 0);
          fact = aran_spherical_seriesd_beta (0) / aran_spherical_seriesd_beta (l);
          for (m = 0; m <= l; m++)
            {
              partial[m][l] = fact / aran_spherical_seriesd_alpha (l, m) *
                conj (y00 * srcterm[m]);
            }
        }

      for (n = l + 1; n <= src->posdeg; n++)
        {
          hterm = aran_spherical_harmonic_multiple_get_term (n - l, 0,
                                                             harmonics);
          srcterm = _spherical_seriesd_get_pos_term (src, n, 0);

          fact = aran_spherical_seriesd_beta (n - l) /
            aran_spherical_seriesd_beta (n);

          /* partial sum: l, 0, n */
          partial[0][n] = 0.;

          for (o = MAX (-n, l - n); o <= MIN (n, n - l); o++)
            {
              guint abs_o = ABS (o);
              gcomplex128 h = hterm[abs_o];     /* h = Y_(l-n)^(o) */
              gcomplex128 s = srcterm[abs_o];

              if (o < 0)
                {
                  h = _sph_sym (h, abs_o);
                  s = _sph_sym (s, abs_o);
                }

              partial[0][n] +=
                (rpow[n - l] * aran_spherical_seriesd_alpha (n - l, abs_o) /
                 aran_spherical_seriesd_alpha (n, abs_o)) * conj (h * s);
            }

          partial[0][n] *= fact;

          if (l > 0)            /* don't overwrite partial[0][n] */
            {
              /* partial sum: l, l, n */
              partial[l][n] = 0.;

              for (o = MAX (-n, l + l - n); o <= MIN (n, n); o++)
                {
                  guint abs_o = ABS (o);
                  gcomplex128 h = hterm[ABS (o - l)];   /* h = Y_(n-l)^(o-l) */
                  gcomplex128 s = srcterm[abs_o];

                  if ((o - l) < 0)
                    h = _sph_sym (h, ABS (o - l));

                  if (o < 0)
                    s = _sph_sym (s, abs_o);

                  partial[l][n] +=
                    (rpow[n - l] *
                     aran_spherical_seriesd_alpha (n - l, ABS (o - l)) /
                     aran_spherical_seriesd_alpha (n, abs_o)) * conj (h * s);
                }

              partial[l][n] *= fact;
            }

          /* recurrence formula on partial sums */
          fact = 0.5 * r / (n - l);
          f_sint_expp = fact * sint * expp;
          fact *= 2. * cost;
          for (m = 1; m < l; m++)
            {
              partial[m][n] = f_sint_expp * partial_p_1[m - 1][n] -
                conj (f_sint_expp) * partial_p_1[m + 1][n] +
                fact * partial_p_1[m][n];
/*               partial[m][n] = 0.; */
            }
        }

      /* set local terms from partial sums */
      for (m = 0; m <= l; m++)
        {
          gcomplex128 sum = 0.;

          for (n = l; n <= src->posdeg; n++)
            sum += partial[m][n];

          dstterm[m] += conj (sum) * (aran_spherical_seriesd_beta (l) *
                                      aran_spherical_seriesd_alpha (l, m));
        }

      flip = partial;
      partial = partial_p_1;
      partial_p_1 = flip;

    }
}

static void aran_multipole_to_local_kkylin (const AranSphericalSeriesd * src,
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
  gcomplex128 buf1[(d + 1) * (d + 1)], buf2[(d + 1) * (d + 1)];
  gcomplex128 *row1[d + 1], *row2[d + 1];
  gcomplex128 **partial, **partial_m_1, **flip;
  gcomplex128 two_cott_expp, exp2p, two_l_p_n_inv_r_sint_expp;

/*   if ((src->negdeg-1) > dst->posdeg) */
/*     g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__); */

  g_return_if_fail (src->negdeg - 1 <= dst->posdeg);

  /* vertical translation */
  if (ABS (sint) < 1.e-5)       /* sint in m2l translation is zero or high values */
    {
      aran_spherical_seriesd_multipole_to_local_vertical (src, dst, r,
                                                          cost, cosp, sinp);
      return;
    }

  /* partial sums buffers */
  partial = row1;
  partial_m_1 = row2;

  for (m = 0; m <= d; m++)
    {
      row1[m] = buf1 + m * (d + 1);
      row2[m] = buf2 + m * (d + 1);
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
      gdouble fact;

      dstterm = _spherical_seriesd_get_pos_term (dst, l, 0);

      for (n = 0; n < src->negdeg; n++)
        {
          hterm = aran_spherical_harmonic_multiple_get_term (l + n, 0,
                                                             harmonics);
          srcterm = _spherical_seriesd_get_neg_term (src, n, 0);

          fact = rpow[l + n + 1] * aran_spherical_seriesd_beta (l + n) /
            aran_spherical_seriesd_beta (n);

          /* partial sum: l, 0, n */
          partial[0][n] = 0.;

          for (o = -n; o <= n; o++)
            {
              guint abs_o = ABS (o);
              gcomplex128 h = hterm[abs_o];     /* h = Y_(l+n)^(o) */
              gcomplex128 s = srcterm[abs_o];

              if (o < 0)
                {
                  h = _sph_sym (h, abs_o);
                  s = _sph_sym (s, abs_o);
                }

              partial[0][n] +=
                (aran_spherical_seriesd_alpha (n, abs_o) /
                 aran_spherical_seriesd_alpha (l + n, abs_o))
                * h * s;
            }

          partial[0][n] *= fact;

          if (l > 0)
            {
              /* partial sum: l, 1, n */
              partial[1][n] = 0.;

              for (o = -n; o <= n; o++)
                {
                  guint abs_o = ABS (o);
                  gcomplex128 h = hterm[ABS (o + 1)];   /* h = Y_(l+n)^(o+1) */
                  gcomplex128 s = srcterm[abs_o];

                  if ((o + 1) < 0)
                    h = _sph_sym (h, ABS (o + 1));

                  if (o < 0)
                    s = _sph_sym (s, abs_o);

                  partial[1][n] +=
                    (aran_spherical_seriesd_alpha (n, abs_o) /
                     aran_spherical_seriesd_alpha (l + n, ABS (o + 1))) * h *
                    s;
                }
              partial[1][n] *= fact;
            }


          /* recurrence formula on partial sums */
          two_cott_expp = 2. * cost / sint * expp;
          exp2p = expp * expp;
          two_l_p_n_inv_r_sint_expp = 2. * (l + n) * inv_r / sint * expp;
          for (m = 2; m <= l; m++)
            {
              partial[m][n] = two_cott_expp * partial[m - 1][n] +
                exp2p * partial[m - 2][n] -
                two_l_p_n_inv_r_sint_expp * partial_m_1[m - 1][n];
            }
        }

      /* set multipole terms from partial sums */
      fact = aran_spherical_seriesd_beta (l) * sign;
      for (m = 0; m <= l; m++)
        {
          gcomplex128 sum = 0.;

          for (n = 0; n < src->negdeg; n++)
            sum += partial[m][n];

          dstterm[m] += conj (sum) *
            (aran_spherical_seriesd_alpha (l, m) * fact);
        }

      flip = partial;
      partial = partial_m_1;
      partial_m_1 = flip;
      sign = -sign;
    }
}

/**
 * aran_spherical_seriesd_translate_kkylin:
 * @src: source expansion series.
 * @xsrc: @src center.
 * @dst: destination expansion series.
 * @xdst: @dst center.
 *
 * Performs the same operation as aran_spherical_seriesd_translate() but
 * with a better algorithm (from K. Kylin). Refer to
 * http://cca.arc.nasa.gov/Interns/kkylin-report.pdf for details.
 */
void aran_spherical_seriesd_translate_kkylin (const AranSphericalSeriesd * src,
                                              const VsgVector3d * xsrc,
                                              AranSphericalSeriesd * dst,
                                              const VsgVector3d * xdst)
{
  VsgVector3d tmp;
  gdouble r, cost, sint, cosp, sinp;

  vsg_vector3d_sub (xdst, xsrc, &tmp);

  vsg_vector3d_to_spherical_internal (&tmp, &r, &cost, &sint, &cosp, &sinp);

  aran_local_translate_kkylin (src, dst, r, cost, sint, cosp, sinp);

  if (src->negdeg > 0)
    {
      aran_multipole_translate_kkylin (src, dst, r, -cost, sint, -cosp, -sinp);
    }
}

/**
 * aran_spherical_seriesd_to_local_kkylin:
 * @src: source expansion series.
 * @xsrc: @src center.
 * @dst: destination expansion series.
 * @xdst: @dst center.
 *
 * Performs the same operation as aran_spherical_seriesd_to_local() but
 * with a better algorithm (from K. Kylin). Refer to
 * http://cca.arc.nasa.gov/Interns/kkylin-report.pdf for details.
 */
void aran_spherical_seriesd_to_local_kkylin (const AranSphericalSeriesd * src,
                                             const VsgVector3d * xsrc,
                                             AranSphericalSeriesd * dst,
                                             const VsgVector3d * xdst)
{
  VsgVector3d tmp;
  gdouble r, cost, sint, cosp, sinp;

  vsg_vector3d_sub (xdst, xsrc, &tmp);

  vsg_vector3d_to_spherical_internal (&tmp, &r, &cost, &sint, &cosp, &sinp);

  aran_local_translate_kkylin (src, dst, r, cost, sint, cosp, sinp);

  if (src->negdeg > 0)
    {
      aran_multipole_to_local_kkylin (src, dst, r, cost, sint, cosp, sinp);
    }
}
