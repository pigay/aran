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

typedef struct _AranSquareBufferd AranSquareBufferd;
typedef gdouble (*AranSquareBufferdGenerator) (guint l, guint m,
                                                  AranSquareBufferd *buf);
struct _AranSquareBufferd {
  gint l;
  gdouble *buffer;
  gdouble **direct;
  AranSquareBufferdGenerator generator;
};

void aran_square_bufferd_require (AranSquareBufferd *buf,
                                  guint max);

AranSquareBufferd *
aran_square_bufferd_new (AranSquareBufferdGenerator generator,
                             guint l)
{
  AranSquareBufferd *ret;

  g_return_val_if_fail (generator != NULL, NULL);

  ret = g_malloc (sizeof (AranSquareBufferd));

  ret->l = -1;
  ret->generator = generator;
  ret->buffer = NULL;
  ret->direct = NULL;

  aran_square_bufferd_require (ret, l);

  return ret;
}
void aran_square_bufferd_free (AranSquareBufferd *buf)
{
  if (buf != NULL)
    {
      if (buf->buffer != NULL)
        g_free (buf->buffer);

      if (buf->direct != NULL)
        g_free (buf->direct);

      buf->buffer = NULL;
      buf->direct = NULL;

      g_free (buf);
    }
}
void aran_square_bufferd_require (AranSquareBufferd *buf,
                                    guint max)
{
  guint l, size;
  gint i, j;

  g_return_if_fail (buf != NULL);

  if (buf->l >= (gint) max) return;

  l = MAX (buf->l, 1);

  while (l < (gint) max) l *= 2;

  size = ((l+1)*(l+1));

  buf->direct = g_realloc (buf->direct, (l+1) * sizeof (gdouble *));
  buf->buffer = g_realloc (buf->buffer, size * sizeof (gdouble));

  for (i=0; i<=l; i ++)
    {
      buf->direct[i] = buf->buffer + (i*(l+1));
    }

  buf->l = l;

  for (i=0; i<=l; i ++)
    {
      for (j=0; j<=l; j ++)
        buf->direct[i][j] = buf->generator (i, j, buf);
    }
}
gdouble *aran_square_bufferd_get_unsafe (AranSquareBufferd *buf,
                                            guint l, guint m)
{
  return &(buf->direct[l][m]);
}


static AranSquareBufferd *_betal_over_betan_buffer = NULL;

static gdouble _betal_over_betan_generator (guint l, guint n, AranSquareBufferd * buf)
{
  return aran_spherical_seriesd_beta (l) /
    aran_spherical_seriesd_beta (n);
}

static void _atexit ()
{
  aran_square_bufferd_free (_betal_over_betan_buffer);
  _betal_over_betan_buffer = NULL;
}

static void _betal_over_betan_require (guint deg)
{
  aran_spherical_seriesd_beta_require (deg);
 
  if (_betal_over_betan_buffer == NULL)
    {
      _betal_over_betan_buffer = aran_square_bufferd_new (_betal_over_betan_generator, deg);

      g_atexit (_atexit);
    }
  else
    {
      aran_square_bufferd_require (_betal_over_betan_buffer, deg);
    }
}

static gdouble _betal_over_betan (guint l, guint n)
{
  return *aran_square_bufferd_get_unsafe (_betal_over_betan_buffer, l, n);
}

typedef struct _AranTranslateBufferd AranTranslateBufferd;
typedef gdouble (*AranTranslateBufferdGenerator) (guint l, guint m, guint n,
                                                  AranTranslateBufferd *buf);
struct _AranTranslateBufferd {
  gint l;
  gdouble *buffer;
  gdouble **direct2;
  gdouble ***direct1;
  AranTranslateBufferdGenerator generator;
};

static void aran_translate_bufferd_require (AranTranslateBufferd *buf,
                                            guint max);

static AranTranslateBufferd *
aran_translate_bufferd_new (AranTranslateBufferdGenerator generator,
                            guint l)
{
  AranTranslateBufferd *ret;

  g_return_val_if_fail (generator != NULL, NULL);

  ret = g_malloc (sizeof (AranTranslateBufferd));

  ret->l = -1;
  ret->generator = generator;
  ret->buffer = NULL;
  ret->direct2 = NULL;
  ret->direct1 = NULL;

  aran_translate_bufferd_require (ret, l);

  return ret;
}

static void aran_translate_bufferd_free (AranTranslateBufferd *buf)
{
  if (buf != NULL)
    {
      if (buf->buffer != NULL)
        g_free (buf->buffer);

      if (buf->direct1 != NULL)
        g_free (buf->direct1);

      if (buf->direct2 != NULL)
        g_free (buf->direct2);

      buf->buffer = NULL;
      buf->direct1 = NULL;
      buf->direct2 = NULL;

      g_free (buf);
    }
}

static void aran_translate_bufferd_require (AranTranslateBufferd *buf,
                                            guint max)
{
  guint l, size;
  gint i, j, k;

  g_return_if_fail (buf != NULL);

  if (buf->l >= (gint) max) return;

  l = MAX (buf->l, 1);

  while (l < (gint) max) l *= 2;

  size = ((l+1)*(l+1)*(l+2))/2;

  buf->direct1 = g_realloc (buf->direct1, (l+1) * sizeof (gdouble *));
  buf->direct2 = g_realloc (buf->direct2, ((l+1)*(l+2))/2 * sizeof (gdouble *));
  buf->buffer = g_realloc (buf->buffer, size * sizeof (gdouble));

  for (i=0; i<=l; i ++)
    {
      for (j=0; j<=i; j++)
        buf->direct2[(i*(i+1))/2 + j] = buf->buffer + (l+1)*((i*(i+1))/2+j);

      buf->direct1[i] = buf->direct2 + (i*(i+1))/2;
    }

  buf->l = l;

  for (i=0; i<=l; i ++)
    {
      for (j=0; j<=i; j ++)
        {
          for (k=0; k<=l; k++)
            {
              /* g_printerr ("%u %u %u %ld\n", i, j, k, &buf->direct1[i][j][k] - buf->buffer); */
              buf->direct1[i][j][k] = buf->generator (i, j, k, buf);
            }
        }
    }
}

static inline gdouble *aran_translate_bufferd_get_unsafe (AranTranslateBufferd *buf,
                                                          guint l, guint m, guint n)
{
  return &(buf->direct1[l][m][n]);
}

static AranTranslateBufferd *_precomputed_translate_vertical_buffer = NULL;

static gdouble _precomputed_translate_vertical_generator (guint l, guint m, guint n, AranTranslateBufferd * buf)
{
  gdouble normaliz = _betal_over_betan (l, n);

  gdouble factor = aran_spherical_seriesd_alpha (l, m) *
    aran_spherical_seriesd_alpha (n, m);

  return normaliz * factor /
    aran_spherical_seriesd_alpha (l + n, 0);
}

static void _atexit2 ()
{
  aran_translate_bufferd_free (_precomputed_translate_vertical_buffer);
  _precomputed_translate_vertical_buffer = NULL;
}

static void _precomputed_translate_vertical_require (guint deg)
{
  aran_spherical_seriesd_beta_require (deg);
  aran_spherical_seriesd_alpha_require (deg+deg);
  _betal_over_betan_require (deg);

  if (_precomputed_translate_vertical_buffer == NULL)
    {
      _precomputed_translate_vertical_buffer = aran_translate_bufferd_new (_precomputed_translate_vertical_generator, deg);

      g_atexit (_atexit2);
    }
  else
    {
      aran_translate_bufferd_require (_precomputed_translate_vertical_buffer, deg);
    }
}

static inline gdouble _precomputed_translate_vertical (guint l, guint m, guint n)
{
  return *aran_translate_bufferd_get_unsafe (_precomputed_translate_vertical_buffer, l, m, n);
}

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
  gcomplex128 *srcterm, *dstterm;
  gint d = MAX (src->posdeg, dst->posdeg);

  if (src->posdeg > dst->posdeg)
    g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__);

  aran_spherical_seriesd_beta_require (d);
  aran_spherical_seriesd_alpha_require (d);
  _betal_over_betan_require (d);

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
              gdouble normaliz = _betal_over_betan (l, n);
              /* gdouble normaliz = aran_spherical_seriesd_beta (l) / */
              /*   aran_spherical_seriesd_beta (n); */
              gdouble factor;
	      gdouble h;

              srcterm = _spherical_seriesd_get_pos_term (src, n, 0);
 
              /* o-m = 0 */
              factor =
                aran_spherical_seriesd_alpha (n - l, 0) *
                aran_spherical_seriesd_alpha (l, m) /
                aran_spherical_seriesd_alpha (n, m);

	      /* h= Y_(n-l)^(o-m) */

              /*
               * in this case, h=Y_(n-l)^0, which simplifies with "normaliz"
               * removing beta(n-l)
               * and then becomes h = P_(n-l)^0 = (cost)^(n-l)
               */
              h = ((n-l)%2 == 0)? 1. : cost;

              *dstterm += h * srcterm[m] * factor * normaliz *
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
  _betal_over_betan_require (d);

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
                _betal_over_betan (l, n);
                /* aran_spherical_seriesd_beta (l) / */
                /* aran_spherical_seriesd_beta (n); */
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
  gcomplex128 *srcterm, *dstterm;
  gint d = MAX (src->negdeg, dst->negdeg) - 1;

  if (src->negdeg > dst->negdeg)
    g_warning ("could loose precision in \"%s\"\n", __PRETTY_FUNCTION__);

  aran_spherical_seriesd_alpha_require (d);
  aran_spherical_seriesd_beta_require (d);
  _betal_over_betan_require (d+1);


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
              gdouble normaliz = _betal_over_betan (l, n);
              /* gdouble normaliz = aran_spherical_seriesd_beta (l) / */
              /*   aran_spherical_seriesd_beta (n); */
              gdouble factor;
	      gdouble h;
              srcterm = _spherical_seriesd_get_neg_term (src, n, 0);

              /* m-o = 0 */
              factor =
                aran_spherical_seriesd_alpha (l - n, 0) *
                aran_spherical_seriesd_alpha (n, m) /
                aran_spherical_seriesd_alpha (l, m);

	      /* h= Y_(l-n)^(m-o) */

              /*
               * in this case, h=Y_(l-n)^0, which simplifies with "normaliz"
               * removing beta(l-n)
               * and then becomes h = P_(l-n)^0 = (cost)^(l-n)
               */
              h = ((l-n)%2 == 0)? 1. : cost;

              *dstterm += h * srcterm[m] * factor * normaliz *
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
  _betal_over_betan_require (d);

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
                _betal_over_betan (l, n);
              /* gdouble normaliz = aran_spherical_seriesd_beta (l - n) * */
              /*   aran_spherical_seriesd_beta (l) / */
              /*   aran_spherical_seriesd_beta (n); */
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

/*
 * Compute M2L translation in a vertical direction, positive or negative,
 * depending on the value of cost (+1. or -1.)
 *
 * Formula:
 * \check{Y}_l^m (\vec{r_0'}) =
 * (-1)^l \sum_{n=m}^{\infty}
 * \frac{1}{|\vec{r_0'}-\vec{r_0}|^{l+n+1}}
 * \frac{\beta_l \beta_{l+n}}{\beta_n}
 * \frac{\alpha_l^m \alpha_n^{-m}}{\alpha_{l+n}^{0}}
 * \overline{Y_{l+n}^{0}(\theta, \phi_{\vec{r_0'}-\vec{r_0}})}
 * \overline{\hat{Y}_{n}^{-m}(\vec{r_0})}
 */
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
  gcomplex128 *dstterm;

  aran_spherical_seriesd_alpha_require (d);
  aran_spherical_seriesd_beta_require (d);
  _betal_over_betan_require (MAX(dst->posdeg, src->negdeg));
  _precomputed_translate_vertical_require(MAX(dst->posdeg, src->negdeg));

  inv_r = 1. / r;
  pow = 1.;
  for (l = 0; l <= d; l++)
    {
      rpow[l] = pow;

      /* o == -m */
      /* h= Y_(l+n)^(m+o) */
      /*
       * in this case, h=Y_(l+n)^0, which simplifies with beta(l+n)
       * and then becomes h = P_(l+n)^0 = (cost)^(l+n)
       * we integrate (cost)^(l+n) within rpow[l+n+1]
       */
      if (l%2 == 0) rpow[l] *= cost;

      pow *= inv_r;
    }

  for (l = 0; l <= dst->posdeg; l++)
    {
      for (m = 0; m <= l; m++)
        {
          gcomplex128 sum = 0.;

          dstterm = _spherical_seriesd_get_pos_term (dst, l, m);

          for (n = m; n < src->negdeg; n++)
            {
              gcomplex128 srcterm;
              gdouble translate_factor;

              srcterm = *_spherical_seriesd_get_neg_term (src, n, m);

              /* translate_factor = beta(l)/beta(n) * alpha(l,m) * alpha(n,m) /
               * alpha (l + n, 0);
               */
              translate_factor = _precomputed_translate_vertical (l, m, n);

              sum += srcterm * (rpow[l + n + 1] * translate_factor);
            }

          /* combination of (-1)^l and Y_n^(-m)*/
          sum = ((l+m)%2 == 0)? sum : -sum;

          *dstterm += sum;
        }
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
  _betal_over_betan_require (MAX(dst->posdeg, src->negdeg));

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
                _betal_over_betan (l, n);
                /* aran_spherical_seriesd_beta (l) / */
                /* aran_spherical_seriesd_beta (n); */
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
