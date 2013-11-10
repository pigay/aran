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

#include "aranwigner.h"

#include "aranwigner-private.h"

#include <math.h>
#include <string.h>

#include <glib/gprintf.h>

#define PHASE(m) (((m)%2 == 0) ? 1. : -1.)

/**
 * AranWigner:
 *
 * Opaque structure. Possesses only private data.
 */

static guint get_l_mprime_m_size (guint l)
{
  return ((l + 1) * (l + 2) * (4 * l + 3)) / 6;
}

static guint get_l_mprime_size (guint l)
{
  return ((l + 1) * (l + 2)) / 2;
}

static guint get_l_size (guint l)
{
  return l + 1;
}

static void aran_wigner_realloc (AranWigner * aw, guint l)
{
  gint i, j, l_mprime;
  guint l_mprime_m_size = get_l_mprime_m_size (l);
  guint l_mprime_size = get_l_mprime_size (l);
  guint l_size = get_l_size (l);
  guint alloc_size = l_mprime_m_size * sizeof (gcomplex128) +
    l_mprime_size * sizeof (gcomplex128 *) + l_size * sizeof (gcomplex128 **);

  aw->l_mprime_m_terms = g_realloc (aw->l_mprime_m_terms, alloc_size);

  aw->l_mprime_terms = (gcomplex128 **) (aw->l_mprime_m_terms + l_mprime_m_size);
  aw->l_terms = (gcomplex128 ***) (aw->l_mprime_terms + l_mprime_size);

  aw->l_mprime_terms[0] = aw->l_mprime_m_terms;
  aw->l_terms[0] = aw->l_mprime_terms;

  l_mprime = 1;
  for (i = 1; i <= l; i++)
    {
      aw->l_mprime_terms[l_mprime] = aw->l_mprime_terms[l_mprime - 1] + (2 * i - 1);
      aw->l_terms[i] = aw->l_mprime_terms + l_mprime;
      l_mprime++;
      for (j = 1; j <= i; j++)
        {
          aw->l_mprime_terms[l_mprime] = aw->l_mprime_terms[l_mprime - 1] + (2 * i + 1);
          l_mprime++;
        }
    }
}

static gboolean aran_wigner_require_d (AranWigner * aw, guint lmax)
{
  gint l, mp, m;
  gdouble cb, sb, cb2, sb2, tb2;
  gcomplex128 d1_0_0, d1_1_1;
  gcomplex128 d1_1_m1;

  if (((gint) lmax) <= aw->lmax)
    return FALSE;

  aran_wigner_realloc (aw, lmax);

  cb = cos (aw->beta);
  sb = sin (aw->beta);
  cb2 = cos (aw->beta * 0.5);
  sb2 = sin (aw->beta * 0.5);
  tb2 = sb2 / cb2;

  aw->l_terms[0][0][0] = 1.; /* d_0^{0,0} */

  if (lmax == 0)
    return TRUE;

  aw->l_terms[1][0][1 + 0] = cb;                         /* d_1^{0,0} */
  aw->l_terms[1][1][1 - 1] = sb2 * sb2;                  /* d_1^{-1,1} */
  aw->l_terms[1][1][1 + 0] = sb / sqrt (2.);             /* d_1^{0,1} */
  aw->l_terms[1][1][1 + 1] = cb2 * cb2;                  /* d_1^{1,1} */
  aw->l_terms[1][0][1 - 1] = aw->l_terms[1][1][1 + 0];   /* d_1^{0,-1} */
  aw->l_terms[1][0][1 + 1] = - aw->l_terms[1][1][1 + 0]; /* d_1^{0,1} */

  if (((gint) lmax) <= 1)
    return TRUE;

  d1_0_0 = aw->l_terms[1][0][1 + 0];
  d1_1_1 = aw->l_terms[1][1][1 + 1];
  d1_1_m1 = aw->l_terms[1][1][1 - 1];

  /* l >= 2 */
  for (l = 2; l <= lmax; l++)
    {
      gdouble two_l_m_1 = l + l - 1.;
      gdouble sq_l = l * l;
      gdouble sq_l_m_1 = (l - 1.) * (l - 1.);

      /* block 1 */
      for (mp = 0; mp <= l - 2; mp++)
        {
          gdouble sq_mp = mp * mp;

          for (m = -mp; m <= mp; m++)
            {
              gdouble sq_m = m * m;
              gdouble a =
                (l * two_l_m_1) / sqrt ((sq_l - sq_mp) * (sq_l - sq_m));
              gcomplex128 b = (d1_0_0 - ((mp * m) / (l * (l - 1.))));
              gdouble c = sqrt ((sq_l_m_1 - sq_mp) * (sq_l_m_1 - sq_m)) /
                ((l - 1.) * two_l_m_1);

              /* d_l^{m,m'} = \frac{l(2l-1)}{\sqrt{(l^2-m^2)(l^2-(m')^2}}
               * \left\{ (d_1^{00}-\frac{mm'}{l(l-1)})d_{l-1}^{mm'} -
               * \frac{\sqrt{[(l-1)^2-m^2][(m-1)^2-(m')^2]}}{(l-1)(2l-1)}d_{l-2}^{mm'} \right\}
               */
              aw->l_terms[l][mp][l + m] =
                a * (b * aw->l_terms[l - 1][mp][(l - 1) + m] -
                     c * aw->l_terms[l - 2][mp][(l - 2) + m]);
            }
        }

      /* block 2 */
      /* last two diagonal terms */
      /* d_l^{l,l}  = ld_1^{1,1} d_{l-1}^{l-1,l-1} */
      aw->l_terms[l][l][l + l] = d1_1_1 *
        aw->l_terms[l - 1][l - 1][(l - 1) + (l - 1)];

      /* d_l^{l-1,l-1}  = (ld_1^{0,0} - l + 1)  d_{l-1}^{l-1,l-1} */
      aw->l_terms[l][l - 1][l + l - 1] = (l * d1_0_0 - l + 1.) *
        aw->l_terms[l - 1][l - 1][(l - 1) + (l - 1)];

      /* block 3 */
      /* d_l^{l,-l} = d_1^{1,-1} d_{l-1}^{l-1,-l+1} */
      aw->l_terms[l][l][l - l] = d1_1_m1 *
        aw->l_terms[l - 1][l - 1][(l - 1) - (l - 1)];

      /* d_l^{l-1,-l+1} = (ld_1^{0,0} + l - 1)  d_{l-1}^{l-1,-l+1} */
      aw->l_terms[l][l - 1][l - (l - 1)] = (l * d1_0_0 + l - 1.) *
        aw->l_terms[l - 1][l - 1][(l - 1) - (l - 1)];

      /* block 4 */
      /* last column terms */
      for (mp = l; mp >= 1; mp--)
        {
          /* 1\le m' \le l : d_l^{l,m'-1} = -\sqrt{\frac{(l+m')}{(l-m'+1)}} \tan\frac{\beta}{2}
           * d_l^{l,m'}
           */
          aw->l_terms[l][mp-1][l + l] = -sqrt ((l + mp) / (l - mp + 1.)) * tb2 *
            aw->l_terms[l][mp][l + l]; /* d_l^{l,mp-1} */
        }

      /* block 5 */
      /* penultimate column */
      for (mp = l - 1; mp >= 1; mp--)
        {
	  gdouble a = sqrt (((gdouble)l+mp)/(((gdouble)l+l)*(l-mp+1.)));

          /* 1\le m' \lt l : d_l^{l,m'-1} = (l\cos\frac{\beta}{2}-m'+1)\sqrt{\frac{(l+m')}{2l(l-m'+1)}}
           * \frac{d_l^{l,m'}}{d_1^{1,1}}
           */
          aw->l_terms[l][mp - 1][l + l - 1] = (l*cb-mp+1.) * a *
	    aw->l_terms[l][mp][l + l] / d1_1_1;
        }

      /* block 6 */
      /* last rows */
      for (mp = l - 1; mp <= l; mp++)
        {
          for (m = 0; m < mp; m++)
            {
              /* l-1\le m'\le l , 0\le m < m' :
               * d_l^{m,m'} = (-1)^{m+m'} d_l^{m',m}
               */
              aw->l_terms[l][mp][l + m] = PHASE (mp + m) *
                aw->l_terms[l][m][l + mp];
            }
        }

      /* block 7 */
      for (m = 0; m < l; m++)
        {
          /* 0\le m < l :
           * d_l^{-m-1,l} = \sqrt{\frac{l-m}{l+m+1}} \tan{\frac{\beta}{2}} d_l^{-m,l}
           */
          aw->l_terms[l][l][l - m - 1] =  sqrt ((l - m) / (l + m + 1.)) * tb2 *
            aw->l_terms[l][l][l - m];
        }

      /* block 8 */
      for (m = 0; m < l; m++)
        {
          gdouble a = sqrt (((gdouble)l-m)/(((gdouble)l+l)*(l+m+1.)));

          /* 0\le m < l :
           * d_l^{-m-1,l-1} = (l\cos{\beta}+m+1) \sqrt{\frac{l-m}{2l(l+m+1)}} \frac{d_l^{-m,l}}{\cos^2\frac{\beta}{2}}
           */
          aw->l_terms[l][l-1][l - m - 1] =  (l*cb+m+1.) * a *
            aw->l_terms[l][l][l - m] / d1_1_1;
        }

      /* block 10 */
      for (mp = 0; mp <= l; mp++)
        {
          for (m = mp+1; m <= l; m++)
            {
              /* 0\le m' \le l , m'+1 \le m \le l : d_l^{m,m'} = (-1)^{m+m'}d_l^{m',m}
               */
              aw->l_terms[l][mp][l + m] = PHASE (m + mp) * aw->l_terms[l][m][l + mp];
              /* 0\le m' \le l , m'+1 \le m \le l : d_l^{-m,m'} = d_l^{-m',m}
               */
              aw->l_terms[l][mp][l - m] = aw->l_terms[l][m][l - mp];
            }
        }
    }

  return TRUE;
}

/**
 * aran_wigner_require:
 * @aw: an #AranWigner structure.
 * @l: requested degree.
 *
 * Ensures that @aw is allocated to @l degree inclusive and that coefficients
 * are correctly computed.
 */
void aran_wigner_require (AranWigner * aw, guint lmax)
{
  /* compute new d_l^{mm'} coefficients or exit */
  if (!aran_wigner_require_d (aw, lmax)) return;

  /* compute e^{-im\alpha} and e^{-im\gamma} coefficients */
  if (ABS (aw->alpha) >= 1.e-5)
    {
      gcomplex128 expa = cos (aw->alpha) - G_I * sin (aw->alpha);
      gcomplex128 *expma = g_alloca ((lmax+1) * sizeof (gcomplex128));
      gint l, mprime, m;

      expma[0] = 1.;
      for (mprime=1; mprime<=lmax; mprime++)
        {
          expma[mprime] = expma[mprime-1] * expa;
        }

      for (l=0; l<=lmax; l++)
        {
          for (mprime = 0; mprime <= l; mprime++)
            {
              for (m = l; m >= 1; m--)
                {
                  aw->l_terms[l][mprime][l - m] *= conj(expma[m]); /* d_l^{-m,m'} <- e^(im\alpha} d_l^{-m,m'} */
                }
              for (m = 0; m <= l; m++)
                {
                  aw->l_terms[l][mprime][l + m] *= expma[m];       /* d_l^{m,m'} <- e^(-im\alpha} d_l^{m,m'} */
                }
            }
        }
    }

  if (ABS (aw->gamma) >= 1.e-5)
    {
      gcomplex128 expg = cos (aw->gamma) - G_I * sin (aw->gamma);
      gcomplex128 *expmg = g_alloca ((lmax+1) * sizeof (gcomplex128));
      gint l, mprime, m;

      expmg[0] = 1.;
      for (mprime=1; mprime<=lmax; mprime++)
        {
          expmg[mprime] = expmg[mprime-1] * expg;
        }

      for (l=0; l<=lmax; l++)
        {
          for (mprime = 0; mprime <= l; mprime++)
            {
              for (m = -l; m <= l; m++)
                aw->l_terms[l][mprime][l + m] *= expmg[mprime]; /* D_l^{m,m'} = e^(-im\alpha} d_l^{m,m'} e^{-im'\gamma} */
            }
        }
    }

  aw->lmax = lmax;
}

/**
 * aran_wigner_term:
 * @aw: an #AranWigner structure.
 * @l: a #guint.
 * @mprime: a #guint. Condition @mprime <= @l must hold.
 * @m: a #gint. Condition -@l <= @m <= @l must hold.
 *
 * Provides access to the @aw term corresponding to Wigner coefficient d_l^{m,m'}.
 *
 * Returns: address of the term.
 */
gcomplex128 *aran_wigner_term (AranWigner * aw, guint l, guint mprime, gint m)
{
  g_return_val_if_fail (aw != NULL, NULL);
  return ARAN_WIGNER_TERM (aw, l, mprime, m);
  /* return aw->l_terms[l][mprime] + l + m; */
}

/**
 * aran_wigner_new:
 * @alpha: Rotation angle (-pi < @alpha < pi).
 * @beta: Rotation angle (-pi < @beta < pi).
 * @gamma: Rotation angle (-pi < @gamma < pi).
 * @l: initial degree preallocation if positive.
 *
 * Creates a new #AranWigner for angle @beta. 
 *
 * Returns: A Wigner D-function coefficients buffer.
 */
AranWigner *aran_wigner_new (gdouble alpha, gdouble beta, gdouble gamma, gint l)
{
  AranWigner *aw = g_malloc (sizeof (AranWigner));

  aw->alpha = alpha;
  aw->beta = beta;
  aw->gamma = gamma;

  aw->lmax = -1;

  aw->l_mprime_m_terms = NULL;

  if (l >= 0)
    aran_wigner_require (aw, l);

  return aw;
}

/**
 * aran_wigner_free:
 * @aw: an #AranWigner structure.
 *
 * Frees all memory allocated for @aw.
 */
void aran_wigner_free (AranWigner * aw)
{
  g_return_if_fail (aw != NULL);

  if (aw->l_mprime_m_terms != NULL)
    g_free (aw->l_mprime_m_terms);

  g_free (aw);
}

/**
 * aran_wigner_copy:
 * @src: an #AranWigner structure.
 * @dst: an #AranWigner structure.
 *
 * Copies @src contents into @dst.
 */
void aran_wigner_copy (AranWigner *src, AranWigner *dst)
{
  dst->beta = src->beta;

  aran_wigner_realloc (dst, src->lmax);

  memcpy(dst->l_terms[0][0], src->l_terms[0][0],
         get_l_mprime_m_size (dst->lmax) * sizeof (gcomplex128));
}

/**
 * aran_wigner_write:
 * @aw: an #AranWigner structure.
 * @file: output file.
 *
 * Writes @aw to @file.
 */
void aran_wigner_write (AranWigner * aw, FILE * file)
{
  gint i, j, k;

  g_return_if_fail (aw != NULL);

  g_fprintf (file, "[");

  for (i = 0; i <= aw->lmax; i++)
    {
      for (k = -i; k <= i; k++)
        {
          for (j = 0; j <= i; j++)
            {
              gcomplex128 term = *aran_wigner_term (aw, i, j, k);
              g_fprintf (file, "[%d,%d,%d]=%g+i%g ", i, k, j,
                         creal(term), cimag(term));
            }
        }
    }

  g_fprintf (file, "]");

}
