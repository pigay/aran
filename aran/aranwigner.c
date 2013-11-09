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

static gboolean aran_wigner_require_d (AranWigner * aw, guint l)
{
  gint i, j, k;
  gdouble cb, sb, cb2, sb2, tb2;
  gcomplex128 d1_0_0, d1_1_1;

  if (((gint) l) <= aw->lmax)
    return FALSE;

  aran_wigner_realloc (aw, l);

  cb = cos (aw->beta);
  sb = sin (aw->beta);
  cb2 = cos (aw->beta * 0.5);
  sb2 = sin (aw->beta * 0.5);
  tb2 = sb2 / cb2;

  aw->l_terms[0][0][0] = 1.; /* d_0^{0,0} */

  if (l == 0)
    return TRUE;

  aw->l_terms[1][0][1 + 0] = cb;                        /* d_1^{0,0} */
  aw->l_terms[1][1][1 - 1] = sb2 * sb2;                 /* d_1^{1,-1} */
  aw->l_terms[1][1][1 + 0] = -sb / sqrt (2.);           /* d_1^{1,0} */
  aw->l_terms[1][1][1 + 1] = cb2 * cb2;                 /* d_1^{1,1} */
  aw->l_terms[1][0][1 + 1] = -aw->l_terms[1][1][1 + 0]; /* d_1^{0,1} */
  aw->l_terms[1][0][1 - 1] = aw->l_terms[1][1][1 + 0];  /* d_1^{0,-1} */

  if (((gint) l) <= 1)
    return TRUE;

  d1_0_0 = aw->l_terms[1][0][1 + 0];
  d1_1_1 = aw->l_terms[1][1][1 + 1];

  /* l >= 2 */
  for (i = 2; i <= l; i++)
    {
      gdouble two_i_m_1 = i + i - 1.;
      gdouble sq_i = i * i;
      gdouble sq_i_m_1 = (i - 1.) * (i - 1.);

      /* j > 0, -j <= k <= j */
      for (j = 0; j <= i - 2; j++)
        {
          gdouble sq_j = j * j;

          for (k = -j; k <= j; k++)
            {
              gdouble sq_k = k * k;
              gdouble a =
                (i * two_i_m_1) / sqrt ((sq_i - sq_j) * (sq_i - sq_k));
              gcomplex128 b = (d1_0_0 - ((j * k) / (i * (i - 1.))));
              gdouble c = sqrt ((sq_i_m_1 - sq_j) * (sq_i_m_1 - sq_k)) /
                ((i - 1.) * two_i_m_1);

              aw->l_terms[i][j][i + k] =
                a * (b * aw->l_terms[i - 1][j][(i - 1) + k] -
                     c * aw->l_terms[i - 2][j][(i - 2) + k]);
            }
        }

      /* last two diagonal terms */
      aw->l_terms[i][i][i + i] = d1_1_1 *
        aw->l_terms[i - 1][i - 1][(i - 1) + (i - 1)]; /* d_i^{i,i} */

      aw->l_terms[i][i - 1][i + i - 1] = (i * d1_0_0 - i + 1.) *
        aw->l_terms[i - 1][i - 1][(i - 1) + (i - 1)]; /* d_i^{i-1,i-1} */

      /* last column terms */
      for (k = i; k >= 1 - i; k--)
        {
          aw->l_terms[i][i][i + k - 1] = -sqrt ((i + k) / (i - k + 1.)) * tb2 *
            aw->l_terms[i][i][i + k];
        }

      /* penultimate column */
      for (k = i - 1; k >= 2 - i; k--)
        {
	  gdouble a = sqrt (((gdouble)i+k)/(((gdouble)i+i)*(i-k+1.)));

          aw->l_terms[i][i - 1][i + k - 1] = (i*cb-k+1.) * a *
	    aw->l_terms[i][i][i + k] / d1_1_1;
        }

      /* extra diagonal terms (|k| > j) */
      for (k = 1; k <= i; k++)
        {
          for (j = 0; j < k; j++)
            {
              aw->l_terms[i][j][i + k] = PHASE (j + k) *
                aw->l_terms[i][k][i + j];
              aw->l_terms[i][j][i - k] = aw->l_terms[i][k][i - j];
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
              for (m = -l; m < 0; m++)
                {
                  aw->l_terms[l][mprime][l + m] *= /* PHASE (mprime + m) * */ conj(expma[ABS(m)]);
                }

              for (m = 0; m <= l; m++)
                aw->l_terms[l][m][l + mprime] *= expma[m];
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
              for (m = -l; m < 0; m++)
                aw->l_terms[l][mprime][l + m] *= expmg[mprime];
              for (m = 0; m <= l; m++)
                aw->l_terms[l][m][l + mprime] *= expmg[mprime];
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
      for (j = 0; j <= i; j++)
        {
          for (k = -i; k <= i; k++)
            {
              gcomplex128 term = *aran_wigner_term (aw, i, j, k);
              g_fprintf (file, "[%d,%d,%d]=%g+i%g ", i, j, k,
                         creal(term), cimag(term));
            }
        }
    }

  g_fprintf (file, "]");

}
