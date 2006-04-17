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

#include "aranwigner.h"

#include "aranwigner-private.h"

#include <math.h>

#include <glib/gprintf.h>

#define PHASE(m) (((m)%2 == 0) ? 1. : -1.)

/**
 * AranWigner:
 *
 * Opaque structure. Possesses only private data.
 */

static guint get_l_m1_m2_size (guint l)
{
  return ((l + 1) * (l + 2) * (4 * l + 3)) / 6;
}

static guint get_l_m1_size (guint l)
{
  return ((l + 1) * (l + 2)) / 2;
}

static guint get_l_size (guint l)
{
  return l + 1;
}

static void aran_wigner_realloc (AranWigner * aw, guint l)
{
  gint i, j, l_m1;
  guint l_m1_m2_size = get_l_m1_m2_size (l);
  guint l_m1_size = get_l_m1_size (l);
  guint l_size = get_l_size (l);
  guint alloc_size = l_m1_m2_size * sizeof (gdouble) +
    l_m1_size * sizeof (gdouble *) + l_size * sizeof (gdouble **);

  aw->l_m1_m2_terms = g_realloc (aw->l_m1_m2_terms, alloc_size);

  aw->l_m1_terms = (gdouble **) (aw->l_m1_m2_terms + l_m1_m2_size);
  aw->l_terms = (gdouble ***) (aw->l_m1_terms + l_m1_size);

  aw->l_m1_terms[0] = aw->l_m1_m2_terms;
  aw->l_terms[0] = aw->l_m1_terms;

  l_m1 = 1;
  for (i = 1; i <= l; i++)
    {
      aw->l_m1_terms[l_m1] = aw->l_m1_terms[l_m1 - 1] + (2 * i - 1);
      aw->l_terms[i] = aw->l_m1_terms + l_m1;
      l_m1++;
      for (j = 1; j <= i; j++)
        {
          aw->l_m1_terms[l_m1] = aw->l_m1_terms[l_m1 - 1] + (2 * i + 1);
          l_m1++;
        }
    }
}

/**
 * aran_wigner_require:
 * @aw: an #AranWigner structure.
 * @l: requested degree.
 *
 * Ensures that @aw is allocated to @l degree inclusive and that coefficients
 * are correctly computed.
 */
void aran_wigner_require (AranWigner * aw, guint l)
{
  gint i, j, k;
  gdouble cb, sb, cb2, sb2, tb2, d1_0_0, d1_1_1;

  if (((gint) l) <= aw->lmax)
    return;

  aran_wigner_realloc (aw, l);

  cb = cos (aw->beta);
  sb = sin (aw->beta);
  cb2 = cos (aw->beta * 0.5);
  sb2 = sin (aw->beta * 0.5);
  tb2 = sb2 / cb2;

  /* l == 0 */
  if (aw->lmax < 0)
    {
      aw->l_terms[0][0][0] = 1.;
      aw->lmax = 0;
    }

  if (l == 0)
    return;

  /* l == 1 */
  if (aw->lmax == 0)
    {
      aw->l_terms[1][0][1 + 0] = cb;
      aw->l_terms[1][1][1 - 1] = sb2 * sb2;
      aw->l_terms[1][1][1 + 0] = -sb / sqrt (2.);
      aw->l_terms[1][1][1 + 1] = cb2 * cb2;
      aw->l_terms[1][0][1 + 1] = -aw->l_terms[1][1][1 + 0];
      aw->l_terms[1][0][1 - 1] = aw->l_terms[1][1][1 + 0];
      aw->lmax = 1;
    }

  if (((gint) l) <= 1)
    return;

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
              gdouble b = (d1_0_0 - ((j * k) / (i * (i - 1.))));
              gdouble c = sqrt ((sq_i_m_1 - sq_j) * (sq_i_m_1 - sq_k)) /
                ((i - 1.) * two_i_m_1);

              aw->l_terms[i][j][i + k] =
                a * (b * aw->l_terms[i - 1][j][(i - 1) + k] -
                     c * aw->l_terms[i - 2][j][(i - 2) + k]);
            }
        }

      /* last two diagonal terms */
      aw->l_terms[i][i][i + i] = d1_1_1 *
        aw->l_terms[i - 1][i - 1][(i - 1) + (i - 1)];

      aw->l_terms[i][i - 1][i + i - 1] = (i * d1_0_0 - i + 1.) *
        aw->l_terms[i - 1][i - 1][(i - 1) + (i - 1)];

      /* last column terms */
      for (k = i; k >= 1 - i; k--)
        {
          aw->l_terms[i][i][i + k - 1] = -sqrt ((i + k) / (i - k + 1.)) * tb2 *
            aw->l_terms[i][i][i + k];
        }

      /* penultimate column */
      for (k = i - 1; k > 0; k--)
        {
          gdouble i_cb = i * cb;
          gdouble a = ((i_cb - k + 1.) / (i_cb - k));

          aw->l_terms[i][i - 1][i + k - 1] =
            -a * sqrt ((i + k) / (i - k + 1.)) *
            tb2 * aw->l_terms[i][i - 1][i + k];
        }

      /* we cut [i][i - 1] k loop in order to avoid problems for k==0 when
       * beta == pi/2.
       */
      aw->l_terms[i][i - 1][i + 0] = -sqrt (2. * (two_i_m_1) / (i - 1.)) *
        cb2 * sb2 * aw->l_terms[i - 1][i - 2][(i - 1) + 0];

      /* remaining of penultimate column */
      for (k = 0; k >= 2 - i; k--)
        {
          gdouble i_cb = i * cb;
          gdouble a = ((i_cb - k + 1.) / (i_cb - k));

          aw->l_terms[i][i - 1][i + k - 1] =
            -a * sqrt ((i + k) / (i - k + 1.)) *
            tb2 * aw->l_terms[i][i - 1][i + k];
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

  aw->lmax = l;
}

/**
 * aran_wigner_term:
 * @aw: an #AranWigner structure.
 * @l: a #guint.
 * @m1: a #guint. Condition @m1 <= @l must hold.
 * @m2: a #gint. Condition -@l <= @m2 <= @l must hold.
 *
 * Provides access to the corresponding term in @aw.
 *
 * Returns: address of the term.
 */
gdouble *aran_wigner_term (AranWigner * aw, guint l, guint m1, gint m2)
{
  g_return_val_if_fail (aw != NULL, NULL);
  return aw->l_terms[l][m1] + l + m2;
}

/**
 * aran_wigner_new:
 * @beta: Rotation angle (-pi < @beta < pi).
 * @l: initial degree preallocation if positive.
 *
 * Creates a new #AranWigner for angle @beta. 
 *
 * Returns: A Wigner D-function coefficients buffer.
 */
AranWigner *aran_wigner_new (gdouble beta, gint l)
{
  AranWigner *aw = g_malloc (sizeof (AranWigner));

  aw->beta = beta;

  aw->lmax = -1;

  aw->l_m1_m2_terms = NULL;

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

  if (aw->l_m1_m2_terms != NULL)
    g_free (aw->l_m1_m2_terms);

  g_free (aw);
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
              g_fprintf (file, "[%d,%d,%d]=%g ", i, j, k,
                         *aran_wigner_term (aw, i, j, k));
            }
        }
    }

  g_fprintf (file, "]");

}
