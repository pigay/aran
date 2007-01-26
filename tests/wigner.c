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

#include "aran-config.h"

#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include "aran/aranwigner.h"
#include "aran/aran.h"
#include "aran/arancoefficientbufferd.h"

#define PHASE(m) (((m)%2 == 0) ? 1. : -1.)

static gdouble epsilon = 1.e-6;
static guint L = 20;
static guint N = 30;

static void parse_args (int argc, char **argv)
{
  int iarg = 1;
  char *arg;

  while (iarg < argc)
    {
      arg = argv[iarg];

      if (g_ascii_strcasecmp (arg, "-err") == 0)
        {
          gdouble tmp = 0;
          iarg++;

          arg = (iarg < argc) ? argv[iarg] : NULL;

          if (sscanf (arg, "%lf", &tmp) == 1 && tmp > 0.)
            epsilon = tmp;
          else
            g_printerr ("Invalid error limit value (-err %s)\n", arg);
        }
      else if (g_ascii_strcasecmp (arg, "--version") == 0)
        {
          g_printerr ("%s version %s\n", argv[0], PACKAGE_VERSION);
          exit (0);
        }
      else
        {
          g_printerr ("Invalid argument \"%s\"\n", arg);
        }

      iarg++;
    }
}

gdouble _lnfact_func (guint n, AranCoefficientBufferd * acb)
{
  gdouble ret;

  if (n <= 1)
    return 0.;

  ret = *aran_coefficient_bufferd_get_unsafe (acb, n - 1) + log (n);

  return ret;
}

static AranCoefficientBufferd *_lnfact = NULL;

static gdouble lnfact (guint n)
{
  return *aran_coefficient_bufferd_get_unsafe (_lnfact, n);
}

static gdouble wigner (guint l, gint m1, gint m2, gdouble beta)
{
  gdouble ret = 0.;
  gint l_m_m1, l_p_m1, l_m_m2, l_p_m2, m1_m_m2;
  gint k, kmin, kmax;
  gdouble cb2, sb2, a;

  g_assert (ABS (m1) <= ((gint) l));
  g_assert (ABS (m2) <= ((gint) l));

  l_m_m1 = l - m1;
  l_p_m1 = l + m1;
  l_m_m2 = l - m2;
  l_p_m2 = l + m2;
  m1_m_m2 = m1 - m2;

  kmin = MAX (0, -m1_m_m2);
  kmax = MIN (l_m_m1, l_p_m2);

  if (kmin > kmax)
    return ret;

  aran_coefficient_bufferd_require (_lnfact, 2*l);

  a = (lnfact (l_p_m1) + lnfact (l_m_m1) + lnfact (l_p_m2) +
       lnfact (l_m_m2)) * 0.5;

  cb2 = cos (beta * 0.5);
  sb2 = sin (beta * 0.5);

  for (k = kmin; k <= kmax; k++)
    {
      gdouble b = lnfact (k) + lnfact (l_m_m1 - k) + lnfact (l_p_m2 - k) +
        lnfact (m1_m_m2 + k);
      gdouble c = pow (cb2, l_m_m1 + l_p_m2 - 2 * k);
      gdouble d = pow (sb2, m1_m_m2 + 2 * k);

      ret += PHASE (m1_m_m2 + k) * exp (a - b) * c * d;
    }

  return ret;
}

static gint check (gdouble beta)
{
  AranWigner *aw;
  gint l, m1;
  gint m2;
  gint faults = 0;

  aw = aran_wigner_new (beta, L);

  for (l = 0; l <= L; l++)
    {
      for (m1 = 0; m1 <= l; m1++)
        {
          for (m2 = -l; m2 <= l; m2++)
            {
              gdouble res = *aran_wigner_term (aw, l, m1, m2);
              gdouble ref = wigner (l, m1, m2, beta);
              gdouble err;

              err = fabs (res - ref);
              if (fabs (ref) > epsilon)
                err /= fabs (ref);

              if (err >= epsilon || !finite (err))
                {
                  faults ++;
                  g_printerr ("Error wigner (l=%d, m1=%d, m2=%d, beta=%g) "\
                              "%g != %g, err=%g\n",
                              l, m1, m2, beta, res, ref, err);
                }
            }
        }
    }

  aran_wigner_free (aw);

  return faults;
}

int main (int argc, char **argv)
{
  gint i;
  gint ret = 0;

  aran_init ();

  parse_args (argc, argv);

  _lnfact = aran_coefficient_bufferd_new (_lnfact_func, L);

  /* general testing */
  for (i=0; i<N; i ++)
    {
      ret += check (i * G_PI / N);
    }

  /* additional values */
  ret += check (G_PI / 2.);
  ret += check (0.84106867056793033);

  aran_coefficient_bufferd_free (_lnfact);

  return ret;
}
