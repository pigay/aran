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

#include "aran/aranlegendre.h"

#define N 10

#include <math.h>

static gdouble epsilon = 1.e-10;

static void check (guint order, guint diff, gdouble (*f) (gdouble x))
{
  guint i;
  gdouble x, yres, yref, err;

  for (i=0; i<N ; i++)
    {
      x = 2.*i/(N-1.) - 1.;

      yres = aran_legendre_associated_evaluate (order, diff, x);
      yref = f (x);

      err = yref-yres;
      if (yref != 0.) err /= yref;

      if (fabs (err) > epsilon || !finite (err))
	{
	  g_printerr ("p_%u^%u (%f) : %e, %e -> %e\n",
		      order, diff, x, yref, yres, err);
	}
    }
}

/* similar recurrence function */
static gdouble Plm (gint l, gint m, gdouble x)
{
  if (l<0 || m>l) return 0.;
  if (l == 0) return 1.;
  if (l == m) return - (l+l-1) * sqrt (1.-x*x) * Plm (l-1, m-1, x);

  return ((l+l-1) * x * Plm (l-1, m, x) - (l+m-1) * Plm (l-2, m, x)) / (l-m);
}

static void multiple_check (guint order)
{
  guint i, j, k;
  gdouble x, yres, yref, err;
  gdouble m[((order+1)*(order+2))/2];
  gdouble *ptr;


  for (i=0; i<N ; i++)
    {
      x = 2.*i/(N-1.) - 1.;

      aran_legendre_associated_evaluate_multiple (order, x, m);

      for (j=0; j<=order; j ++)
	{
	  for (k=0; k<=j; k ++)
	    {
/* 	      yref = aran_legendre_associated_evaluate (j, k, x); */
	      yref = Plm (j, k, x);
	      ptr = aran_legendre_associated_multiple_get_term (j, k, m);
	      yres = *ptr;

	      err = yref-yres;
	      if (yref != 0.) err /= yref;

	      if ((fabs (err) > epsilon && yref > epsilon) || !finite (err))
		{
		  g_printerr ("p_%u^%u (%f) : multiple %e, %e -> %e\n",
			      j, k, x, yref, yres, err);
		}

	    }
	}
    }
}

/*
 * Associated Legendre polynomials formulas taken from:
 * http://mathworld.wolfram.com/LegendrePolynomial.html
 */

static gdouble p_0_0 (gdouble x)
{
  return 1;
}

static gdouble p_1_0 (gdouble x)
{
  return x;
}

static gdouble p_1_1 (gdouble x)
{
  return - sqrt (1. - x*x);
}

static gdouble p_2_0 (gdouble x)
{
  return 0.5 * (3.*x*x - 1.);
}

static gdouble p_2_1 (gdouble x)
{
  return - 3. * x * sqrt (1. - x*x);
}

static gdouble p_2_2 (gdouble x)
{
  return 3. * (1. - x*x);
}

static gdouble p_3_0 (gdouble x)
{
  return 0.5 * x * (5.*x*x - 3.);
}

static gdouble p_3_1 (gdouble x)
{
  return 1.5 * (1. - 5.*x*x) * sqrt (1. - x*x);
}

static gdouble p_3_2 (gdouble x)
{
  return 15. * x * (1. - x*x);
}

static gdouble p_3_3 (gdouble x)
{
  return - 15. * (1. - x*x) * sqrt (1. - x*x);
}

static gdouble p_4_0 (gdouble x)
{
  gdouble x2 = x*x;

  return 1./8. * ((35.*x2- 30.)*x2 + 3.);
}

static gdouble p_4_1 (gdouble x)
{
  gdouble x2 = x*x;

  return 2.5 * x * (3. - 7.*x2) * sqrt (1. - x2);
}

static gdouble p_4_2 (gdouble x)
{
  gdouble x2 = x*x;

  return 7.5 * (7.*x2 - 1.) * (1. - x2);
}

static gdouble p_4_3 (gdouble x)
{
  gdouble x2 = x*x;

  return - 105. * x * (1. - x2) * sqrt (1. - x2);
}

static gdouble p_4_4 (gdouble x)
{
  gdouble x2 = x*x;

  return 105. * (1. - x2) * (1. - x2);
}

static gdouble p_5_0 (gdouble x)
{
  gdouble x2 = x*x;

  return 1./8. * x * (63.*x2*x2 - 70.*x2 + 15.);
}

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
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

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

      iarg ++;
    }
}
int main (int argc, char **argv)
{
  int ret = 0;

  parse_args (argc, argv);

  check (0, 0, p_0_0);

  check (1, 0, p_1_0);
  check (1, 1, p_1_1);

  check (2, 0, p_2_0);
  check (2, 1, p_2_1);
  check (2, 2, p_2_2);

  check (3, 0, p_3_0);
  check (3, 1, p_3_1);
  check (3, 2, p_3_2);
  check (3, 3, p_3_3);

  check (4, 0, p_4_0);
  check (4, 1, p_4_1);
  check (4, 2, p_4_2);
  check (4, 3, p_4_3);
  check (4, 4, p_4_4);

  check (5, 0, p_5_0);

  /* check two similar formulas. maximum degree is low since higher values 
   * lead to numerical artifacts (that aren't truly bugs :). 
   */
  multiple_check (10);

  return ret;
}

