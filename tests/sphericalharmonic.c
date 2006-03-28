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

#include "aran-config.h"

#include <stdio.h>

#include <stdlib.h>

#include <math.h>

#include "aran/aransphericalharmonic.h"

#define N 30

static gdouble epsilon = 1.e-10;

static void check (guint order, guint diff,
		   gcomplex128 (*f) (gdouble t, gdouble p))
{
  guint i, j;
  gdouble t, p;
  gcomplex128 yres, yref, err;

  for (i=0; i<N ; i ++)
    {
      p = 2.*G_PI * i/(N-1.);

      for (j=0; j<N; j ++)
	{
	  t = G_PI * j/(N-1.);

	  yres = aran_spherical_harmonic_evaluate (order, diff, t, p);
	  yref = f (t, p);

	  err = yref-yres;

	  if (yref != 0.) err /= yref;

	  if (cabs (err) > epsilon && cabs(yref) > epsilon)
	    {
	      g_printerr ("Y_%u^%d (%f,%f) : (%e,%e), (%e,%e) -> %e\n",
			  order, diff,
			  t, p,
			  creal (yref), cimag (yref),
			  creal (yres), cimag (yres),
			  cabs (err));
	    }
	}
    }
}

static void multiple_check (guint order)
{
  guint i, j, k, n;
  gdouble t, p;
  gcomplex128 yres, yref, err;
  gcomplex128 m[((order+1)*(order+2))/2];
  gcomplex128 *ptr;

  for (i=0; i<N ; i ++)
    {
      p = 2.*G_PI * i/(N-1.);

      for (j=0; j<N; j ++)
	{
	  t = G_PI * j/(N-1.);

	  aran_spherical_harmonic_evaluate_multiple (order, t, p, m);

	  for (k=0; k<=order; k ++)
	    {
	      for (n=0; n<=k; n ++)
		{
		  yref = aran_spherical_harmonic_evaluate (k, n, t, p);
		  ptr = aran_spherical_harmonic_multiple_get_term (k, n, m);
		  yres = *ptr;

		  err = yref-yres;

		  if (yref != 0.) err /= yref;

		  if (cabs (err) > epsilon && cabs(yref) > epsilon)
		    {
		      g_printerr ("Y_%u^%d (%f,%f) multiple: "
				  "(%e,%e), (%e,%e) -> %e\n",
				  k, n,
				  t, p,
				  creal (yref), cimag (yref),
				  creal (yres), cimag (yres),
				  cabs (err));
		    }
		}
	    }
	}
    }
}
/*
 * Spherical harmonics formulas taken from:
 * http://mathworld.wolfram.com/SphericalHarmonic.html
 */

static gcomplex128 Y_0_0 (gdouble t, gdouble p)
{
  return 0.5 / sqrt (G_PI);
}

static gcomplex128 Y_1_m1 (gdouble t, gdouble p)
{
  return 0.5 * sqrt (1.5 / G_PI) * sin (t) * cexp (- G_I * p);
}

static gcomplex128 Y_1_0 (gdouble t, gdouble p)
{
  return 0.5 * sqrt (3. / G_PI) * cos (t);
}

static gcomplex128 Y_1_1 (gdouble t, gdouble p)
{
  return - 0.5 * sqrt (1.5 / G_PI) * sin (t) * cexp (G_I * p);
}

static gcomplex128 Y_2_m2 (gdouble t, gdouble p)
{
  gdouble sint = sin (t);
  return 0.25 * sqrt (7.5 / G_PI) * sint*sint * cexp (- 2. * G_I * p);
}

static gcomplex128 Y_2_m1 (gdouble t, gdouble p)
{
  return 0.5 * sqrt (7.5 / G_PI) * sin (t) * cos (t) * cexp (- G_I * p);
}

static gcomplex128 Y_2_0 (gdouble t, gdouble p)
{
  gdouble cost = cos (t);
  return 0.25 * sqrt (5. / G_PI) * (3.*cost*cost - 1.);
}

static gcomplex128 Y_2_1 (gdouble t, gdouble p)
{
  return - 0.5 * sqrt (7.5 / G_PI) * sin (t) * cos (t) * cexp (G_I * p);
}

static gcomplex128 Y_2_2 (gdouble t, gdouble p)
{
  gdouble sint = sin (t);
  return 0.25 * sqrt (7.5 / G_PI) * sint*sint * cexp (2. * G_I * p);
}

static gcomplex128 Y_3_m3 (gdouble t, gdouble p)
{
  gdouble sint = sin (t);
  return 0.125 * sqrt (35. / G_PI) * sint*sint*sint * cexp (- 3. * G_I * p);
}

static gcomplex128 Y_3_m2 (gdouble t, gdouble p)
{
  gdouble sint = sin (t);
  return 0.25 * sqrt (52.5 / G_PI) * sint*sint * cos (t) *
    cexp (- 2. * G_I * p);
}

static gcomplex128 Y_3_m1 (gdouble t, gdouble p)
{
  gdouble cost = cos (t);
  return 0.125 * sqrt (21. / G_PI) *sin (t) * (5. * cost*cost - 1.) *
    cexp (- G_I * p);
}

static gcomplex128 Y_3_0 (gdouble t, gdouble p)
{
  gdouble cost = cos (t);
  return 0.25 * sqrt (7. / G_PI) * (5. * cost*cost - 3.) * cost;
}

static gcomplex128 Y_3_1 (gdouble t, gdouble p)
{
  gdouble cost = cos (t);
  return - 0.125 * sqrt (21. / G_PI) *sin (t) * (5. * cost*cost - 1.) *
    cexp (G_I * p);
}

static gcomplex128 Y_3_2 (gdouble t, gdouble p)
{
  gdouble sint = sin (t);
  return 0.25 * sqrt (52.5 / G_PI) * sint*sint * cos (t) *
    cexp (2. * G_I * p);
}

static gcomplex128 Y_3_3 (gdouble t, gdouble p)
{
  gdouble sint = sin (t);
  return - 0.125 * sqrt (35. / G_PI) * sint*sint*sint * cexp (3. * G_I * p);
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

  check (0, 0, Y_0_0);

  check (1, -1, Y_1_m1);
  check (1, 0, Y_1_0);
  check (1, 1, Y_1_1);

  check (2, -2, Y_2_m2);
  check (2, -1, Y_2_m1);
  check (2, 0, Y_2_0);
  check (2, 1, Y_2_1);
  check (2, 2, Y_2_2);

  check (3, -3, Y_3_m3);
  check (3, -2, Y_3_m2);
  check (3, -1, Y_3_m1);
  check (3, 0, Y_3_0);
  check (3, 1, Y_3_1);
  check (3, 2, Y_3_2);
  check (3, 3, Y_3_3);

  multiple_check (50);

  return ret;
}

