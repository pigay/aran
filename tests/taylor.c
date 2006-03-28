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

#include <stdlib.h>

#include <complex.h>

#include <math.h>

#include "aran/aran.h"
#include "aran/aranlaurentseriesd.h"

#define N 10

static const guint order = 20;

static gdouble tolerance = 1.E-7;
static gboolean verbose = FALSE;

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
	    tolerance = tmp;
	  else
	    g_printerr ("Invalid error limit value (-err %s)\n", arg);
	}
      else if (g_ascii_strcasecmp (arg, "-v") == 0)
	{
	  verbose = TRUE;
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

static void err (gcomplex128 z, gcomplex128 ref, gcomplex128 res)
{
  gcomplex128 diff = (ref-res);

  if (ref != 0.) diff /= ref;

  if (verbose)
    g_printerr ("testing (%f,%f) : %e\n", creal (z), cimag (z), cabs (diff));

  if (cabs (diff) > tolerance)
    g_printerr (" Error (%f,%f) : (%e,%e) != (%e,%e) %e\n",
		creal (z), cimag (z),
		creal (ref), cimag (ref),
		creal (res), cimag (res),
		cabs (diff));

}

static void check (AranLaurentSeriesd *als,
		   gcomplex128 (*func) (gcomplex128 z),
		   gcomplex128 center, gdouble radius)
{
  gcomplex128 z, res, ref;
  guint i, j;

  z = 0.;
  ref = func (z+center);
  res = aran_laurent_seriesd_evaluate (als, z);
  err (z+center, ref, res);

  for (i=0; i<=N; i ++)
    {
      gdouble r = i * radius / N;

      for (j=0; j<i; j ++)
	{
	  gdouble arg = 2. * j * G_PI / i;

	  z = r*cexp (I*arg);

	  ref = func (z+center);
	  res = aran_laurent_seriesd_evaluate (als, z);

	  err (z+center, ref, res);
	}
    }
}

/* FIXME: check overflows if required order becomes too high. */
/* or sould we use gdouble variable? */
static guint64 facto (guint64 n)
{
  guint64 i, facto = 1;

  for (i=2; i<=n; i ++)
    {
      facto *= i;
    }

  return facto;
}

static void build_cexp_series (AranLaurentSeriesd *als)
{
  guint8 i;
  for (i=0; i<=aran_laurent_seriesd_get_posdeg (als); i ++)
    {
      gcomplex128 *term = aran_laurent_seriesd_get_term (als, i);

      *term = 1./facto (i);
    }
}

int main (int argc, char **argv)
{
  int ret = 0;
  AranLaurentSeriesd *als, *alst;

  aran_init();

  parse_args (argc, argv);

  /* create a taylor series. */
  als = aran_laurent_seriesd_new (order, 0);

  build_cexp_series (als);
  check (als, cexp, 0., 1.);

  /* translate taylor series inside its convergence radius */
  alst = aran_laurent_seriesd_new (order, 0);

  aran_laurent_seriesd_translate (als, 0., alst, I*.5);
  check (alst, cexp, I*.5, .5); /* convergence radius halved */

  aran_laurent_seriesd_free (als);
  aran_laurent_seriesd_free (alst);

  return ret;
}

