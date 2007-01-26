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

#include "aran/aransphericalharmonic.h"

#define N 30

static gdouble epsilon = 1.e-10;

static void multiple_check (guint order)
{
  guint i, j, k, n;
  gdouble t, p;
  gcomplex128 yres, yref, err;
  gcomplex128 m[((order+1)*(order+2))/2];
  gcomplex128 spec[((order+1)*(order+2))/2];

  for (i=0; i<N ; i ++)
    {
      p = 2.*G_PI * i/(N-1.);

      for (j=1; j<N-1; j ++)
	{
	  t = G_PI * j/(N-1.);

	  aran_spherical_harmonic_pre_gradient_multiple (order, t, p, m, spec);

	  for (k=0; k<=order; k ++)
	    {
	      for (n=0; n<=k; n ++)
		{
		  yref = *aran_spherical_harmonic_multiple_get_term (k, n, m) /
                    sin (t);
		  yres = *aran_spherical_harmonic_multiple_get_term (k, n,
                                                                     spec);

		  err = yref-yres;

		  if (yref != 0.) err /= yref;

		  if ((cabs (err) > epsilon && cabs(yref) > epsilon) ||
		      !finite (cabs (err)))
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


  multiple_check (50);

  return ret;
}

