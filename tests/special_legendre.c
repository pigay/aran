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

#include "aran/aranlegendre.h"

static const guint N = 100;

static gdouble epsilon = 1.e-5;

static void multiple_check (guint order)
{
  guint i, j, k;
  gdouble x, yres, yref, err;
  gdouble m[((order+1)*(order+2))/2];
  gdouble spec[((order+1)*(order+2))/2];

  for (i=1; i<N ; i++)
    {
      x = i*2./N - 1.;

      aran_legendre_associated_evaluate_multiple (order, x, m);
      aran_legendre_associated_evaluate_special_internal (order, x,
                                                          sqrt (1.-x*x),
                                                          m, spec);
      for (j=1; j<=order; j ++)
	{
	  for (k=1; k<=j; k ++)
	    {
	      yres = *aran_legendre_associated_multiple_get_term (j, k, spec);
	      yref = *aran_legendre_associated_multiple_get_term (j, k, m) /
                sqrt (1.-x*x);

	      err = yref-yres;

	      /* don't use relative error on too small values. */
	      if (fabs (yref) > 1.e-4) err /= yref;

	      /* avoid testing on high values because yref evaluation is
		 inaccurate by definition. */
	      if ((fabs (err) > epsilon && fabs (yref) < 1.e+3) ||
		!finite (err))
		{
		  g_printerr ("p_%u^%u (%f) : multiple %e, %e -> %e\n",
			      j, k, x, yref, yres, err);
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

  multiple_check (40);

  return ret;
}

