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

#include "aran/aranlaurentseriesd.h"
#include "aran/arandevelopment2d.h"

#define N 10

static gdouble epsilon = 1.e-9;

typedef struct _Particle Particle;

struct _Particle
{
  VsgVector2d vector;

  gdouble charge;

  gcomplex128 field;
};

/* check Local expansion against Newton gradient for several points inside a ball */
static void check (const gchar *log,
                   AranDevelopment2d *dev, gdouble radius,
		   VsgVector2d *center,
		   gcomplex128 (*f) (gcomplex128 x))
{
  guint j, k;
  gdouble t;
  gdouble err;
  gdouble r, cost, sint;

  for (j=0; j<N; j ++)
    {
      /* iterate through Teta interval [0, pi) */
      t = G_PI * j/(N-1.);
      cost = cos (t);
      sint = sin (t);

      for (k=0; k<N; k ++)
        {
          /* iterate through r interval [0, radius] */
          gcomplex128 vec;
          gcomplex128 vref;
          gcomplex128 vres;
          gcomplex128 verr;

          r = radius * k/(N-1.);

          vec = r * (cost + G_I * sint);
          vres = aran_laurent_seriesd_evaluate (dev->local, vec);

          vres = -vres;

          vec += (center->x + G_I * center->y);

          vref = f (vec);

          verr = vref - vres;

          err = cabs(verr) / cabs(vref);

          if (fabs (err) > epsilon || !finite (err))
            {
              g_printerr ("Error %s (r=%f, t=%f) : "              \
                          "ref(%f,%f) != res(%f,%f) -> %e\n",
                          log,
                          r, t,
                          creal(vref), cimag(vref),
                          creal(vres), cimag(vres),
                          fabs (err));
            }
        }
    }
}



Particle particle;

/* Newton force field evaluation at position vec */
gcomplex128 newtongrad (gcomplex128 vec)
{
  return 1. / (vec - (particle.vector.x + G_I * particle.vector.y));
}

static guint deg = 29;

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
      else if (g_ascii_strcasecmp (arg, "-pr") == 0)
	{
	  guint tmp = 0;
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (sscanf (arg, "%u", &tmp) == 1 && tmp > 0)
            deg = tmp;
	  else
	    g_printerr ("Invalid precision order value (-pr %s)\n", arg);
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
  AranDevelopment2d *dev;
  VsgPRTree2dNodeInfo node_info;

  parse_args (argc, argv);

  particle.vector.x = 3.;
  particle.vector.y = -1.;

  particle.charge = 1.;

  node_info.center.x = 5.;
  node_info.center.y = 0.;

  dev = aran_development2d_new (deg, 0);
  aran_development2d_p2l (&particle.vector, particle.charge, &node_info, dev);

  aran_laurent_seriesd_write (dev->local, stderr);
  g_printerr ("\n\n");
  check ("devel", dev, 1., &node_info.center, newtongrad);

  aran_development2d_free (dev);

  return ret;
}

