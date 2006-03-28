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
#include "aran/aransphericalseriesd.h"

#define N 10

static gdouble epsilon = 1.e-7;

static void check (const gchar *log,
                   AranSphericalSeriesd *ass, gdouble radius,
		   VsgVector3d *center,
		   gcomplex128 (*f) (VsgVector3d *x))
{
  guint i, j, k;
  gdouble t, p;
  gcomplex128 yres, yref, err;
  gdouble r, cost, sint, cosp, sinp;

  for (i=0; i<N ; i ++)
    {
      p = 2.*G_PI * i/(N-1.);
      cosp = cos (p);
      sinp = sin (p);

      for (j=0; j<N; j ++)
	{
	  t = G_PI * j/(N-1.);
	  cost = cos (t);
	  sint = sin (t);

	  for (k=0; k<N; k ++)
	    {
	      VsgVector3d vec;
	      r = radius + radius * k/(N-1.);

	      yres = aran_spherical_seriesd_evaluate_internal (ass, r,
							       cost, sint,
							       cosp, sinp);

	      vsg_vector3d_from_spherical_internal (&vec, r,
						    cost, sint,
						    cosp, sinp);
	      vsg_vector3d_add (&vec, center, &vec);

	      yref = f (&vec);

	      err = yref-yres;

	      if (yref != 0.) err /= yref;

	      if (cabs (err) > epsilon && cabs(yref) > epsilon)
		{
		  g_printerr ("Error %s (%f,%f,%f) : " \
                              "(%e+%ej), (%e+%ej) -> %e\n",
                              log,
			      r, t, p,
			      creal (yref), cimag (yref),
			      creal (yres), cimag (yres),
			      cabs (err));
		}
	    }
	}
    }
}



static VsgVector3d p;

static AranSphericalSeriesd *create_newton (guint deg)
{
  AranSphericalSeriesd *ass = aran_spherical_seriesd_new (0, deg);
  gint l, m;
  gcomplex128 harmonics[((deg+1)*(deg+2))/2];
  gdouble r, cost, sint, cosp, sinp;
  gcomplex128 expp;
  gdouble fact;

  vsg_vector3d_to_spherical_internal (&p, &r, &cost, &sint, &cosp, &sinp);
  expp = cosp + G_I * sinp;

  aran_spherical_harmonic_evaluate_multiple_internal (deg, cost, sint, expp,
						      harmonics);


  *aran_spherical_seriesd_get_term (ass, 0, 0) = 0.;

  fact = 1.;

  for (l=0; l<deg; l ++)
    {
      gcomplex128 term;

      term = fact * (4.*G_PI / (l+l+1.)) *
	*aran_spherical_harmonic_multiple_get_term (l, 0, harmonics);

      *aran_spherical_seriesd_get_term (ass, -l-1, 0) = conj (term);

      for (m=1; m<=l; m ++)
	{
	  term = fact * (4.*G_PI / (l+l+1.)) *
	    *aran_spherical_harmonic_multiple_get_term (l, m, harmonics);

	    *aran_spherical_seriesd_get_term (ass, -l-1, m) = conj (term);

	}
      fact *= r;
    }

  return ass;
}

gcomplex128 newtonpot (VsgVector3d *vec)
{
  VsgVector3d tmp;

  vsg_vector3d_sub (vec, &p, &tmp);

  return 1. / vsg_vector3d_norm (&tmp);
}

static VsgVector3d zero = {0., 0., 0.};
static VsgVector3d tr = {1., 1.1, 1.2};

static guint deg = 30;

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
  AranSphericalSeriesd *ass;
  AranSphericalSeriesd *ast;
  AranSphericalSeriesd *ast2;

  parse_args (argc, argv);

  p.x = 0.1;
  p.y = 0.2;
  p.z = 0.3;

  ass = create_newton (deg);

/*   aran_spherical_seriesd_write (ass, stderr); */
/*   g_printerr ("\n"); */
  check ("devel", ass, 5., &zero, newtonpot);

  ast = aran_spherical_seriesd_clone (ass);

  aran_spherical_seriesd_set_zero (ast);
  aran_spherical_seriesd_translate (ass, &zero, ast, &tr);

/*   aran_spherical_seriesd_write (ast, stderr); */
/*   g_printerr ("\n"); */
  check ("translated", ast, 5., &tr, newtonpot);

  ast2 = aran_spherical_seriesd_clone (ass);

  aran_spherical_seriesd_set_zero (ast2);
  aran_spherical_seriesd_translate_kkylin (ass, &zero, ast2, &tr);

/*   aran_spherical_seriesd_write (ast2, stderr); */
/*   g_printerr ("\n"); */
  check ("translated kkylin", ast2, 5., &tr, newtonpot);

  aran_spherical_seriesd_set_zero (ass);
  aran_spherical_seriesd_translate (ast, &tr, ass, &zero);

/*   aran_spherical_seriesd_write (ass, stderr); */
/*   g_printerr ("\n"); */
  check ("translated back", ass, 5., &zero, newtonpot);

  aran_spherical_seriesd_free (ass);
  aran_spherical_seriesd_free (ast);
  aran_spherical_seriesd_free (ast2);

  return ret;
}

