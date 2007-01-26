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

#include <stdlib.h>

#include <complex.h>

#include <math.h>

#include "aran/aran.h"
#include "aran/arandevelopment3d.h"


typedef struct _PointAccum PointAccum;

struct _PointAccum
{
  VsgVector3d vector;

  gdouble density;

  gcomplex128 accum;

  guint id;
};



/* static void p2p (PointAccum *one, PointAccum *other) */
/* { */
/*   if (one != other) */
/*     { */
/*       VsgVector3d tmp; */
/*       gdouble inv_r; */

/*       /\* destination - source *\/ */
/*       vsg_vector3d_sub (&one->vector, &other->vector, &tmp); */

/*       inv_r = 1. / vsg_vector3d_norm (&tmp); */

/*       /\* G(xd, xs) = 1./(xd-xs) * src_density *\/ */

/*       one->accum += inv_r * other->density; */
/*       other->accum += inv_r * one->density; */
/*     } */
/* } */

static gdouble kernel (PointAccum *one, PointAccum *other)
{
  VsgVector3d tmp;
  gdouble inv_r;

  /* destination - source */
  vsg_vector3d_sub (&one->vector, &other->vector, &tmp);

  inv_r = 1. / vsg_vector3d_norm (&tmp);

  return inv_r * one->density;
}

static void p2m (PointAccum *particle, const VsgVector3d *center,
		 AranDevelopment3d *devel)
{
  VsgVector3d tmp;
  guint deg = aran_spherical_seriesd_get_negdeg (devel->multipole);
  gint l, m;
  gcomplex128 harmonics[((deg+1)*(deg+2))/2];
  gdouble r, cost, sint, cosp, sinp;
  gcomplex128 expp;
  gdouble fact;

  vsg_vector3d_sub (&particle->vector, center, &tmp);

  vsg_vector3d_to_spherical_internal (&tmp, &r, &cost, &sint, &cosp, &sinp);
  expp = cosp + G_I * sinp;

  aran_spherical_harmonic_evaluate_multiple_internal (deg, cost, sint, expp,
						      harmonics);


  *aran_spherical_seriesd_get_term (devel->multipole, 0, 0) += 0.;

  fact = particle->density;

  for (l=0; l<deg; l ++)
    {
      gcomplex128 *ptr;
      gcomplex128 term;

      term = fact * (4.*G_PI / (l+l+1.)) *
	*aran_spherical_harmonic_multiple_get_term (l, 0, harmonics);

      ptr = aran_spherical_seriesd_get_term (devel->multipole, -l-1, 0);
      ptr[0] += conj (term);

      for (m=1; m<=l; m ++)
	{

	  term = fact * (4.*G_PI / (l+l+1.)) *
	    *aran_spherical_harmonic_multiple_get_term (l, m, harmonics);

          ptr[m] += conj (term);
	}
      fact *= r;
    }

}

static gdouble epsilon = 1.E-3;
static guint order = 24;

void parse_args (int argc, char **argv)
{
  int iarg = 1;
  char *arg;

  while (iarg < argc)
    {
      arg = argv[iarg];

      if (g_ascii_strcasecmp (arg, "-pr") == 0)
	{
	  guint tmp = 0;
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (sscanf (arg, "%u", &tmp) == 1 && tmp > 0)
	      order = tmp;
	  else
	    g_printerr ("Invalid precision order value (-pr %s)\n", arg);
	}
      else if (g_ascii_strcasecmp (arg, "-err") == 0)
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

PointAccum points[] = {
  {{1., 0., 0.},0.05, 0., 0},
  {{1., -.5, 0.},0.0, 0., 0},
};

VsgVector3d centers[] = {
  {0.75, 0.25, -0.25},
  {0.75, -0.75, -0.25},
};

int main (int argc, char **argv)
{
  AranDevelopment3d *one, *other;
  int ret = 0;
  gcomplex128 zres, zref;

  aran_init();

  parse_args (argc, argv);

  one = aran_development3d_new (0, order);
  other = aran_development3d_new (order, 0);
 
  aran_development3d_set_zero (one);
  aran_development3d_set_zero (other);

  p2m (&points[0], &centers[0], one);

  zres = aran_development3d_multipole_evaluate (&centers[0],
                                               one,
                                               &points[1].vector);

  zref = kernel (&points[0], &points[1]);

  if (cabs (zref-zres) > epsilon)
    {
      g_printerr ("Direct multipole error (%e + i%e) (%e + i%e)\n",
                  creal(zref), cimag (zref),
                  creal(zres), cimag (zres));
    }

  aran_development3d_m2l (&centers[0], one, &centers[1], other);

  zres = aran_development3d_local_evaluate (&centers[1],
                                            other,
                                            &points[1].vector);

  if (cabs (zref-zres) > epsilon)
    {
      g_printerr ("Multipole to local error (%e + i%e) (%e + i%e)\n",
                  creal(zref), cimag (zref),
                  creal(zres), cimag (zres));
    }

  aran_development3d_free (one);
  aran_development3d_free (other);

  return ret;
}

