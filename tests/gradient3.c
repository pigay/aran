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

static gdouble epsilon = 1.e-9;

static void local_to_cartesian (gdouble r,
                                gdouble cost, gdouble sint,
                                gdouble cosp, gdouble sinp,
                                gdouble dr, gdouble dt, gdouble dp,
                                VsgVector3d *grad)
{
  grad->x = sint*cosp*dr + cost*cosp*dt - sinp*dp;
  grad->y = sint*sinp*dr + cost*sinp*dt + cosp*dp;
  grad->z = cost*dr - sint*dt;
}

static void check (const gchar *log,
                   AranSphericalSeriesd *ass, gdouble radius,
		   VsgVector3d *center,
		   void (*f) (VsgVector3d *x, VsgVector3d *grad))
{
  guint i, j, k;
  gdouble t, p;
  gdouble err;
  gdouble r, cost, sint, cosp, sinp;
  gcomplex128 dr, dt, dp;

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
	      VsgVector3d vref;
	      VsgVector3d vres;
	      VsgVector3d verr;

	      r = radius * k/(N-1.);

	      vsg_vector3d_from_spherical_internal (&vec, r,
						    cost, sint,
						    cosp, sinp);

/*               aran_spherical_seriesd_local_gradient_evaluate (ass, &vec, */
/*                                                               &vres); */

              aran_spherical_seriesd_local_gradient_evaluate_internal (ass,
                                                                       r,
                                                                       cost, sint,
                                                                       cosp, sinp,
                                                                       &dr, &dt, &dp);
              local_to_cartesian (r, cost, sint, cosp, sinp,
                                  creal (dr), creal (dt), creal (dp),
                                  &vres);

              
              vsg_vector3d_scalp (&vres, -1., &vres);

	      vsg_vector3d_add (&vec, center, &vec);

	      f (&vec, &vref);

	      vsg_vector3d_sub (&vref, &vres, &verr);

              err = vsg_vector3d_norm (&verr) / vsg_vector3d_norm (&vref);

	      if (fabs (err) > epsilon)
		{
		  g_printerr ("Error %s (%f,%f,%f) : " \
                              "(%f,%f,%f) != (%f,%f,%f) -> %e\n",
                              log,
			      r, t, p,
/* 			      vec.x, vec.y, vec.z, */
			      vref.x, vref.y, vref.z,
			      vres.x, vres.y, vres.z,
			      fabs (err));
		}
	    }
	}
    }
}



static VsgVector3d p;

static AranSphericalSeriesd *create_newton (guint deg, VsgVector3d *center)
{
  AranSphericalSeriesd *ass = aran_spherical_seriesd_new (deg, 0);
  gint l, m;
  gcomplex128 harmonics[((deg+1)*(deg+2))/2];
  gdouble r, cost, sint, cosp, sinp;
  gcomplex128 expp;
  gdouble fact, inv_r;
  VsgVector3d tmp;

  vsg_vector3d_sub (&p, center, &tmp);

  vsg_vector3d_to_spherical_internal (&tmp, &r, &cost, &sint, &cosp, &sinp);
  expp = cosp + G_I * sinp;

  aran_spherical_harmonic_evaluate_multiple_internal (deg, cost, sint, expp,
						      harmonics);


  *aran_spherical_seriesd_get_term (ass, 0, 0) = 0.;

  inv_r = 1./ r;
  fact = inv_r;

  for (l=0; l<=deg; l ++)
    {
      gcomplex128 term;

      term = fact * (4.*G_PI / (l+l+1.)) *
	*aran_spherical_harmonic_multiple_get_term (l, 0, harmonics);

      *aran_spherical_seriesd_get_term (ass, l, 0) = conj (term);

      for (m=1; m<=l; m ++)
	{
	  term = fact * (4.*G_PI / (l+l+1.)) *
	    *aran_spherical_harmonic_multiple_get_term (l, m, harmonics);

	    *aran_spherical_seriesd_get_term (ass, l, m) = conj (term);

	}
      fact *= inv_r;
    }

  return ass;
}

void newtongrad (VsgVector3d *vec, VsgVector3d *grad)
{
  gdouble r;

  vsg_vector3d_sub (vec, &p, grad);

  r = vsg_vector3d_norm (grad);

  vsg_vector3d_scalp (grad, 1./ (r*r*r), grad);
}

//static VsgVector3d zero = {0., 0., 0.};
static VsgVector3d tr = {4.5, 0., 0.};
static VsgVector3d center = {5., 0., 0.};

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
  AranSphericalSeriesd *ass;
  AranSphericalSeriesd *ast;
  AranSphericalSeriesd *ast2;

  parse_args (argc, argv);

  p.x = 3.;
  p.y = -1.;
  p.z = -1.;

  ass = create_newton (deg, &center);

/*   aran_spherical_seriesd_write (ass, stderr); */
/*   g_printerr ("\n\n"); */
  check ("devel", ass, 1., &center, newtongrad);

  ast = aran_spherical_seriesd_clone (ass);
  aran_spherical_seriesd_set_zero (ast);

  aran_spherical_seriesd_translate (ass, &center, ast, &tr);

/*   aran_spherical_seriesd_write (ast, stderr); */
/*   g_printerr ("\n\n"); */
  check ("translated", ast, 0.5, &tr, newtongrad);

  ast2 = aran_spherical_seriesd_clone (ass);
  aran_spherical_seriesd_set_zero (ast2);

  aran_spherical_seriesd_translate_kkylin (ass, &center, ast2, &tr);

/*   aran_spherical_seriesd_write (ast2, stderr); */
/*   g_printerr ("\n\n"); */
  check ("translated kkylin", ast2, 0.5, &tr, newtongrad);

  aran_spherical_seriesd_free (ass);
  aran_spherical_seriesd_free (ast);
  aran_spherical_seriesd_free (ast2);

  return ret;
}

