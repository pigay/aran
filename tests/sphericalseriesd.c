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

#include "aran/aran.h"
#include "aran/aranlegendre.h"
#include "aran/aransphericalharmonic.h"
#include "aran/aransphericalseriesd.h"

#define N 10

static gdouble epsilon = 1.e-10;

static void check (AranSphericalSeriesd *ass, gdouble radius,
		   gcomplex128 (*f) (gdouble r, gdouble cost, gdouble sint,
				     gdouble cosp, gdouble sinp))
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
	      r = radius + radius * k/(N-1.);

	      yres = aran_spherical_seriesd_evaluate_internal (ass, r,
							       cost, sint,
							       cosp, sinp);
	      yref = f (r, cost, sint, cosp, sinp);

	      err = yref-yres;

	      if (yref != 0.) err /= yref;

	      if ((cabs (err) > epsilon && cabs(yref) > epsilon) ||
		  !finite (cabs (err)))
		{
		  g_printerr ("Error (%f,%f,%f) : (%e,%e), (%e,%e) -> %e\n",
			      r, t, p,
			      creal (yref), cimag (yref),
			      creal (yres), cimag (yres),
			      cabs (err));
		}
	    }
	}
    }
}

/* static AranSphericalSeriesd *create_series1 () */
/* { */
/*   AranSphericalSeriesd *ass = aran_spherical_seriesd_new (2, 1); */

/*   *aran_spherical_seriesd_get_term (ass, 0, 0) = 1.; */

/*   return ass; */
/* } */

/* static gcomplex128 function1 (gdouble r, gdouble cost, gdouble sint, */
/* 			      gdouble cosp, gdouble sinp) */
/* { */
/*   gcomplex128 res = 0.; */
/*   gcomplex128 expp; */

/*   expp = cosp + G_I * sinp; */

/*   res += aran_spherical_harmonic_evaluate_internal (0, 0, cost, sint, 1.); */

/*   return res; */
/* } */

/* static AranSphericalSeriesd *create_series2 () */
/* { */
/*   AranSphericalSeriesd *ass = aran_spherical_seriesd_new (2, 1); */

/*   *aran_spherical_seriesd_get_term (ass, 0, 0) = 1.; */
/*   *aran_spherical_seriesd_get_term (ass, 1, 0) = 1.1; */

/*   return ass; */
/* } */

/* static gcomplex128 function2 (gdouble r, gdouble cost, gdouble sint, */
/* 			      gdouble cosp, gdouble sinp) */
/* { */
/*   gcomplex128 res = 0.; */
/*   gcomplex128 expp; */

/*   expp = cosp + G_I * sinp; */

/*   res += aran_spherical_harmonic_evaluate_internal (0, 0, cost, sint, 1.); */

/*   res += 1.1 * r * */
/*     aran_spherical_harmonic_evaluate_internal (1, 0, cost, sint, 1.); */

/*   return res; */
/* } */

/* static AranSphericalSeriesd *create_series3 () */
/* { */
/*   AranSphericalSeriesd *ass = aran_spherical_seriesd_new (2, 1); */

/*   *aran_spherical_seriesd_get_term (ass, 1, 1) = 1.2; */

/*   return ass; */
/* } */

/* static gcomplex128 function3 (gdouble r, gdouble cost, gdouble sint, */
/* 			      gdouble cosp, gdouble sinp) */
/* { */
/*   gcomplex128 res = 0.; */
/*   gcomplex128 expp; */

/*   expp = cosp + G_I * sinp; */

/*   res += 1.2 * r * */
/*     aran_spherical_harmonic_evaluate_internal (1, 1, cost, sint, expp); */

/*   return res; */
/* } */

/* static AranSphericalSeriesd *create_series4 () */
/* { */
/*   AranSphericalSeriesd *ass = aran_spherical_seriesd_new (2, 1); */

/*   *aran_spherical_seriesd_get_term (ass, 1, -1) = 1.3; */

/*   return ass; */
/* } */

/* static gcomplex128 function4 (gdouble r, gdouble cost, gdouble sint, */
/* 			      gdouble cosp, gdouble sinp) */
/* { */
/*   gcomplex128 res = 0.; */
/*   gcomplex128 expp; */

/*   expp = cosp + G_I * sinp; */

/*   res += 1.3 * r * */
/*     aran_spherical_harmonic_evaluate_internal (1, -1, cost, sint, expp); */

/*   return res; */
/* } */

/* static AranSphericalSeriesd *create_series5 () */
/* { */
/*   AranSphericalSeriesd *ass = aran_spherical_seriesd_new (0, 2); */

/*   *aran_spherical_seriesd_get_term (ass, -1, 0) = 1.3; */
/*   *aran_spherical_seriesd_get_term (ass, -2, -1) = 1.4; */
/*   *aran_spherical_seriesd_get_term (ass, -2, 1) = 1.5; */

/*   return ass; */
/* } */

/* static gcomplex128 function5 (gdouble r, gdouble cost, gdouble sint, */
/* 			      gdouble cosp, gdouble sinp) */
/* { */
/*   gcomplex128 res = 0.; */
/*   gcomplex128 expp; */

/*   expp = cosp + G_I * sinp; */

/*   res += 1.3 * 1./r * */
/*     aran_spherical_harmonic_evaluate_internal (0, 0, cost, sint, 1.); */

/*   res += 1.4 * 1./(r*r) * */
/*     aran_spherical_harmonic_evaluate_internal (1, -1, cost, sint, expp); */

/*   res += 1.5 * 1./(r*r) * */
/*     aran_spherical_harmonic_evaluate_internal (1, 1, cost, sint, expp); */

/*   return res; */
/* } */

static VsgVector3d p = {1., 0., 0.};

static AranSphericalSeriesd *legendre_addition (guint l)
{
  AranSphericalSeriesd *ass = aran_spherical_seriesd_new (l, 0);
  gcomplex128 harmonics[((l+1)*(l+2))/2];
  VsgVector3d tmp;
  gdouble r, cost, sint, cosp, sinp;
  gcomplex128 expp;
  gint m;

  r = vsg_vector3d_norm (&p);
  vsg_vector3d_scalp (&p, 1./r, &tmp);

  cost = tmp.z;
  sint = sqrt (tmp.x*tmp.x + tmp.y*tmp.y);

  cosp = tmp.x / sint;
  sinp = tmp.y / sint;

  expp = cosp + G_I * sinp;

  aran_spherical_harmonic_evaluate_multiple_internal (l, cost, sint, expp,
						      harmonics);

  aran_spherical_seriesd_set_zero (ass);

  for (m=0; m<=l; m ++)
    {
      gcomplex128 term = (4.*G_PI / (l+l+1)) *
	*aran_spherical_harmonic_multiple_get_term (l, m, harmonics);

      *aran_spherical_seriesd_get_term (ass, l, m) = conj (term);
     }

  return ass;
}

static gcomplex128 P_l (guint l, gdouble cost2, gdouble sint2,
			gdouble cosp2, gdouble sinp2)
{
  VsgVector3d tmp;
  gdouble r1, cost1, sint1, cosp1, sinp1;
  gdouble cosp1_m_p2;

  r1 = vsg_vector3d_norm (&p);
  vsg_vector3d_scalp (&p, 1./r1, &tmp);

  cost1 = tmp.z;
  sint1 = sqrt (tmp.x*tmp.x + tmp.y*tmp.y);

  cosp1 = tmp.x / sint1;
  sinp1 = tmp.y / sint1;

  cosp1_m_p2 = cosp1*cosp2 + sinp1*sinp2;

  return aran_legendre_evaluate (l, cost1*cost2 + sint1*sint2*cosp1_m_p2);
}

static gcomplex128 P_0 (gdouble r, gdouble cost, gdouble sint,
			gdouble cosp, gdouble sinp)
{
  return P_l (0, cost, sint, cosp, sinp);
}

static gcomplex128 P_1 (gdouble r, gdouble cost, gdouble sint,
			gdouble cosp, gdouble sinp)
{
  return P_l (1, cost, sint, cosp, sinp) * r;
}

static gcomplex128 P_2 (gdouble r, gdouble cost, gdouble sint,
			gdouble cosp, gdouble sinp)
{
  return P_l (2, cost, sint, cosp, sinp) * r*r;
}

static gcomplex128 P_3 (gdouble r, gdouble cost, gdouble sint,
			gdouble cosp, gdouble sinp)
{
  return P_l (3, cost, sint, cosp, sinp) * r*r*r;
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
  AranSphericalSeriesd *ass;
  AranSphericalSeriesd *ass2;

  aran_init ();

  parse_args (argc, argv);

  ass = aran_spherical_seriesd_new (2, 1);

  /* check basic allocation copy */
  *aran_spherical_seriesd_get_term (ass, -1, 0) = -1.;
  *aran_spherical_seriesd_get_term (ass, 0, 0) = 1.;
  *aran_spherical_seriesd_get_term (ass, 1, -1) = 1. - 1 * G_I;
  *aran_spherical_seriesd_get_term (ass, 1, 0) = 1.;
  *aran_spherical_seriesd_get_term (ass, 1, 1) = 1. + 1 * G_I;

  ass2 = aran_spherical_seriesd_clone (ass);

  aran_spherical_seriesd_free (ass);
  aran_spherical_seriesd_free (ass2);

/*   /\* check single harmonics and small sums *\/ */
/*   ass = create_series1 (); */
/*   check (ass, 1., function1); */
/*   aran_spherical_seriesd_free (ass); */

/*   ass = create_series2 (); */
/*   check (ass, 1., function2); */
/*   aran_spherical_seriesd_free (ass); */

/*   ass = create_series3 (); */
/*   check (ass, 1., function3); */
/*   aran_spherical_seriesd_free (ass); */

/*   ass = create_series4 (); */
/*   check (ass, 1., function4); */
/*   aran_spherical_seriesd_free (ass); */

/*   ass = create_series5 (); */
/*   check (ass, 1., function5); */
/*   aran_spherical_seriesd_free (ass); */

  /* check legendre addition theorem formulas */
  ass = legendre_addition (0);
  check (ass, 1., P_0);
  aran_spherical_seriesd_free (ass);

  ass = legendre_addition (1);
  check (ass, 1., P_1);
  aran_spherical_seriesd_free (ass);

  ass = legendre_addition (2);
  check (ass, 1., P_2);
  aran_spherical_seriesd_free (ass);

  ass = legendre_addition (3);
  check (ass, 1., P_3);
  aran_spherical_seriesd_free (ass);

  p.x = 0.8660254037844386;
  p.y = 0.8660254037844386;
  p.z = 0.8660254037844386;

  ass = legendre_addition (2);
  check (ass, 2., P_2);
  aran_spherical_seriesd_free (ass);

  p.x = -0.5;
  p.y = 0.5;
  p.z = 2.5;

  ass = legendre_addition (3);
  check (ass, 4., P_3);
  aran_spherical_seriesd_free (ass);


  return ret;
}

