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

#include "aran/aranwigner.h"
#include "aran/aran.h"
#include "aran/aransphericalharmonic.h"
#include "aran/aransphericalseriesd.h"

#define PHASE(m) (((m)%2 == 0) ? 1. : -1.)

static gdouble epsilon = 1.e-7;
static guint N = 10;
static guint L = 29;

static VsgVector3d zero = {0., 0., 0.};

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
          iarg++;

          arg = (iarg < argc) ? argv[iarg] : NULL;

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
            L = tmp;
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

      iarg++;
    }
}

static AranSphericalSeriesd *create_taylor (guint deg, VsgVector3d *center,
                                            VsgVector3d *p)
{
  AranSphericalSeriesd *ass;
  gint l, m;
  const guint size = ((deg+1)*(deg+2))/2;
  gcomplex128 harmonics[size];
  gdouble r, cost, sint, cosp, sinp;
  gcomplex128 expp;
  gdouble fact, inv_r;
  VsgVector3d tmp;

  ass = aran_spherical_seriesd_new (deg, 0);
  aran_spherical_seriesd_set_zero (ass);

  vsg_vector3d_sub (p, center, &tmp);

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

static gint buffer_diff (gint L, gcomplex128 *ref, gcomplex128 *res)
{
  gint faults = 0;
  gint l, m;

  for (l = 0; l <= L; l++)
    {
      for (m = 0; m <= l; m++)
        {
          gcomplex128 err;
          gint i = (l*(l+1))/2 + m;

          err = ref[i] - res[i];

          if (cabs (ref[i]) > epsilon)
            err /= ref[i];

          if (cabs (err) > epsilon || !finite (cabs (err)))
            {
              g_printerr ("AranSphericalSeriesd error (%g,%g) != (%g,%g)\n",
                          creal (ref[i]), cimag (ref[i]),
                          creal (res[i]), cimag (res[i]));
              faults ++;
            }
        }
    }

  return faults;
}

static gint ass_diff (AranSphericalSeriesd *res, AranSphericalSeriesd *ref)
{
  gint faults = 0;
  gint L = aran_spherical_seriesd_get_posdeg (ref);

  if (L != aran_spherical_seriesd_get_posdeg (res))
    {
      g_printerr ("AranSphericalSeriesd degree error %d != %d\n", L,
                  aran_spherical_seriesd_get_posdeg (res));
      faults ++;
    }
  else
    {
      faults += buffer_diff (L,
                             aran_spherical_seriesd_get_term (ref, 0, 0),
                             aran_spherical_seriesd_get_term (res, 0, 0));
    }

  L = aran_spherical_seriesd_get_negdeg (ref);

  if (L != aran_spherical_seriesd_get_negdeg (res))
    {
      g_printerr ("AranSphericalSeriesd degree error %d != %d\n", L,
                  aran_spherical_seriesd_get_negdeg (res));
      faults ++;
    }
  else if (L > 0)
    {
      faults += buffer_diff (L-1,
                             aran_spherical_seriesd_get_term (ref, -1, 0),
                             aran_spherical_seriesd_get_term (res, -1, 0));
    }

  return faults;
}

gint check (gdouble alpha, gdouble beta, gdouble gamma)
{
  gint faults = 0;
  VsgVector3d p = {3., 3., 3.};
  VsgVector3d prot = VSG_V3D_ZERO;
  VsgMatrix4d mat = VSG_M4D_ID;
  AranSphericalSeriesd *ass, *rot, *chk, *back;

  vsg_matrix4d_rotate_euler (&mat, alpha, beta, gamma);

  vsg_matrix4d_vecmult (&mat, &p, &prot);

/*   vsg_vector3d_write (&p, stderr); */
/*   g_printerr ("\n"); */
/*   vsg_vector3d_write (&prot, stderr); */
/*   g_printerr ("\n"); */

  ass = create_taylor (L, &zero, &p);
  chk = create_taylor (L, &zero, &prot);

  rot = aran_spherical_seriesd_new (L, 0);
  back = aran_spherical_seriesd_new (L, 0);

  aran_spherical_seriesd_rotate (ass, alpha, beta, gamma, rot);
  aran_spherical_seriesd_rotate_inverse (rot, alpha, beta, gamma, back);

  faults += ass_diff (chk, rot);

  faults += ass_diff (back, ass);

  aran_spherical_seriesd_free (ass);
  aran_spherical_seriesd_free (rot);
  aran_spherical_seriesd_free (chk);
  aran_spherical_seriesd_free (back);

  return faults;
}

int main (int argc, char **argv)
{
  gint i, j;
  gint ret = 0;

  aran_init ();

  parse_args (argc, argv);

  for (i=0; i<N; i ++)
    {
      gdouble theta = i*G_PI/N;

      for (j=0; j<N; j ++)
        {
          gdouble phi = j*2.*G_PI/N;

          ret += check (phi, theta, 0.);
          ret += check (phi, theta, phi);
        }
    }
  return ret;
}
