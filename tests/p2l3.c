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
#include "aran/aransphericalseriesd.h"
#include "aran/arandevelopment3d.h"

#define N 10

static gdouble epsilon = 1.e-9;

typedef struct _Particle Particle;

struct _Particle
{
  VsgVector3d vector;

  gdouble charge;

  VsgVector3d field;
  gdouble potential;
};

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

/* check Local expansion against Newton gradient for several points inside a ball */
static void check (const gchar *log,
                   AranDevelopment3d *dev, gdouble radius,
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
      /* iterate through Phi interval [0, 2pi) */
      p = 2.*G_PI * i/(N-1.);
      cosp = cos (p);
      sinp = sin (p);

      for (j=0; j<N; j ++)
	{
          /* iterate through Teta interval [0, pi) */
	  t = G_PI * j/(N-1.);
	  cost = cos (t);
	  sint = sin (t);

	  for (k=0; k<N; k ++)
	    {
              /* iterate through r interval [0, radius] */
	      VsgVector3d vec;
	      VsgVector3d vref;
	      VsgVector3d vres;
	      VsgVector3d verr;

	      r = radius * k/(N-1.);

	      vsg_vector3d_from_spherical_internal (&vec, r,
						    cost, sint,
						    cosp, sinp);

/*               aran_spherical_seriesd_gradient_evaluate (ass, &vec, */
/*                                                         &vres); */

              aran_spherical_seriesd_gradient_evaluate_internal (dev->local,
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

	      if (fabs (err) > epsilon || !finite (err))
		{
		  g_printerr ("Error %s (r=%f, t=%f, p=%f) : " \
                              "ref(%f,%f,%f) != res(%f,%f,%f) -> %e\n",
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


/* compute some particle contribution to a Local expansion */
void aran_development3d_p2l (const VsgVector3d *position, const gdouble charge,
                             const VsgPRTree3dNodeInfo *dst_node,
                             AranDevelopment3d *devel)
{
  VsgVector3d tmp;
  guint deg = aran_spherical_seriesd_get_posdeg (devel->local);
  gint l, m;
  gcomplex128 harmonics[((deg+1)*(deg+2))/2];
  gdouble r, cost, sint, cosp, sinp;
  gcomplex128 expp;
  gdouble fact, inv_r;

  vsg_vector3d_sub (position, &dst_node->center, &tmp);

  vsg_vector3d_to_spherical_internal (&tmp, &r, &cost, &sint, &cosp, &sinp);
  expp = cosp + G_I * sinp;

  aran_spherical_harmonic_evaluate_multiple_internal (deg, cost, sint, expp,
                                                      harmonics);


  *aran_spherical_seriesd_get_term (devel->local, 0, 0) += 0.;

  inv_r = 1. / r;
  fact = charge * inv_r;

  for (l=0; l<=deg; l ++)
    {
      gcomplex128 *ptr;
      gcomplex128 term;

      term = fact * (4.*G_PI / (l+l+1.)) *
        *aran_spherical_harmonic_multiple_get_term (l, 0, harmonics);

      /* *aran_spherical_seriesd_get_term (devel->local, l, 0) = conj (term); */
      ptr = aran_spherical_seriesd_get_term (devel->local, l, 0);
      ptr[0] += conj (term);

      for (m=1; m<=l; m ++)
        {

          term = fact * (4.*G_PI / (l+l+1.)) *
            *aran_spherical_harmonic_multiple_get_term (l, m, harmonics);

          /* *aran_spherical_seriesd_get_term (devel->local, l, m) = conj (term); */
          ptr[m] += conj (term);
        }
      fact *= inv_r;
    }

}

Particle particle;

/* Newton force field evaluation at position vec */
void newtongrad (VsgVector3d *vec, VsgVector3d *grad)
{
  gdouble r;

  vsg_vector3d_sub (vec, &particle.vector, grad);

  r = vsg_vector3d_norm (grad);

  vsg_vector3d_scalp (grad, particle.charge / (r*r*r), grad);
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
  AranDevelopment3d *dev;
  VsgPRTree3dNodeInfo node_info;

  parse_args (argc, argv);

  particle.vector.x = 3.;
  particle.vector.y = -1.;
  particle.vector.z = -1.;

  particle.charge = 1.;

  node_info.center.x = 5.;
  node_info.center.y = 0.;
  node_info.center.z = 0.;

  dev = aran_development3d_new (deg, 0);
  aran_development3d_p2l (&particle.vector, particle.charge, &node_info, dev);

  /* aran_spherical_seriesd_write (dev->local, stderr); */
  /* g_printerr ("\n\n"); */
  check ("devel", dev, 1., &node_info.center, newtongrad);

  aran_development3d_free (dev);

  return ret;
}

