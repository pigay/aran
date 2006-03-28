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
#include "aran/aransolver3d.h"
#include "aran/aranbinomial.h"


/* tree bbox size */
#define TR (1.)
/* circle radius */
#define R (0.999 * TR)



typedef struct _PointAccum PointAccum;

struct _PointAccum
{
  VsgVector3d vector;

  gdouble density;

  gcomplex128 accum;

  guint id;
};



static void p2p (PointAccum *one, PointAccum *other)
{
  if (one != other)
    {
      VsgVector3d tmp;
      gdouble inv_r;

      /* destination - source */
      vsg_vector3d_sub (&one->vector, &other->vector, &tmp);

      inv_r = 1. / vsg_vector3d_norm (&tmp);

      /* G(xd, xs) = 1./(xd-xs) * src_density */

      one->accum += inv_r * other->density;
      other->accum += inv_r * one->density;
    }
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

static void l2p (const VsgVector3d *center, AranDevelopment3d *devel,
		 PointAccum *particle)
{
  particle->accum += aran_development3d_local_evaluate (center, devel,
							&particle->vector);
}

static void _direct (PointAccum **points, guint np)
{
  guint i, j;

  for (i=0; i<np; i ++)
    {
      for (j=i+1; j<np; j ++)
        {
          p2p (points[i], points[j]);
        }
    }
}

static void _one_circle_distribution (PointAccum **points,
				      AranSolver3d *solver);
static void _random_distribution (PointAccum **points,
				  AranSolver3d *solver);

static void _grid_distribution (PointAccum **points,
                                AranSolver3d *solver);

static gdouble err_lim = 1.E-3;
static guint np = 12;
static guint order = 20;
static gboolean check = TRUE;
static gboolean direct = FALSE;
static guint maxbox = 1;

static AranMultipole2MultipoleFunc3d m2m =
(AranMultipole2MultipoleFunc3d) aran_development3d_m2m;

static AranMultipole2LocalFunc3d m2l =
(AranMultipole2LocalFunc3d) aran_development3d_m2l;

static AranLocal2LocalFunc3d l2l =
(AranLocal2LocalFunc3d) aran_development3d_l2l;

static void (*_distribution) (PointAccum **, AranSolver3d *solver) =
_one_circle_distribution;


static
void parse_args (int argc, char **argv)
{
  int iarg = 1;
  char *arg;

  while (iarg < argc)
    {
      arg = argv[iarg];

      if (g_ascii_strcasecmp (arg, "-np") == 0)
	{
	  guint tmp = 0;
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (sscanf (arg, "%u", &tmp) == 1)
	      np = tmp;
	  else
	    g_printerr ("Invalid particles number (-np %s)\n", arg);
	}
      else if (g_ascii_strcasecmp (arg, "-pr") == 0)
	{
	  guint tmp = 0;
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (sscanf (arg, "%u", &tmp) == 1 && tmp > 0)
	      order = tmp;
	  else
	    g_printerr ("Invalid precision order value (-pr %s)\n", arg);
	}
      else if (g_ascii_strcasecmp (arg, "-s") == 0)
	{
	  guint tmp = 0;
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (sscanf (arg, "%u", &tmp) == 1 && tmp > 0)
	      maxbox = tmp;
	  else
	    g_printerr ("Invalid maximum box size value (-s %s)\n", arg);
	}
      else if (g_ascii_strcasecmp (arg, "-err") == 0)
	{
	  gdouble tmp = 0;
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (sscanf (arg, "%lf", &tmp) == 1 && tmp > 0.)
	      err_lim = tmp;
	  else
	    g_printerr ("Invalid error limit value (-err %s)\n", arg);
	}
      else if (g_ascii_strcasecmp (arg, "-nocheck") == 0)
	{
	  check = FALSE;
	}

      else if (g_ascii_strcasecmp (arg, "-check") == 0)
	{
	  check = TRUE;
	}
      else if (g_ascii_strcasecmp (arg, "-direct") == 0)
	{
	  direct = TRUE;
	}
      else if (g_ascii_strcasecmp (arg, "-dist") == 0)
	{
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (g_ascii_strcasecmp (arg, "circle") == 0)
	    {
	      _distribution = _one_circle_distribution;
	    }
	  else if (g_ascii_strcasecmp (arg, "random") == 0)
	    {
	      _distribution = _random_distribution;
	    }
	  else if (g_ascii_strcasecmp (arg, "grid") == 0)
	    {
	      _distribution = _grid_distribution;
	    }
	  else 
	    {
	      g_printerr ("Invalid distribution name (-dist %s)\n", arg);
	    }
	}
      else if (g_ascii_strcasecmp (arg, "-translation") == 0)
	{
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (g_ascii_strcasecmp (arg, "normal") == 0)
	    {
              m2m = (AranMultipole2MultipoleFunc3d) aran_development3d_m2m;

              m2l = (AranMultipole2LocalFunc3d) aran_development3d_m2l;

              l2l = (AranLocal2LocalFunc3d) aran_development3d_l2l;
	    }
	  else if (g_ascii_strcasecmp (arg, "kkylin") == 0)
	    {
              m2m =
                (AranMultipole2MultipoleFunc3d) aran_development3d_m2m_kkylin;

              m2l = (AranMultipole2LocalFunc3d) aran_development3d_m2l_kkylin;

              l2l = (AranLocal2LocalFunc3d) aran_development3d_l2l_kkylin;
	    }
	  else 
	    {
	      g_printerr ("Invalid translation name (-translation %s)\n", arg);
	    }
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

static void _one_circle_distribution (PointAccum **points, AranSolver3d *solver)
{
  guint i;

  for (i=0; i<np; i++)
    {
      PointAccum *point = g_malloc0 (sizeof (PointAccum));

      point->vector.x = R * cos(2.*G_PI*i/np);
      point->vector.y = R * sin(2.*G_PI*i/np);
      point->vector.z = 0.0;
/*       point->vector.z = point->vector.x; */
      point->density = 1./np;
      point->accum = 0.;
      point->id = i;

      points[i] = point;

      if (!direct) aran_solver3d_insert_point (solver, point);
    }
}

static void _random_distribution (PointAccum **points,
				  AranSolver3d *solver)
{
  guint i, cptr = 0;
  gdouble tol = 1.e-2;

  aran_solver3d_set_tolerance (solver, tol);

  for (i=0; i<np; i++)
    {
      PointAccum *point = g_malloc0 (sizeof (PointAccum));
      PointAccum *found;

      point->vector.x = g_random_double_range (-R,R);
      point->vector.y = g_random_double_range (-R,R);
      point->vector.z = g_random_double_range (-R,R);
      point->density = g_random_double_range (0., 2./np);
      point->accum = 0.;
      point->id = cptr;

      if (direct)
        {
          points[cptr] = point;
          cptr ++;
        }
      else
        {
          found = aran_solver3d_find_point (solver, point);

          if (found != NULL)
            {
              found->density += point->density;

              g_free (point);
            }
          else
            {
              points[cptr] = point;
              aran_solver3d_insert_point (solver, point);
              cptr ++;
            }
        }
    }

  np = cptr;
}

static void _grid_distribution (PointAccum **points,
                                AranSolver3d *solver)
{
  guint i, j, k, size, cptr = 0;
  gdouble x, y, z, density;

  size = floor (pow (np, 1./3.));

  np = size*size*size;
  density = 1./np;

  for (i=0; i<size; i ++)
    {
      x = R * (2.*i/(size-1.)-1.);
      for (j=0; j<size; j ++)
        {
          y = R * (2.*j/(size-1.)-1.);
          for (k=0; k<size; k ++)
            {
              PointAccum *point = g_malloc0 (sizeof (PointAccum));

              z = R * (2.*k/(size-1.)-1.);

              point->vector.x = x;
              point->vector.y = y;
              point->vector.z = z;
              point->density = density;
              point->accum = 0.;
              point->id = cptr;

              if (!direct)
                aran_solver3d_insert_point (solver, point);

              points[cptr] = point;

              cptr ++;
            }
        }
    }
}

int main (int argc, char **argv)
{
  VsgVector3d lbound = {-TR, -TR, -TR};
  VsgVector3d ubound = {TR, TR, TR};
  PointAccum **points;
  VsgPRTree3d *prtree;
  AranSolver3d *solver;
  int ret = 0;
  guint i;

  aran_init();

  parse_args (argc, argv);

  points = g_malloc0 (np * sizeof (PointAccum *));

  prtree =
    vsg_prtree3d_new_full (&lbound, &ubound,
			    (VsgPoint3dLocFunc) vsg_vector3d_vector3d_locfunc,
			    (VsgPoint3dDistFunc) vsg_vector3d_dist,
			    (VsgRegion3dLocFunc) NULL, maxbox);

  solver = aran_solver3d_new (prtree, ARAN_TYPE_DEVELOPMENT3D,
			      aran_development3d_new (0, order),
			      (AranZeroFunc) aran_development3d_set_zero);

  aran_solver3d_set_functions (solver,
			       (AranParticle2ParticleFunc3d) p2p,
			       (AranParticle2MultipoleFunc3d) p2m,
                               m2m,
			       m2l,
			       l2l,
			       (AranLocal2ParticleFunc3d)l2p);

  _distribution (points, solver);

/*   g_printerr ("ok depth = %d size = %d\n", */
/*               aran_solver3d_depth (solver), */
/*               aran_solver3d_point_count (solver)); */

  if (direct) _direct (points, np);
  else aran_solver3d_solve (solver);

/*   vsg_prtree3d_write (prtree, stderr); */

  aran_solver3d_free (solver);

  if (check)
    {
      for (i=0; i<np; i++)
	{
	  guint j;
	  gcomplex128 sum = 0.;
	  gcomplex128 err;

	  for (j=0; j<np; j++)
	    {
	      if (i != j)
		{
		  VsgVector3d tmp;
		  gdouble inv_r;

		  vsg_vector3d_sub (&points[i]->vector, &points[j]->vector,
				    &tmp);

		  inv_r = 1. / vsg_vector3d_norm (&tmp);
		  sum += inv_r * points[j]->density;
		}
	    }

	  err = (points[i]->accum - sum) /
	    MAX(cabs (points[i]->accum), cabs (sum));

	  if (cabs (err) > err_lim)
	    {
	      g_printerr ("Error: pt%u (%e,%e): (%e,%e)!=(%e,%e) -> ",
			  i,
			  creal (err), cimag (err),
			  creal (points[i]->accum), cimag (points[i]->accum),
			  creal (sum), cimag (sum));
	      vsg_vector3d_write (&points[i]->vector, stderr);
	      g_printerr ("\n");
	      ret ++;
	    }
	}
    }

  for (i=0; i<np; i++)
    {
      g_free (points[i]);
    }

  g_free (points);

  return ret;
}

