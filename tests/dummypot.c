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
#include "aran/aransolver2d.h"
#include "aran/aranbinomial.h"
#include "aran/aranprofile.h"
#include "aran/aranprofiledb.h"

#include "glib/gprintf.h"

/* approximation degree */
#define K 20

/* circle radius */
#define R 0.9



typedef struct _PointAccum PointAccum;

struct _PointAccum
{
  VsgVector2d vector;

  gdouble density;

  gcomplex128 accum;

  guint id;
};

static void point_accum_clear_accum (PointAccum *pa)
{
  pa->accum = 0.;
}

void p2p (PointAccum *one, PointAccum *other)
{
  if (one != other)
    {
      /* destination - source */
      gcomplex128 zd_m_zs = (one->vector.x + G_I*one->vector.y) -
	(other->vector.x + G_I*other->vector.y);
      gcomplex128 inv_zd_m_zs = 1. / zd_m_zs;

      /* G(zd, zs) = 1./(zd-zs) * src_density */

      one->accum += inv_zd_m_zs * other->density;
      other->accum += - inv_zd_m_zs * one->density;
    }
}

void p2m (PointAccum *particle, const VsgPRTree2dNodeInfo *devel_node,
          AranDevelopment2d *devel)
{
  guint i;
  gcomplex128 *multipole = aran_laurent_seriesd_get_term (devel->multipole, 0);
  gcomplex128 tmp, zp_m_zm;

  tmp = particle->density;

  zp_m_zm = (particle->vector.x + G_I*particle->vector.y) -
    (devel_node->center.x + G_I*devel_node->center.y);

  /* a_0 = 0 */
  multipole[0] += 0.;

  for (i=1; i<aran_laurent_seriesd_get_negdeg (devel->multipole); i++)
    {
      /* a_i = (zp-zm)^(k-1) */
      multipole[i] += tmp;
      tmp *= zp_m_zm;
    }
}

void l2p (const VsgPRTree2dNodeInfo *devel_node, AranDevelopment2d *devel,
          PointAccum *particle)
{

  particle->accum += aran_development2d_local_evaluate (devel_node,
                                                        devel,
							&particle->vector);
}

static void _profile_operators (const gchar *db_filename)
{
  GKeyFile *profiles_file = g_key_file_new ();
  gchar *profiles_data;
  gchar comment[1024] = {'\0', };
  PointAccum p1 = {{0.1, 0.1}, 0.1, 0., 0};
  PointAccum p2 = {{0.5, 0.5}, 0.1, 0., 1};
  gdouble chisq, t;
  AranPoly1d *ap1d = aran_poly1d_new (2);
  gint maxbox = 100;
  const gchar *profiles_group;
  FILE *f = stdout;

  profiles_group = g_getenv ("ARAN_PROFILE_GROUP");
  if (profiles_group == NULL) profiles_group = ARAN_PROFILE_DB_DEFAULT_GROUP;

  if (db_filename != NULL)
    {
      g_key_file_load_from_file (profiles_file, db_filename,
                                 G_KEY_FILE_KEEP_COMMENTS |
                                 G_KEY_FILE_KEEP_TRANSLATIONS, NULL);
      f = fopen (db_filename, "w");
    }

  ap1d->degree = 0;
  t = aran_profile_nearfunc_2d ((AranParticle2ParticleFunc2d) p2p,
                                (AranParticleInitFunc2d) point_accum_clear_accum,
                                &p1, &p2, maxbox);
  ap1d->terms[0] = t / (maxbox * maxbox);
  aran_poly1d_write_key_file (ap1d, profiles_file, profiles_group, "p2p");

  ap1d->degree = 1;
  chisq =
    aran_poly1d_profile_p2m_2d ((AranParticle2MultipoleFunc2d) p2m,
                                (AranZeroFunc) aran_development2d_set_zero,
                                &p1,
                                (AranDevelopmentNewFunc)
                                aran_development2d_new,
                                (GDestroyNotify) aran_development2d_free,
                                ap1d, 20);
  aran_poly1d_write_key_file (ap1d, profiles_file, profiles_group, "p2m");
  g_sprintf (comment, " \"%s\" function fitting: chisq=%g", "p2m", chisq);
  g_key_file_set_comment (profiles_file, profiles_group, "p2m", comment,NULL);

  ap1d->degree = 1;
  t =
    aran_poly1d_profile_l2p_2d ((AranLocal2ParticleFunc2d) l2p,
                                NULL,
                                &p1,
                                (AranDevelopmentNewFunc) aran_development2d_new,
                                (GDestroyNotify) aran_development2d_free,
                                ap1d, 20);
  aran_poly1d_write_key_file (ap1d, profiles_file, profiles_group, "l2p");
  g_sprintf (comment, " \"%s\" function fitting: chisq=%g", "l2p", chisq);
  g_key_file_set_comment (profiles_file, profiles_group, "l2p", comment,NULL);
  aran_poly1d_free (ap1d);

  profiles_data = g_key_file_to_data (profiles_file, NULL, NULL);
  g_key_file_free (profiles_file);

  g_free (profiles_data);
  g_fprintf (f, "%s", profiles_data);

  if (db_filename != NULL) fclose (f);
}

static void _one_circle_distribution (PointAccum **points,
				      AranSolver2d *solver);
static void _random_distribution (PointAccum **points,
				  AranSolver2d *solver);

static gdouble err_lim = 1.E-6;
static guint np = 12;
static guint order = K;
static gboolean check = TRUE;
static guint maxbox = 1;

static void (*_distribution) (PointAccum **, AranSolver2d *solver) =
_one_circle_distribution;

static
void parse_args (int argc, char **argv)
{
  int iarg = 1;
  char *arg;

  while (iarg < argc)
    {
      arg = argv[iarg];

      if (g_ascii_strcasecmp (arg, "-profile") == 0)
	{
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

          _profile_operators (arg);

          exit (0);
	}
      else if (g_ascii_strcasecmp (arg, "-np") == 0)
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
	  else 
	    {
	      g_printerr ("Invalid distribution name (-dist %s)\n", arg);
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

static void _one_circle_distribution (PointAccum **points, AranSolver2d *solver)
{
  guint i;

  for (i=0; i<np; i++)
    {
      PointAccum *point = g_malloc0 (sizeof (PointAccum));

      point->vector.x = R * cos(2.*G_PI*i/np);
      point->vector.y = R * sin(2.*G_PI*i/np);
      point->density = 1.;
      point->accum = 0.;
      point->id = i;

      points[i] = point;

      aran_solver2d_insert_point (solver, point);
    }
}

static void _random_distribution (PointAccum **points,
				  AranSolver2d *solver)
{
  guint i;

  for (i=0; i<np; i++)
    {
      PointAccum *point = g_malloc0 (sizeof (PointAccum));

      point->vector.x = g_random_double_range (-R,R);
      point->vector.y = g_random_double_range (-R,R);
      point->density = 1.;
      point->accum = 0.;
      point->id = i;

      points[i] = point;

      aran_solver2d_insert_point (solver, point);
    }
}

int main (int argc, char **argv)
{
  VsgVector2d lbound = {-1., -1.};
  VsgVector2d ubound = {1., 1.};
  PointAccum **points;
  VsgPRTree2d *prtree;
  AranSolver2d *solver;
  int ret = 0;
  guint i;

  aran_init();

  parse_args (argc, argv);

  points = g_malloc (np * sizeof (PointAccum*));

  prtree = 
    vsg_prtree2d_new_full (&lbound, &ubound,
			    (VsgPoint2dLocFunc) vsg_vector2d_vector2d_locfunc,
			    (VsgPoint2dDistFunc) vsg_vector2d_dist,
			    (VsgRegion2dLocFunc) NULL, maxbox);

  aran_binomial_require (2*order);

  solver = aran_solver2d_new (prtree, ARAN_TYPE_DEVELOPMENT2D,
			      aran_development2d_new (0, order),
			      (AranZeroFunc) aran_development2d_set_zero);

  aran_solver2d_set_functions (solver,
			       (AranParticle2ParticleFunc2d) p2p,
			       (AranParticle2MultipoleFunc2d) p2m,
			       (AranMultipole2MultipoleFunc2d) aran_development2d_m2m,
			       (AranMultipole2LocalFunc2d) aran_development2d_m2l,
			       (AranLocal2LocalFunc2d) aran_development2d_l2l,
			       (AranLocal2ParticleFunc2d)l2p);

  _distribution (points, solver);

  aran_solver2d_solve (solver);

/*   vsg_prtree2d_write (prtree, stderr); */

  aran_solver2d_free (solver);

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
		  gcomplex128 zd_m_zs =
		    (points[i]->vector.x + I*points[i]->vector.y) -
		    (points[j]->vector.x + I*points[j]->vector.y);

		  sum += 1./zd_m_zs * points[j]->density;
		}
	    }

	  err = (points[i]->accum - sum) /
	    MAX(cabs (points[i]->accum), cabs (sum));

	  if (cabs (err) > err_lim || !finite (cabs (err)))
	    {
	      g_printerr ("Error: pt%u (%e,%e) -> ",
			  i,
			  creal (err), cimag (err));
	      vsg_vector2d_write (&points[i]->vector, stderr);
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

