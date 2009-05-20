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
#include <string.h>

#include "aran/aran.h"
#include "aran/aransolver2d.h"
#include "aran/aranbinomial.h"

/* Particle type */
typedef struct _PointAccum PointAccum;

struct _PointAccum
{
  VsgVector2d vector;

  gdouble density;

  gcomplex128 accum;

  gint id;
};

/* Points number */
#define N 3600000

/* approximation degree */
#define K 20

/* circle radius */
#define R 0.9

static gint rk = 0;
static gint sz = 1;

static GPtrArray *points = NULL;
static PointAccum *check_points = NULL;

static guint32 _random_seed = 0;
static gint _flush_interval = 1000;



static void p2p (PointAccum *one, PointAccum *other)
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

/*       g_printerr ("%d : p2p one=", rk); */
/*       vsg_vector2d_write (&one->vector, stderr); */
/*       g_printerr (" other="); */
/*       vsg_vector2d_write (&other->vector, stderr); */
/*       g_printerr ("\n"); */
    }
}

static void p2m (PointAccum *particle, const VsgPRTree2dNodeInfo *devel_node,
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

static void l2p (const VsgPRTree2dNodeInfo *devel_node,
                 AranDevelopment2d *devel,
		 PointAccum *particle)
{
  particle->accum += aran_development2d_local_evaluate (devel_node,
                                                        devel,
							&particle->vector);
/*   g_printerr ("%d : l2p v=", rk); */
/*   vsg_vector2d_write (&particle->vector, stderr); */
/*   g_printerr (" accum=(%e,%e)\n", creal (particle->accum), */
/*               cimag (particle->accum)); */

}

PointAccum *point_accum_alloc (gboolean resident, gpointer user_data)
{
  PointAccum *ret;
  ret = g_malloc0 (sizeof (PointAccum));

  if (resident)
    g_ptr_array_add (points, ret);

  return ret;
}

void point_accum_destroy (PointAccum *data, gboolean resident,
                 gpointer user_data)
{
  if (resident)
    g_ptr_array_remove (points, data);

  g_free (data);
}

#ifdef VSG_HAVE_MPI

/* migration pack/unpack functions */
void point_accum_migrate_pack (PointAccum *pt, VsgPackedMsg *pm,
                               gpointer user_data)
{
  vsg_packed_msg_send_append (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR2D);
  vsg_packed_msg_send_append (pm, &pt->density, 1, MPI_DOUBLE);
  vsg_packed_msg_send_append (pm, &pt->accum, 1, ARAN_MPI_TYPE_GCOMPLEX128);
  vsg_packed_msg_send_append (pm, &pt->id, 1, MPI_INT);
}

void point_accum_migrate_unpack (PointAccum *pt, VsgPackedMsg *pm,
                                 gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR2D);
  vsg_packed_msg_recv_read (pm, &pt->density, 1, MPI_DOUBLE);
  vsg_packed_msg_recv_read (pm, &pt->accum, 1, ARAN_MPI_TYPE_GCOMPLEX128);
  vsg_packed_msg_recv_read (pm, &pt->id, 1, MPI_INT);
}

/* visit forward pack/unpack functions */
/* only the coordinates part is transmitted */
void point_accum_visit_fw_pack (PointAccum *pt, VsgPackedMsg *pm,
                                gpointer user_data)
{
/*   g_printerr ("%d : fw pack d=%e v=", */
/*               rk, */
/*               pt->density); */
/*   vsg_vector2d_write (&pt->vector, stderr); */
/*   g_printerr ("\n"); */

  vsg_packed_msg_send_append (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR2D);
  vsg_packed_msg_send_append (pm, &pt->density, 1, MPI_DOUBLE);
}

void point_accum_visit_fw_unpack (PointAccum *pt, VsgPackedMsg *pm,
                                  gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR2D);
  vsg_packed_msg_recv_read (pm, &pt->density, 1, MPI_DOUBLE);
  pt->accum = 0.;


/*   g_printerr ("%d : fw unpack d=%e v=", */
/*               rk, */
/*               pt->density); */
/*   vsg_vector2d_write (&pt->vector, stderr); */
/*   g_printerr ("\n"); */
}

/* visit forward pack/unpack functions */
/* only the accum part is transmitted */
void point_accum_visit_bw_pack (PointAccum *pt, VsgPackedMsg *pm,
                                gpointer user_data)
{
  vsg_packed_msg_send_append (pm, &pt->accum, 1, ARAN_MPI_TYPE_GCOMPLEX128);

/*   g_printerr ("%d : bw pack accum=(%e,%e)\n", */
/*               rk, */
/*               creal (pt->accum), cimag (pt->accum)); */
}

void point_accum_visit_bw_unpack (PointAccum *pt, VsgPackedMsg *pm,
                                  gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, &pt->accum, 1, ARAN_MPI_TYPE_GCOMPLEX128);
}

void point_accum_visit_bw_reduce (PointAccum *a, PointAccum *b,
                                  gpointer user_data)
{
/*   g_printerr ("%d : bw red a=(%e,%e) b=(%e,%e) a+b=(%e,%e)\n", */
/*               rk, */
/*               creal (a->accum), cimag (a->accum), */
/*               creal (b->accum), cimag (b->accum), */
/*               creal (a->accum+b->accum), cimag (a->accum+b->accum)); */

  b->accum += a->accum;
}

VsgParallelVTable point_accum_vtable = {
  (VsgMigrableAllocDataFunc) point_accum_alloc, NULL,
  (VsgMigrableDestroyDataFunc) point_accum_destroy, NULL,
  {(VsgMigrablePackDataFunc) point_accum_migrate_pack, NULL,
   (VsgMigrablePackDataFunc) point_accum_migrate_unpack, NULL,
  },
  {(VsgMigrablePackDataFunc) point_accum_visit_fw_pack, NULL,
   (VsgMigrablePackDataFunc) point_accum_visit_fw_unpack, NULL,
  },
  {(VsgMigrablePackDataFunc) point_accum_visit_bw_pack, NULL,
   (VsgMigrablePackDataFunc) point_accum_visit_bw_unpack, NULL,
   (VsgMigrableReductionDataFunc) point_accum_visit_bw_reduce, NULL,
  },

};

#endif /* VSG_HAVE_MPI */

static void _one_circle_distribution (GPtrArray *points,
				      AranSolver2d *solver);
static void _random_distribution (GPtrArray *points,
				  AranSolver2d *solver);

static gdouble err_lim = 1.E-6;
static guint np = 12;
static guint order = K;
static gboolean check = TRUE;
static guint maxbox = 1;
static gboolean _verbose = FALSE;
static gboolean _hilbert = FALSE;

static void (*_distribution) (GPtrArray *, AranSolver2d *solver) =
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

	  if (sscanf (arg, "%u", &tmp) == 1 && tmp < N)
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
      else if (g_ascii_strcasecmp (arg, "-hilbert") == 0)
	{
	  _hilbert = TRUE;
	}
      else if (g_strncasecmp (arg, "-v", 2) == 0 ||
               g_strncasecmp (arg, "--verbose", 9) == 0)
        {
          _verbose = TRUE;
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

static void _one_circle_distribution (GPtrArray *points, AranSolver2d *solver)
{
  gint i;
  PointAccum *point;
  gdouble dtheta = 2. * G_PI / (np);
  gdouble theta0 = 0.;

  for (i=0; i<np; i++)
    {
      PointAccum tmp;

#ifdef VSG_HAVE_MPI
      if (i%_flush_interval == 0)
        {
          aran_solver2d_migrate_flush (solver);
          if (i%(_flush_interval*10) == 0)
            {
              if (_verbose && rk == 0)
                g_printerr ("%d: contiguous dist before %dth point\n", rk, i);

              aran_solver2d_distribute_contiguous_leaves (solver);
            }
        }
#endif /* VSG_HAVE_MPI */

      tmp.vector.x = R * cos (theta0 + i * dtheta);
      tmp.vector.y = R * sin (theta0 + i * dtheta);
      tmp.density = 1.;
      tmp.accum = 0.;
      tmp.id = i;

      if (check) memcpy (&check_points[i], &tmp, sizeof (PointAccum));

      if (i%sz != rk) continue;

      if (i % 10000 == 0 && _verbose)
        g_printerr ("%d: insert %dth point\n", rk, i);

      point = point_accum_alloc (TRUE, NULL);

      memcpy (point, &tmp, sizeof (PointAccum));

      aran_solver2d_insert_point (solver, point);
    }

#ifdef VSG_HAVE_MPI
  aran_solver2d_migrate_flush (solver);
  aran_solver2d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */
}

static void _random_distribution (GPtrArray *points,
				  AranSolver2d *solver)
{
  guint i;
  PointAccum *point;
  GRand *rand = g_rand_new_with_seed (_random_seed);

  for (i=0; i<np; i++)
    {
      PointAccum tmp;

#ifdef VSG_HAVE_MPI
      if (i%_flush_interval == 0)
        {
          aran_solver2d_migrate_flush (solver);
          if (i%(_flush_interval*10) == 0)
            {
              if (_verbose && rk == 0)
                g_printerr ("%d: contiguous dist before %dth point\n", rk, i);

              aran_solver2d_distribute_contiguous_leaves (solver);
            }
        }
#endif /* VSG_HAVE_MPI */

      tmp.vector.x = g_rand_double_range (rand, -R, R);;
      tmp.vector.y = g_rand_double_range (rand, -R, R);;
      tmp.density = 1.;
      tmp.accum = 0.;
      tmp.id = i;

      if (check) memcpy (&check_points[i], &tmp, sizeof (PointAccum));

      if (i%sz != rk) continue;

      if (i % 10000 == 0 && _verbose)
        g_printerr ("%d: insert %dth point\n", rk, i);

      point = point_accum_alloc (TRUE, NULL);

      memcpy (point, &tmp, sizeof (PointAccum));

      aran_solver2d_insert_point (solver, point);
    }

#ifdef VSG_HAVE_MPI
  aran_solver2d_migrate_flush (solver);
  aran_solver2d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */

  g_rand_free (rand);
}

void empty_array (gpointer var, gpointer data)
{
  g_free (var);
}

static gdouble maxerr = 0.;

void check_point_accum (PointAccum *point, gint *ret)
{
  gint i;
  gcomplex128 err;
  gdouble abserr;

  i = point->id;

  err = (point->accum - check_points[i].accum) /
    MAX(cabs (point->accum), cabs (check_points[i].accum));

  abserr = cabs (err);
  if (maxerr < abserr) maxerr = abserr;
 
  if (abserr > err_lim || !finite (abserr))
    {
      g_printerr ("%d : Error on pt%d res=(%e,%e) ref(%e,%e) err=(%e,%e) pt=",
                  rk, i,
                  creal (point->accum), cimag (point->accum),
                  creal (check_points[i].accum), cimag (check_points[i].accum),
                  creal (err), cimag (err));
      vsg_vector2d_write (&point->vector, stderr);
      g_printerr ("\n");
      (*ret) ++;
    }
}

int main (int argc, char **argv)
{
#ifdef VSG_HAVE_MPI
  VsgPRTreeParallelConfig pconfig = {NULL,};
#endif

  VsgVector2d lbound = {-1., -1.};
  VsgVector2d ubound = {1., 1.};
  VsgPRTree2d *prtree;
  AranSolver2d *solver;
  int ret = 0;
  guint i;
  GTimer *timer = NULL;

#ifdef VSG_HAVE_MPI
  MPI_Init (&argc, &argv);

  MPI_Comm_size (MPI_COMM_WORLD, &sz);
  MPI_Comm_rank (MPI_COMM_WORLD, &rk);
#endif

  aran_init();

  parse_args (argc, argv);

#ifdef VSG_HAVE_MPI
  pconfig.communicator = MPI_COMM_WORLD;

  pconfig.point = point_accum_vtable;

  aran_development2d_vtable_init (&pconfig.node_data, 0, order);
#endif

  points = g_ptr_array_new ();

  if (check)
    check_points = g_malloc0 (np * sizeof (PointAccum));

  prtree = 
    vsg_prtree2d_new_full (&lbound, &ubound,
			    (VsgPoint2dLocFunc) vsg_vector2d_vector2d_locfunc,
			    (VsgPoint2dDistFunc) vsg_vector2d_dist,
			    (VsgRegion2dLocFunc) NULL, maxbox);

  aran_binomial_require (2*order);

  solver = aran_solver2d_new (prtree, ARAN_TYPE_DEVELOPMENT2D,
			      aran_development2d_new (0, order),
			      (AranZeroFunc) aran_development2d_set_zero);

#ifdef VSG_HAVE_MPI
  aran_solver2d_set_parallel (solver, &pconfig);
#endif

  aran_solver2d_set_functions (solver,
			       (AranParticle2ParticleFunc2d) p2p,
			       (AranParticle2MultipoleFunc2d) p2m,
			       (AranMultipole2MultipoleFunc2d) aran_development2d_m2m,
			       (AranMultipole2LocalFunc2d) aran_development2d_m2l,
			       (AranLocal2LocalFunc2d) aran_development2d_l2l,
			       (AranLocal2ParticleFunc2d)l2p);

  if (_hilbert)
    {
      /* configure for hilbert curve order traversal */
      aran_solver2d_set_children_order_hilbert (solver);
    }

  _distribution (points, solver);

  if (_verbose)
    {
      g_printerr ("%d : solve begin\n", rk);

#ifdef VSG_HAVE_MPI
      MPI_Barrier (MPI_COMM_WORLD);
#endif

      timer = g_timer_new ();
    }

  aran_solver2d_solve (solver);

  if (_verbose)
    {
#ifdef VSG_HAVE_MPI
      MPI_Barrier (MPI_COMM_WORLD);
#endif

      g_printerr ("%d : solve ok elapsed=%f seconds\n", rk,
                  g_timer_elapsed (timer, NULL));

      g_timer_destroy (timer);
    }

  if (check)
    {
      for (i=0; i<np; i++)
	{
          PointAccum *pi = &check_points[i];
	  guint j;

	  for (j=0; j<np; j++)
	    {
	      if (i != j)
		{
                  PointAccum *pj = &check_points[j];
		  gcomplex128 zd_m_zs =
		    (pi->vector.x + I*pi->vector.y) -
		    (pj->vector.x + I*pj->vector.y);

		  pi->accum += 1./zd_m_zs * pj->density;
		}
	    }
        }

      aran_solver2d_foreach_point (solver, (GFunc) check_point_accum, &ret);

      if (_verbose)
        g_printerr ("%d : max err = %e\n", rk, maxerr);

      g_free (check_points);
    }

  aran_solver2d_free (solver);

#ifdef VSG_HAVE_MPI
  aran_development2d_vtable_clear (&pconfig.node_data);
#endif

  /* destroy the points */
  g_ptr_array_foreach (points, empty_array, NULL);
  g_ptr_array_free (points, TRUE);

#ifdef VSG_HAVE_MPI
  MPI_Finalize ();
#endif

  return ret;
}

