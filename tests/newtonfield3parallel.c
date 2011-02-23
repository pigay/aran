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
#include <string.h>

#include <complex.h>

#include <math.h>

#include <glib/gprintf.h>

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

  VsgVector3d field;
/*   gcomplex128 accum; */

  guint id;
};

static gint rk = 0;
static gint sz = 1;

static GPtrArray *points = NULL;
static PointAccum *check_points = NULL;
static guint _cp_size = 0;

static guint32 _random_seed = 0;
static gint _flush_interval = 1000;


static const gdouble epsilon = 1.e-5;

static void p2p_one_way (PointAccum *one, const PointAccum *other)
{
  if (one != other)
    {
      VsgVector3d tmp;
      VsgVector3d tmp1;
      gdouble r, inv_r, inv_r3;

      /* destination - source */
      vsg_vector3d_sub (&one->vector, &other->vector, &tmp);

      r = vsg_vector3d_norm (&tmp);

      if (r > epsilon)
        {
          inv_r = 1. / r;
          inv_r3 = inv_r*inv_r*inv_r;

/*           one->accum += inv_r * other->density; */

          vsg_vector3d_scalp (&tmp, - inv_r3*other->density, &tmp1);
          vsg_vector3d_add (&one->field, &tmp1, &one->field);
        }
    }
}

static void p2p (PointAccum *one, PointAccum *other)
{
  if (one != other)
    {
      VsgVector3d tmp;
      VsgVector3d tmp1;
      VsgVector3d tmp2;
      gdouble r, inv_r, inv_r3;

      /* destination - source */
      vsg_vector3d_sub (&one->vector, &other->vector, &tmp);

      r = vsg_vector3d_norm (&tmp);

      if (r > epsilon)
        {
          inv_r = 1. / r;
          inv_r3 = inv_r*inv_r*inv_r;

/*           one->accum += inv_r * other->density; */
          vsg_vector3d_scalp (&tmp, - inv_r3*other->density, &tmp1);
          vsg_vector3d_add (&one->field, &tmp1, &one->field);

/*           other->accum += inv_r * one->density; */
          vsg_vector3d_scalp (&tmp, inv_r3*one->density, &tmp2);
          vsg_vector3d_add (&other->field, &tmp2, &other->field);

/*           g_printerr ("%d : p2p %d %d ", rk, one->id, other->id); */
/*           vsg_vector3d_write (&one->vector, stderr); */
/*           g_printerr (" "); */
/*           vsg_vector3d_write (&one->field, stderr); */
/*           g_printerr ("\n"); */
        }
    }
}

static void p2m (PointAccum *particle, const VsgPRTree3dNodeInfo *dst_node,
                 AranDevelopment3d *devel)
{
  VsgVector3d tmp;
  guint deg = aran_spherical_seriesd_get_negdeg (devel->multipole);
  gint l, m;
  gcomplex128 harmonics[((deg+1)*(deg+2))/2];
  gdouble r, cost, sint, cosp, sinp;
  gcomplex128 expp;
  gdouble fact;

  vsg_vector3d_sub (&particle->vector, &dst_node->center, &tmp);

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

static void l2p (const VsgPRTree3dNodeInfo *devel_node, AranDevelopment3d *devel,
                 PointAccum *particle)
{
  VsgVector3d tmp;

  aran_development3d_local_gradient_evaluate (devel_node, devel,
                                              &particle->vector,
                                              &tmp);

  vsg_vector3d_add (&particle->field, &tmp, &particle->field);
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
  vsg_packed_msg_send_append (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_send_append (pm, &pt->density, 1, MPI_DOUBLE);
  vsg_packed_msg_send_append (pm, &pt->field, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_send_append (pm, &pt->id, 1, MPI_INT);
}

void point_accum_migrate_unpack (PointAccum *pt, VsgPackedMsg *pm,
                                 gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_recv_read (pm, &pt->density, 1, MPI_DOUBLE);
  vsg_packed_msg_recv_read (pm, &pt->field, 1, VSG_MPI_TYPE_VECTOR3D);
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
/*   vsg_vector3d_write (&pt->vector, stderr); */
/*   g_printerr ("\n"); */

  vsg_packed_msg_send_append (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_send_append (pm, &pt->density, 1, MPI_DOUBLE);
}

void point_accum_visit_fw_unpack (PointAccum *pt, VsgPackedMsg *pm,
                                  gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_recv_read (pm, &pt->density, 1, MPI_DOUBLE);
  pt->field = VSG_V3D_ZERO;
  pt->id = -1;

/*   g_printerr ("%d : fw unpack d=%e v=", */
/*               rk, */
/*               pt->density); */
/*   vsg_vector3d_write (&pt->vector, stderr); */
/*   g_printerr ("\n"); */
}

/* visit forward pack/unpack functions */
/* only the accum part is transmitted */
void point_accum_visit_bw_pack (PointAccum *pt, VsgPackedMsg *pm,
                                gpointer user_data)
{
/*   g_printerr ("%d : bw pack v=",rk); */
/*   vsg_vector3d_write (&pt->field, stderr); */
/*   g_printerr ("\n"); */
/*   fflush (stderr); */

  vsg_packed_msg_send_append (pm, &pt->field, 1, VSG_MPI_TYPE_VECTOR3D);

}

void point_accum_visit_bw_unpack (PointAccum *pt, VsgPackedMsg *pm,
                                  gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, &pt->field, 1, VSG_MPI_TYPE_VECTOR3D);

/*   g_printerr ("%d : bw unpack v=", rk); */
/*   vsg_vector3d_write (&pt->field, stderr); */
/*   g_printerr ("\n"); */
/*   fflush (stderr); */
}

void point_accum_visit_bw_reduce (PointAccum *a, PointAccum *b,
                                  gpointer user_data)
{
/*   g_printerr ("%d : bw red id=%d a=", */
/*               rk, b->id); */
/*   vsg_vector3d_write (&a->field, stderr); */
/*   g_printerr (" b="); */
/*   vsg_vector3d_write (&b->field, stderr); */

  vsg_vector3d_add (&a->field, &b->field, &b->field);

/*   g_printerr (" a+b="); */
/*   vsg_vector3d_write (&b->field, stderr); */
/*   g_printerr ("\n"); */
/*   fflush (stderr); */
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

static void _traverse_count_local_nodes (VsgPRTree3dNodeInfo *node_info,
                                         gint *count)
{
  if (! VSG_PRTREE3D_NODE_INFO_IS_REMOTE (node_info))
    {
      (*count) ++;
    }
}

static void _pt_write (PointAccum *pt, FILE *file)
{
  fprintf (file, "i=%d v=", pt->id);
  vsg_vector3d_write (&pt->vector, file);
  fprintf (file, " d=%e ", pt->density);
  vsg_vector3d_write (&pt->field, file);
  fprintf (file, "\n");
}

static void _traverse_fg_write (VsgPRTree3dNodeInfo *node_info, FILE *file)
{
  fprintf (file, "node c=");
  vsg_vector3d_write (&node_info->center, file);
  fprintf (file, " dev=");

  if (!VSG_PRTREE3D_NODE_INFO_IS_REMOTE (node_info))
    {
      aran_development3d_write ((AranDevelopment3d *) node_info->user_data,
                                file);
      fprintf (file, "\n");

      g_slist_foreach (node_info->point_list, (GFunc) _pt_write, file);
    }

  fprintf (file, "\n");
}

static void _tree_write (VsgPRTree3d *tree, gchar *prefix)
{
  gchar fn[1024];
  FILE *f;

  g_sprintf (fn, "%s%03d.txt", prefix, rk);
  f = fopen (fn, "w");

  vsg_prtree3d_traverse (tree, G_PRE_ORDER,
                         (VsgPRTree3dFunc) _traverse_fg_write,
                         f);
  fclose (f);

}

static void _vtp_pt_write (PointAccum *pt, FILE *file)
{
  fprintf (file, "%g %g %g\n",
           pt->vector.x, pt->vector.y, pt->vector.z);
}

static void _vtp_traverse_bg_write (VsgPRTree3dNodeInfo *node_info, FILE *file)
{
  if (! VSG_PRTREE3D_NODE_INFO_IS_REMOTE (node_info))
    {
      gdouble x = node_info->center.x;
      gdouble y = node_info->center.y;
      gdouble z = node_info->center.z;
      gdouble dx = node_info->ubound.x - node_info->center.x;
      gdouble dy = node_info->ubound.y - node_info->center.y;
      gdouble dz = node_info->ubound.z - node_info->center.z;

      fprintf (file, "%g %g %g\n", x - dx, y, z);
      fprintf (file, "%g %g %g\n", x + dx, y, z);

      fprintf (file, "%g %g %g\n", x, y - dy, z);
      fprintf (file, "%g %g %g\n", x, y + dy, z);

      fprintf (file, "%g %g %g\n", x, y, z - dz);
      fprintf (file, "%g %g %g\n", x, y, z + dz);
    }
}

static void _vtp_traverse_fg_write (VsgPRTree3dNodeInfo *node_info, FILE *file)
{
  g_slist_foreach (node_info->point_list, (GFunc) _vtp_pt_write, file);
}

static void _vtp_tree_write (AranSolver3d *solver, gchar *prefix)
{
  gchar fn[1024];
  FILE *f;
  gint np = aran_solver3d_point_count (solver);
  gint nn = 0;
  gint i;

  g_sprintf (fn, "%s%03d.vtp", prefix, rk);
  f = fopen (fn, "w");

  aran_solver3d_traverse (solver, G_PRE_ORDER,
                          (VsgPRTree3dFunc) _traverse_count_local_nodes,
                          &nn);

  fprintf (f, "<?xml version=\"1.0\" ?>\n"                       \
           "<VTKFile type=\"PolyData\" version=\"0.1\">\n"       \
           "<PolyData>\n"                                        \
           "<Piece NumberOfPoints=\"%d\" NumberOfVerts=\"%d\" "  \
           "NumberOfLines=\"%d\" NumberOfStrips=\"0\" "          \
           "NumberOfPolys=\"0\">\n",
           np + 6*nn, np, 3*nn);

  fprintf (f, "<Points>\n");
  fprintf (f, "<DataArray type=\"Float64\" NumberOfComponents=\"3\" "   \
           "format=\"ascii\">\n");
  aran_solver3d_traverse (solver, G_PRE_ORDER,
                          (VsgPRTree3dFunc) _vtp_traverse_fg_write,
                          f);
  aran_solver3d_traverse (solver, G_PRE_ORDER,
                          (VsgPRTree3dFunc) _vtp_traverse_bg_write,
                          f);
  fprintf (f, "</DataArray>\n</Points>\n");

  fprintf (f, "<Verts>\n");
  fprintf (f, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (i=0; i<np; i++) fprintf (f, "%d ", i);
  fprintf (f, "\n</DataArray>\n");
  fprintf (f, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" >");
  for (i=0; i<np; i++) fprintf (f, "%d ", i+1);
  fprintf (f, "\n</DataArray>\n</Verts>\n");

  fprintf (f, "<Lines>\n");
  fprintf (f, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (i=0; i<3*nn; i++) fprintf (f, "%d %d ", np + 2*i, np + 2*i+1);
  fprintf (f, "\n</DataArray>\n");
  fprintf (f, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" >");
  for (i=0; i<3*nn; i++) fprintf (f, "%d ", 2*(i+1));
  fprintf (f, "\n</DataArray>\n</Lines>\n");

  fprintf (f, "\n</Piece>\n");
  fprintf (f, "</PolyData>\n</VTKFile>\n");
  fclose (f);

}

static void _one_circle_fill (AranSolver3d *solver);
static void _random_fill (AranSolver3d *solver);

static void _random2_fill (AranSolver3d *solver);

static void _plummer_fill (AranSolver3d *solver);

static void _uvsphere_fill (AranSolver3d *solver);

static gchar *_load_file = NULL;
static void _load_fill (AranSolver3d *solver);

static void _grid_fill (AranSolver3d *solver);

static gdouble err_lim = 1.E-3;
static guint np = 12;
static guint order = 20;
static gboolean check = TRUE;
static gboolean direct = FALSE;
static guint maxbox = 1;
static gboolean _verbose = FALSE;
static gboolean _write = FALSE;
static gboolean _hilbert = FALSE;

static AranMultipole2MultipoleFunc3d m2m =
(AranMultipole2MultipoleFunc3d) aran_development3d_m2m;

static AranMultipole2LocalFunc3d m2l =
(AranMultipole2LocalFunc3d) aran_development3d_m2l;

static AranLocal2LocalFunc3d l2l =
(AranLocal2LocalFunc3d) aran_development3d_l2l;

static void (*_fill) (AranSolver3d *solver) =
_one_circle_fill;


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
              _fill = _one_circle_fill;
            }
          else if (g_ascii_strcasecmp (arg, "random2") == 0)
            {
              _fill = _random2_fill;
            }
          else if (g_ascii_strcasecmp (arg, "random") == 0)
            {
              _fill = _random_fill;
            }
          else if (g_ascii_strcasecmp (arg, "grid") == 0)
            {
              _fill = _grid_fill;
            }
          else if (g_ascii_strcasecmp (arg, "uvsphere") == 0)
            {
              _fill = _uvsphere_fill;
            }
          else if (g_ascii_strcasecmp (arg, "plummer") == 0)
            {
              _fill = _plummer_fill;
            }
          else if (g_ascii_strcasecmp (arg, "load") == 0)
            {
              _fill = _load_fill;

              iarg ++;
              arg = (iarg<argc) ? argv[iarg] : NULL;

              _load_file = g_malloc (1024*sizeof (gchar));
              sscanf (arg, "%s", _load_file);
            }
          else
            {
              g_printerr ("Invalid fill name (-dist %s)\n", arg);
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
          else if (g_ascii_strcasecmp (arg, "rotate") == 0)
            {
              m2m =
                (AranMultipole2MultipoleFunc3d) aran_development3d_m2m_rotate;

              m2l = (AranMultipole2LocalFunc3d) aran_development3d_m2l_rotate;

              l2l = (AranLocal2LocalFunc3d) aran_development3d_l2l_rotate;
            }
          else
            {
              g_printerr ("Invalid translation name (-translation %s)\n", arg);
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
      else if (g_strncasecmp (arg, "--write", 9) == 0)
        {
          _write = TRUE;
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

static void _one_circle_fill (AranSolver3d *solver)
{
  guint i;
  PointAccum *point;

  for (i=0; i<np; i++)
    {
      PointAccum tmp;

#ifdef VSG_HAVE_MPI
      if (i%_flush_interval == 0)
        {
          aran_solver3d_migrate_flush (solver);
          if (i%(_flush_interval*10) == 0)
            {
              if (_verbose && rk == 0)
                g_printerr ("%d: contiguous dist before %dth point\n", rk, i);

              aran_solver3d_distribute_contiguous_leaves (solver);
              _flush_interval *=2;
            }
        }
#endif /* VSG_HAVE_MPI */

      tmp.vector.x = R * cos(2.*G_PI*i/np);
      tmp.vector.y = R * sin(2.*G_PI*i/np);
      tmp.vector.z = 0.0;
/*       tmp.vector.z = tmp.vector.x; */
      tmp.density = 1./np;
      tmp.field = VSG_V3D_ZERO;
      tmp.id = i;

      if (check) memcpy (&check_points[i], &tmp, sizeof (PointAccum));

      if (i%sz != rk) continue;

      if (i % 10000 == 0 && _verbose)
        g_printerr ("%d: insert %dth point\n", rk, i);

      point = point_accum_alloc (TRUE, NULL);

      memcpy (point, &tmp, sizeof (PointAccum));

      if (!direct) aran_solver3d_insert_point (solver, point);
    }
#ifdef VSG_HAVE_MPI
  aran_solver3d_migrate_flush (solver);
  aran_solver3d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */
}

static void _random_fill (AranSolver3d *solver)
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
          aran_solver3d_migrate_flush (solver);
          if (i%(_flush_interval*10) == 0)
            {
              if (_verbose && rk == 0)
                g_printerr ("%d: contiguous dist before %dth point\n", rk, i);

              aran_solver3d_distribute_contiguous_leaves (solver);
              _flush_interval *=2;
            }
        }
#endif /* VSG_HAVE_MPI */

      tmp.vector.x = g_rand_double_range (rand, -R, R);;
      tmp.vector.y = g_rand_double_range (rand, -R, R);;
      tmp.vector.z = g_rand_double_range (rand, -R, R);;
      tmp.density = 1.;
      tmp.field = VSG_V3D_ZERO;
      tmp.id = i;

      if (check) memcpy (&check_points[i], &tmp, sizeof (PointAccum));

      if (i%sz != rk) continue;

      if (i % 10000 == 0 && _verbose)
        g_printerr ("%d: insert %dth point\n", rk, i);

      point = point_accum_alloc (TRUE, NULL);

      memcpy (point, &tmp, sizeof (PointAccum));

      aran_solver3d_insert_point (solver, point);
    }

#ifdef VSG_HAVE_MPI
  aran_solver3d_migrate_flush (solver);
  aran_solver3d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */

  g_rand_free (rand);
}

static void _random2_fill (AranSolver3d *solver)
{
  guint i;
  PointAccum *point;
  GRand *rand = g_rand_new_with_seed (_random_seed);

  point = point_accum_alloc (TRUE, NULL);

  for (i=0; i<np; i++)
    {

#ifdef VSG_HAVE_MPI
      if (i%_flush_interval == 0)
        {
          aran_solver3d_migrate_flush (solver);
          if (i%(_flush_interval*10) == 0)
            {
              if (_verbose && rk == 0)
                g_printerr ("%d: contiguous dist before %dth point\n", rk, i);

              aran_solver3d_distribute_contiguous_leaves (solver);
              _flush_interval *=2;
            }
        }
#endif /* VSG_HAVE_MPI */

      point->vector.x = g_rand_double_range (rand, -R, R);;
      point->vector.y = g_rand_double_range (rand, -R, R);;
      point->vector.z = g_rand_double_range (rand, -R, R);;
      point->density = 1.;
      point->field = VSG_V3D_ZERO;
      point->id = i;

      if (check) memcpy (&check_points[i], point, sizeof (PointAccum));

      if (aran_solver3d_insert_point_local (solver, point))
        {
          if (i % 10000 == 0 && _verbose)
            g_printerr ("%d: insert %dth point\n", rk, i);

          point = point_accum_alloc (TRUE, NULL);
        }

    }

  point_accum_destroy (point, TRUE, NULL);

#ifdef VSG_HAVE_MPI
  aran_solver3d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */

  g_rand_free (rand);
}

static void _plummer_fill (AranSolver3d *solver)
{
  gdouble rmax = 300.;
  guint i, real_np = 0;
  PointAccum *point;
  GRand *rand = g_rand_new_with_seed (_random_seed);

  point = point_accum_alloc (TRUE, NULL);

  if (rk == 0)
    {
      point->vector.x = rmax;
      point->vector.y = rmax;
      point->vector.z = rmax;

      aran_solver3d_insert_point (solver, point);
      aran_solver3d_remove_point (solver, point);

      point->vector.x = - rmax;
      point->vector.y = - rmax;
      point->vector.z = - rmax;

      aran_solver3d_insert_point (solver, point);
      aran_solver3d_remove_point (solver, point);

    }

  aran_solver3d_migrate_flush (solver);

  for (i=0; i<np; i++)
    {
      gdouble x1, x2, x3, r, tmp;

#ifdef VSG_HAVE_MPI
      if (i%_flush_interval == 0)
        {
          aran_solver3d_migrate_flush (solver);
          if (i%(_flush_interval*10) == 0)
            {
              if (_verbose && rk == 0)
                g_printerr ("%d: contiguous dist before %dth point\n", rk, i);

              aran_solver3d_distribute_contiguous_leaves (solver);
              _flush_interval *=2;
            }
        }
#endif /* VSG_HAVE_MPI */

      x1 = g_rand_double_range (rand, 0., 1.);
      x2 = g_rand_double_range (rand, 0., 1.);
      x3 = g_rand_double_range (rand, 0., 1.);

      r = 1./ sqrt (pow (x1, -2./3.) -1.);

      point->vector.z = R * (1. - (x2+x2)) * r;

      tmp = sqrt (r*r - point->vector.z*point->vector.z);

      point->vector.x = R * tmp * cos (2.*G_PI * x3);
      point->vector.y = R * tmp * sin (2.*G_PI * x3);

      if (fabs (point->vector.x) > rmax) continue;
      if (fabs (point->vector.y) > rmax) continue;
      if (fabs (point->vector.z) > rmax) continue;

      real_np ++;

      point->density = 1./np;
      point->field = VSG_V3D_ZERO;
      point->id = i;

      if (check) memcpy (&check_points[i], point, sizeof (PointAccum));

      if (aran_solver3d_insert_point_local (solver, point))
        {
          if (i % 10000 == 0 && _verbose)
            g_printerr ("%d: insert %dth point\n", rk, i);

          point = point_accum_alloc (TRUE, NULL);
        }
    }

  g_printerr ("%d : real_np %d\n", rk, real_np);

  point_accum_destroy (point, TRUE, NULL);

#ifdef VSG_HAVE_MPI
  aran_solver3d_migrate_flush (solver);
  aran_solver3d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */

  g_rand_free (rand);
}

static void _load_fill (AranSolver3d *solver)
{
  gdouble x, y, z, d;
  guint i;
  PointAccum *point;
  FILE *file = fopen (_load_file, "r");

  point = point_accum_alloc (TRUE, NULL);

  while (fscanf (file, "%d %lf %lf %lf %lf", &i, &x, &y, &z, &d) != EOF)
    {

/*       g_printerr ("%d : loaded point %d\n", rk, i); */

#ifdef VSG_HAVE_MPI
      if (i%_flush_interval == 0)
        {
          aran_solver3d_migrate_flush (solver);
          if (i%(_flush_interval*10) == 0)
            {
              if (_verbose && rk == 0)
                g_printerr ("%d: contiguous dist before %dth point\n", rk, i);

              aran_solver3d_distribute_contiguous_leaves (solver);
              _flush_interval *=2;
            }
        }
#endif /* VSG_HAVE_MPI */

      point->vector.x = x;
      point->vector.y = y;
      point->vector.z = z;
      point->density = d;
      point->field = VSG_V3D_ZERO;
      point->id = i;

      np ++;

      if (check)
        {
          if (np > _cp_size)
            {
              _cp_size *= 2;
              check_points = g_realloc (check_points, _cp_size);
            }

          memcpy (&check_points[i], point, sizeof (PointAccum));
        }

      if (aran_solver3d_insert_point_local (solver, point))
        {
          if (i % 10000 == 0 && _verbose)
            g_printerr ("%d: insert %dth point\n", rk, i);

          point = point_accum_alloc (TRUE, NULL);
        }

    }

  point_accum_destroy (point, TRUE, NULL);

#ifdef VSG_HAVE_MPI
  aran_solver3d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */

}

/* UV sphere fill */
static void _uvsphere_fill (AranSolver3d *solver)
{
  gint i;
  PointAccum *point;
  VsgVector3d lb, ub;
  GRand *rand = g_rand_new_with_seed (_random_seed);
  gdouble r;

  aran_solver3d_get_bounds (solver, &lb, &ub);
  r = MIN (lb.x, ub.x);

  point = point_accum_alloc (TRUE, NULL);

  for (i=0; i< np; i++)
    {
      gdouble theta, phi;

      theta = g_rand_double_range (rand, 0.01 * G_PI, 0.99 * G_PI);
      phi = g_rand_double_range (rand, 0., 2.*G_PI);

      vsg_vector3d_from_spherical (&point->vector, r, theta, phi);
      point->density = 1.;
      point->field = VSG_V3D_ZERO;
      point->id = i;

      if (check) memcpy (&check_points[i], point, sizeof (PointAccum));

      if (aran_solver3d_insert_point_local (solver, point))
        {
          if (i % 10000 == 0 && _verbose)
            g_printerr ("%d: insert %dth point\n", rk, i);

          point = point_accum_alloc (TRUE, NULL);
        }

#ifdef VSG_HAVE_MPI
      if (i%(_flush_interval*10) == 0)
        {
          if (_verbose && rk == 0)
            g_printerr ("%d: contiguous dist before %dth point\n", rk, i);
          aran_solver3d_distribute_contiguous_leaves (solver);
              _flush_interval *=2;
        }
#endif /* VSG_HAVE_MPI */
    }

  point_accum_destroy (point, TRUE, NULL);

#ifdef VSG_HAVE_MPI
  aran_solver3d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */

  g_rand_free (rand);
}

static void _grid_fill (AranSolver3d *solver)
{
/*   guint i, j, k, size, cptr = 0; */
/*   gdouble x, y, z, density; */

/*   size = floor (pow (np, 1./3.)); */

/*   np = size*size*size; */
/*   density = 1./np; */

/*   for (i=0; i<size; i ++) */
/*     { */
/*       x = R * (2.*i/(size-1.)-1.); */
/*       for (j=0; j<size; j ++) */
/*         { */
/*           y = R * (2.*j/(size-1.)-1.); */
/*           for (k=0; k<size; k ++) */
/*             { */
/*               PointAccum *point = g_malloc0 (sizeof (PointAccum)); */

/*               z = R * (2.*k/(size-1.)-1.); */

/*               point->vector.x = x; */
/*               point->vector.y = y; */
/*               point->vector.z = z; */
/*               point->density = density; */
/*               point->field = VSG_V3D_ZERO; */
/* /\*               point->accum = 0.; *\/ */
/*               point->id = cptr; */

/*               if (!direct) */
/*                 aran_solver3d_insert_point (solver, point); */

/*               points[cptr] = point; */

/*               cptr ++; */
/*             } */
/*         } */
/*     } */
}

gdouble maxerr = 0.;

void check_point_field (PointAccum *point, gint *ret)
{
  gint i;
  gdouble err;
  gdouble denom;
  PointAccum *check;
  VsgVector3d tmp;

  i = point->id;

  check = &check_points[i];

  denom = vsg_vector3d_norm (&check->field);
  vsg_vector3d_sub (&point->field, &check->field, &tmp);
  err = vsg_vector3d_norm (&tmp);
  if (denom > 0.) err /= denom;

  maxerr = MAX (maxerr, fabs (err));

  if (fabs (err) > err_lim || !finite (err))
    {
      g_printerr ("%d : Field simulation error: %d relative=(%e) "
                  "pos=(%f,%f,%f)\n"
                  "                        computed=(%f,%f,%f) "
                  "exact=(%f,%f,%f)\n",
                  rk, i, fabs(err),
                  point->vector.x, point->vector.y, point->vector.z,
                  point->field.x, point->field.y, point->field.z,
                  check->field.x, check->field.y, check->field.z);
      (*ret) ++;
    }
/*   else */
/*     { */
/*       g_printerr ("%d : Field simulation ok: %d relative=(%e) exact=(%f,%f,%f)\n", */
/*                   rk, i, fabs(err), check->field.x, check->field.y, check->field.z); */

/*     } */
}

void check_parallel_points ()
{
  VsgVector3d lbound = {-TR, -TR, -TR};
  VsgVector3d ubound = {TR, TR, TR};
  VsgPRTree3d *prtree;
  AranSolver3d *solver;
  int i;

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

  if (_hilbert)
    {
      /* configure for hilbert curve order traversal */
      aran_solver3d_set_children_order_hilbert (solver);
    }

  for (i=0; i<np; i++)
    {
      aran_solver3d_insert_point (solver, &check_points[i]);
    }

  aran_solver3d_solve (solver);
 
  aran_solver3d_free (solver);
}

int main (int argc, char **argv)
{
#ifdef VSG_HAVE_MPI
  VsgPRTreeParallelConfig pconfig = {{NULL,}};
#endif

  VsgVector3d lbound = {-TR, -TR, -TR};
  VsgVector3d ubound = {TR, TR, TR};
  VsgPRTree3d *prtree;
  AranSolver3d *solver;
  int ret = 0;
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

  aran_development3d_vtable_init (&pconfig.node_data, 0, order);
#endif

  points = g_ptr_array_new ();

  if (check)
    {
      _cp_size = MAX (np, 128);
      check_points = g_malloc0 (_cp_size * sizeof (PointAccum));
    }

  prtree =
    vsg_prtree3d_new_full (&lbound, &ubound,
                            (VsgPoint3dLocFunc) vsg_vector3d_vector3d_locfunc,
                            (VsgPoint3dDistFunc) vsg_vector3d_dist,
                            (VsgRegion3dLocFunc) NULL, maxbox);

  solver = aran_solver3d_new (prtree, ARAN_TYPE_DEVELOPMENT3D,
                              aran_development3d_new (0, order),
                              (AranZeroFunc) aran_development3d_set_zero);

#ifdef VSG_HAVE_MPI
  aran_solver3d_set_parallel (solver, &pconfig);
#endif

  aran_solver3d_set_functions (solver,
                               (AranParticle2ParticleFunc3d) p2p,
                               (AranParticle2MultipoleFunc3d) p2m,
                               m2m,
                               m2l,
                               l2l,
                               (AranLocal2ParticleFunc3d)l2p);

  if (_hilbert)
    {
      /* configure for hilbert curve order traversal */
      aran_solver3d_set_children_order_hilbert (solver);
    }

  _fill (solver);

/*   g_printerr ("ok depth = %d size = %d\n", */
/*               aran_solver3d_depth (solver), */
/*               aran_solver3d_point_count (solver)); */

  if (_verbose)
    {
      g_printerr ("%d : solve begin\n", rk);

#ifdef VSG_HAVE_MPI
      MPI_Barrier (MPI_COMM_WORLD);
#endif

      timer = g_timer_new ();
    }

  aran_solver3d_solve (solver);

  if (_verbose)
    {
#ifdef VSG_HAVE_MPI
      MPI_Barrier (MPI_COMM_WORLD);
#endif

      g_printerr ("%d : solve ok elapsed=%f seconds\n", rk,
                  g_timer_elapsed (timer, NULL));

      g_timer_destroy (timer);
    }

  if (_write)
    {
      gchar fn[1024];
      FILE *f;

      g_sprintf (fn, "tree%03d.txt", rk);
      f = fopen (fn, "w");
      vsg_prtree3d_write (prtree, f);
      fclose (f);

      _tree_write (prtree, "solv");
      _vtp_tree_write (solver, "solv");
    }

  if (check)
    {
      gint i, j;

      if (sz == 1)
        {
          for (i=0; i<np; i ++)
            {
              PointAccum *pi = &check_points[i];

              for (j=0; j<np; j ++)
                {
                  if (j != i)
                    {
                      PointAccum *pj = &check_points[j];
                      p2p_one_way (pi, pj);
                    }
                }

            }
        }
      else
        check_parallel_points ();

      aran_solver3d_foreach_point (solver, (GFunc) check_point_field, &ret);

      if (_verbose)
        g_printerr ("%d : max err = %e\n", rk, maxerr);

      g_free (check_points);
    }

  aran_solver3d_free (solver);

#ifdef VSG_HAVE_MPI
  aran_development3d_vtable_clear (&pconfig.node_data);
#endif

  g_ptr_array_free (points, TRUE);

  if (_load_file != NULL) g_free (_load_file);

#ifdef VSG_HAVE_MPI
  MPI_Finalize ();
#endif

  return ret;
}

