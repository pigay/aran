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
#include <unistd.h>

#include <complex.h>

#include <math.h>

#include <glib/gprintf.h>

#include "aran/aran.h"
#include "aran/aransolver3d.h"
#include "aran/aranbinomial.h"
#include "aran/aranprofiledb.h"


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

  guint64 id;
};

static gint rk = 0;
static gint sz = 1;

/* static GPtrArray *points = NULL; */
static PointAccum *check_points = NULL;
static guint64 _cp_size = 0;

static guint32 _random_seed = 0;
static gint _flush_interval = 1000;


static const gdouble epsilon = 1.e-5;

static void point_accum_clear_accum (PointAccum *pa)
{
  vsg_vector3d_set (&pa->field, 0., 0., 0.);
}

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

void p2p (PointAccum *one, PointAccum *other)
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

          /* g_printerr ("%d : p2p %lu %lu  \n", rk, one->id, other->id); */
/*           vsg_vector3d_write (&one->vector, stderr); */
/*           g_printerr (" "); */
/*           vsg_vector3d_write (&one->field, stderr); */
/*           g_printerr ("\n"); */
        }
    }
}

void p2m (PointAccum *particle, const VsgPRTree3dNodeInfo *dst_node,
          AranDevelopment3d *dst)
{
  aran_development3d_p2m(&particle->vector, particle->density, dst_node, dst);
}

void l2p (const VsgPRTree3dNodeInfo *devel_node, AranDevelopment3d *devel,
          PointAccum *particle)
{
  VsgVector3d tmp;

  aran_development3d_local_gradient_evaluate (devel_node, devel,
                                              &particle->vector,
                                              &tmp);

  vsg_vector3d_add (&particle->field, &tmp, &particle->field);
}

void p2l (PointAccum *particle, const VsgPRTree3dNodeInfo *dst_node,
          AranDevelopment3d *dst)
{
  /* g_printerr ("%d : p2l %lu / %#lx %#lx %#lx %d\n", rk, particle->id, dst_node->id.x, dst_node->id.y, dst_node->id.z, dst_node->id.depth); */

  aran_development3d_p2l(&particle->vector, particle->density, dst_node, dst);
}

void m2p (const VsgPRTree3dNodeInfo *devel_node, AranDevelopment3d *devel,
          PointAccum *particle)
{
  VsgVector3d tmp;

  /* g_printerr ("%d : m2p %lu \n", rk, particle->id); */

  aran_development3d_m2pv (devel_node, devel, &particle->vector, &tmp);

  vsg_vector3d_add (&particle->field, &tmp, &particle->field);

  /* g_printerr (" f="); */
  /* vsg_vector3d_write (&particle->field, stderr); */
  /* g_printerr (" \n"); */
}

#define _USE_G_SLICES GLIB_CHECK_VERSION (2, 10, 0)

PointAccum *point_accum_alloc (gboolean resident, gpointer user_data)
{
  PointAccum *ret;
#if _USE_G_SLICES
  ret = g_slice_new0 (PointAccum);
#else
  ret = g_malloc0 (sizeof (PointAccum));
#endif

  /* if (resident) */
  /*   g_ptr_array_add (points, ret); */

  return ret;
}

void point_accum_destroy (PointAccum *data, gboolean resident,
                 gpointer user_data)
{
  /* if (resident) */
  /*   g_ptr_array_remove (points, data); */

#if _USE_G_SLICES
  g_slice_free (PointAccum, data);
#else
  g_free (data);
#endif
}

#ifdef VSG_HAVE_MPI

/* migration pack/unpack functions */
void point_accum_migrate_pack (PointAccum *pt, VsgPackedMsg *pm,
                               gpointer user_data)
{
  vsg_packed_msg_send_append (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_send_append (pm, &pt->density, 1, MPI_DOUBLE);
  vsg_packed_msg_send_append (pm, &pt->field, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_send_append (pm, &pt->id, 1, MPI_UNSIGNED_LONG);
}

void point_accum_migrate_unpack (PointAccum *pt, VsgPackedMsg *pm,
                                 gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_recv_read (pm, &pt->density, 1, MPI_DOUBLE);
  vsg_packed_msg_recv_read (pm, &pt->field, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_recv_read (pm, &pt->id, 1, MPI_UNSIGNED_LONG);
}

/* visit forward pack/unpack functions */
/* only the coordinates part is transmitted */
void point_accum_visit_fw_pack (PointAccum *pt, VsgPackedMsg *pm,
                                gpointer user_data)
{
  /* g_printerr ("%d : fw pack %lu %p\n", rk, pt->id, pt); */
/*   g_printerr ("%d : fw pack d=%e v=", */
/*               rk, */
/*               pt->density); */
/*   vsg_vector3d_write (&pt->vector, stderr); */
/*   g_printerr ("\n"); */

  vsg_packed_msg_send_append (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_send_append (pm, &pt->density, 1, MPI_DOUBLE);
  vsg_packed_msg_send_append (pm, &pt->id, 1, MPI_UNSIGNED_LONG);
}

void point_accum_visit_fw_unpack (PointAccum *pt, VsgPackedMsg *pm,
                                  gpointer user_data)
{
  vsg_packed_msg_recv_read (pm, &pt->vector, 1, VSG_MPI_TYPE_VECTOR3D);
  vsg_packed_msg_recv_read (pm, &pt->density, 1, MPI_DOUBLE);
  pt->field = VSG_V3D_ZERO;
  vsg_packed_msg_recv_read (pm, &pt->id, 1, MPI_UNSIGNED_LONG);
  pt->id = -pt->id;

  /* g_printerr ("%d : fw unpack %lu\n", rk, pt->id); */
}

/* visit forward pack/unpack functions */
/* only the accum part is transmitted */
void point_accum_visit_bw_pack (PointAccum *pt, VsgPackedMsg *pm,
                                gpointer user_data)
{
  /* g_printerr ("%d : bw pack %lu\n",rk, pt->id); */

  vsg_packed_msg_send_append (pm, &pt->field, 1, VSG_MPI_TYPE_VECTOR3D);

}

void point_accum_visit_bw_unpack (PointAccum *pt, VsgPackedMsg *pm,
                                  gpointer user_data)
{

  /* g_printerr ("%d : bw unpack %lu %p\n",rk, pt->id, pt); */
  vsg_packed_msg_recv_read (pm, &pt->field, 1, VSG_MPI_TYPE_VECTOR3D);
}

void point_accum_visit_bw_reduce (PointAccum *a, PointAccum *b,
                                  gpointer user_data)
{
/*   g_printerr ("%d : bw red id=%lu a=", */
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

typedef struct _FileAndIndent FileAndIndent;
struct _FileAndIndent {
  FILE *file;
  gchar *indent;
};

static void _pt_write (PointAccum *pt, FileAndIndent *fai)
{
  fprintf (fai->file, "%si=%lu v=", fai->indent, pt->id);
  vsg_vector3d_write (&pt->vector, fai->file);
  fprintf (fai->file, " d=%e ", pt->density);
  vsg_vector3d_write (&pt->field, fai->file);
  fprintf (fai->file, "\n");
}

static void _traverse_fg_write (VsgPRTree3dNodeInfo *node_info, FILE *file)
{
  gchar *indent = alloca (node_info->depth * 2 * sizeof (gchar)+1);
  memset (indent, ' ', node_info->depth * 2 * sizeof (gchar));
  indent[node_info->depth * 2 * sizeof (gchar)] = '\0';

  fprintf (file, "%snode c=", indent);
  vsg_prtree_key3d_write (&node_info->id, file);

  if (!VSG_PRTREE3D_NODE_INFO_IS_REMOTE (node_info))
    {
      FileAndIndent fai = {file, indent};
      /* fprintf (file, " dev="); */
      /* aran_development3d_write ((AranDevelopment3d *) node_info->user_data, */
      /*                           file); */
      /* fprintf (file, "\n"); */

      if (node_info->point_count >0)
        {
          fprintf (file, "\n");
          g_slist_foreach (node_info->point_list, (GFunc) _pt_write, &fai);
        }
    }
  else
    fprintf (file, "remote");

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
  guint64 np = aran_solver3d_point_count (solver);
  gint nn = 0;
  guint64 i;

  g_sprintf (fn, "%s%03d.vtp", prefix, rk);
  f = fopen (fn, "w");

  aran_solver3d_traverse (solver, G_PRE_ORDER,
                          (VsgPRTree3dFunc) _traverse_count_local_nodes,
                          &nn);

  fprintf (f, "<?xml version=\"1.0\" ?>\n"                       \
           "<VTKFile type=\"PolyData\" version=\"0.1\">\n"       \
           "<PolyData>\n"                                        \
           "<Piece NumberOfPoints=\"%lu\" NumberOfVerts=\"%lu\" "  \
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
  for (i=0; i<np; i++) fprintf (f, "%lu ", i);
  fprintf (f, "\n</DataArray>\n");
  fprintf (f, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" >");
  for (i=0; i<np; i++) fprintf (f, "%lu ", i+1);
  fprintf (f, "\n</DataArray>\n</Verts>\n");

  fprintf (f, "<Lines>\n");
  fprintf (f, "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  for (i=0; i<3*nn; i++) fprintf (f, "%lu %lu ", np + 2*i, np + 2*i+1);
  fprintf (f, "\n</DataArray>\n");
  fprintf (f, "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\" >");
  for (i=0; i<3*nn; i++) fprintf (f, "%lu ", 2*(i+1));
  fprintf (f, "\n</DataArray>\n</Lines>\n");

  fprintf (f, "\n</Piece>\n");
  fprintf (f, "</PolyData>\n</VTKFile>\n");
  fclose (f);

}

static void _one_circle_fill (AranSolver3d *solver);
static void _random_fill (AranSolver3d *solver);

static void _random2_fill (AranSolver3d *solver);
static void _unbalanced_fill (AranSolver3d *solver);

static void _plummer_fill (AranSolver3d *solver);

static void _uvsphere_fill (AranSolver3d *solver);

static gchar *_load_file = NULL;
static void _load_fill (AranSolver3d *solver);

static void _grid_fill (AranSolver3d *solver);

static gdouble err_lim = 1.E-3;
static guint64 np = 12;
static guint order = 20;
static gboolean check = TRUE;
static gboolean direct = FALSE;
static guint maxbox = 1;
static guint virtual_maxbox = 0;
static gboolean _verbose = FALSE;
static gboolean _write = FALSE;
static gboolean _hilbert = FALSE;
static guint semifar_threshold = G_MAXUINT;

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
          guint64 tmp = 0;
          iarg ++;

          arg = (iarg<argc) ? argv[iarg] : NULL;

          if (sscanf (arg, "%lu", &tmp) == 1)
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
      else if (g_ascii_strncasecmp (arg, "-virtual-maxbox", 15) == 0)
	{
	  guint tmp = 0;
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (sscanf (arg, "%u", &tmp) == 1)
            virtual_maxbox = tmp;
	  else
	    g_printerr ("Invalid virtual maxbox (-virtual-maxbox %s)\n", arg);
	}
      else if (g_ascii_strcasecmp (arg, "-err") == 0)
        {
          gdouble tmp = 0;
          iarg ++;

          arg = (iarg<argc) ? argv[iarg] : NULL;

          if (sscanf (arg, "%lf", &tmp) == 1)
              err_lim = tmp;
          else
            g_printerr ("Invalid error limit value (-err %s)\n", arg);
        }
      else if (g_ascii_strcasecmp (arg, "-semifar") == 0)
	{
	  guint tmp = 0;
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (sscanf (arg, "%u", &tmp) == 1)
	      semifar_threshold = tmp;
	  else
	    g_printerr ("Invalid semifar threshold (-semifar %s)\n", arg);
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
          else if (g_ascii_strcasecmp (arg, "unbalanced") == 0)
            {
              _fill = _unbalanced_fill;
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
      else if (g_ascii_strncasecmp (arg, "-v", 2) == 0 ||
               g_ascii_strncasecmp (arg, "--verbose", 9) == 0)
        {
          _verbose = TRUE;
        }
      else if (g_ascii_strncasecmp (arg, "--write", 9) == 0)
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
  guint64 i;
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
                g_printerr ("%d: contiguous dist before %luth point\n", rk, i);

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
        g_printerr ("%d: insert %luth point\n", rk, i);

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
  guint64 i;
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
                g_printerr ("%d: contiguous dist before %luth point\n", rk, i);

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
        g_printerr ("%d: insert %luth point\n", rk, i);

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
  guint64 i;
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
                g_printerr ("%d: contiguous dist before %luth point\n", rk, i);

              aran_solver3d_distribute_contiguous_leaves (solver);
              _flush_interval *=2;
            }
        }
#endif /* VSG_HAVE_MPI */

      point->vector.x = g_rand_double_range (rand, -R, R);;
      point->vector.y = g_rand_double_range (rand, -R, R);;
      point->vector.z = g_rand_double_range (rand, -R, R);;
      point->density = 1./np;
      point->field = VSG_V3D_ZERO;
      point->id = i;

      /* if (i == 3) point->density = 0.; */
      /* if (i == 4) point->density = 0.; */
      /* if (i == 5) point->density = 0.; */

      if (check) memcpy (&check_points[i], point, sizeof (PointAccum));

      if (aran_solver3d_insert_point_local (solver, point))
        {
          if (i % 10000 == 0 && _verbose)
            g_printerr ("%d: insert %luth point\n", rk, i);

          point = point_accum_alloc (TRUE, NULL);
        }

    }

  point_accum_destroy (point, TRUE, NULL);

#ifdef VSG_HAVE_MPI
  aran_solver3d_migrate_flush (solver);
  aran_solver3d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */

  g_rand_free (rand);
}


static void _unbalanced_fill (AranSolver3d *solver)
{
  PointAccum *point;

  np --;

  _random2_fill (solver);

  point = point_accum_alloc (TRUE, NULL);
  point->vector.x = -2. * R;
  point->vector.y = 0.;
  point->vector.z = 0.;

  point->density = 1.;
  point->field = VSG_V3D_ZERO;
  point->id = np;

  if (check) memcpy (&check_points[np], point, sizeof (PointAccum));
  if (!aran_solver3d_insert_point_local (solver, point))
      point_accum_destroy (point, TRUE, NULL);

  np ++;

#ifdef VSG_HAVE_MPI
  aran_solver3d_migrate_flush (solver);
  aran_solver3d_distribute_contiguous_leaves (solver);
#endif /* VSG_HAVE_MPI */
}

static void _plummer_fill (AranSolver3d *solver)
{
  gdouble rmax = 300.;
  guint64 i, real_np = 0;
  PointAccum lb, ub;
  PointAccum *point;
  GRand *rand = g_rand_new_with_seed (_random_seed);

  point = point_accum_alloc (TRUE, NULL);

  ub.vector.x = rmax;
  ub.vector.y = rmax;
  ub.vector.z = rmax;

  lb.vector.x = - rmax;
  lb.vector.y = - rmax;
  lb.vector.z = - rmax;

  if (rk == 0)
    {
      aran_solver3d_insert_point (solver, &ub);
      aran_solver3d_insert_point (solver, &lb);
    }

  aran_solver3d_migrate_flush (solver);

  aran_solver3d_remove_point (solver, &lb);
  aran_solver3d_remove_point (solver, &ub);

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
                g_printerr ("%d: contiguous dist before %luth point\n", rk, i);

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

      point->density = 1./np;
      point->field = VSG_V3D_ZERO;
      point->id = real_np;

      if (check) memcpy (&check_points[real_np], point, sizeof (PointAccum));

      real_np ++;

      if (aran_solver3d_insert_point_local (solver, point))
        {
          if (i % 10000 == 0 && _verbose)
            g_printerr ("%d: insert %luth point\n", rk, i);

          point = point_accum_alloc (TRUE, NULL);
        }
    }

  g_printerr ("%d : real_np %lu\n", rk, real_np);

  np = real_np;

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
  guint64 i;
  PointAccum *point;
  FILE *file = fopen (_load_file, "r");

  point = point_accum_alloc (TRUE, NULL);

  while (fscanf (file, "%lu %lf %lf %lf %lf", &i, &x, &y, &z, &d) != EOF)
    {

/*       g_printerr ("%d : loaded point %d\n", rk, i); */

#ifdef VSG_HAVE_MPI
      if (i%_flush_interval == 0)
        {
          aran_solver3d_migrate_flush (solver);
          if (i%(_flush_interval*10) == 0)
            {
              if (_verbose && rk == 0)
                g_printerr ("%d: contiguous dist before %luth point\n", rk, i);

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
            g_printerr ("%d: insert %luth point\n", rk, i);

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
  guint64 i;
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
            g_printerr ("%d: insert %luth point\n", rk, i);

          point = point_accum_alloc (TRUE, NULL);
        }

#ifdef VSG_HAVE_MPI
      if (i%(_flush_interval*10) == 0)
        {
          if (_verbose && rk == 0)
            g_printerr ("%d: contiguous dist before %luth point\n", rk, i);
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
/*   guint64 i, j, k, size, cptr = 0; */
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
  guint64 i;
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
      g_printerr ("%d : Field simulation error: %lu relative=(%e) "
                  "pos=(%f,%f,%f)\n"
                  "                        computed=(%f,%f,%f) "
                  "exact=(%f,%f,%f)\n",
                  rk, i, fabs(err),
                  point->vector.x, point->vector.y, point->vector.z,
                  point->field.x, point->field.y, point->field.z,
                  check->field.x, check->field.y, check->field.z);
      (*ret) ++;
    }
  /* else */
  /*   { */
  /*     g_printerr ("%d : Field simulation ok: %lu relative=(%e) exact=(%f,%f,%f)\n", */
  /*                 rk, i, fabs(err), check->field.x, check->field.y, check->field.z); */

  /*   } */
}

void check_parallel_points (AranSolver3d *solver)
{
  VsgVector3d lbound;
  VsgVector3d ubound;
  VsgPRTree3d *prtree;
  AranSolver3d *solver2;
  guint64 i;

  aran_solver3d_get_bounds (solver, &lbound, &ubound);

  prtree =
    vsg_prtree3d_new_full (&lbound, &ubound,
                            (VsgPoint3dLocFunc) vsg_vector3d_vector3d_locfunc,
                            (VsgPoint3dDistFunc) vsg_vector3d_dist,
                            (VsgRegion3dLocFunc) NULL, maxbox);

  solver2 = aran_solver3d_new (prtree, ARAN_TYPE_DEVELOPMENT3D,
                              aran_development3d_new (0, order),
                              (AranZeroFunc) aran_development3d_set_zero);

  aran_solver3d_set_functions (solver2,
                               (AranParticle2ParticleFunc3d) p2p,
                               (AranParticle2MultipoleFunc3d) p2m,
                               m2m,
                               m2l,
                               l2l,
                               (AranLocal2ParticleFunc3d) l2p);


  if (semifar_threshold < G_MAXUINT)
    {
      /* if optimal threshold was requested, we need to compare with the same value */
      if (semifar_threshold == 0)
        aran_solver3d_get_functions_full (solver, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                                          &semifar_threshold);

      aran_solver3d_set_functions_full (solver2,
                                        (AranParticle2ParticleFunc3d) p2p,
                                        (AranParticle2MultipoleFunc3d) p2m,
                                        m2m,
                                        m2l,
                                        l2l,
                                        (AranLocal2ParticleFunc3d) l2p,
                                        (AranParticle2LocalFunc3d) p2l,
                                        (AranMultipole2ParticleFunc3d) m2p,
                                        semifar_threshold);
    }

  if (_hilbert)
    {
      /* configure for hilbert curve order traversal */
      aran_solver3d_set_children_order_hilbert (solver2);
    }

  for (i=0; i<np; i++)
    {
      aran_solver3d_insert_point (solver2, &check_points[i]);
    }

  aran_solver3d_solve (solver2);

  aran_solver3d_free (solver2);
}

gboolean _nf_isleaf_virtual_maxbox (const VsgPRTree3dNodeInfo *node_info,
                                    gpointer virtual_maxbox)
{
  /* shared nodes shouldn't be considered as virtual leaves in this case
   * because point_count is only a local count. For example, a shared node
   * without any local child would be considered as a virtual leaf whatever is
   * its global point_count */
  if (VSG_PRTREE3D_NODE_INFO_IS_SHARED (node_info)) return FALSE;

  return node_info->point_count <= * ((guint *) virtual_maxbox);
}

guint getpeak(pid_t pid)
{
  char statusfile[1024], line[1024], unit[12];
  FILE *f;
  guint peak;

  if (pid == 0) pid = getpid();

  sprintf(statusfile, "/proc/%d/status", pid);
  f = fopen(statusfile, "r");

  while (fscanf(f, "%[^\n]\n", line) != EOF)
    {
      if (sscanf(line, "VmPeak: %u %s", &peak, unit) > 1)
        break;
    }
  fclose (f);

  return peak;
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

  /* points = g_ptr_array_new (); */

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

  if (virtual_maxbox != 0)
    aran_solver3d_set_nf_isleaf (solver, _nf_isleaf_virtual_maxbox,
                                 &virtual_maxbox);

  aran_solver3d_set_functions (solver,
                               (AranParticle2ParticleFunc3d) p2p,
                               (AranParticle2MultipoleFunc3d) p2m,
                               m2m,
                               m2l,
                               l2l,
                               (AranLocal2ParticleFunc3d) l2p);

  if (semifar_threshold < G_MAXUINT)
    {
      aran_solver3d_set_functions_full (solver,
                                        (AranParticle2ParticleFunc3d) p2p,
                                        (AranParticle2MultipoleFunc3d) p2m,
                                        m2m,
                                        m2l,
                                        l2l,
                                        (AranLocal2ParticleFunc3d) l2p,
                                        (AranParticle2LocalFunc3d) p2l,
                                        (AranMultipole2ParticleFunc3d) m2p,
                                        semifar_threshold);

      if (semifar_threshold == 0)
        {
          PointAccum p1 = {{0.1, 0.1, 0.1}, 0.1, {0., 0., 0.}, 0};
          PointAccum p2 = {{-0.1, -0.1, -0.1}, 0.1, {0., 0., 0.}, 1};

          /* compute operators timings to be able to compute optimal solver parameters */
          aran_solver3d_profile_operators (solver, (AranParticleInitFunc3d) point_accum_clear_accum,
                                           &p1, &p2);

          /* alternatively, we could get timings from profile databases */
          /* aran_profile_db_read_file ("./profiledb-newtonfield3.ini", NULL); */
          /* aran_solver3d_db_profile_operators (solver, (gdouble) order); */

        }
      
    }

  if (_hilbert)
    {
      /* configure for hilbert curve order traversal */
      aran_solver3d_set_children_order_hilbert (solver);
    }

  if (_verbose)
    {
      g_printerr ("%d : fill begin\n", rk);
      g_printerr ("%d : memory peak1 count = %u\n", rk, getpeak(0));


#ifdef VSG_HAVE_MPI
      MPI_Barrier (MPI_COMM_WORLD);
#endif

      timer = g_timer_new ();
    }

  _fill (solver);

  if (_verbose)
    {
      g_printerr ("%d : fill elapsed=%f seconds\n", rk,
                  g_timer_elapsed (timer, NULL));

      g_printerr ("%d : tree depth = %d\n", rk,
                  aran_solver3d_depth (solver));

      g_printerr ("%d : particle count=%d\n", rk,
                  aran_solver3d_point_count (solver));

      g_timer_destroy (timer);
  /* g_mem_profile(); */
    }

  if (_verbose)
    {
      g_printerr ("%d : solve begin\n", rk);
      g_printerr ("%d : memory peak2 count = %u\n", rk, getpeak(0));


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
      g_printerr ("%d : memory peak3 count = %u\n", rk, getpeak(0));


      g_timer_destroy (timer);

      {
        glong zero_count, p2p_count, p2m_count, m2m_count;
        glong m2l_count, l2l_count, l2p_count, p2l_count, m2p_count;
        glong p2p_remote_count, m2l_remote_count;

        aran_solver3d_get_stats (solver, &zero_count,
                                 &p2p_count, &p2m_count,
                                 &m2m_count, &m2l_count,
                                 &l2l_count, &l2p_count,
                                 &p2l_count, &m2p_count,
                                 &p2p_remote_count,
                                 &m2l_remote_count);

        g_printerr ("%d : zero count=%ld\n", rk, zero_count);
        g_printerr ("%d : p2p count=%ld\n", rk, p2p_count);
        g_printerr ("%d : p2p remote count=%ld\n", rk, p2p_remote_count);
        g_printerr ("%d : p2m count=%ld\n", rk, p2m_count);
        g_printerr ("%d : m2m count=%ld\n", rk, m2m_count);
        g_printerr ("%d : m2l count=%ld\n", rk, m2l_count);
        g_printerr ("%d : m2l remote count=%ld\n", rk, m2l_remote_count);
        g_printerr ("%d : l2l count=%ld\n", rk, l2l_count);
        g_printerr ("%d : l2p count=%ld\n", rk, l2p_count);
        g_printerr ("%d : p2l count=%ld\n", rk, p2l_count);
        g_printerr ("%d : m2p count=%ld\n", rk, m2p_count);
      }
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
      guint64 i, j;

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
        check_parallel_points (solver);

      aran_solver3d_foreach_point (solver, (GFunc) check_point_field, &ret);

      if (_verbose)
        g_printerr ("%d : max err = %e\n", rk, maxerr);

      g_free (check_points);
    }

  aran_solver3d_free (solver);

#ifdef VSG_HAVE_MPI
  aran_development3d_vtable_clear (&pconfig.node_data);
#endif

  /* g_ptr_array_free (points, TRUE); */

  if (_load_file != NULL) g_free (_load_file);

#ifdef VSG_HAVE_MPI
  MPI_Finalize ();
#endif

  return ret;
}

