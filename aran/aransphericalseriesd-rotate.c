/* LIBARAN - Fast Multipole Method library
 * Copyright (C) 2006-2007 Pierre Gay
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
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

#include "aransphericalseriesd.h"
#include "aransphericalseriesd-private.h"

#include "aransphericalharmonic.h"

#include "aranwignerrepo.h"

#include <string.h>
#include <math.h>

/* rotate a Series buffer (multipole or local) */
static void _buffer_rotate (AranWigner * aw, guint deg,
                            gcomplex128 *expma, gcomplex128 *expmg,
                            gcomplex128 * src, gcomplex128 * dst)
{
  gint l, m1, m2;

  for (l = 0; l <= deg; l++)
    {
      for (m1 = 0; m1 <= l; m1++)
        {
          gcomplex128 sum = 0.;

          for (m2 = -l; m2 < 0; m2++)
            {
              gdouble wigner_term = PHASE (m1 + m2) *
                *aran_wigner_term (aw, l, m1, m2);
              gcomplex128 src_l_m2 = src[(l * (l + 1)) / 2 + ABS (m2)];

              sum += expmg[m1] / expma[ABS (m2)] * wigner_term *
                _sph_sym (src_l_m2, ABS (m2));
            }

          for (m2 = 0; m2 <= l; m2++)
            {
              gdouble wigner_term = *aran_wigner_term (aw, l, m2, m1);
              gcomplex128 src_l_m2 = src[(l * (l + 1)) / 2 + m2];

              sum += expma[m2] * expmg[m1] * wigner_term * src_l_m2;
            }

          dst[(l * (l + 1)) / 2 + m1] += sum;
        }
    }
}

/* inverse rotate a Series buffer (multipole or local) */
static void _buffer_rotate_inverse (AranWigner * aw, guint deg,
                                    gcomplex128 *expma, gcomplex128 *expmg,
                                    gcomplex128 * src, gcomplex128 * dst)
{
  gint l, m1, m2;

  for (l = 0; l <= deg; l++)
    {
      for (m1 = 0; m1 <= l; m1++)
        {
          gcomplex128 sum = 0.;

          for (m2 = -l; m2 < 0; m2++)
            {
              gdouble wigner_term = *aran_wigner_term (aw, l, m1, m2);
              gcomplex128 src_l_m2 = src[(l * (l + 1)) / 2 + ABS (m2)];

              sum += conj (expma[m1] / expmg[ABS (m2)]) * wigner_term *
                _sph_sym (src_l_m2, ABS (m2));
            }

          for (m2 = 0; m2 <= l; m2++)
            {
              gdouble wigner_term = PHASE (m1 + m2) *
                *aran_wigner_term (aw, l, m2, m1);
              gcomplex128 src_l_m2 = src[(l * (l + 1)) / 2 + m2];

              sum += conj (expmg[m2] * expma[m1]) * wigner_term * src_l_m2;
            }

          dst[(l * (l + 1)) / 2 + m1] += sum;
        }
    }
}

/* rotate a Series buffer (multipole or local) only by alpha + beta angles */
static void _buffer_rotate_ab (AranWigner * aw, guint deg,
                               gcomplex128 *expma,
                               gcomplex128 * src, gcomplex128 * dst)
{
  gint l, m1, m2;

  for (l = 0; l <= deg; l++)
    {
      for (m1 = 0; m1 <= l; m1++)
        {
          gcomplex128 sum = 0.;

          for (m2 = -l; m2 < 0; m2++)
            {
              gdouble wigner_term = PHASE (m1 + m2) *
                *aran_wigner_term (aw, l, m1, m2);
              gcomplex128 src_l_m2 = src[(l * (l + 1)) / 2 + ABS (m2)];

              sum += wigner_term / expma[ABS (m2)] *
                _sph_sym (src_l_m2, ABS (m2));
            }

          for (m2 = 0; m2 <= l; m2++)
            {
              gdouble wigner_term = *aran_wigner_term (aw, l, m2, m1);
              gcomplex128 src_l_m2 = src[(l * (l + 1)) / 2 + m2];

              sum += expma[m2] * wigner_term * src_l_m2;
            }

          dst[(l * (l + 1)) / 2 + m1] += sum;
        }
    }
}

/* inverse rotate a Series buffer (multipole or local) only by alpha + beta
 * angles
 */
static void _buffer_rotate_ab_inverse (AranWigner * aw, guint deg,
                                    gcomplex128 *expma,
                                    gcomplex128 * src, gcomplex128 * dst)
{
  gint l, m1, m2;

  for (l = 0; l <= deg; l++)
    {
      for (m1 = 0; m1 <= l; m1++)
        {
          gcomplex128 sum = 0.;

          for (m2 = -l; m2 < 0; m2++)
            {
              gdouble wigner_term = *aran_wigner_term (aw, l, m1, m2);
              gcomplex128 src_l_m2 = src[(l * (l + 1)) / 2 + ABS (m2)];

              sum += conj (expma[m1]) * wigner_term *
                _sph_sym (src_l_m2, ABS (m2));
            }

          for (m2 = 0; m2 <= l; m2++)
            {
              gdouble wigner_term = PHASE (m1 + m2) *
                *aran_wigner_term (aw, l, m2, m1);
              gcomplex128 src_l_m2 = src[(l * (l + 1)) / 2 + m2];

              sum += conj (expma[m1]) * wigner_term * src_l_m2;
            }

          dst[(l * (l + 1)) / 2 + m1] += sum;
        }
    }
}

/**
 * aran_spherical_seriesd_rotate:
 * @src: source #AranSphericalSeriesd.
 * @alpha: 1st Euler angle.
 * @beta: 2nd Euler angle. Condition 0 <= /theta < pi must hold.
 * @gamma: 3rd Euler angle.
 * @dst: destination #AranSphericalSeriesd.
 *
 * Computes the rotation of @src and accumulates the result into @dst. Angles
 * are in the Euler ZYZ convention.
 */
void aran_spherical_seriesd_rotate (const AranSphericalSeriesd * src,
                                    gdouble alpha, gdouble beta,
                                    gdouble gamma, AranSphericalSeriesd * dst)
{
  gint l;
  guint8 nd = MIN (dst->negdeg, src->negdeg);
  guint8 pd = MIN (dst->posdeg, src->posdeg);
  guint8 lmax = MAX (pd+1, nd);

  gcomplex128 expma[lmax];
  gcomplex128 expmg[lmax];
  gcomplex128 expa = cos (alpha) - G_I * sin (alpha);
  gcomplex128 expg = cos (gamma) - G_I * sin (gamma);

  AranWigner *aw = aran_wigner_repo_lookup (beta);

  aran_wigner_require (aw, lmax);

  expma[0] = 1.;
  expmg[0] = 1.;
  for (l = 1; l < lmax; l++)
    {
      expma[l] = expma[l - 1] * expa;
      expmg[l] = expmg[l - 1] * expg;
    }

  _buffer_rotate (aw, pd, expma, expmg,
                  _spherical_seriesd_get_pos_term (src, 0, 0),
                  _spherical_seriesd_get_pos_term (dst, 0, 0));

  if (nd > 0)
    {
      _buffer_rotate (aw, nd - 1, expma, expmg,
                      _spherical_seriesd_get_neg_term (src, 0, 0),
                      _spherical_seriesd_get_neg_term (dst, 0, 0));
    }
}

/**
 * aran_spherical_seriesd_rotate_inverse:
 * @src: source #AranSphericalSeriesd.
 * @alpha: 1st Euler angle.
 * @beta: 2nd Euler angle. Condition 0 <= /theta < pi must hold.
 * @gamma: 3rd Euler angle.
 * @dst: destination #AranSphericalSeriesd.
 *
 * Computes the inverse rotation of @src and accumulates the result into
 * @dst. Angles are in the Euler ZYZ convention. The term inverse means that
 * performing aran_spherical_seriesd_rotate_inverse to the result of a call to
 * aran_spherical_seriesd_rotate() is equivalent to identity.
 */
void aran_spherical_seriesd_rotate_inverse (const AranSphericalSeriesd * src,
                                            gdouble alpha, gdouble beta,
                                            gdouble gamma,
                                            AranSphericalSeriesd * dst)
{
  gint l;
  guint8 nd = MIN (dst->negdeg, src->negdeg);
  guint8 pd = MIN (dst->posdeg, src->posdeg);
  guint8 lmax = MAX (pd+1, nd);

  gcomplex128 expma[lmax];
  gcomplex128 expmg[lmax];
  gcomplex128 expa = cos (alpha) - G_I * sin (alpha);
  gcomplex128 expg = cos (gamma) - G_I * sin (gamma);

  AranWigner *aw = aran_wigner_repo_lookup (beta);

  aran_wigner_require (aw, lmax);

  expma[0] = 1.;
  expmg[0] = 1.;
  for (l = 1; l < lmax; l++)
    {
      expma[l] = expma[l - 1] * expa;
      expmg[l] = expmg[l - 1] * expg;
    }

  _buffer_rotate_inverse (aw, pd, expma, expmg,
                  _spherical_seriesd_get_pos_term (src, 0, 0),
                  _spherical_seriesd_get_pos_term (dst, 0, 0));

  if (nd > 0)
    {
      _buffer_rotate_inverse (aw, nd - 1, expma, expmg,
                      _spherical_seriesd_get_neg_term (src, 0, 0),
                      _spherical_seriesd_get_neg_term (dst, 0, 0));
    }
}

/**
 * aran_spherical_seriesd_to_local_rotate:
 * @src: source expansion series.
 * @xsrc: @src center.
 * @dst: destination expansion series.
 * @xdst: @dst center.
 *
 * Performs the same operation as aran_spherical_seriesd_to_local() but
 * with the "Point and Shoot" algorithm which achieves O(p^3) complexity. 
 */
void aran_spherical_seriesd_to_local_rotate (const AranSphericalSeriesd *src,
                                             const VsgVector3d *xsrc,
                                             AranSphericalSeriesd *dst,
                                             const VsgVector3d *xdst)
{
  gint l;
  guint8 nd = MAX (dst->negdeg, src->negdeg);
  guint8 pd = MAX (dst->posdeg, src->posdeg);
  guint8 lmax = MAX (pd+1, nd);

  gcomplex128 expma[lmax];
  gcomplex128 expa;

  AranWigner *aw;
  VsgVector3d dir;
  gdouble r, theta, phi;
  gdouble cost = 1.;
  AranSphericalSeriesd *rot = (AranSphericalSeriesd *)
    g_alloca (ARAN_SPHERICAL_SERIESD_SIZE (src->posdeg, src->negdeg));

  AranSphericalSeriesd *trans = (AranSphericalSeriesd *)
    g_alloca (ARAN_SPHERICAL_SERIESD_SIZE (dst->posdeg, dst->negdeg));

  /* create work AranSphericalSeriesd */
  memcpy (rot, src, sizeof (AranSphericalSeriesd));
  memcpy (trans, dst, sizeof (AranSphericalSeriesd));

  aran_spherical_seriesd_set_zero (rot);
  aran_spherical_seriesd_set_zero (trans);

  /* get translation vector */
  vsg_vector3d_sub (xdst, xsrc, &dir);

  if (dir.z < 0.)
    {
      cost = -1.;
      vsg_vector3d_scalp (&dir, -1., &dir);
    }

  /* get rotation angles */
  vsg_vector3d_to_spherical (&dir, &r, &theta, &phi);

  /* prepare rotation tools: AranWigner + complex exponentials arrays */
  aw = aran_wigner_repo_lookup (theta);

  aran_wigner_require (aw, lmax);

  expa = cos (-phi) - G_I * sin (-phi);

  expma[0] = 1.;
  for (l = 1; l < lmax; l++)
    {
      expma[l] = expma[l - 1] * expa;
    }

  /* rotate src */
  _buffer_rotate_ab (aw, src->posdeg, expma,
                     _spherical_seriesd_get_pos_term (src, 0, 0),
                     _spherical_seriesd_get_pos_term (rot, 0, 0));

  if (src->negdeg > 0)
    {
      _buffer_rotate_ab (aw, src->negdeg - 1, expma,
                         _spherical_seriesd_get_neg_term (src, 0, 0),
                         _spherical_seriesd_get_neg_term (rot, 0, 0));
    }
/*   aran_spherical_seriesd_rotate (src, -phi, theta, 0., rot); */

  /* translate src to dst center */
  aran_spherical_seriesd_multipole_to_local_vertical (rot, trans,
                                                      r, cost,
                                                      1., 0.);

  /* rotate back and accumulate into dst */
  _buffer_rotate_ab_inverse (aw, dst->posdeg, expma,
                             _spherical_seriesd_get_pos_term (trans, 0, 0),
                             _spherical_seriesd_get_pos_term (dst, 0, 0));

  if (dst->negdeg > 0)
    {
      _buffer_rotate_ab_inverse (aw, dst->negdeg - 1, expma,
                                 _spherical_seriesd_get_neg_term (trans, 0, 0),
                                 _spherical_seriesd_get_neg_term (dst, 0, 0));
    }
/*   aran_spherical_seriesd_rotate_inverse (trans, -phi, theta, 0., dst); */
}

/**
 * aran_spherical_seriesd_translate_rotate:
 * @src: source expansion series.
 * @xsrc: @src center.
 * @dst: destination expansion series.
 * @xdst: @dst center.
 *
 * Performs the same operation as aran_spherical_seriesd_translate() but
 * with the "Point and Shoot" algorithm which achieves O(p^3) complexity. 
 */
void aran_spherical_seriesd_translate_rotate (const AranSphericalSeriesd *src,
                                              const VsgVector3d *xsrc,
                                              AranSphericalSeriesd *dst,
                                              const VsgVector3d *xdst)
{
  gint l;
  guint8 nd = MAX (dst->negdeg, src->negdeg);
  guint8 pd = MAX (dst->posdeg, src->posdeg);
  guint8 lmax = MAX (pd+1, nd);

  gcomplex128 expma[lmax];
  gcomplex128 expa;

  AranWigner *aw;
  VsgVector3d dir;
  gdouble r, theta, phi;
  gdouble cost = 1.;
  AranSphericalSeriesd *rot = (AranSphericalSeriesd *)
    g_alloca (ARAN_SPHERICAL_SERIESD_SIZE (src->posdeg, src->negdeg));

  AranSphericalSeriesd *trans = (AranSphericalSeriesd *)
    g_alloca (ARAN_SPHERICAL_SERIESD_SIZE (dst->posdeg, dst->negdeg));

  /* create work AranSphericalSeriesd */
  memcpy (rot, src, sizeof (AranSphericalSeriesd));
  memcpy (trans, dst, sizeof (AranSphericalSeriesd));

  aran_spherical_seriesd_set_zero (rot);
  aran_spherical_seriesd_set_zero (trans);

  /* get translation vector */
  vsg_vector3d_sub (xdst, xsrc, &dir);

  if (dir.z < 0.)
    {
      cost = -1.;
      vsg_vector3d_scalp (&dir, -1., &dir);
    }

  /* get rotation angles */
  vsg_vector3d_to_spherical (&dir, &r, &theta, &phi);

  /* prepare rotation tools: AranWigner + complex exponentials arrays */
  aw = aran_wigner_repo_lookup (theta);

  aran_wigner_require (aw, lmax);

  expa = cos (-phi) - G_I * sin (-phi);

  expma[0] = 1.;
  for (l = 1; l < lmax; l++)
    {
      expma[l] = expma[l - 1] * expa;
    }

  /* rotate src */
  _buffer_rotate_ab (aw, src->posdeg, expma,
                     _spherical_seriesd_get_pos_term (src, 0, 0),
                     _spherical_seriesd_get_pos_term (rot, 0, 0));

  if (src->negdeg > 0)
    {
      _buffer_rotate_ab (aw, src->negdeg - 1, expma,
                      _spherical_seriesd_get_neg_term (src, 0, 0),
                      _spherical_seriesd_get_neg_term (rot, 0, 0));
    }
/*   aran_spherical_seriesd_rotate (src, -phi, theta, 0., rot); */

  /* translate src to dst center */
  aran_spherical_seriesd_translate_vertical (rot, trans,
                                             r, cost,
                                             1., 0.);

  /* rotate back and accumulate into dst */
  _buffer_rotate_ab_inverse (aw, dst->posdeg, expma,
                  _spherical_seriesd_get_pos_term (trans, 0, 0),
                  _spherical_seriesd_get_pos_term (dst, 0, 0));

  if (dst->negdeg > 0)
    {
      _buffer_rotate_ab_inverse (aw, dst->negdeg - 1, expma,
                      _spherical_seriesd_get_neg_term (trans, 0, 0),
                      _spherical_seriesd_get_neg_term (dst, 0, 0));
    }
/*   aran_spherical_seriesd_rotate_inverse (trans, -phi, theta, 0., dst); */
}
