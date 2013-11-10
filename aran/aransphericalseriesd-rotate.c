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
#include "aranwigner-private.h"

/* #include <vsg/vsgd-inline.h> */

#include <string.h>
#include <math.h>

/* rotate a Series buffer (multipole or local) */
static void _buffer_rotate (AranWigner * aw, guint deg,
                            gcomplex128 * src, gcomplex128 * dst)
{
  gcomplex128 src_l[deg];
  gcomplex128 src_l_neg[deg];
  gint l, mprime, m;
  gint l_lp1_over_2 = 0;

  for (l = 0; l <= deg; l++)
    {
      l_lp1_over_2 += l;

      for (m=0; m<=l; m++)
        {
          src_l[m] = src[l_lp1_over_2 + m];
          src_l_neg[m] = _sph_sym (src_l[m], m);
        }

      for (mprime = 0; mprime <= l; mprime++)
        {
          gcomplex128 sum = 0.;

          for (m = -l; m < 0; m++)
            {
              gcomplex128 wigner_term =
                *ARAN_WIGNER_TERM (aw, l, mprime, m);

              sum += wigner_term * src_l_neg[-m];
              /* m<0
               * RY_l^{m'} += D_l^{m,m'} (-1)^m\overline{Y}_l^{-m}
               */
            }

          for (m = 0; m <= l; m++)
            {
              gcomplex128 wigner_term = *ARAN_WIGNER_TERM (aw, l, mprime, m);

              sum += wigner_term * src_l[m];
              /* m>=0
               * RY_l^{m'} += D_l^{m,m'} Y_l^{m}
               */
            }

          dst[l_lp1_over_2 + mprime] += sum;
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
  guint8 nd = MIN (dst->negdeg, src->negdeg);
  guint8 pd = MIN (dst->posdeg, src->posdeg);
  guint8 lmax = MAX (pd+1, nd);

  AranWigner *aw = aran_wigner_repo_lookup (alpha, beta, gamma);

  aran_wigner_require (aw, lmax);

  _buffer_rotate (aw, src->posdeg,
                  _spherical_seriesd_get_pos_term (src, 0, 0),
                  _spherical_seriesd_get_pos_term (dst, 0, 0));

  if (src->negdeg > 0)
    {
      _buffer_rotate (aw, src->negdeg - 1,
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
  guint8 nd = MIN (dst->negdeg, src->negdeg);
  guint8 pd = MIN (dst->posdeg, src->posdeg);
  guint8 lmax = MAX (pd+1, nd);

  AranWigner *aw = aran_wigner_repo_lookup (-gamma, -beta, -alpha);

  aran_wigner_require (aw, lmax);

  _buffer_rotate (aw, pd,
                  _spherical_seriesd_get_pos_term (src, 0, 0),
                  _spherical_seriesd_get_pos_term (dst, 0, 0));

  if (nd > 0)
    {
      _buffer_rotate (aw, nd - 1,
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
  /* AranWigner *aw; */
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

  aran_spherical_seriesd_rotate (src, -phi, theta, 0., rot);

  /* translate src to dst center */
  aran_spherical_seriesd_multipole_to_local_vertical (rot, trans,
                                                      r, cost,
                                                      1., 0.);

  aran_spherical_seriesd_rotate_inverse (trans, -phi, theta, 0., dst);
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

  /* rotate src */
  aran_spherical_seriesd_rotate (src, -phi, theta, 0., rot);
  /* aran_spherical_seriesd_rotate_ab (src, -phi, theta, 0., rot); */

  /* translate src to dst center */
  aran_spherical_seriesd_translate_vertical (rot, trans,
                                             r, cost,
                                             1., 0.);

  /* rotate back and accumulate into dst */
  aran_spherical_seriesd_rotate_inverse (trans, -phi, theta, 0., dst);
}
