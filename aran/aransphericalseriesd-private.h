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

#ifndef __ARAN_SPHERICAL_SERIESD_PRIVATE_H__
#define __ARAN_SPHERICAL_SERIESD_PRIVATE_H__

#include "aransphericalseriesd.h"

#define PHASE(m) (((m)%2 == 0) ? 1. : -1.)
#define ARAN_SPHERICAL_SERIESD_SIZE(pd, nd) ( \
sizeof (AranSphericalSeriesd) + \
_spherical_seriesd_size ((pd), (nd)) * \
sizeof (gcomplex128) \
)

struct _AranSphericalSeriesd
{
  guint8 posdeg;
  guint8 negdeg;
};

static inline size_t _spherical_seriesd_size (guint8 posdeg, guint8 negdeg)
{
  return ((posdeg + 1) * (posdeg + 2) + (negdeg) * (negdeg + 1)) / 2;
}

static inline
  gcomplex128 * _spherical_seriesd_get_pos_term (const AranSphericalSeriesd *
                                                 ass, guint l, gint m)
{
  gcomplex128 *res;

  g_return_val_if_fail (l <= ass->posdeg, NULL);
  g_return_val_if_fail (ABS (m) <= l, NULL);

  res = ((gcomplex128 *) (ass + 1)) + ((l * (l + 1)) / 2 + m);

  return res;
}

static inline
  gcomplex128 * _spherical_seriesd_get_neg_term (const AranSphericalSeriesd *
                                                 ass, guint l, gint m)
{
  gcomplex128 *res;

  g_return_val_if_fail (l < ass->negdeg, NULL);
  g_return_val_if_fail (ABS (m) <= l, NULL);

  res = ((gcomplex128 *) (ass + 1)) + ((ass->posdeg + 1) * (ass->posdeg + 2) +
                                       l * (l + 1)) / 2 + m;

  return res;
}

static inline
  gcomplex128 * _spherical_seriesd_get_term (const AranSphericalSeriesd * ass,
                                             gint l, gint m)
{
  guint absl;
  g_return_val_if_fail (ass != NULL, NULL);

  if (l < 0)
    {
      absl = -l - 1;

      return _spherical_seriesd_get_neg_term (ass, absl, m);
    }

  absl = l;
  return _spherical_seriesd_get_pos_term (ass, absl, m);

}

static inline gcomplex128 _sph_sym (gcomplex128 z, gint m)
{
  z = conj (z);

  if (m % 2 != 0)
    return -z;

  return z;
}

void
aran_spherical_seriesd_translate_vertical (const AranSphericalSeriesd * src,
                                           AranSphericalSeriesd * dst,
                                           gdouble r,
                                           gdouble cost,
                                           gdouble cosp, gdouble sinp);

void
aran_spherical_seriesd_multipole_to_local_vertical
(const AranSphericalSeriesd * src,
 AranSphericalSeriesd * dst,
 gdouble r,
 gdouble cost,
 gdouble cosp, gdouble sinp);

void aran_spherical_seriesd_beta_require (guint deg);
void aran_spherical_seriesd_alpha_require (guint deg);

gdouble aran_spherical_seriesd_beta (gint l);
gdouble aran_spherical_seriesd_alpha (guint n, guint p);


#endif /* __ARAN_SPHERICAL_SERIESD_PRIVATE_H__ */
