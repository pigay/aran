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

#ifndef __ARAN_WIGNER_PRIVATE_H__
#define __ARAN_WIGNER_PRIVATE_H__

#include <aran/aranwigner.h>

struct _AranWigner
{
  gint lmax;
  gdouble alpha;
  gdouble beta;
  gdouble gamma;
  gcomplex128 *l_mprime_m_terms;
  gcomplex128 **l_mprime_terms;
  gcomplex128 ***l_terms;
};

#define ARAN_WIGNER_TERM(aw,l,mprime,m) (  \
  (aw)->l_terms[(l)][(mprime)] + (l) + (m) \
)

#endif /* __ARAN_WIGNER_PRIVATE_H__ */
