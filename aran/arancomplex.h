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

#ifndef __ARAN_COMPLEX_H__
#define __ARAN_COMPLEX_H__

#include <aran/aranccomplex.h>

#include <glib.h>

G_BEGIN_DECLS;

/* typedefs */

/* Glib lacks some things:) */
typedef float complex gcomplex64;
typedef double complex gcomplex128;

#define G_I (I)

G_END_DECLS;

#endif /* __ARAN_COMPLEX_H__ */
