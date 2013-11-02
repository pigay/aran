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

#ifndef __ARAN_WIGNER_REPO_H__
#define __ARAN_WIGNER_REPO_H__

#include <glib.h>

#include <aran/aranwigner.h>

AranWigner *aran_wigner_repo_lookup (gdouble alpha, gdouble beta,
                                     gdouble gamma);
AranWigner *aran_wigner_repo_steal (gdouble alpha, gdouble beta,
                                    gdouble gamma);
void aran_wigner_repo_forget (gdouble alpha, gdouble beta,
                              gdouble gamma);
void aran_wigner_repo_forget_all ();

#endif /* __ARAN_WIGNER_REPO_H__ */
