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

#ifndef __ARAN_PROFILE_DB_H__
#define __ARAN_PROFILE_DB_H__

#include <vsg/vsgd.h>
#include <aran/aran.h>
#include <aran/aranpoly1d.h>

#define ARAN_PROFILE_DB_DEFAULT_GROUP "default"


void aran_profile_db_set_verbose (gboolean verbose);
gboolean aran_profile_db_read_file (const gchar *filename, const gchar *group);

AranPoly1d *aran_profile_db_get_by_name (gchar *name);
AranPoly1d *aran_profile_db_get_by_address (gpointer address);

void aran_profile_db_set_value (gchar *name, gpointer address,
                                AranPoly1d *value);

gdouble aran_profile_db_name_eval (gchar *name, gdouble x);
gdouble aran_profile_db_address_eval (gpointer address, gdouble x);


#endif /* __ARAN_PROFILE_DB_H__ */
