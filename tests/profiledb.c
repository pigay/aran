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

#include "aran/aran.h"
#include "aran/aranprofiledb.h"
#include "aran/arandevelopment2d.h"
#include "aran/arandevelopment3d.h"

static gchar _filename[1024] = "profiledb-test.ini";
static gchar _group[1024] = ARAN_PROFILE_DB_DEFAULT_GROUP;
static gboolean _verbose = FALSE;

void parse_args (int argc, char **argv)
{
  int iarg = 1;
  char *arg;

  while (iarg < argc)
    {
      arg = argv[iarg];

      if (g_ascii_strcasecmp (arg, "--file") == 0)
	{
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (arg == NULL || sscanf (arg, "%1024s", _filename) != 1)
            g_printerr ("Error: invalid argument for option \"--file\"\n");
	}
      else if (g_ascii_strcasecmp (arg, "--group") == 0)
	{
	  iarg ++;

	  arg = (iarg<argc) ? argv[iarg] : NULL;

	  if (arg == NULL || sscanf (arg, "%1024s", _group) != 1)
            g_printerr ("Error: invalid argument for option \"--group\"\n");
	}
      else if (g_ascii_strcasecmp (arg, "--verbose") == 0 ||
               g_ascii_strcasecmp (arg, "-v") == 0)
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

static void _check_operator_profile (gpointer op, gchar *name)
{
  AranPoly1d *p1, *p2;

  p1 = aran_profile_db_get_by_name (name);
  p2 = aran_profile_db_get_by_address (op);
  if (p1 != p2)
    g_printerr ("Error: aran_profile_db_get_by_name != " \
                "aran_profile_db_get_by_address for \"%s\"\n", name);
  else
    {
      gdouble x = aran_profile_db_name_eval (name, 1.);

      if (x != x)
        g_printerr ("Error evaluating \"%s\" by name: %g\n", name, x);
      else
        {
          x = aran_profile_db_address_eval (op, 1.);
          if (x != x)
            g_printerr ("Error evaluating \"%s\" by address: %g\n", name, x);
          else if (_verbose)
            g_printerr ("Evaluating \"%s\" ok (f(%g)=%g)\n", name, 1., x);
        }
    }
}

#define CHECK_OPERATOR_PROFILE(op) _check_operator_profile (op, #op)

int main (int argc, char **argv)
{
  int ret = 0;

  aran_init();

  parse_args (argc, argv);

  /* avoid all informational messages when running in silent mode */
  aran_profile_db_set_verbose (_verbose);

  /* load specified profiles db */
  if (! aran_profile_db_read_file (_filename, _group))
    g_printerr ("Failed to load profile db \"%s\" with group \"%s\"\n",
                _filename, _group);

  /* check some operators' rpofiles existence */
  CHECK_OPERATOR_PROFILE (aran_development2d_m2m);
  CHECK_OPERATOR_PROFILE (aran_development2d_m2l);
  CHECK_OPERATOR_PROFILE (aran_development2d_l2l);

  CHECK_OPERATOR_PROFILE (aran_development3d_m2m);
  CHECK_OPERATOR_PROFILE (aran_development3d_m2l);
  CHECK_OPERATOR_PROFILE (aran_development3d_l2l);

  /* try to evaluate a bad symbol, check if respnse is correct */
  {
    gdouble x = aran_profile_db_address_eval (&ret, 1.);

    if (x == x)
      g_printerr ("Error evaluating bad profile symbol (%g should be NaN)\n",
                  x);
    else if (_verbose)
      g_printerr ("Evaluating bad profile symbol ok (f(%g)=%g)\n", 1., x);
  }

  return ret;
}

