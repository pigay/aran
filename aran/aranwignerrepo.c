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

#include "aranwignerrepo.h"

#include "aranwigner-private.h"

#include "arancomplex.h"

#include <math.h>

#define PHASE(m) (((m)%2 == 0) ? 1. : -1.)

static GTree *repo = NULL;
static gdouble epsilon = 1.e-3;

static gint _angle_compare (gdouble *a, gdouble *b)
{
  gdouble adiff = a[0]-b[0];
  gdouble bdiff = a[1]-b[1];
  gdouble gdiff = a[2]-b[2];

  if (adiff < -epsilon || bdiff < -epsilon || gdiff < -epsilon) return -1;
  if (adiff > epsilon || bdiff > epsilon || gdiff > epsilon) return 1;

  return 0;
}


/**
 * aran_wigner_repo_lookup:
 * @alpha: an angle in radians.
 * @beta: an angle in radians.
 * @gamma: an angle in radians.
 *
 * Looks for corresponding #AranWigner structure in the repository. Creates
 * it if necessary.
 *
 * Returns: the #AranWigner corresponding to a rotation of angle @beta.
 */
AranWigner *aran_wigner_repo_lookup (gdouble alpha, gdouble beta,
                                     gdouble gamma)
{
  AranWigner *ret;
  gdouble abg[3] = {alpha, beta, gamma};

  if (repo == NULL)
    {
      repo = g_tree_new_full ((GCompareDataFunc) _angle_compare,
                              NULL, NULL,
                              (GDestroyNotify) aran_wigner_free);
      g_atexit (aran_wigner_repo_forget_all);
    }

  ret = g_tree_lookup (repo, abg);

  if (ret == NULL)
    {
      ret = aran_wigner_new (alpha, beta, gamma, 0);
      g_tree_insert (repo, &(ret->alpha), ret);
    }

  /* g_printerr ("wigner repo lookup %g,%g,%g %g,%g,%g\n", alpha, beta, gamma, ret->alpha, ret->beta, ret->gamma); */

  return ret;
}

/**
 * aran_wigner_repo_steal:
 * @alpha: an angle in radians.
 * @beta: an angle in radians.
 * @gamma: an angle in radians.
 *
 * Removes the corresponding #AranWigner structure from the repository.
 *
 * Returns: the #AranWigner corresponding to a rotation of angle @beta.
 */
AranWigner *aran_wigner_repo_steal (gdouble alpha, gdouble beta,
                                    gdouble gamma)
{
  AranWigner *ret;

  g_return_val_if_fail (repo != NULL, NULL);

  ret = aran_wigner_repo_lookup (alpha, beta, gamma);

  g_tree_steal (repo, &beta);

  return ret;
}

/**
 * aran_wigner_repo_forget:
 * @alpha: an angle in radians.
 * @beta: an angle in radians.
 * @gamma: an angle in radians.
 *
 * Destroys the corresponding #AranWigner structure in the repository if it is
 * present. Returns quietly otherwise.
 */
void aran_wigner_repo_forget (gdouble alpha, gdouble beta,
                              gdouble gamma)
{
  gdouble abg[3] = {alpha, beta, gamma};

  g_return_if_fail (repo != NULL);

  g_tree_remove (repo, abg);
}

/**
 * aran_wigner_repo_forget_all:
 *
 * Destroys all #AranWigner structure of the repository. Call this if you are
 * done with Wigner rotation and you want to free some memory.
 */
void aran_wigner_repo_forget_all ()
{
  if (repo != NULL)
    {
      g_tree_destroy (repo);
      repo = NULL;
    }  
}

