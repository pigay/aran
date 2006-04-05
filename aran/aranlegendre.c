/* LIBARAN - Fast Multipole Method library
 * Copyright (C) 2006 Pierre Gay
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

#include "aranlegendre.h"

#include <math.h>

#include "aranbinomial.h"

/**
 * aran_legendre_associated_multiple_get_term:
 * @l: a #guint.
 * @m: a #guint. Condition  0 <= @m <= @l must hold.
 * @result: coefficient array.
 *
 * Localizes a specific term in an associated Legendre array (typically the
 * result of a call to aran_legendre_associated_evaluate_multiple()).
 *
 * Returns: address of the P_@l^@m term in @result.
 */
gdouble *aran_legendre_associated_multiple_get_term (gint l, gint m,
						     gdouble *result)
{
  return result + (((l)*(l+1))/2+m);
}


/**
 * aran_legendre_associated_evaluate_internal:
 * @l: a #guint.
 * @m: a #guint. Condition  0 <= @m <= @l must hold.
 * @x: polynomial variable.
 * @sqrt_1_m_x2: must be equivalent to sqrt (1.-@x*@x).
 *
 * Computes the associated Legendre polynomial P_@l^@m at @x.
 *
 * Returns: value of the polynomial.
 */
gdouble aran_legendre_associated_evaluate_internal (guint l, guint m,
						    gdouble x,
						    gdouble sqrt_1_m_x2)
{
  gdouble result[((l+1)*(l+2))/2];

  aran_legendre_associated_evaluate_multiple_internal (l, x, sqrt_1_m_x2,
                                                       result);

  return *aran_legendre_associated_multiple_get_term (l, m, result);
}

/**
 * aran_legendre_associated_evaluate:
 * @l: a #guint.
 * @m: a #guint. Condition  0 <= @m <= @l must hold.
 * @x: polynomial variable.
 *
 * Computes the associated Legendre polynomial P_@l^@m at @x.
 *
 * Returns: value of the polynomial.
 */
gdouble aran_legendre_associated_evaluate (guint l, guint m, gdouble x)
{
  return aran_legendre_associated_evaluate_internal (l, m, x, sqrt (1.-x*x));
}

/**
 * aran_legendre_evaluate_multiple:
 * @l: a #guint.
 * @x: polynomial variable.
 * @result: result array. Must be at least of size @l+1.
 *
 * Computes the Legendre polynomials from P_0 to P_@l at @x.
 */
void aran_legendre_evaluate_multiple (guint l, gdouble x, gdouble *result)
{
  gint i;

  result[0] = 1.;

  if (l > 0)
    {
      result[1] = x;

      for (i=2; i<=l; i ++)
        {
          result[i] = ((i+i-1.) * x * result[i-1] - (i-1) * result[i-2]) / i;
        }
    }
}

/**
 * aran_legendre_evaluate:
 * @l: a #guint.
 * @x: polynomial variable.
 *
 * Computes the Legendre polynomial P_@l at @x.
 *
 * Returns: value of the polynomial.
 */
gdouble aran_legendre_evaluate (guint l, gdouble x)
{
  gdouble P[l+1];

  aran_legendre_evaluate_multiple (l, x, P);
  return P[l];
}

/**
 * aran_legendre_associated_evaluate_multiple_internal:
 * @l: a #guint.
 * @x: polynomial variable.
 * @sqrt_1_m_x2: must be equivalent to sqrt (1.-@x*@x).
 * @result: Array where to store the result. Must be at least of size:
 *          ((@l+1)(@l+2))/2.
 *
 * Computes the associated Legendre polynomial P_@l^m at @x for all values of
 * m from 0 to @l.
 */
void aran_legendre_associated_evaluate_multiple_internal (guint l,
							  gdouble x,
							  gdouble sqrt_1_m_x2,
							  gdouble *result)
{
  gint i, j;
  gdouble *ptr, *ptr1, *ptr2;

  ptr = aran_legendre_associated_multiple_get_term (0, 0, result);
  ptr[0] = 1.;

  if (l>0)
    {
      ptr = aran_legendre_associated_multiple_get_term (1, 0, result);
      
      ptr[0] = x;
      ptr[1] = - sqrt_1_m_x2;

      for (i=2; i<=l; i ++)
        {
          ptr = aran_legendre_associated_multiple_get_term (i, 0, result);
          ptr1 = aran_legendre_associated_multiple_get_term (i-1, 0, result);
          ptr2 = aran_legendre_associated_multiple_get_term (i-2, 0, result);

          for (j=0; j<i-1; j ++)
            {
              ptr[j] = ((i+i-1) * x * ptr1[j] - (i+j-1) * ptr2[j]) / (i-j);
            }

          ptr[i-1] = (i+i-1) * x * ptr1[i-1];
          ptr[i] =  - (i+i-1) * sqrt_1_m_x2 * ptr1[i-1];
        }
    }

}

/**
 * aran_legendre_associated_evaluate_multiple:
 * @l: a #guint.
 * @x: polynomial variable.
 * @result: Array where to store the result. Must be at least of size:
 *          ((@l+1)(@l+2))/2.
 *
 * Computes the associated Legendre polynomial P_@l^m at @x for all values of
 * m from 0 to @l.
 */
void aran_legendre_associated_evaluate_multiple (guint l,
						 gdouble x,
						 gdouble *result)
{
  aran_legendre_associated_evaluate_multiple_internal (l, x, sqrt (1.-x*x),
						       result);
}

/**
 * aran_legendre_associated_evaluate_special_internal:
 * @l: a #guint.
 * @x: polynomial variable.
 * @sqrt_1_m_x2: must be equivalent to sqrt (1.-@x*@x).
 * @legendre: Array where to store the result. Must be at least of size:
 *            ((@l+1)(@l+2))/2.
 * @special_result: special polynomial values.Must be at least of size:
 *                  ((@l+1)(@l+2))/2.
 *
 * Computes the associated Legendre polynomial P_@l^m at @x for all values of
 * m from 0 to @l. A second result table is filled with special polynomials
 * value useful in some derivatives computation.
 */
void
aran_legendre_associated_evaluate_special_internal (guint l,
                                                    gdouble x,
                                                    gdouble sqrt_1_m_x2,
                                                    gdouble *legendre,
                                                    gdouble *special_result)
{
  guint i, j;
  gdouble *ptr;
  gdouble *special, *special1, *special2;
  gdouble factor = G_MAXDOUBLE; /* infinity? */

  if (sqrt_1_m_x2 != 0.) factor = 1. / sqrt_1_m_x2;

  ptr = aran_legendre_associated_multiple_get_term (0, 0, legendre);
  special = aran_legendre_associated_multiple_get_term (0, 0, special_result);

  special[0] = ptr[0] * factor; /* spec[0][0] = leg[0][0]*factor */

  if (l >= 1)
    {
      special1 = special;
      special = aran_legendre_associated_multiple_get_term (1, 0,
                                                            special_result);

      /* spec[1][1] = leg[0][0]*factor */
      special[1] = - ptr[0];

      ptr = aran_legendre_associated_multiple_get_term (1, 0, legendre);

      /* spec[1][0] = leg[1][0]*factor */
      special[0] = ptr[0] * factor;

      for (i=2; i<=l; i ++)
        {
          special2 = special1;
          special1 = special;
          special =
            aran_legendre_associated_multiple_get_term (i, 0, special_result);

          /* spec[i][i] = -(2i-1)*leg[i-1][i-1]*fact */
          special[i] = - (i+i-1.) * ptr[i-1];

          ptr = aran_legendre_associated_multiple_get_term (i, 0, legendre);

          /* spec[i][0] = leg[i][0]*fact */
          special[0] = ptr[0] * factor;

          for (j=1; j<i-1; j ++)
            {
              /* spec[i][j] = */
              /*   ((2i-1)*x*spec[i-1][j]-(i+j-1)*spec[i-2][j])/(i-j) */
              special[j] =
                ((i+i-1.) * x * special1[j] - (i+j-1.) * special2[j]) / (i-j);
            }

          /* spec[i][i-1] = (2i-1)*x*spec[i-1][j] */
          special[i-1] = (i+i-1.)* x * special1[i-1];
        }
    }

}

