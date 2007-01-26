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

#include "aransphericalharmonic.h"

#include <math.h>

#include "aranbinomial.h"
#include "aranlegendre.h"

#include "aranbinomialbufferd.h"

static AranBinomialBufferd *_spherical = NULL;

static gdouble _spherical_generator (guint l, guint m,
                                     AranBinomialBufferd *buf)
{
  if ((m == 0)) return sqrt ((l+l+1.)/(4.*G_PI));

  return *aran_binomial_bufferd_get_unsafe (buf, l, m-1) /
    sqrt ((l-m+1)*(l+m));
}

static inline gdouble *_get_term (guint l, guint m)
{
  return aran_binomial_bufferd_get_unsafe (_spherical, l, m);
}

static void _spherical_atexit ()
{
  aran_binomial_bufferd_free (_spherical);
  _spherical = NULL;
}

/**
 * aran_spherical_harmonic_require:
 * @degree: a #guint.
 *
 * Configures the software to allow computations of spherical harmonics up to
 * degree @degree.
 *
 */
void aran_spherical_harmonic_require (guint degree)
{
  if (_spherical == NULL)
    {
      _spherical = aran_binomial_bufferd_new (_spherical_generator, degree);
      g_atexit (_spherical_atexit);
    }
  else
    {
      aran_binomial_bufferd_require (_spherical, degree);
    }
}

/**
 * aran_spherical_harmonic_evaluate_internal:
 * @l: a #guint.
 * @m: a #guint.
 * @cost: cos (theta).
 * @sint: sin (theta).
 * @expmp: exp (i*phi), a #gcomplex128.
 *
 * Evaluates Y_l^m (theta, phi), given the parameters @cost, @sint, @expmp.
 *
 * Returns: value of Y_l^m (theta, phi).
 */
gcomplex128 aran_spherical_harmonic_evaluate_internal (guint l, gint m,
						       gdouble cost,
						       gdouble sint,
						       gcomplex128 expmp)
{
  gcomplex128 res;
  gdouble tmp;
  gboolean mneg = FALSE;

  aran_spherical_harmonic_require (l);

  if (m < 0)
    {
      mneg = TRUE;
      m = -m;
    }

  tmp = aran_legendre_associated_evaluate_internal (l, m, cost, sint);
  tmp *= *_get_term (l, m);

  res = tmp * expmp;

  if (mneg)
    {
      res = conj (res);
      if (m%2 == 1) res = -res;
    }

  return res;
}

/**
 * aran_spherical_harmonic_evaluate:
 * @l: a #guint.
 * @m: a #guint.
 * @theta: a #gdouble.
 * @phi: a #gdouble.
 *
 * Evaluates Y_l^m (@theta, @phi).
 *
 * Returns: value of Y_l^m (@theta, @phi).
 */
gcomplex128 aran_spherical_harmonic_evaluate (guint l, gint m,
					      gdouble theta, gdouble phi)
{
  gdouble mp = fabs (phi * m);
  gcomplex128 expmp = cos (mp) + G_I * sin (mp);

  return aran_spherical_harmonic_evaluate_internal (l, m,
						    cos (theta), sin (theta),
						    expmp);
}

/**
 * aran_spherical_harmonic_evaluate_multiple_internal:
 * @l: a #guint.
 * @cost: cos (theta).
 * @sint: sin (theta).
 * @expp: exp (i*phi).
 * @result: resulting array of #gcomplex128. Must be at least ((@l+1)*(@l+2))/2
 *          long.
 *
 * Computes the value of Y_i^j (theta, phi) for 0 <= i <= @l and -i <= j <= i.
 * Values are stored in the array @result.
 */
void aran_spherical_harmonic_evaluate_multiple_internal (guint l,
							 gdouble cost,
							 gdouble sint,
							 gcomplex128 expp,
							 gcomplex128 *result)
{
  gint i, j;
  gdouble legendre[((l+1)*(l+2))/2];
  gdouble *lptr;
  gdouble *term;
  gcomplex128 expppow[l+1];
  gcomplex128 pow = 1.;

  aran_spherical_harmonic_require (l);

  aran_legendre_associated_evaluate_multiple_internal (l, cost, sint,
						       legendre);

  lptr = legendre;
  term = _get_term (0, 0);

  for (i=0; i<=l; i ++)
    {
      expppow[i] = pow;

      for (j=0; j<=i; j ++)
	{
	  *result = expppow[j] * (*lptr) * (*term);

	  lptr ++;
	  result ++;
	  term ++;
	}
      pow *= expp;
    }
}

/**
 * aran_spherical_harmonic_evaluate_multiple:
 * @l: a #guint.
 * @theta: a #gdouble.
 * @phi: a #gdouble.
 * @result: resulting array of #gcomplex128. Must be at least ((@l+1)*(@l+2))/2
 *          long.

 *
 * Computes the value of Y_i^j (@theta, @phi) for 0 <= i <= @l and
 * -i <= j <= i. Values are stored in the array @result .
 */
void aran_spherical_harmonic_evaluate_multiple (guint l,
						gdouble theta,
						gdouble phi,
						gcomplex128 *result)
{
  gcomplex128 expp = cos (phi) + G_I * sin (phi);

  return aran_spherical_harmonic_evaluate_multiple_internal (l,
							     cos (theta),
							     sin (theta),
							     expp,
							     result);
}

/**
 * aran_spherical_harmonic_pre_gradient_multiple_internal:
 * @l: a #guint.
 * @cost: cos (theta).
 * @sint: sin (theta).
 * @expp: exp (i*phi).
 * @harmonics:  array of #gcomplex128 holding the harmonics. Must be at least
 *              ((@l+1)*(@l+2))/2 long.
 * @special: resulting array of #gcomplex128. Must be at least
 *           ((@l+1)*(@l+2))/2 long.
 *
 * Computes values for Y_i^j (theta, phi) as in
 * aran_spherical_harmonic_evaluate_multiple() along with some special values
 * returned in @special that are useful in the evaluation of the gradient of
 * #AranSphericalSeriesd
 */
void
aran_spherical_harmonic_pre_gradient_multiple_internal (guint l,
                                                        gdouble cost,
                                                        gdouble sint,
                                                        gcomplex128 expp,
                                                        gcomplex128 *harmonics,
                                                        gcomplex128 *special)
{
  gint i, j;
  gdouble legendre[((l+1)*(l+2))/2];
  gdouble special_legendre[((l+1)*(l+2))/2];
  gdouble *lptr;
  gdouble *slptr;
  gdouble *term;
  gcomplex128 expppow[l+1];
  gcomplex128 pow = 1.;

  aran_spherical_harmonic_require (l);

  aran_legendre_associated_evaluate_multiple_internal (l, cost, sint,
						       legendre);

  aran_legendre_associated_evaluate_special_internal (l, cost, sint,
                                                      legendre,
                                                      special_legendre);

  lptr = legendre;
  slptr = special_legendre;
  term = _get_term (0, 0);

  for (i=0; i<=l; i ++)
    {
      expppow[i] = pow;

      for (j=0; j<=i; j ++)
	{
          gcomplex128 tmp = expppow[j] * (*term);
	  *harmonics = tmp * (*lptr);
	  *special = tmp * (*slptr);

	  lptr ++;
	  slptr ++;
	  harmonics ++;
	  special ++;
	  term ++;
	}
      pow *= expp;
    }
}

/**
 * aran_spherical_harmonic_pre_gradient_multiple:
 * @l: a #guint.
 * @theta: a #gdouble.
 * @phi: a #gdouble.
 * @harmonics:  array of #gcomplex128 holding the harmonics. Must be at least
 *              ((@l+1)*(@l+2))/2 long.
 * @special: resulting array of #gcomplex128. Must be at least
 *           ((@l+1)*(@l+2))/2 long.
 *
 * Computes values for Y_i^j (@theta, @phi) as in
 * aran_spherical_harmonic_evaluate_multiple() along with some special values
 * returned in @special that are useful in the evaluation of the gradient of
 * #AranSphericalSeriesd
 */
void
aran_spherical_harmonic_pre_gradient_multiple (guint l,
                                               gdouble theta,
                                               gdouble phi,
                                               gcomplex128 *harmonics,
                                               gcomplex128 *special)
{
  gcomplex128 expp = cos (phi) + G_I * sin (phi);

  return aran_spherical_harmonic_pre_gradient_multiple_internal (l,
                                                                 cos (theta),
                                                                 sin (theta),
                                                                 expp,
                                                                 harmonics,
                                                                 special);
}

/**
 * aran_spherical_harmonic_multiple_get_term:
 * @l: a #guint.
 * @m: a #gint. Condition 0 <= @m <= @l must hold.
 * @result: array of #gcomplex128. Typically the previous result of a previous
 *          call to aran_spherical_harmonic_evaluate_multiple().
 *
 * Gives the address of term Y_l^m in the @result array.
 *
 * Returns: address of requested term.
 */
gcomplex128 *aran_spherical_harmonic_multiple_get_term (gint l, gint m,
						        gcomplex128 *result)
{
  return result + (((l)*(l+1))/2+m);
}
