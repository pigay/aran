#include "aranlinear.h"
#include "aranfit.h"
#include "aranpolynomialfit.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include <glib/gprintf.h>

/* evaluation of successive Legendre polynomials */
static void legendrev (gint n, gdouble x, gdouble *f)
{
  gint i;

  f[0] = 1.;

  if (n>=1)
    {
      f[1] = x;

      for (i=2; i<n; i++)
        {
          f[i] = ((i+i-1.) * x * f[i-1] - (i-1.) * f[i-2]) / i;
        }
    }
}

/*
 * evaluation of Legendre polynomials roots
 * from Numerical Recipes in C. p152.
 */
static void aran_gauleg (gint n, gdouble *nodes, gdouble *weights,
                         gdouble epsilon)
{
  gint m, j, i;
  gdouble x, x1, P, P_m_1, P_m_2, Pp;

  m = (n+1)/2; /* roots are symmetric: only compute positive ones */

  for (i=1; i<=m; i ++)
    {
      x = cos (M_PI * (i-0.5)/(n+0.5)); /* ith root approximation */

      do {
        P = 1.;
        P_m_1 = 0.;

        for (j=1; j<=n; j ++)
          {
            P_m_2 = P_m_1;
            P_m_1 = P;
            P = ((j+j-1.) * x * P_m_1 - (j-1.) * P_m_2) / j;
          }

        Pp = n * (x * P - P_m_1) / (x*x - 1.); /* P derivative */

        x1 = x;

        x = x - P/Pp; /* Newton method */

      } while (fabs (x-x1) > epsilon);

      nodes[i-1] = -x;
      nodes[n-i] = x;

      weights[i-1] = 2. / ((1. - x*x) * Pp*Pp);
      weights[n-i] = weights[i-1];
    }
}

/* Linear combination of successive Legendre polynomials */

#define ARAN_LEGENDRE_COMBINATION_TERMS(xcc) \
((gdouble *) (((AranLegendreCombination *)(xcc))+1))

#define ARAN_LEGENDRE_COMBINATION_OFFSET(degree) \
  (sizeof (AranLegendreCombination) + (degree) * sizeof (gdouble))

typedef struct _AranLegendreCombination AranLegendreCombination;
struct _AranLegendreCombination
{
  guint degree;

  gdouble xmin, xmax;
  gdouble ymin, ymax;
};

#define aran_legendre_combination_new(deg) \
  aran_legendre_combination_new_full (deg, -1., 1., -1., 1.)

#define xcc_terms(xcc) ARAN_LEGENDRE_COMBINATION_TERMS (xcc)

static AranLegendreCombination *
aran_legendre_combination_new_full (guint degree,
                                    gdouble xmin, gdouble xmax,
                                    gdouble ymin, gdouble ymax)
{
  AranLegendreCombination *ret =
    g_malloc (ARAN_LEGENDRE_COMBINATION_OFFSET (degree));

  ret->degree = degree;

  ret->xmin = xmin;
  ret->xmax = xmax;
  ret->ymin = ymin;
  ret->ymax = ymax;

  return ret;
}

static void aran_legendre_combination_free (AranLegendreCombination *xcc)
{
  g_free (xcc);
}

static gdouble
aran_legendre_combination_eval (AranLegendreCombination *xcc,
                                gdouble x)
{
  gint i;
  gdouble ret = 0.;

  if (xcc == NULL) return 0.;

  x = 2. * (x - xcc->xmin) / (xcc->xmax - xcc->xmin) - 1.;

#ifndef _ARAN_NO_CLENSHAW
  {
    gdouble q0, q1, q;
    gdouble *terms = xcc_terms (xcc);

    q0 = 0.;
    q1 = 0.;

    for (i=xcc->degree-2; i >= 0; i--)
      {
        q = ((i+i+3.) * x * q0 - (i+2.) * q1 - terms[i+1]) / (i+1.);
        q1 = q0;
        q0 = q;
      }

    ret = terms[0] + q1 - x * q0;
  }

#else

  {
    gdouble leg[xcc->degree];
    gdouble *terms = xcc_terms (xcc);

    legendrev (xcc->degree, x, leg);

/*     for (i=0; i<xcc->degree; i++) */
/*       g_fprintf (stdout, "%g ", */
/*                    leg[i]); */
/*     g_fprintf (stdout, "\n"); */

    for (i=0; i<xcc->degree; i++)
          ret += terms[i] * leg[i];
  }
#endif

  ret = (ret+1.)*(xcc->ymax-xcc->ymin)*0.5+xcc->ymin;

  return ret;
}

static void
aran_legendre_combination_to_poly1d (AranLegendreCombination *xcc,
                                     AranPoly1d *ap1d)
{
  gint degree = ap1d->degree+1;
  gdouble *terms = ap1d->terms;
  guint i, j;
  gdouble lagrange_terms[degree][degree];
  gdouble *lagrange[degree];
  gdouble x[degree], w[degree];

  for (i=0; i<degree; i++)
    lagrange[i] = &lagrange_terms[i][0];

  /* Lagrange interpolation on Legendre polynomial roots */
  aran_gauleg (degree, x, w, DBL_EPSILON);

  for (i=0; i<degree; i++)
    {
      terms[i] = aran_legendre_combination_eval (xcc, x[i]);

      lagrange[i][0] = 1.;
      for (j=1; j<degree; j++)
        {
          lagrange[i][j] = lagrange[i][j-1] * x[i];
        }
    }

  aran_linear_solve (degree, lagrange, terms);
}

/**
 * aran_poly1d_fit:
 * @nsamples: number of points to fit
 * @x: abscissas
 * @f: ordinates
 * @ap1d: output polynomial
 *
 * Fits (@x,@f) point to a @ap1d polynomial.
 */
gdouble aran_poly1d_fit (gint nsamples, gdouble *x, gdouble *f,
                         AranPoly1d *ap1d)
{
  return aran_poly1d_fit_sigma (nsamples, x, f, NULL, ap1d);
}

/**
 * aran_poly1d_fit_sigma:
 * @nsamples: number of points to fit
 * @x: abscissas
 * @f: ordinates
 * @sigma: standard deviations or %NULL
 * @ap1d: output polynomial
 *
 * Fits (@x,@f) point to a @ap1d polynomial given @sigma standard deviations.
 */
gdouble aran_poly1d_fit_sigma (gint nsamples, gdouble *x, gdouble *f,
                               gdouble *sigma, AranPoly1d *ap1d)
{
  gint degree = ap1d->degree;
  gdouble u_terms[nsamples][degree+1];
  gdouble v_terms[degree+1][degree+1];
  gdouble *u[nsamples];
  gdouble *v[degree+1];
  gdouble w[degree+1];
  gdouble chisq;
  gdouble xmin = G_MAXDOUBLE, xmax = - G_MAXDOUBLE;
  gdouble ymin = G_MAXDOUBLE, ymax = - G_MAXDOUBLE;
  gint i;

  AranLegendreCombination *leg;

  /* init matrices */
  for (i=0; i<nsamples; i++)
    {
      u[i] = u_terms[i];
    }
  for (i=0; i<=degree; i++)
    {
      v[i] = v_terms[i];
    }

  /* compute bounds */
  for (i=0; i<nsamples; i++)
    {
      if (x[i] < xmin) xmin = x[i];
      if (x[i] > xmax) xmax = x[i];
      if (f[i] < ymin) ymin = f[i];
      if (f[i] > ymax) ymax = f[i];
    }

  /* scale */
  for (i=0; i<nsamples; i++)
    {
      x[i] = 2. * (x[i] - xmin) / (xmax - xmin) - 1.;
      f[i] = 2. * (f[i] - ymin) / (ymax - ymin) - 1.;
    }

  /* create a polynomial in scaled legendre form */
  leg = aran_legendre_combination_new_full (degree+1, xmin, xmax, ymin, ymax);

  /* fit the legendre polynomial to the function */
  aran_svd_fit (nsamples, degree+1, x, f, sigma,
                ARAN_LEGENDRE_COMBINATION_TERMS (leg),
                u, v, w, &chisq, legendrev);

  /* convert the polynomial to canonical form */
  aran_legendre_combination_to_poly1d (leg, ap1d);

  aran_legendre_combination_free (leg);

  return chisq;

}
