#include "aranfit.h"

#include <string.h>

#include <glib/gprintf.h>

#define TOLERANCE DBL_EPSILON * 100.

/**
 * aran_svd_fit:
 * @ndata: number of data to fit
 * @nfunc: number of basis functions
 * @x: abscissas
 * @y: ordinates
 * @sigma: individual standard deviations for data
 * @a: coefficients of the fitting functions (output)
 * @u: @ndata*@nfunc workspace array
 * @v: @nfunc*@nfunc workspace array
 * @w: @nfunc workspace array
 * @chisq: Chi square of the fitted functions (output)
 * @func: User supplied fitting functions
 *
 * Computes the fitting of (@x,@y) points with @func fitting functions.
 * See Numerical recipes in C (chapter 15.4).
 */
void aran_svd_fit (gint ndata, gint nfunc, gdouble *x, gdouble *y,
                   gdouble *sigma, gdouble *a, gdouble **u, gdouble **v,
                   gdouble *w, gdouble *chisq, AranFittingFunc func)
{
/* from Numerical recipes in C (chapter 15.4) */

  gint i, j;
  gdouble wmax, tmp, threshold, sum;
  gdouble b[ndata];
  gdouble afunc[nfunc];
  gdouble sigma2[ndata];

  if (sigma != NULL)
    for (i=0; i<ndata; i++)
      sigma2[i] = sigma[i];
  else
    for (i=0; i<ndata; i++)
      sigma2[i] = 1.;

  for (i=0; i<ndata; i++)
    {
      func (nfunc, x[i], afunc);
      tmp = 1. / sigma2[i];
      for (j=0; j<nfunc; j++) u[i][j] = afunc[j]*tmp;
      b[i] = y[i]*tmp;
    }

  aran_svd (ndata, nfunc, u, w, v);

  wmax = 0.;

  for (j=0; j<nfunc; j++)
    if (w[j] > wmax) wmax = w[j];

  threshold = wmax * TOLERANCE;

  for (j=0; j<nfunc; j++)
    if (w[j] <= threshold) w[j] = 0.;

  aran_svd_solve (ndata, nfunc, u, w, v, a, b);

  *chisq = 0.;

  for (i=0; i<ndata; i++)
    {
      func (nfunc, x[i], afunc);
      sum = 0.;
      for (j=0; j<nfunc; j++) sum += afunc[j]*a[j];
      tmp = (y[i]-sum) / sigma2[i];
      *chisq += tmp*tmp;
    }
}
