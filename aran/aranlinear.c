#include "aranlinear.h"

#include <stdio.h>
#include <math.h>

#include <glib/gprintf.h>

#define SIGN(a, b) ((b<0)? -(a) : (a))

void aran_pivot_init (guint n, gint *piv)
{
  gint i;

  for (i=0; i<n; i ++)
    piv[i] = i;
}

static void aran_pivot_swap (gint *piv, gint i, gint j)
{
  gint tmp = piv[i];

  piv[i] = piv[j];
  piv[j] = tmp;
}

static void aran_pivot_search (guint n, gdouble **mat, gint *piv, gint j)
{
  gint i = j;
  gint imax = i;
  gdouble max = 0.;

  while (i<n)
    {
      gdouble term = fabs (mat[piv[i]][j]);
      if (term > max)
	{
	  max = term;
	  imax = i;
	}
      i ++;
    }

  if (j != imax)
    {
      aran_pivot_swap (piv, j, imax);
/*       g_fprintf (stderr, "pivi %d %d\n", j, imax); */
    }
}

static gdouble aran_partial_sum (gint nmax, gdouble **mat, gint *piv, gint i, 
                                 gint j)
{
  gdouble res = 0.;
  gint k;

  for (k=0; k<=nmax; k ++)
    {
      res += mat[piv[i]][k] * mat[piv[k]][j];
    }

  return res;
}

gint aran_lu_factorize (guint n, gdouble **mat, gint *piv)
{
  gint ret = TRUE;
  gint i, j;

  for (j=0; j<n; j ++)
    {
      for (i=0; i<j; i ++)
	{
	  mat[piv[i]][j] -= aran_partial_sum (i-1, mat, piv, i, j);
	}

      for (i=j; i<n; i ++)
        {
          mat[piv[i]][j] -= aran_partial_sum (j-1, mat, piv, i, j);
        }

      aran_pivot_search (n, mat, piv, j);

      if (fabs (mat[piv[j]][j]) < 1.e-10)
	{

          ret = FALSE;
	}

      for (i=j+1; i<n; i ++)
        {
          mat[piv[i]][j] /= mat[piv[j]][j];
        }

    }

  return ret;
}

void aran_lu_updown (guint n, gdouble **mat, gint *piv, gdouble *rhs)
{
  gint i, j;
  gdouble sol[n];

  for (i=0; i<n; i ++)
    {
      gdouble tmp = rhs[piv[i]];

      for (j=0; j<i; j ++)
	{
/* 	  tmp -= mat[piv[i]][j] * rhs[piv[j]]; */
	  tmp -= mat[piv[i]][j] * sol[j];
	}

      sol[i] = tmp;
    }

  for (i=n-1; i>=0; i --)
    {
      gdouble tmp = sol[i];

      for (j=i+1; j<n; j ++)
	{
/* 	  tmp -= mat[piv[i]][j] * sol[j]; */
	  tmp -= mat[piv[i]][j] * rhs[j];
	}

      rhs[i] = tmp / mat[piv[i]][i];
    }
}
gint aran_linear_solve (guint n, gdouble **mat, gdouble *vec)
{
  gint piv[n];
  
  aran_pivot_init (n, piv);

  if (!aran_lu_factorize (n, mat, piv))
    {
/*       gint i, j; */
/*       for (i=0; i<n; i++) */
/*         { */
/*           for (j=0; j<n; j++) */
/*             g_fprintf (stderr, "%g ", mat[i][j]); */
/*           g_fprintf (stderr, "\n"); */
/*         } */
      return FALSE;
    }

  aran_lu_updown (n, mat, piv, vec);

  return TRUE;
}

gdouble aran_vec_norm (guint n, gdouble *vec)
{
  gdouble result = 0.;
  gint i;

  for (i=0; i<n; i ++)
    result += fabs (vec[i]);

  return (result);
}

void aran_vec_sub (guint n, gdouble *left, gdouble *right, 
                   gdouble *result)
{
  gint i;

  for (i=0; i<n; i ++)
    result[i] = left[i] - right[i];

}

void aran_mat_vec_mult (guint n, gdouble **mat, gdouble *vec, 
                        gdouble *result)
{
  gint i, j;

  for (i=0; i<n; i ++)
    {
      gdouble tmp = 0.;
      for (j=0; j<n; j ++)
        tmp += mat[i][j] * vec[j];
      result[i] = tmp;
    }
}

gdouble aran_lu_determinant (guint n, gdouble **mat, gint *piv)
{
  gint i;
  gdouble ret = 1.;

  for (i=0; i<n; i++)
    {
      ret *= mat[piv[i]][i];
    }

  return ret;
}

void aran_mat_mat_mult (guint m, guint n, guint o, 
                        gdouble **left, gdouble **right, 
                        gdouble **result)
{
  gint i, j, k;

  for (i=0; i<m; i ++)
    {
      for (j=0; j<o; j ++)
        {
          gdouble tmp = 0.;
          for (k=0; k<n; k ++)
            tmp += left[i][k] * right[k][j];
          result[i][j] = tmp;
        }
    }
}

void aran_mat_transpose (guint n, gdouble **mat, gdouble **result)
{
  gint i, j;

  for (i=0; i<n; i ++)
    {
      for (j=i+1; j<n; j ++)
        { /* avoid problems if arguments are aliased */
          gdouble tmp = mat[j][i];
          result[j][i] = mat[i][j];
          result[i][j] = tmp;
        }
    }
}

/* SVD inspired from Numerical Recipes in C www.nr.com */

static gdouble pythag (gdouble a, gdouble b)
{
  gdouble absa, absb, frac;

  absa = fabs (a);
  absb = fabs (b);

  if (absa > absb)
    {
      frac = absb/absa;
  
      return absa*sqrt (1. + frac*frac);
    }

  frac = absa/absb;
  return (absb == 0. ? 0. : absb*sqrt (1. + frac*frac));
}

#define SVD_MAXITER 30

void aran_svd (guint m, guint n, gdouble **a, gdouble *w, 
               gdouble **v)
{
  gint flag, i, its, j, jj, k, l, nm;
  gdouble anorm, c, f, g, h, s, scale, x, y, z, tmp;
  gdouble rv1[n];

  g = 0.;
  scale = 0.;
  anorm = 0.;

  for (i=0; i<n; i++)
    {
      rv1[i] = scale*g;

      g = 0.;
      s = 0.;
      scale = 0.;

      if (i<m)
        {
          for (k=i; k<m; k++) scale += fabs (a[k][i]);

          if (scale)
            {
              for (k=i; k<m; k++)
                {
                  a[k][i] /= scale;
                  s += a[k][i]*a[k][i];
                }
              f = a[i][i];
              tmp = sqrt (s);
              g = - SIGN (tmp, f);
              h = f*g-s;
              a[i][i] = f-g;
              for (j=i+1; j<n; j++)
                {
                  s = 0.;
                  for (k=i; k<m; k++) s += a[k][i]*a[k][j];
                  f = s/h;
                  for (k=i; k<m; k++) a[k][j] += f*a[k][i];
                }
              for (k=i; k<m; k++) a[k][i] *= scale;
            }
        }

      w[i] = scale * g;

      g = 0.;
      s = 0.;
      scale = 0.;

      if (i<m && i != (n-1)) 
        {
          for (k=i+1; k<n; k++) scale += fabs (a[i][k]);
          if (scale)
            {
              for (k=i+1; k<n; k++)
                {
                  a[i][k] /= scale;
                  s += a[i][k]*a[i][k];
                }
              f = a[i][i+1];
              tmp = sqrt (s);
              g = -SIGN (tmp, f);
              h = f*g-s;
              a[i][i+1] = f-g;
              for (k=i+1; k<n; k++) rv1[k] = a[i][k]/h;
              for (j=i+1; j<m; j++)
                {
                  s = 0.;
                  for (k=i+1; k<n; k++) s += a[j][k]*a[i][k];
                  for (k=i+1; k<n; k++) a[j][k] += s*rv1[k];
                }
              for (k=i+1; k<n; k++) a[i][k] *= scale;
            }
        }
      tmp = fabs (w[i]) + fabs (rv1[i]);
      anorm = MAX (anorm, tmp);
    }

  for (i=n-1; i>=0; i--)
    {
      if (i<(n-1))
        {
          if (g)
            {
              for (j=i+1; j<n; j++)
                v[j][i] = (a[i][j]/a[i][i+1])/g;

              for (j=i+1; j<n; j++)
                {
                  s = 0.;
                  for (k=i+1; k<n; k++) s += a[i][k]*v[k][j];
                  for (k=i+1; k<n; k++) v[k][j] += s*v[k][i];
                }
            }

          for (j=i+1; j<n; j++) v[i][j] = v[j][i] = 0.;
        }
      v[i][i] = 1.;
      g = rv1[i];
    }

  for (i=MIN (m, n)-1; i>=0; i--)
    {
      g = w[i];
      for (j=i+1; j<n; j++) a[i][j] = 0.;
      if (g)
        {
          g = 1./g;
          for (j=i+1; j<n; j++)
            {
              s = 0.;
              for (k=i+1; k<m; k++) s += a[k][i]*a[k][j];
              f = (s/a[i][i])*g;
              for (k=i; k<m; k++) a[k][j] += f*a[k][i];
            }
          for (j=i; j<m; j++) a[j][i] *= g;
        } else for (j=i; j<m; j++) a[j][i] = 0.;
      a[i][i]++;
    }

  for (k=n-1; k>=0; k--)
    {
      for (its=0; its<SVD_MAXITER; its++)
        {
          flag = 1;
          for (l=k; l>=0; l--)
            {
              nm = l-1;
              if ((gdouble) (fabs (rv1[l])+anorm) == anorm)
                {
                  flag = 0;
                  break;
                }
              if ((gdouble) (fabs (w[nm])+anorm) == anorm) break;
            }

          if (flag)
            {
              c = 0.;
              s = 1.;
              for (i=l; i<=k; i++)
                {
                  f = s*rv1[i];
                  rv1[i] = c*rv1[i];
                  if ((gdouble) (fabs (f)+anorm) == anorm) break;
                  g = w[i];
                  h = pythag (f, g);
                  w[i] = h;
                  h = 1./h;
                  c = g*h;
                  s = -f*h;
                  for (j=0; j<m; j++)
                    {
                      y = a[j][nm];
                      z = a[j][i];
                      a[j][nm] = y*c+z*s;
                      a[j][i] = z*c-y*s;
                    }
                }
            }

          z = w[k];

          if (l == k)
            {
/*               g_fprintf (stderr, "convergence\n"); */
              if (z < 0.)
                {
                  w[k] = -z;
                  for (j=0;j<n;j++) v[j][k] = -v[j][k];
                }
              break;
            }

          if (its == SVD_MAXITER-1)
            {
              g_fprintf (stderr,
                         "no convergence in %d svdcmp iterations\n",
                         SVD_MAXITER);
            }
          x = w[l];
          nm = k-1;
          y = w[nm];
          g = rv1[nm];
          h = rv1[k];
          f = ((y-z)* (y+z)+ (g-h)* (g+h))/ (2.0*h*y);
          g = pythag (f, 1.);
          f = ((x-z)* (x+z)+h* ((y/ (f+SIGN (g, f)))-h))/x;
          c = 1.;
          s = 1.;
          for (j=l; j<=nm; j++)
            {
              i = j+1;
              g = rv1[i];
              y = w[i];
              h = s*g;
              g = c*g;
              z = pythag (f, h);
              rv1[j] = z;
              c = f/z;
              s = h/z;
              f = x*c+g*s;
              g = g*c-x*s;
              h = y*s;
              y *= c;
              for (jj=0; jj<n; jj++)
                {
                  x = v[jj][j];
                  z = v[jj][i];
                  v[jj][j] = x*c+z*s;
                  v[jj][i] = z*c-x*s;
                }
              z = pythag (f, h);
              w[j] = z;
              if (z)
                {
                  z = 1./z;
                  c = f*z;
                  s = h*z;
                }
              f = c*g+s*y;
              x = c*y-s*g;
              for (jj=0; jj<m;jj++)
                {
                  y = a[jj][j];
                  z = a[jj][i];
                  a[jj][j] = y*c+z*s;
                  a[jj][i] = z*c-y*s;
                }
            }
          rv1[l] = 0.;
          rv1[k] = f;
          w[k] = x;
        }
    }
}

void aran_svd_solve (guint m, guint n, gdouble **u, gdouble *w,
                     gdouble **v, gdouble *sol, gdouble *rhs)
{
  gdouble tmp[n];
  gint i, j, k;
  gdouble sum;

  for (j=0; j<n; j++)
    {
      sum = 0.;

      if (w[j] != 0.)
        {
          for (i=0; i<m; i++) sum += u[i][j] * rhs[i];
          sum /= w[j];
        }

      tmp[j] = sum;
    }

  for (j=0; j<n; j++)
    {
      sum =0.;

      for (k=0; k<n; k++) sum += v[j][k] * tmp[k];

      sol[j] = sum;
    }
}

