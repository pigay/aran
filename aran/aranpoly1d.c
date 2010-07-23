#include "aranpoly1d.h"

#include "glib/gprintf.h"
#include "string.h"

AranPoly1d *aran_poly1d_new (guint degree)
{
  AranPoly1d *ret = g_malloc (sizeof (AranPoly1d));

  ret->degree = degree;
  ret->terms = g_malloc0 ((degree+1) * sizeof (gdouble));

  return ret;
}

AranPoly1d *aran_poly1d_new_with_terms (int degree, const double *terms)
{
  AranPoly1d *ap1d = aran_poly1d_new (degree);

  if (degree >= 0) memcpy (ap1d->terms, terms, (degree+1) * sizeof (double));

  return ap1d;
}

void aran_poly1d_free (AranPoly1d *ap1d)
{
  if (ap1d->terms != NULL)
    {
      g_free (ap1d->terms);
      ap1d->terms = NULL;
    }
  g_free (ap1d);
}

void aran_poly1d_write (AranPoly1d *ap1d, FILE *file)
{
  gint i;
  for (i=0; i<=ap1d->degree; i++)
    {
      gdouble term = aran_poly1d_term (ap1d, i);

      g_fprintf (file, "%+g*x**%d", term, i);
    }
}

void aran_poly1d_add (AranPoly1d *one, AranPoly1d *other, AranPoly1d *ret)
{
  gint i;
  gdouble l, r;

  g_assert (ret->degree >= MAX (one->degree, other->degree));

  for (i=0; i<=ret->degree; i++)
    {
      if (i<=one->degree) l = aran_poly1d_term (one, i);
      else l = 0.;

      if (i<=other->degree) r = aran_poly1d_term (other, i);
      else r = 0.;

      aran_poly1d_term (ret, i) = l + r;
    }
}

void aran_poly1d_sub (AranPoly1d *one, AranPoly1d *other, AranPoly1d *ret)
{
  gint i;
  gdouble l, r;

  g_assert (ret->degree >= MAX (one->degree, other->degree));

  for (i=0; i<=ret->degree; i++)
    {
      if (i<=one->degree) l = aran_poly1d_term (one, i);
      else l = 0.;

      if (i<=other->degree) r = aran_poly1d_term (other, i);
      else r = 0.;

      aran_poly1d_term (ret, i) = l - r;
    }
}

void aran_poly1d_scalp (AranPoly1d *ap1d, gdouble factor, AranPoly1d *ret)
{
  gint i;

  g_assert (ret->degree >= ap1d->degree);

  for (i=0; i<=ap1d->degree; i++)
    {
      aran_poly1d_term (ret, i) = factor * aran_poly1d_term (ap1d, i);
    }

  for (i=ap1d->degree+1; i<=ret->degree; i++)
    {
      aran_poly1d_term (ret, i) = 0.;
    }
}

gdouble aran_poly1d_eval (AranPoly1d *ap1d, gdouble x)
{
  gint i;
  gdouble ret = aran_poly1d_term (ap1d, ap1d->degree);

  /* Horner scheme */
  for (i=ap1d->degree-1; i>=0; i--)
    {
      ret = ret * x +  aran_poly1d_term (ap1d, i);
    }

  return ret;
}

void aran_poly1d_write_key_file (AranPoly1d *ap1d, GKeyFile *kf,
                                 const gchar *group, const gchar *key)
{
  g_key_file_set_double_list (kf, group, key, ap1d->terms, ap1d->degree+1);
}

