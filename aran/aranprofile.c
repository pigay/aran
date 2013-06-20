#include "aranprofile.h"

#include "aranpolynomialfit.h"

#include "aranprofiledb.h"

#include <glib/gprintf.h>

#define _APVAR(var) _aran_profile_ ## var

#define ARAN_PROFILE_CODE_NUMBASE(init,code,var,timebase,numbase)       \
  {                                                                     \
    GTimer *_APVAR (timer);                                             \
    gdouble _APVAR (time) = 0.;                                         \
    gdouble _APVAR (t1);                                                \
    glong _APVAR (num) = 0;                                             \
    glong _APVAR (i);                                                   \
                                                                        \
    {init}                                                              \
    {code}                                                              \
                                                                        \
    _APVAR (timer) = g_timer_new ();                                    \
                                                                        \
    for (_APVAR (i)=0;                                                  \
         _APVAR (i)<(numbase);                                          \
         _APVAR (i) ++)                                                 \
      {                                                                 \
        {init}                                                          \
        {code}                                                          \
      }                                                                 \
    _APVAR (t1) = g_timer_elapsed (_APVAR (timer), NULL);               \
    _APVAR(num) = ((timebase) * (numbase)) / _APVAR (t1);               \
    _APVAR(num) = MAX((numbase), _APVAR (num));                         \
    ARAN_RUSAGE_PROFILE_BEGIN;                                          \
                                                                        \
    for (_APVAR (i) = 0;                                                \
         _APVAR (i) < _APVAR (num);                                     \
         _APVAR (i) ++)                                                 \
      {                                                                 \
        {init}                                                          \
        {code}                                                          \
      }                                                                 \
                                                                        \
    ARAN_RUSAGE_PROFILE_END ({_APVAR (time) =                           \
          ARAN_RUSAGE_PROFILE_VAR (utime) +                             \
          ARAN_RUSAGE_PROFILE_VAR (stime);});                           \
    (var) = _APVAR(time) / (_APVAR (num));                              \
    /* g_printerr("%ld %g %g\n", _APVAR(num), _APVAR(time), g_timer_elapsed (_APVAR (timer), NULL)); */ \
    g_timer_destroy (_APVAR (timer));                                   \
  }


#define ARAN_PROFILE_CODE(init,code,var,timebase)       \
  ARAN_PROFILE_CODE_NUMBASE(init,code,var,timebase,1)


#define _TIMEBASE1 (.05)
#define _TIMEBASE1_INIT (.01)
#define _NUMBASE1 (100)

static void _nop_init (gpointer arg)
{

}

gdouble aran_profile_p2p_2d (AranParticle2ParticleFunc2d p2p,
                             AranParticleInitFunc2d init,
                             VsgPoint2 p1,
                             VsgPoint2 p2)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = (AranParticleInitFunc2d) _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (p1); init (p2);}, {
                               p2p (p1, p2);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (p1); init (p2);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_nearfunc_2d (AranParticle2ParticleFunc2d p2p,
                                  AranParticleInitFunc2d init,
                                  VsgPoint2 p1,
                                  VsgPoint2 p2,
                                  gint maxbox)
{
  GSList *l1 = g_slist_append (NULL, p1);
  GSList *l2 = g_slist_append (NULL, p2);
  gdouble t = 0., tinit = 0.;
  gint i, j;

  if (init == NULL) init = (AranParticleInitFunc2d) _nop_init;

  l1->next = l1;
  l2->next = l2;

  ARAN_PROFILE_CODE_NUMBASE ({init (p1); init (p2);}, {
      for (i=0; i<maxbox; i ++)
        {
          VsgPoint2 pp1 = l1->data;

          for (j=0; j<maxbox; j++)
            {
              VsgPoint2 pp2 = l2->data;

              p2p (pp1, pp2);

              l2 = l2->next;
            }
          l1 = l1->next;
        }
    },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (p1); init (p2);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  g_slist_free_1 (l1);
  g_slist_free_1 (l2);

  return t - tinit;
}

gdouble aran_profile_p2m_2d (AranParticle2MultipoleFunc2d p2m,
                             AranZeroFunc init,
                             VsgPoint2 p,
                             VsgPRTree2dNodeInfo *nodeinfo,
                             gpointer dst)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {
                               p2m (p, nodeinfo, dst);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_m2m_2d (AranMultipole2MultipoleFunc2d m2m,
                             AranZeroFunc init,
                             VsgPRTree2dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree2dNodeInfo *dstinfo,
                             gpointer dst)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {
                               m2m (srcinfo, src, dstinfo, dst);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_m2l_2d (AranMultipole2LocalFunc2d m2l,
                             AranZeroFunc init,
                             VsgPRTree2dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree2dNodeInfo *dstinfo,
                             gpointer dst)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {
                               m2l (srcinfo, src, dstinfo, dst);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_l2l_2d (AranLocal2LocalFunc2d l2l,
                             AranZeroFunc init,
                             VsgPRTree2dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree2dNodeInfo *dstinfo,
                             gpointer dst)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {
                               l2l (srcinfo, src, dstinfo, dst);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_l2p_2d (AranLocal2ParticleFunc2d l2p,
                             AranParticleInitFunc3d init,
                             VsgPRTree2dNodeInfo *nodeinfo,
                             gpointer dst,
                             VsgPoint2 p)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {
                               l2p (nodeinfo, dst, p);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_p2p_3d (AranParticle2ParticleFunc3d p2p,
                             AranParticleInitFunc3d init,
                             VsgPoint3 p1,
                             VsgPoint3 p2)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = (AranParticleInitFunc3d) _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (p1); init (p2);}, {
                               p2p (p1, p2);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (p1); init (p2);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  g_printerr ("p2p3d detail tinit=%g t=%g\n", tinit, t);

  return t - tinit;
}

gdouble aran_profile_nearfunc_3d (AranParticle2ParticleFunc3d p2p,
                                  AranParticleInitFunc3d init,
                                  VsgPoint2 p1,
                                  VsgPoint2 p2,
                                  gint maxbox)
{
  GSList *l1 = g_slist_append (NULL, p1);
  GSList *l2 = g_slist_append (NULL, p2);
  gdouble t = 0., tinit = 0.;
  gint i, j;

  if (init == NULL) init = (AranParticleInitFunc3d) _nop_init;

  l1->next = l1;
  l2->next = l2;

  ARAN_PROFILE_CODE_NUMBASE ({init (p1); init (p2);}, {
      for (i=0; i<maxbox; i ++)
        {
          VsgPoint2 pp1 = l1->data;

          for (j=0; j<maxbox; j++)
            {
              VsgPoint2 pp2 = l2->data;

              p2p (pp1, pp2);

              l2 = l2->next;
            }
          l1 = l1->next;
        }
    },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (p1); init (p2);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  g_slist_free_1 (l1);
  g_slist_free_1 (l2);

  return t - tinit;
}

gdouble aran_profile_p2m_3d (AranParticle2MultipoleFunc3d p2m,
                             AranZeroFunc init,
                             VsgPoint3 p,
                             VsgPRTree3dNodeInfo *nodeinfo,
                             gpointer dst)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {
                               p2m (p, nodeinfo, dst);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_m2m_3d (AranMultipole2MultipoleFunc3d m2m,
                             AranZeroFunc init,
                             VsgPRTree3dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree3dNodeInfo *dstinfo,
                             gpointer dst)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {
                               m2m (srcinfo, src, dstinfo, dst);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_m2l_3d (AranMultipole2LocalFunc3d m2l,
                             AranZeroFunc init,
                             VsgPRTree3dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree3dNodeInfo *dstinfo,
                             gpointer dst)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {
                               m2l (srcinfo, src, dstinfo, dst);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_l2l_3d (AranLocal2LocalFunc3d l2l,
                             AranZeroFunc init,
                             VsgPRTree3dNodeInfo *srcinfo,
                             gpointer src,
                             VsgPRTree3dNodeInfo *dstinfo,
                             gpointer dst)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {
                               l2l (srcinfo, src, dstinfo, dst);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

gdouble aran_profile_l2p_3d (AranLocal2ParticleFunc3d l2p,
                             AranParticleInitFunc3d init,
                             VsgPRTree3dNodeInfo *nodeinfo,
                             gpointer dst,
                             VsgPoint3 p)
{
  gdouble t = 0., tinit = 0.;

  if (init == NULL) init = _nop_init;

  ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {
                               l2p (nodeinfo, dst, p);
                             },
                             t, _TIMEBASE1, _NUMBASE1);
  ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {},
                             tinit, _TIMEBASE1_INIT, _NUMBASE1);

  return t - tinit;
}

static void _write_samples (GKeyFile *profiles_file, const gchar *profiles_group,
                            const gchar *key, gint nsamples, gdouble *abscissas,
                            gdouble *samples)
{
  gint len = 30 * nsamples;
  gchar *str = alloca (len * sizeof (gchar));
  gchar *str2 = str;
  gint i;

  for (i=0; i<nsamples; i++)
    {
      gint ret;
      ret = g_snprintf (str2, len, "%g:%g;", abscissas[i], samples[i]);
      len -= ret;
      str2 += ret;
    }

  g_key_file_set_string (profiles_file, profiles_group, key, str);
}

void aran_profile_key_file_poly1d_write (GKeyFile *profiles_file,
                                         const gchar *profiles_group,
                                         const gchar *key, AranPoly1d *ap1d,
                                         gdouble chisq, gint nsamples,
                                         gdouble *abscissas, gdouble *samples)
{
  gchar comment[1024] = {'\0', };
  gchar samples_key[1024] = {'\0', };

  aran_poly1d_write_key_file (ap1d, profiles_file, profiles_group, key);

  g_sprintf (comment, " \"%s\" function fitting: chisq=%g", key, chisq);
  g_key_file_set_comment (profiles_file, profiles_group, key, comment, NULL);

  g_sprintf(samples_key, ARAN_PROFILE_DB_SAMPLES_PREFIX "%s", key);
  _write_samples (profiles_file, profiles_group, samples_key, nsamples, abscissas, samples);
}

#define MAXDEG (50)
#define _TIMEBASE2 (.05)
#define _TIMEBASE2_INIT (.01)
#define _NUMBASE2 (100)

gdouble aran_poly1d_profile_p2m_2d (AranParticle2MultipoleFunc2d p2m,
                                    AranZeroFunc init,
                                    VsgPoint2 p,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_p2m_2d_samples (p2m, init, p, _new, _free, ap1d,
                                             nsamples, abscissas, samples);
}

gdouble aran_poly1d_profile_p2m_2d_samples (AranParticle2MultipoleFunc2d p2m,
                                            AranZeroFunc init,
                                            VsgPoint2 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree2dNodeInfo dstinfo = {.center={0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, { p2m (p, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_p2mi_2d_samples (AranParticle2MultipoleInternalFunc2d p2mi,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree2dNodeInfo dstinfo = {.center={0., 0.},};
  VsgVector2d p = {0.1, 0.1};
  gdouble charge = 1.;
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, { p2mi (&p, charge, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_p2l_2d_samples (AranParticle2LocalFunc2d p2l,
                                            AranZeroFunc init,
                                            VsgPoint2 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree2dNodeInfo dstinfo = {.center={0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new ((guint8) abscissas[i], 0);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, { p2l (p, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_p2li_2d_samples (AranParticle2LocalInternalFunc2d p2li,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree2dNodeInfo dstinfo = {.center={0., 0.},};
  VsgVector2d p = {0.1, 0.1};
  gdouble charge = 1.;
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new ((guint8) abscissas[i], 0);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, { p2li (&p, charge, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_m2m_2d (AranMultipole2MultipoleFunc2d m2m,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_m2m_2d_samples (m2m, init, _new, _free, ap1d,
                                             nsamples, abscissas, samples);
}

gdouble aran_poly1d_profile_m2m_2d_samples (AranMultipole2MultipoleFunc2d m2m,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree2dNodeInfo srcinfo = {.center={0.1, 0.1},};
  VsgPRTree2dNodeInfo dstinfo = {.center={0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer src, dst;

      src = _new (0, (guint8) abscissas[i]);
      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);},
                                 {m2m (&srcinfo, src, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (src);
      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}


gdouble aran_poly1d_profile_m2l_2d (AranMultipole2LocalFunc2d m2l,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_m2l_2d_samples (m2l, init, _new, _free, ap1d,
                                             nsamples, abscissas, samples);
}

gdouble aran_poly1d_profile_m2l_2d_samples (AranMultipole2LocalFunc2d m2l,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree2dNodeInfo srcinfo = {.center={0.1, 0.1},};
  VsgPRTree2dNodeInfo dstinfo = {.center={0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer src, dst;

      src = _new (0, (guint8) abscissas[i]);
      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);},
                                 {m2l (&srcinfo, src, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (src);
      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}


gdouble aran_poly1d_profile_l2l_2d (AranLocal2LocalFunc2d l2l,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_l2l_2d_samples (l2l, init, _new, _free, ap1d,
                                             nsamples, abscissas, samples);
}

gdouble aran_poly1d_profile_l2l_2d_samples (AranLocal2LocalFunc2d l2l,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree2dNodeInfo srcinfo = {.center={0.1, 0.1},};
  VsgPRTree2dNodeInfo dstinfo = {.center={0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer src, dst;

      src = _new (0, (guint8) abscissas[i]);
      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);},
                                 {l2l (&srcinfo, src, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (src);
      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_l2p_2d (AranLocal2ParticleFunc2d l2p,
                                    AranParticleInitFunc2d init,
                                    VsgPoint2 p,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_l2p_2d_samples (l2p, init, p, _new, _free, ap1d,
                                             nsamples, abscissas, samples);
}

gdouble aran_poly1d_profile_l2p_2d_samples (AranLocal2ParticleFunc2d l2p,
                                            AranParticleInitFunc2d init,
                                            VsgPoint2 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree2dNodeInfo dstinfo = {.center={0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {l2p (&dstinfo, dst, p);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}


gdouble aran_poly1d_profile_l2pi_2d_samples (AranLocal2ParticleInternalFunc2d l2pi,
                                            AranParticleInitFunc2d init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree2dNodeInfo dstinfo = {.center={0., 0.},};
  VsgVector2d p = {0.1, 0.1};
  gcomplex128 accum = 0.;
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit = 0, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({/* init (p); */}, {accum = l2pi (&dstinfo, dst, &p);},
                                 t, _TIMEBASE2, _NUMBASE2);

      /* ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {}, */
      /*                            tinit, _TIMEBASE2_INIT, _NUMBASE2); */

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  accum = cabs(accum);
  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}


gdouble aran_poly1d_profile_p2m_3d (AranParticle2MultipoleFunc3d p2m,
                                    AranZeroFunc init,
                                    VsgPoint3 p,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_p2m_3d_samples (p2m, init, p, _new, _free, ap1d,
                                             nsamples, abscissas, samples);

}

gdouble aran_poly1d_profile_p2m_3d_samples (AranParticle2MultipoleFunc3d p2m,
                                            AranZeroFunc init,
                                            VsgPoint3 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq = 0.;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, { p2m (p, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_p2mi_3d_samples (AranParticle2MultipoleInternalFunc3d p2mi,
                                             AranZeroFunc init,
                                             AranDevelopmentNewFunc _new,
                                             GDestroyNotify _free,
                                             AranPoly1d *ap1d,
                                             gint nsamples,
                                             gdouble *abscissas,
                                             gdouble *samples)
{
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  VsgVector3d p = {0.1, 0.1, 0.1};
  gdouble charge = 1.;
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq = 0.;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, { p2mi (&p, charge, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_p2l_3d_samples (AranParticle2LocalFunc3d p2l,
                                            AranZeroFunc init,
                                            VsgPoint3 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq = 0.;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new ((guint8) abscissas[i], 0);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, { p2l (p, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_p2li_3d_samples (AranParticle2LocalInternalFunc3d p2li,
                                             AranZeroFunc init,
                                             AranDevelopmentNewFunc _new,
                                             GDestroyNotify _free,
                                             AranPoly1d *ap1d,
                                             gint nsamples,
                                             gdouble *abscissas,
                                             gdouble *samples)
{
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  VsgVector3d p = {0.1, 0.1, 0.1};
  gdouble charge = 1.;
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq = 0.;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new ((guint8) abscissas[i], 0);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, { p2li (&p, charge, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_m2m_3d (AranMultipole2MultipoleFunc3d m2m,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_m2m_3d_samples (m2m, init, _new, _free, ap1d,
                                             nsamples, abscissas, samples);
}

gdouble aran_poly1d_profile_m2m_3d_samples (AranMultipole2MultipoleFunc3d m2m,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree3dNodeInfo srcinfo = {.center={0.1, 0.1, 0.1},};
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer src, dst;

      src = _new (0, (guint8) abscissas[i]);
      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);},
                                 {m2m (&srcinfo, src, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (src);
      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_m2l_3d (AranMultipole2LocalFunc3d m2l,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_m2l_3d_samples (m2l, init, _new, _free, ap1d,
                                             nsamples, abscissas, samples);
}

gdouble aran_poly1d_profile_m2l_3d_samples (AranMultipole2LocalFunc3d m2l,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree3dNodeInfo srcinfo = {.center={0.1, 0.1, 0.1},};
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer src, dst;

      src = _new (0, (guint8) abscissas[i]);
      dst = _new (0, (guint8) abscissas[i]);

      m2l (&srcinfo, src, &dstinfo, dst);
      ARAN_PROFILE_CODE_NUMBASE ({init (dst);},
                                 {m2l (&srcinfo, src, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (src);
      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_l2l_3d (AranLocal2LocalFunc3d l2l,
                                    AranZeroFunc init,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_l2l_3d_samples (l2l, init, _new, _free, ap1d,
                                             nsamples, abscissas, samples);

}

gdouble aran_poly1d_profile_l2l_3d_samples (AranLocal2LocalFunc3d l2l,
                                            AranZeroFunc init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree3dNodeInfo srcinfo = {.center={0.1, 0.1, 0.1},};
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  if (init == NULL) init = _nop_init;

  for (i=0; i<nsamples; i ++)
    {
      gpointer src, dst;

      src = _new (0, (guint8) abscissas[i]);
      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);},
                                 {l2l (&srcinfo, src, &dstinfo, dst);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (dst);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (src);
      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_l2p_3d (AranLocal2ParticleFunc3d l2p,
                                    AranParticleInitFunc3d init,
                                    VsgPoint3 p,
                                    AranDevelopmentNewFunc _new,
                                    GDestroyNotify _free,
                                    AranPoly1d *ap1d,
                                    gint nsamples)
{
  gint i;
  gdouble abscissas[nsamples];
  gdouble samples[nsamples];

  for (i=0; i<nsamples; i ++)
    {
      gint deg = (MAXDEG * i) / (nsamples-1);

      abscissas[i] = (gdouble) deg;
    }

  return aran_poly1d_profile_l2p_3d_samples (l2p, init, p, _new, _free, ap1d,
                                             nsamples, abscissas, samples);
}

gdouble aran_poly1d_profile_l2p_3d_samples (AranLocal2ParticleFunc3d l2p,
                                            AranParticleInitFunc3d init,
                                            VsgPoint3 p,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit, chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {l2p (&dstinfo, dst, p);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

gdouble aran_poly1d_profile_l2pi_3d_samples (AranLocal2ParticleInternalFunc3d l2pi,
                                            AranParticleInitFunc3d init,
                                            AranDevelopmentNewFunc _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples,
                                            gdouble *abscissas,
                                            gdouble *samples)
{
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  VsgVector3d p = {0.1, 0.1, 0.1};
  gcomplex128 accum = 0.;
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit = 0., chisq;
  gint i;

  if (init == NULL) init = _nop_init;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({/* init (p); */}, {accum = l2pi (&dstinfo, dst, &p);},
                                 t, _TIMEBASE2, _NUMBASE2);

      /* ARAN_PROFILE_CODE_NUMBASE ({init (p);}, {}, */
      /*                            tinit, _TIMEBASE2_INIT, _NUMBASE2); */

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      accum = cabs (accum);

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

void vsg_vector3d_set_zero (VsgVector3d *v)
{
  vsg_vector3d_set (v, 0., 0., 0.);
}

gdouble aran_poly1d_profile_l2pvi_3d_samples (AranLocal2ParticleGradInternalFunc3d l2pvi,
                                              AranParticleInitFunc3d init,
                                              AranDevelopmentNewFunc _new,
                                              GDestroyNotify _free,
                                              AranPoly1d *ap1d,
                                              gint nsamples,
                                              gdouble *abscissas,
                                              gdouble *samples)
{
  VsgPRTree3dNodeInfo dstinfo = {.center={0., 0., 0.},};
  VsgVector3d p = {0.1, 0.1, 0.1};
  VsgVector3d g;
  gdouble work_abscissas[nsamples];
  gdouble work_samples[nsamples];
  gdouble t, tinit = 0., chisq;
  gint i;

  g_return_val_if_fail (abscissas != NULL, -1.);
  if (samples == NULL) samples = work_samples;

  for (i=0; i<nsamples; i ++)
    {
      gpointer dst;

      dst = _new (0, (guint8) abscissas[i]);

      ARAN_PROFILE_CODE_NUMBASE ({vsg_vector3d_set_zero (&g);}, {l2pvi (&dstinfo, dst, &p, &g);},
                                 t, _TIMEBASE2, _NUMBASE2);

      ARAN_PROFILE_CODE_NUMBASE ({vsg_vector3d_set_zero (&g);}, {},
                                 tinit, _TIMEBASE2_INIT, _NUMBASE2);

      samples[i] = t - tinit;

      work_abscissas[i] = abscissas[i];
      work_samples[i] = samples[i];

      _free (dst);
    }

  chisq = aran_poly1d_fit (nsamples, work_abscissas, work_samples, ap1d);

  return chisq;
}

