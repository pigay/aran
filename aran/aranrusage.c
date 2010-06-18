#include "aranrusage.h"

void aran_rusage_set (AranRusage *rusage)
{
  getrusage (rusage->who, &rusage->rusage);
}

void aran_rusage_set_rel (AranRusage *rusage)
{
  struct rusage tmp;
  
  getrusage (rusage->who, &tmp);

#define _SUB(one,other,result,field) \
(result)->field = (one)->field - (other)->field

  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_utime.tv_sec);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_utime.tv_usec);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_stime.tv_sec);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_stime.tv_usec);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_maxrss);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_ixrss);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_idrss);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_minflt);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_majflt);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_nswap);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_inblock);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_oublock);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_msgsnd);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_msgrcv);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_nsignals);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_nvcsw);
  _SUB (&tmp, &rusage->rusage, &rusage->rusage, ru_nivcsw);

#undef _SUB
}

time_t aran_rusage_get_usec (AranRusage *rusage)
{
  return rusage->rusage.ru_utime.tv_sec;
}

time_t aran_rusage_get_uusec (AranRusage *rusage)
{
  return rusage->rusage.ru_utime.tv_usec;
}

time_t aran_rusage_get_ssec (AranRusage *rusage)
{
  return rusage->rusage.ru_stime.tv_sec;
}

time_t aran_rusage_get_susec (AranRusage *rusage)
{
  return rusage->rusage.ru_stime.tv_usec;
}

gdouble aran_rusage_get_utime (AranRusage *rusage)
{
  return rusage->rusage.ru_utime.tv_sec + 1.e-6 * rusage->rusage.ru_utime.tv_usec;
}
gdouble aran_rusage_get_stime (AranRusage *rusage)
{
  return rusage->rusage.ru_stime.tv_sec + 1.e-6 * rusage->rusage.ru_stime.tv_usec;
}
glong aran_rusage_get_maxrss (AranRusage *rusage)
{
  return rusage->rusage.ru_maxrss;
}

glong aran_rusage_get_ixrss (AranRusage *rusage)
{
  return rusage->rusage.ru_ixrss;
}

glong aran_rusage_get_idrss (AranRusage *rusage)
{
  return rusage->rusage.ru_idrss;
}

glong aran_rusage_get_isrss (AranRusage *rusage)
{
  return rusage->rusage.ru_isrss;
}

glong aran_rusage_get_minflt (AranRusage *rusage)
{
  return rusage->rusage.ru_minflt;
}

glong aran_rusage_get_majflt (AranRusage *rusage)
{
  return rusage->rusage.ru_majflt;
}

glong aran_rusage_get_nswap (AranRusage *rusage)
{
  return rusage->rusage.ru_nswap;
}

glong aran_rusage_get_inblock (AranRusage *rusage)
{
  return rusage->rusage.ru_inblock;
}

glong aran_rusage_get_oublock (AranRusage *rusage)
{
  return rusage->rusage.ru_oublock;
}

glong aran_rusage_get_msgsnd (AranRusage *rusage)
{
  return rusage->rusage.ru_msgsnd;
}

glong aran_rusage_get_msgrcv (AranRusage *rusage)
{
  return rusage->rusage.ru_msgrcv;
}

glong aran_rusage_get_nsignals (AranRusage *rusage)
{
  return rusage->rusage.ru_nsignals;
}
