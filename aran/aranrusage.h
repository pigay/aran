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

#ifndef __ARAN_RUSAGE_H__
#define __ARAN_RUSAGE_H__

#include <sys/time.h>
#include <sys/resource.h>

#include <glib.h>

#define ARAN_RUSAGE_PROFILE_BEGIN                     \
  {                                                   \
    AranRusage _prof_ar1 = {ARAN_RUSAGE_SELF};        \
    aran_rusage_set (&_prof_ar1);

#define ARAN_RUSAGE_PROFILE_VAR(comp)           \
  aran_rusage_get_##comp (&_prof_ar1);

#define ARAN_RUSAGE_PROFILE_END(gets)                \
  aran_rusage_set_rel (&_prof_ar1);                  \
  gets                                               \
  }


typedef enum _AranRusageWho AranRusageWho;

typedef struct _AranRusage AranRusage;

enum _AranRusageWho
{
  ARAN_RUSAGE_DEFAULT = 0,
  ARAN_RUSAGE_SELF = RUSAGE_SELF,
  ARAN_RUSAGE_CHILDREN = RUSAGE_CHILDREN,
};

struct _AranRusage
{
  int who;

  struct rusage rusage;
	
};

void aran_rusage_set (AranRusage *rusage);
void aran_rusage_set_rel (AranRusage *rusage);

time_t aran_rusage_get_usec (AranRusage *rusage);
time_t aran_rusage_get_uusec (AranRusage *rusage);

time_t aran_rusage_get_ssec (AranRusage *rusage);
time_t aran_rusage_get_susec (AranRusage *rusage);

gdouble aran_rusage_get_utime (AranRusage *rusage);
gdouble aran_rusage_get_stime (AranRusage *rusage);

glong aran_rusage_get_maxrss (AranRusage *rusage);
glong aran_rusage_get_ixrss (AranRusage *rusage);
glong aran_rusage_get_idrss (AranRusage *rusage);
glong aran_rusage_get_isrss (AranRusage *rusage);
glong aran_rusage_get_minflt (AranRusage *rusage);
glong aran_rusage_get_majflt (AranRusage *rusage);
glong aran_rusage_get_nswap (AranRusage *rusage);
glong aran_rusage_get_inblock (AranRusage *rusage);
glong aran_rusage_get_oublock (AranRusage *rusage);
glong aran_rusage_get_msgsnd (AranRusage *rusage);
glong aran_rusage_get_msgrcv (AranRusage *rusage);
glong aran_rusage_get_nsignals (AranRusage *rusage);


#endif /*  __ARAN_RUSAGE_H__*/

