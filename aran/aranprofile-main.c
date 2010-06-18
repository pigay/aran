#include "aran-config.h"

#include "aranprofile.h"
#include "aranprofiledb.h"
#include "glib.h"
#include "stdio.h"
#include "glib/gprintf.h"

static GKeyFile *_profiles_file = NULL;
static gchar *_profiles_group = ARAN_PROFILE_DB_DEFAULT_GROUP;
/* static gboolean _verbose = FALSE; */


typedef gdouble (*AranPoly1dProfileFunc) (gpointer operator,
                                            AranZeroFunc init,
                                            gpointer _new,
                                            GDestroyNotify _free,
                                            AranPoly1d *ap1d,
                                            gint nsamples);

typedef struct _AranDevelopmentVTable AranDevelopmentVTable;
struct _AranDevelopmentVTable {
  AranZeroFunc zero_func;
  AranDevelopmentNewFunc new_func;
  GDestroyNotify free_func;
};

static const AranDevelopmentVTable vtable_2d = {
  (AranZeroFunc) aran_development2d_set_zero,
  (AranDevelopmentNewFunc) aran_development2d_new,
  (GDestroyNotify) aran_development2d_free,
};

static const AranDevelopmentVTable vtable_3d = {
  (AranZeroFunc) aran_development3d_set_zero,
  (AranDevelopmentNewFunc) aran_development3d_new,
  (GDestroyNotify) aran_development3d_free,
};

typedef struct _ProfileData ProfileData;

struct _ProfileData {
  gchar *name;
  gpointer operator;
  AranPoly1dProfileFunc profile_func;
  const AranDevelopmentVTable *vtable;
  gint degree;
};

static ProfileData _profiles[] = {
  /* 2D operators */
  {"aran_development2d_m2m", aran_development2d_m2m,
   (AranPoly1dProfileFunc) aran_poly1d_profile_m2m_2d,
   &vtable_2d, 2},
  {"aran_development2d_m2l", aran_development2d_m2l,
   (AranPoly1dProfileFunc) aran_poly1d_profile_m2l_2d,
   &vtable_2d, 2},
  {"aran_development2d_l2l", aran_development2d_l2l,
   (AranPoly1dProfileFunc) aran_poly1d_profile_l2l_2d,
   &vtable_2d, 2},

  /* 3D operators */
  {"aran_development3d_m2m", aran_development3d_m2m,
   (AranPoly1dProfileFunc) aran_poly1d_profile_m2m_3d,
   &vtable_3d, 4},
  {"aran_development3d_m2l", aran_development3d_m2l,
   (AranPoly1dProfileFunc) aran_poly1d_profile_m2l_3d,
   &vtable_3d, 4},
  {"aran_development3d_l2l", aran_development3d_l2l,
   (AranPoly1dProfileFunc) aran_poly1d_profile_l2l_3d,
   &vtable_3d, 4},
  {"aran_development3d_m2m_kkylin", aran_development3d_m2m_kkylin,
   (AranPoly1dProfileFunc) aran_poly1d_profile_m2m_3d,
   &vtable_3d, 3},
  {"aran_development3d_m2l_kkylin", aran_development3d_m2l_kkylin,
   (AranPoly1dProfileFunc) aran_poly1d_profile_m2l_3d,
   &vtable_3d, 3},
  {"aran_development3d_l2l_kkylin", aran_development3d_l2l_kkylin,
   (AranPoly1dProfileFunc) aran_poly1d_profile_l2l_3d,
   &vtable_3d, 3},
  {"aran_development3d_m2m_rotate", aran_development3d_m2m_rotate,
   (AranPoly1dProfileFunc) aran_poly1d_profile_m2m_3d,
   &vtable_3d, 3},
  {"aran_development3d_m2l_rotate", aran_development3d_m2l_rotate,
   (AranPoly1dProfileFunc) aran_poly1d_profile_m2l_3d,
   &vtable_3d, 3},
  {"aran_development3d_l2l_rotate", aran_development3d_l2l_rotate,
   (AranPoly1dProfileFunc) aran_poly1d_profile_l2l_3d,
   &vtable_3d, 3},
  {NULL,},
};

static void _profile (ProfileData pd)
{
  gchar comment[256];
  GError *error;
  gdouble ap1d_terms[pd.degree];
  AranPoly1d ap1d = {
    pd.degree,
    ap1d_terms,
  };
  gdouble chisq;

  chisq =
    pd.profile_func (pd.operator, pd.vtable->zero_func, pd.vtable->new_func,
                     pd.vtable->free_func, &ap1d, 20);

  aran_poly1d_write_key_file (&ap1d, _profiles_file, _profiles_group, pd.name);
  g_sprintf (comment, " \"%s\" function fitting: chisq=%g", pd.name, chisq);
  g_key_file_set_comment (_profiles_file, _profiles_group, pd.name, comment,
                          &error);
}

static GOptionEntry entries[] =
{
/*   { "verbose", 'v', 0, G_OPTION_ARG_NONE, &_verbose, "Be verbose", NULL }, */
  { "group", 'g', 0, G_OPTION_ARG_STRING, &_profiles_group, "Group name",
    NULL },
  { NULL }
};
int main (gint argc, gchar **argv)
{
  gint i;
  gchar *profiles_data;
  GError *error = NULL;
  GOptionContext *context;

  context = g_option_context_new ("- aranprofile");
  g_option_context_add_main_entries (context, entries,
                                     PACKAGE_NAME "-" PACKAGE_VERSION);

  if (!g_option_context_parse (context, &argc, &argv, &error))
    {
      g_print ("option parsing failed: %s\n", error->message);
      return 1;
    }

  aran_init ();

  _profiles_file = g_key_file_new ();

  g_key_file_set_comment (_profiles_file, NULL, NULL,
                          " Aran profiling configuration file",
                          &error);

  i = 0;
  while (_profiles[i].name != NULL)
    {
      _profile (_profiles[i]);

      i++;
    }

  profiles_data = g_key_file_to_data (_profiles_file, NULL, NULL);

  g_key_file_free (_profiles_file);

  g_printf ("%s", profiles_data);

  g_free (profiles_data);

  return 0;
}
