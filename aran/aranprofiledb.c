#include "aranprofiledb.h"

#include "gmodule.h"

#include "string.h"

static GHashTable *_profile_db_names = NULL;
static GHashTable *_profile_db_addresses = NULL;
static GHashTable *_samples_db_names = NULL;
static GHashTable *_samples_db_addresses = NULL;

static GModule *_self_module = NULL;

static gboolean _verbose = TRUE;

static void _atexit ()
{
  if (_profile_db_names != NULL)
    g_hash_table_destroy (_profile_db_names);
  if (_profile_db_addresses != NULL)
    g_hash_table_destroy (_profile_db_addresses);
  if (_samples_db_names != NULL)
    g_hash_table_destroy (_samples_db_names);
  if (_samples_db_addresses != NULL)
    g_hash_table_destroy (_samples_db_addresses);
  if (_self_module != NULL)
    g_module_close (_self_module);
}

static void _require_hashes ()
{
  if (_profile_db_names == NULL)
    {
      _profile_db_addresses = g_hash_table_new (g_direct_hash, g_direct_equal);
      _samples_db_addresses = g_hash_table_new (g_direct_hash, g_direct_equal);

      /* only names hash table has responsibility to destroy stored values */
      _profile_db_names =
        g_hash_table_new_full (g_str_hash, g_str_equal, g_free,
                               (GDestroyNotify) aran_poly1d_free);
      _samples_db_names =
        g_hash_table_new_full (g_str_hash, g_str_equal, g_free,
                               (GDestroyNotify) g_hash_table_destroy);

      /* open this executable for symbol matching */
      _self_module = g_module_open (NULL, G_MODULE_BIND_LAZY);

      g_atexit (_atexit);

      aran_profile_db_read_file (NULL, NULL);
    }
}

static void aran_profile_db_set_value_unsafe (gchar *name, gpointer address,
                                              AranPoly1d *value)
{
  g_hash_table_replace (_profile_db_names, g_strdup (name), value);

  if (address != NULL)
    g_hash_table_replace (_profile_db_addresses, address, value);

}

static void aran_samples_db_set_value_unsafe (gchar *name, gpointer address,
                                              GHashTable *value)
{
  g_hash_table_replace (_samples_db_names, g_strdup (name), value);

  if (address != NULL)
    g_hash_table_replace (_samples_db_addresses, address, value);

}

void aran_profile_db_set_verbose (gboolean verbose)
{
  _verbose = verbose;
}

AranPoly1d *aran_profile_db_get_by_name (gchar *name)
{
  _require_hashes ();

  return g_hash_table_lookup (_profile_db_names, name);
}

AranPoly1d *aran_profile_db_get_by_address (gpointer address)
{
  _require_hashes ();

  return g_hash_table_lookup (_profile_db_addresses, address);
}

GHashTable *aran_samples_db_get_by_name (gchar *name)
{
  _require_hashes ();

  return g_hash_table_lookup (_samples_db_names, name);
}

GHashTable *aran_samples_db_get_by_address (gpointer address)
{
  _require_hashes ();

  return g_hash_table_lookup (_samples_db_addresses, address);
}

static gpointer _get_symbol_by_name (gchar *symbol_name)
{
  gpointer symbol;

  /* try to look for a symbol in this executable */
  if (!g_module_symbol (_self_module, symbol_name, &symbol))
    return NULL;

  return symbol;
}

/* parse array of string of the form "<some_integer>:<some_double>" */
static GHashTable *_read_samples (gchar **strv)
{
  GHashTable *ret = g_hash_table_new_full (g_double_hash, g_double_equal, g_free,
                                           g_free);

  while (*strv != NULL)
    {
      gchar **kvstrv = g_strsplit (*strv, ":", 2);

      if (kvstrv[0] != NULL && kvstrv[1] != NULL)
        {
          gchar *endptr;
          gdouble k;
          gdouble *kp;
          gdouble v;
          gdouble *vp;

          k = (gdouble) g_ascii_strtod (kvstrv[0], &endptr);
          if (endptr == kvstrv[0]) continue;

          v = g_ascii_strtod (kvstrv[1], &endptr);
          if (endptr == kvstrv[1]) continue;

          kp = g_malloc (sizeof (gdouble));
          *kp = k;

          vp = g_malloc (sizeof (gdouble));
          *vp = v;

          g_hash_table_replace (ret, kp, vp);
        }

      g_strfreev (kvstrv);
      strv ++;
    }

  return ret;
}

gboolean aran_profile_db_read_file (const gchar *filename, const gchar *group)
{
  gboolean ok;
  GKeyFile *kf = g_key_file_new ();
  GError *error = NULL;
  gint i;
  gchar **keys;
  gsize nkeys;
  const gchar *orig_filename = filename;

  _require_hashes ();

  if (filename == NULL)
    {
      filename = g_getenv ("ARAN_PROFILE_DB");
      if (filename == NULL)
        {
          filename = ARAN_DATA_DIR G_DIR_SEPARATOR_S "profiledb.ini";
        }
    }

  ok = g_key_file_load_from_file (kf, filename, G_KEY_FILE_NONE, NULL);
  if (! ok)
    {
      if (_verbose || orig_filename != NULL)
        g_warning ("Unable to open profile_db file \"%s\"", filename);
      goto error;
    }

  if (group == NULL) group = g_getenv ("ARAN_PROFILE_GROUP");
  if (group == NULL) group = ARAN_PROFILE_DB_DEFAULT_GROUP;

  ok = g_key_file_has_group (kf, group);
  if (! ok)
    {
      if (_verbose)
        g_warning ("Unable to locate group \"%s\" in profile_db file \"%s\"",
                   group, filename);
      goto error;
    }

  if (_verbose)
    g_debug ("Loading profile db file \"%s\" with group \"%s\"", filename,
             group);

  keys = g_key_file_get_keys (kf, group, &nkeys, &error);

  if (error != NULL)
        {
          g_error ("Failed to load profile \"[%s]\": error=\"%s\"",
                   group, error->message);

          ok = FALSE;
          goto error;
        }
 
  for (i=0; i<nkeys; i++)
    {
      gpointer address = _get_symbol_by_name (keys[i]);

      gdouble *terms;
      gsize degree;

      /* key holds samples insted of polynom */
      if (g_str_has_prefix (keys[i], ARAN_PROFILE_DB_SAMPLES_PREFIX))
        {
          gchar **strv = g_key_file_get_string_list (kf, group, keys[i], NULL, &error);
          GHashTable *samples;

          if (strv == NULL) continue;

          if (_verbose)
            g_debug ("read samples %s\n", keys[i] + strlen (ARAN_PROFILE_DB_SAMPLES_PREFIX));

          samples = _read_samples (strv);

          address = _get_symbol_by_name (keys[i] + strlen (ARAN_PROFILE_DB_SAMPLES_PREFIX));

          aran_samples_db_set_value_unsafe (keys[i] + strlen (ARAN_PROFILE_DB_SAMPLES_PREFIX), address, samples);

          g_strfreev (strv);

          continue;
        }

      terms = g_key_file_get_double_list (kf, group, keys[i], &degree, &error);

      if (error != NULL)
        {
          g_error ("Failed to load profile \"[%s]%s\": error=\"%s\"",
                   group, keys[i], error->message);
          error = NULL;
        }
      else
        {
          AranPoly1d *ap1d = aran_poly1d_new_with_terms (degree-1, terms);

          address = _get_symbol_by_name (keys[i]);

          if (_verbose)
            g_debug ("loaded profile symbol \"%s\" (addr=%p)", keys[i],
                     address);

          aran_profile_db_set_value_unsafe (keys[i], address, ap1d);
        }
    }

  g_strfreev (keys);

 error:
  g_key_file_free (kf);

  return ok;
}

/**
 * aran_profile_db_set_value:
 * @address: storage key for string lookup
 * @name: storage key for pointer lookup
 * @value: value to store
 *
 * Stores @value polynomial into profile database at @name (mandatory) and
 * @address (optional)
 */
void aran_profile_db_set_value (gchar *name, gpointer address,
                                AranPoly1d *value)
{
  g_return_if_fail (name != NULL);

  _require_hashes ();

  aran_profile_db_set_value_unsafe (name, address, value);
}


gdouble aran_profile_db_name_eval (gchar *name, gdouble x)
{
  AranPoly1d *ap1d;
  GHashTable *samples = aran_samples_db_get_by_name (name);

  if (samples)
    {
      gdouble *retval = g_hash_table_lookup (samples, &x);

      if (retval != NULL)
        {
          return *retval;
        }
    }

  ap1d = aran_profile_db_get_by_name (name);

  /* return invalid value in case the symbol wasn't found */
  if (ap1d == NULL) return 0./0.;

  return aran_poly1d_eval (ap1d, x);
}

gdouble aran_profile_db_address_eval (gpointer address, gdouble x)
{
  AranPoly1d *ap1d;
  GHashTable *samples = aran_samples_db_get_by_address (address);

  if (samples)
    {
      gdouble *retval = g_hash_table_lookup (samples, &x);

      if (retval != NULL)
        {
          return *retval;
        }
    }

  ap1d = aran_profile_db_get_by_address (address);

  /* return invalid value in case the symbol wasn't found */
  if (ap1d == NULL) return 0./0.;

  return aran_poly1d_eval (ap1d, x);
}


