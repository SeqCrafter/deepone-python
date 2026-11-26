/*
 *    ViennaRNA/utils/strings.c
 *
 *               c  Ivo L Hofacker and Walter Fontana
 *                        Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

#ifndef HAVE_STRDUP
char *
strdup(const char *s)
{
  char *dup;

  dup = vrna_alloc(strlen(s) + 1);
  strcpy(dup, s);
  return dup;
}

#endif

PUBLIC char *
vrna_strdup_printf(const char *format,
                   ...)
{
  char    *result;
  va_list argp;

  va_start(argp, format);
  result = vrna_strdup_vprintf(format, argp);
  va_end(argp); /* Each va_start() or va_copy() needs a va_end() */

  return result;
}

PUBLIC char *
vrna_strdup_vprintf(const char  *format,
                    va_list     argp)
{
  char    *result;
  int     r;

  result = NULL;

#ifndef HAVE_VASPRINTF
  int     count;
  va_list copy;
  va_copy(copy, argp);

  r = -1;

  /* retrieve the number of characters that the string requires */
#ifdef _WIN32
  /*
   * vsnprintf() in Windows is not ANSI compliant, although it's
   * "...included for compliance to the ANSI standard"
   * Thus, we use _vscprintf() that explicitly counts characters
   */
  count = _vscprintf(format, argp);
#else
  count = vsnprintf(NULL, 0, format, argp);
#endif

  if ((count >= 0) && (count < INT_MAX)) {
    char *buf = (char *)vrna_alloc(sizeof(char) * (count + 1));
    if (buf == NULL)
      r = -1;
    else if ((r = vsnprintf(buf, count + 1, format, copy)) < 0)
      free(buf);
    else
      result = buf;
  }

  va_end(copy);  /* Each va_start() or va_copy() needs a va_end() */
#else
  va_list copy;
  va_copy(copy, argp);

  /* the default is to use vasprintf() if available */
  r = vasprintf(&result, format, copy);

  va_end(copy); /* Each va_start() or va_copy() needs a va_end() */
#endif

  /* check for any memory allocation error indicated by r == -1 */
  if (r == -1) {
    vrna_log_warning("vrna_strdup_printf: memory allocation failure!");
    result = NULL;
  }

  return result;
}

/*-----------------------------------------------------------------*/

PUBLIC void
vrna_seq_toupper(char *sequence)
{
  unsigned int i;

  if (sequence)
    for (i = 0; sequence[i]; i++)
      sequence[i] = toupper(sequence[i]);
}

PUBLIC char **
vrna_strsplit(const char  *string,
              const char  *delimiter)
{
  char          delim[2], *ptr, *ptr2, *token, *save, **split;
  unsigned int  n;

  split = NULL;
  n     = 0;

  if (string) {
    if ((delimiter) && (*delimiter))
      delim[0] = *delimiter;
    else
      delim[0] = '&';

    delim[1] = '\0';

    /* copy string such that we can alter it via strtok() */
    ptr2 = strdup(string);

    /* count how many elements we'll extract */
    ptr = ptr2;

    while (*ptr++)
      if (*ptr == *delim)
        n++;

    /*
     * allocate (n + 1) + 1 elements in split list
     * n + 1 elements plus 1 additional element to indicate
     * the last element in split
     */
    split = (char **)vrna_alloc(sizeof(char *) * (n + 2));

    n     = 0;
#ifdef _WIN32
# ifndef __MINGW32__
    token = strtok_s(ptr2, delim, &save);
# else
    token = strtok_r(ptr2, delim, &save);
# endif
#else
    token = strtok_r(ptr2, delim, &save);
#endif

    while (token != NULL) {
      split[n++]  = vrna_strdup_printf("%s", token);
#ifdef _WIN32
# ifndef __MINGW32__
      token = strtok_s(NULL, delim, &save);
# else
      token = strtok_r(NULL, delim, &save);
# endif
#else
      token = strtok_r(NULL, delim, &save);
#endif
    }

    split[n] = NULL;

    free(ptr2);
  }

  return split;
}

#if 0
PUBLIC char *
vrna_strsplice(const char   *string,
               const char   *delimiter,
               unsigned int **positions,
               unsigned int options)
{
  char *result = NULL;
  
  if (string) {
    if (delimiter) {
      if (options & VRNA_STRSPLICE_IN){
        if (positions) {
          /* count how many more characters we require for the fully spliced string */
          for (size_t n = 0; positions[n] != 0; n++);

          size_t dl = strlen(delimiter);
          size_t l  = strlen(string);

          result = (char *)vrna_alloc(sizeof(char) * (l + dl * n + 1));

          /* finally, construct the spliced sequence */
          size_t start = 0;
          size_t end   = 0;
          size_t last_pos = 0;
          /* handle first case separately */
          memcpy(result, string, sizeof(char) * ((*positions)[0] - 1));
          memcpy(result + (*positions)[0] - 1, delimiter, sizeof(char) * dl);
          start += (*positions)[0] - 1;
          end   += (*positions)[0] - 1 + dl;

          for (size_t i = 1; i < n; i++) {
            memcpy(result + end, string + start, sizeof(char) * positions
          }

        } else {
          result = strdup(string);
        }
      } else if (options & VRNA_STRSPLICE_OUT) {
        
      }
    } else {
      /* no delimiter specified, so we don't need to do anything */
      result = strdup(string);
      if ((options & VRNA_STRSPLICE_OUT) &&
          (positions)) {
        *positions = (unsigned int *)vrna_alloc(sizeof(unsigned int));
        (*positions)[0] = 0;
      }
    }
  }

  return result;
}

#endif

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 * ###########################################
 */

#endif
