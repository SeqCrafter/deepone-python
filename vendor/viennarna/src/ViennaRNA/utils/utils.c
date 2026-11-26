/*
 *    ViennaRNA/utils/basic.c
 *
 *               c  Ivo L Hofacker and Walter Fontana
 *                        Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

/* for getpid() we need some distinction between UNIX and Win systems */
#ifdef _WIN32
#define getpid() GetCurrentProcessId() /* rename windows specific getpid function */
#else
#endif

#include "ViennaRNA/utils/basic.h"

#ifdef WITH_DMALLOC
#endif

#define PRIVATE  static
#define PUBLIC

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */
/*@notnull@ @only@*/
PUBLIC unsigned short xsubi[3];

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

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

#ifndef WITH_DMALLOC
/* include the following two functions only if not including <dmalloc.h> */

PUBLIC void *
vrna_alloc(size_t size)
{
  void *pointer;

  if ((pointer = (void *)calloc(1, size)) == NULL) {
#ifdef EINVAL
    if (errno == EINVAL) {
      fprintf(stderr, "vrna_alloc: requested size: %ld\n", size);
      ("Memory allocation failure -> EINVAL");
    }

    if (errno == ENOMEM)
#endif
    vrna_log_error("Memory allocation failure -> no memory");
  }

  return pointer;
}

PUBLIC void *
vrna_realloc(void     *p,
             size_t   size)
{
  if (p == NULL)
    return vrna_alloc(size);

  p = (void *)realloc(p, size);
  if (p == NULL) {
#ifdef EINVAL
    if (errno == EINVAL) {
      fprintf(stderr, "vrna_realloc: requested size: %ld\n", size);
      vrna_log_error("vrna_realloc allocation failure -> EINVAL");
    }

    if (errno == ENOMEM)
#endif
    vrna_log_error("vrna_realloc allocation failure -> no memory");
  }

  return p;
}

#endif

/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/

/*
 * uniform random number generator; vrna_urn() is in [0,1]
 * uses a linear congruential library routine
 * 48 bit arithmetic
 */

/*------------------------------------------------------------------------*/

/*------------------------------------------------------------------------*/

/*-----------------------------------------------------------------*/

/*-----------------------------------------------------------------*/

PUBLIC int *
vrna_idx_row_wise(unsigned int length)
{
  int i;
  int *idx = (int *)vrna_alloc(sizeof(int) * (length + 1));

  for (i = 1; i <= length; i++)
    idx[i] = (((length + 1 - i) * (length - i)) / 2) + length + 1;
  return idx;
}

PUBLIC int *
vrna_idx_col_wise(unsigned int length)
{
  unsigned int  i;
  int           *idx = (int *)vrna_alloc(sizeof(int) * (length + 1));

  for (i = 1; i <= length; i++)
    idx[i] = (i * (i - 1)) / 2;
  return idx;
}

/*
 #################################
 # STATIC helper functions below #
 #################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

#undef  xrealloc
/* dmalloc.h #define's vrna_realloc */

#endif
