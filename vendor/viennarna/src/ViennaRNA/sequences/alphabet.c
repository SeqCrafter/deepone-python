/*
 *  alphabet.c
 *
 *  Code for handling nucleotide and base pair alphabet
 *
 *  Part of the ViennaRNA Package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <ctype.h>
#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/sequences/alphabet.h"

/*
 * For now, we neglect all non-standard nucleotides in an input sequence, i.e. only
 * ACGTUN is allowed.
 *
 * However, the standard nucleotide ambiguity code table would allow for many more:
 *
 * A = Adenylic acid
 * C = Cytidylic acid
 * G = Guanylic acid
 * T = Thymidylic acid
 * U = Uridylic acid
 * I = Inosylic acid
 * R = A or G = puRine
 * Y = C or T = pYrimidine
 * K = G or T = Keto
 * M = A or C = aMino
 * S = G or C = Strong base pair
 * W = A or T = Weak base pair
 * B = not A (G or C or T)
 * D = not C (A or G or T)
 * H = not G (A or C or T)
 * V = not T/U (A or C or G)
 * N = aNy base  (by convention, X is used for unknown amino acids, N for unknown nucleotides)
 *
 * For the future, we aim to accept all of the above codes.
 */

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */
PRIVATE const char Law_and_Order[] = "_ACGUTXKI";

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */
PRIVATE char *
wrap_get_ptypes(const short *S,
                vrna_md_t   *md);                               /* provides backward compatibility for old ptypes array in pf computations */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC unsigned int
vrna_sequence_length_max(unsigned int options)
{
  ;;

  /*
   * return (unsigned int)sqrt((double)INT_MAX);
   */
  /*
   *  many functions in RNAlib still rely on the sequence length
   *  at pos 0 in the integer encoded sequence array S. Since this
   *  encoding is stored in a short * array, the maximum length
   *  of any sequence is SHRT_MAX
   */
  return (unsigned int)SHRT_MAX;
}

PUBLIC void
vrna_ptypes_prepare(vrna_fold_compound_t  *fc,
                    unsigned int          options)
{
  if (!fc)
    return;

  if (options & VRNA_OPTION_MFE) {
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        {
          if (!fc->ptype) {
            /* temporary hack for multi-strand case */
            if (fc->strands > 1) {
              int min_loop_size = fc->params->model_details.min_loop_size;
              fc->params->model_details.min_loop_size = 0;
              fc->ptype = vrna_ptypes(fc->sequence_encoding2,
                                      &(fc->params->model_details));
              fc->params->model_details.min_loop_size = min_loop_size;
            } else {
              fc->ptype = vrna_ptypes(fc->sequence_encoding2,
                                      &(fc->params->model_details));
            }
          }
        }

        break;

      default:
        break;
    }
  }

  ;
}

PUBLIC char *
vrna_ptypes(const short *S,
            vrna_md_t   *md)
{
  char  *ptype;
  int   n, i, j, k, l, *idx;
  int   min_loop_size = md->min_loop_size;

  n = S[0];

  if ((unsigned int)n > vrna_sequence_length_max(VRNA_OPTION_DEFAULT)) {
    vrna_log_warning("vrna_ptypes@alphabet.c: sequence length of %d exceeds addressable range",
                         n);
    return NULL;
  }

  ptype = (char *)vrna_alloc(sizeof(char) * ((n * (n + 1)) / 2 + 2));
  idx   = vrna_idx_col_wise(n);

  for (k = 1; k < n - min_loop_size; k++)
    for (l = 1; l <= 2; l++) {
      int type, ntype = 0, otype = 0;
      i = k;
      j = i + min_loop_size + l;
      if (j > n)
        continue;

      type = md->pair[S[i]][S[j]];
      while ((i >= 1) && (j <= n)) {
        if ((i > 1) && (j < n))
          ntype = md->pair[S[i - 1]][S[j + 1]];

        if (md->noLP && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */

        ptype[idx[j] + i] = (char)type;
        otype             = type;
        type              = ntype;
        i--;
        j++;
      }
    }
  free(idx);
  return ptype;
}

PUBLIC short *
vrna_seq_encode(const char  *sequence,
                vrna_md_t   *md)
{
  unsigned int  i, l;
  short         *S = NULL;

  if (sequence && md) {
    S = vrna_seq_encode_simple(sequence, md);

    l = (unsigned int)strlen(sequence);

    for (i = 1; i <= l; i++)
      S[i] = md->alias[S[i]];

    S[l + 1]  = S[1];
    S[0]      = S[l];
  }

  return S;
}

PUBLIC short *
vrna_seq_encode_simple(const char *sequence,
                       vrna_md_t  *md)
{
  unsigned int  i, l;
  short         *S = NULL;

  if (sequence && md) {
    l = (unsigned int)strlen(sequence);
    S = (short *)vrna_alloc(sizeof(short) * (l + 2));

    for (i = 1; i <= l; i++) /* make numerical encoding of sequence */
      S[i] = (short)vrna_nucleotide_encode(sequence[i - 1], md);

    S[l + 1]  = S[1];
    S[0]      = (short)l;
  }

  return S;
}

PUBLIC int
vrna_nucleotide_encode(char       c,
                       vrna_md_t  *md)
{
  /* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */
  int code = -1;

  c = toupper(c);

  if (md) {
    if (md->energy_set > 0) {
      code = (int)(c - 'A') + 1;
    } else {
      const char *pos;
      pos = strchr(Law_and_Order, c);
      if (pos == NULL)
        code = 0;
      else
        code = (int)(pos - Law_and_Order);

      if (code > 5)
        code = 0;

      if (code > 4)
        code--;           /* make T and U equivalent */
    }
  }

  return code;
}

PUBLIC unsigned int
vrna_get_ptype_md(int       i,
                  int       j,
                  vrna_md_t *md)
{
  unsigned int tt = (unsigned int)md->pair[i][j];

  return (tt == 0) ? 7 : tt;
}

PUBLIC unsigned int
vrna_get_ptype(int  ij,
               char *ptype)
{
  unsigned int tt = (unsigned int)ptype[ij];

  return (tt == 0) ? 7 : tt;
}

PUBLIC unsigned int
vrna_get_ptype_window(int   i,
                      int   j,
                      char  **ptype)
{
  unsigned int tt = (unsigned int)ptype[i][j - i];

  return (tt == 0) ? 7 : tt;
}

PRIVATE char *
wrap_get_ptypes(const short *S,
                vrna_md_t   *md)
{
  char  *ptype;
  int   n, i, j, k, l, *idx;

  n     = S[0];
  ptype = (char *)vrna_alloc(sizeof(char) * ((n * (n + 1)) / 2 + 2));
  idx   = vrna_idx_row_wise(n);
  int   min_loop_size = md->min_loop_size;

  for (k = 1; k < n - min_loop_size; k++)
    for (l = 1; l <= 2; l++) {
      int type, ntype = 0, otype = 0;
      i = k;
      j = i + min_loop_size + l;
      if (j > n)
        continue;

      type = md->pair[S[i]][S[j]];
      while ((i >= 1) && (j <= n)) {
        if ((i > 1) && (j < n))
          ntype = md->pair[S[i - 1]][S[j + 1]];

        if (md->noLP && (!otype) && (!ntype))
          type = 0; /* i.j can only form isolated pairs */

        ptype[idx[i] - j] = (char)type;
        otype             = type;
        type              = ntype;
        i--;
        j++;
      }
    }
  free(idx);
  return ptype;
}

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

PUBLIC char *
get_ptypes(const short  *S,
           vrna_md_t    *md,
           unsigned int idx_type)
{
  if (S) {
    if ((unsigned int)S[0] > vrna_sequence_length_max(VRNA_OPTION_DEFAULT)) {
      vrna_log_warning("get_ptypes@alphabet.c: sequence length of %d exceeds addressable range",
                           (int)S[0]);
      return NULL;
    }

    if (idx_type)
      return wrap_get_ptypes(S, md);
    else
      return vrna_ptypes(S, md);
  } else {
    return NULL;
  }
}

#endif
