/** \file fold_compound.c **/

/*
 *                Fold Compound API
 *
 *                This file contains everything necessary to create
 *                and destroy data structures of type vrna_fold_compound_s,
 *                the most fundamental data structure throughout RNAlib
 *
 *                c Ronny Lorenz
 *
 *                ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/sequences/alphabet.h"

#ifndef INLINE
#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif
#endif

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

#define WITH_PTYPE          1L    /* passed to set_fold_compound() to indicate that we need to set fc->ptype */
#define WITH_PTYPE_COMPAT   2L    /* passed to set_fold_compound() to indicate that we need to set fc->ptype_compat */

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

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
PRIVATE void
set_fold_compound(vrna_fold_compound_t  *fc,
                  unsigned int          options,
                  unsigned int          aux);

PRIVATE void
sanitize_bp_span(vrna_fold_compound_t *fc,
                 unsigned int         options);

PRIVATE void
add_params(vrna_fold_compound_t *fc,
           vrna_md_t            *md_p,
           unsigned int         options);

PRIVATE vrna_fold_compound_t *
init_fc_single(void);

PRIVATE INLINE void
nullify(vrna_fold_compound_t *fc);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_fold_compound_free(vrna_fold_compound_t *fc)
{
  unsigned int s;

  if (fc) {
    /* first destroy common attributes */
    vrna_mx_mfe_free(fc);
    /* No partition-function matrices exist in this MFE-only implementation. */
    free(fc->iindx);
    free(fc->jindx);
    free(fc->params);
    free(fc->exp_params);

    vrna_hc_free(fc->hc);
    vrna_ud_remove(fc);
    vrna_gr_reset(fc);
    vrna_sequence_remove_all(fc);

    /* now distinguish the fc type */
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        free(fc->sequence);
        free(fc->sequence_encoding);
        free(fc->sequence_encoding2);
        free(fc->ptype);
        free(fc->ptype_pf_compat);
        vrna_sc_free(fc->sc);
        break;
      default:                      /* do nothing */
        break;
    }

    /* free Distance Class Partitioning stuff (should be NULL if not used) */
    free(fc->reference_pt1);
    free(fc->reference_pt2);
    free(fc->referenceBPs1);
    free(fc->referenceBPs2);
    free(fc->bpdist);
    free(fc->mm1);
    free(fc->mm2);

    /* free local folding related stuff (should be NULL if not used) */
    free(fc->ptype_local);
    free(fc->pscore_local);
#ifdef VRNA_WITH_SVM
    vrna_zsc_filter_free(fc);
#endif

    if (fc->free_auxdata)
      fc->free_auxdata(fc->auxdata);

    free(fc);
  }
}

PUBLIC vrna_fold_compound_t *
vrna_fold_compound(const char       *sequence,
                   const vrna_md_t  *md_p,
                   unsigned int     options)
{
  unsigned int          length, aux_options;
  vrna_fold_compound_t  *fc;
  vrna_md_t             md;

  if (sequence == NULL)
    return NULL;

  /* sanity check */
  length = strlen(sequence);
  if (length == 0) {
    vrna_log_warning("vrna_fold_compound@data_structures.c: "
                         "sequence length must be greater 0");
    return NULL;
  }

  if (length > vrna_sequence_length_max(options)) {
    vrna_log_warning("vrna_fold_compound@data_structures.c: "
                         "sequence length of %d exceeds addressable range",
                         length);
    return NULL;
  }

  fc = init_fc_single();

  fc->length    = length;
  fc->sequence  = strdup(sequence);

  aux_options = 0L;

  /* get a copy of the model details */
  if (md_p)
    md = *md_p;
  else
    vrna_md_set_default(&md);

  /* now for the energy parameters */
  add_params(fc, &md, options);

  sanitize_bp_span(fc, options);

  {
    /* regular global structure prediction */
    aux_options |= WITH_PTYPE;

    if (options & VRNA_OPTION_PF)
      aux_options |= WITH_PTYPE_COMPAT;

    set_fold_compound(fc, options, aux_options);

    if (!(options & VRNA_OPTION_EVAL_ONLY)) {
      /* add default hard constraints */
      vrna_hc_init(fc);

      /* add DP matrices (if required) */
      vrna_mx_add(fc, VRNA_MX_DEFAULT, options);
    }
  }

  return fc;
}

PUBLIC int
vrna_fold_compound_prepare(vrna_fold_compound_t *fc,
                           unsigned int         options)
{
  int ret = 1; /* success */

  /* check maximum sequence length restrictions */
  if (fc->length > vrna_sequence_length_max(options)) {
    vrna_log_warning(
      "vrna_fold_compound_prepare@data_structures.c: sequence length of %d exceeds addressable range",
      fc->length);
    return 0;
  }

  /* make sure to always provide sane bp-span settings */
  sanitize_bp_span(fc, options);

  /* prepare Boltzmann factors if required */
  vrna_params_prepare(fc, options);

  /* prepare ptype array(s) */
  vrna_ptypes_prepare(fc, options);

  ;

  /* prepare hard constraints */
  vrna_hc_prepare(fc, options);

  /* prepare soft constraints data structure, if required */
  vrna_sc_prepare(fc, options);

  /* Add DP matrices, if not they are not present or do not fit current settings */
  vrna_mx_prepare(fc, options);

  /* prepare auxiliary grammar rules data structure, if required */
  ret &= vrna_gr_prepare(fc, options);

  return ret;
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE void
sanitize_bp_span(vrna_fold_compound_t *fc,
                 unsigned int         options)
{
  vrna_md_t *md;

  md = &(fc->params->model_details);

  /* make sure that min_loop_size, max_bp_span, and window_size are sane */
  {
    /* non-local fold mode */
    md->window_size = (int)fc->length;
  }

  if ((md->max_bp_span <= 0) || (md->max_bp_span > md->window_size))
    md->max_bp_span = md->window_size;
}

PRIVATE void
add_params(vrna_fold_compound_t *fc,
           vrna_md_t            *md_p,
           unsigned int         options)
{
  /*
   * ALWAYS provide regular energy parameters
   * remove previous parameters if present and they differ from current model
   */
  if (fc->params) {
    if (memcmp(md_p, &(fc->params->model_details), sizeof(vrna_md_t)) != 0) {
      free(fc->params);
      fc->params = NULL;
    }
  }

  if (!fc->params)
    fc->params = vrna_params(md_p);

  vrna_params_prepare(fc, options);
}

PRIVATE void
set_fold_compound(vrna_fold_compound_t  *fc,
                  unsigned int          options,
                  unsigned int          aux)
{
  char          *sequence, **sequences, **ptr;
  unsigned int  length, s;
  vrna_md_t     *md_p;

  sequence  = NULL;
  sequences = NULL;

  md_p = &(fc->params->model_details);

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      sequence = fc->sequence;

      fc->sequence  = NULL;
      fc->length    = 0;

      /* split input sequences at default delimiter '&' */
      sequences = vrna_strsplit(sequence, NULL);

      /* add individual sequences to fold compound */
      for (ptr = sequences; *ptr; ptr++) {
        vrna_sequence_add(fc, *ptr, VRNA_SEQUENCE_RNA);
        free(*ptr);
      }

      free(sequences);
      free(sequence);

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
      if (fc->strands > 1) {
        fc->cutpoint = fc->nucleotides[0].length + 1;

#if 0
        if (md_p->min_loop_size == TURN)
          md_p->min_loop_size = 0;                              /* is it safe to set this here? */
#endif
      }

#endif

      if (!(options & VRNA_OPTION_EVAL_ONLY)) {
        /* temporary hack for multi-strand case */
        if (fc->strands > 1) {
          int min_loop_size = md_p->min_loop_size;
          md_p->min_loop_size = 0;
          fc->ptype = (aux & WITH_PTYPE) ? vrna_ptypes(fc->sequence_encoding2, md_p) : NULL;
          md_p->min_loop_size = min_loop_size;
        } else {
          fc->ptype = (aux & WITH_PTYPE) ? vrna_ptypes(fc->sequence_encoding2, md_p) : NULL;
        }
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY
        /* backward compatibility ptypes */
        fc->ptype_pf_compat =
          (aux & WITH_PTYPE_COMPAT) ? get_ptypes(fc->sequence_encoding2, md_p, 1) : NULL;
#endif
      }

      break;

    default:                      /* do nothing ? */
      break;
  }

  vrna_sequence_prepare(fc);

  if (!(options & VRNA_OPTION_WINDOW) && (fc->length <= vrna_sequence_length_max(options))) {
    fc->iindx = vrna_idx_row_wise(fc->length);
    fc->jindx = vrna_idx_col_wise(fc->length);
  }
}

PRIVATE vrna_fold_compound_t *
init_fc_single(void)
{
  vrna_fold_compound_t  init = {
    .type = VRNA_FC_TYPE_SINGLE
  };
  vrna_fold_compound_t  *fc = vrna_alloc(sizeof(vrna_fold_compound_t));

  if (fc) {
    memcpy(fc, &init, sizeof(vrna_fold_compound_t));
    nullify(fc);
  }

  return fc;
}

PRIVATE INLINE void
nullify(vrna_fold_compound_t *fc)
{
  if (fc) {
    fc->length            = 0;
    fc->strands           = 0;
    fc->cutpoint          = -1;
    fc->strand_number     = NULL;
    fc->strand_order      = NULL;
    fc->strand_order_uniq = NULL;
    fc->strand_start      = NULL;
    fc->strand_end        = NULL;
    fc->nucleotides       = NULL;
    fc->alignment         = NULL;

    fc->hc            = NULL;
    fc->matrices      = NULL;
    fc->exp_matrices  = NULL;
    fc->params        = NULL;
    fc->exp_params    = NULL;
    fc->iindx         = NULL;
    fc->jindx         = NULL;

    fc->stat_cb       = NULL;
    fc->auxdata       = NULL;
    fc->free_auxdata  = NULL;

    fc->domains_struc = NULL;
    fc->domains_up    = NULL;
    fc->aux_grammar   = NULL;

    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        fc->sequence            = NULL;
        fc->sequence_encoding   = NULL;
        fc->encoding5           = NULL;
        fc->encoding3           = NULL;
        fc->sequence_encoding2  = NULL;
        fc->ptype               = NULL;
        fc->ptype_pf_compat     = NULL;
        fc->sc                  = NULL;

        break;

      }

    fc->maxD1         = 0;
    fc->maxD2         = 0;
    fc->reference_pt1 = NULL;
    fc->reference_pt2 = NULL;
    fc->referenceBPs1 = NULL;
    fc->referenceBPs2 = NULL;
    fc->bpdist        = NULL;
    fc->mm1           = NULL;
    fc->mm2           = NULL;

    fc->window_size = -1;
    fc->ptype_local = NULL;
#ifdef VRNA_WITH_SVM
    fc->zscore_data = NULL;
#endif
  }
}
