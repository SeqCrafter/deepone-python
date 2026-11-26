/* constraints handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/utils/basic.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define STATE_CLEAN         (unsigned char)0
#define STATE_DIRTY_UP      (unsigned char)1
#define STATE_DIRTY_BP      (unsigned char)2
#define STATE_UNINITIALIZED (unsigned char)4

#include "hc_depot.inc"

#include "wrap_exterior_hc.inc"
#include "wrap_hairpin_hc.inc"
#include "wrap_internal_hc.inc"
#include "wrap_multibranch_hc.inc"

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
PRIVATE unsigned char
default_pair_constraint(vrna_fold_compound_t  *fc,
                        unsigned int          i,
                        unsigned int          j);

PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc);

PRIVATE void
hc_update_up(vrna_fold_compound_t *vc);

PRIVATE void
default_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options);

PRIVATE void
default_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options);

PRIVATE void
prepare_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options);

PRIVATE void
prepare_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options);

PRIVATE void
hc_prepare_eval_f(vrna_hc_t *hc);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC void
vrna_hc_init(vrna_fold_compound_t *vc)
{
  unsigned int  n;
  vrna_hc_t     *hc;

  n = vc->length;

  /* free previous hard constraints */
  vrna_hc_free(vc->hc);

  /* allocate memory new hard constraints data structure */
  hc          = (vrna_hc_t *)vrna_alloc(sizeof(vrna_hc_t));
  hc->type    = VRNA_HC_DEFAULT;
  hc->n       = n;
  hc->strands = vc->strands;
  hc->sn      = vc->strand_number;
  hc->mx      = (unsigned char *)vrna_alloc(sizeof(unsigned char) * ((n + 1) * (n + 1) + 1));
  hc->up_ext  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  hc->up_hp   = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  hc->up_int  = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  hc->up_ml   = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (n + 2));
  hc->depot   = NULL;

  /* set new hard constraints */
  vc->hc = hc;

  /* prefill default values  */
  hc_reset_to_default(vc);

  /* add null pointers for the generalized hard constraint feature */
  hc->f         = NULL;
  hc->data      = NULL;
  hc->free_data = NULL;

  /* update */
  hc_update_up(vc);

  hc->state = STATE_CLEAN;
}

PUBLIC int
vrna_hc_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options)
{
  int ret = 0;

  if (fc) {
    {
      if (fc->hc->state & STATE_UNINITIALIZED) {
        default_hc_up(fc, options);
        default_hc_bp(fc, options);
      }

      if (fc->hc->state & STATE_DIRTY_UP)
        prepare_hc_up(fc, options);

      if (fc->hc->state & STATE_DIRTY_BP)
        prepare_hc_bp(fc, options);

      if (fc->hc->state & ~STATE_CLEAN)
        hc_update_up(fc);

    }

    fc->hc->sn = fc->strand_number;

    hc_prepare_eval_f(fc->hc);

    fc->hc->state = STATE_CLEAN;
    ret           = 1;
  }

  return ret;
}

PUBLIC void
vrna_hc_free(vrna_hc_t *hc)
{
  if (hc) {
    if (hc->type == VRNA_HC_DEFAULT)
      free(hc->mx);
    else if (hc->type == VRNA_HC_WINDOW)
      free(hc->matrix_local);

    hc_depot_free(hc);

    free(hc->up_ext);
    free(hc->up_hp);
    free(hc->up_int);
    free(hc->up_ml);

    if (hc->free_data)
      hc->free_data(hc->data);

    free(hc);
  }
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE unsigned char
default_pair_constraint(vrna_fold_compound_t  *fc,
                        unsigned int          i,
                        unsigned int          j)
{
  unsigned char constraint, can_stack;
  short         *S;
  unsigned int  *sn;
  unsigned int  type;
  vrna_md_t     *md;

  sn          = fc->strand_number;
  md          = &(fc->params->model_details);
  constraint  = VRNA_CONSTRAINT_CONTEXT_NONE;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      S = fc->sequence_encoding2;
      if (((j - i) < (unsigned int)md->max_bp_span) &&
          ((sn[i] != sn[j]) || ((j - i) > (unsigned int)md->min_loop_size))) {
        type = md->pair[S[i]][S[j]];
        switch (type) {
          case 0:
            break;
          case 3:
          /* fallthrough */
          case 4:
            if (md->noGU) {
              break;
            } else if (md->noGUclosure) {
              constraint  = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
              constraint  &= ~(VRNA_CONSTRAINT_CONTEXT_HP_LOOP | VRNA_CONSTRAINT_CONTEXT_MB_LOOP);
              break;
            }

          /* else fallthrough */
          default:
            constraint = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            break;
        }

        if (md->noLP) {
          /* check, whether this nucleotide can actually stack with anything or only forms isolated pairs */
          can_stack = VRNA_CONSTRAINT_CONTEXT_NONE;

          /* can it be enclosed by another base pair? */
          if ((i > 1) &&
              (j < fc->length) &&
              (((j - i + 2) < (unsigned int)md->max_bp_span) || (sn[i - 1] != sn[j + 1])) &&
              (md->pair[S[i - 1]][S[j + 1]]))
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          /* can it enclose another base pair? */
          if ((i + 2 < j) &&
              (((j - i - 2) > (unsigned int)md->min_loop_size) || (sn[i + 1] != sn[j - 1])) &&
              (md->pair[S[i + 1]][S[j - 1]]))
            can_stack = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

          constraint &= can_stack;
        }
      }

      break;

    }

  return constraint;
}

PRIVATE void
default_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned int  i, n;
  vrna_hc_t     *hc;

  hc = fc->hc;

  {
    n = fc->length;

    for (i = 1; i <= n; i++)
      hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
  }
}

PRIVATE void
prepare_hc_up(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned char   option, type, t1, t2;
  unsigned int    i, j, k, n, s, *ss;
  vrna_hc_t       *hc;
  vrna_hc_depot_t *depot;

  hc = fc->hc;

  {
    n     = fc->length;
    ss    = fc->strand_start;
    depot = hc->depot;

    /* 2. apply constraints as stored in depot */
    if ((depot) && (depot->up)) {
      for (s = 0; s < depot->strands; s++) {
        for (k = 1; k <= depot->up_size[s]; k++) {
          /* process nucleotide-specific constraint */
          option  = depot->up[s][k].context;
          i       = ss[s] + k - 1; /* constraint position in current strand order */

          if (depot->up[s][k].nonspec) {
            /* this is actually a must-pair constraint */
            type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            /* acknowledge pairing direction */
            t1  = (depot->up[s][k].direction <= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;
            t2  = (depot->up[s][k].direction >= 0) ? type : VRNA_CONSTRAINT_CONTEXT_NONE;

            if (option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE) {
              /* only allow for possibly non-canonical pairs, do not enforce them */
              for (j = 1; j < i; j++) {
                hc->mx[n * i + j] |= t1;
                hc->mx[n * j + i] |= t1;
              }
              for (j = i + 1; j <= n; j++) {
                hc->mx[n * i + j] |= t2;
                hc->mx[n * j + i] |= t2;
              }
            } else {
              /* force pairing direction */
              for (j = 1; j < i; j++) {
                hc->mx[n * i + j] &= t1;
                hc->mx[n * j + i] &= t1;
              }
              for (j = i + 1; j <= n; j++) {
                hc->mx[n * i + j] &= t2;
                hc->mx[n * j + i] &= t2;
              }
            }

            /* nucleotide mustn't be unpaired */
            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE)
              hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
          } else {
            /* 'regular' nucleotide-specific constraint */
            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              /*
               * force nucleotide to appear unpaired within a certain type of loop
               * do not allow i to be paired with any other nucleotide
               */
              if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                for (j = 1; j < i; j++) {
                  hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
                for (j = i + 1; j <= n; j++) {
                  hc->mx[n * i + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * j + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }

              type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              hc->mx[n * i + i] = type;
            } else {
              type = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

              /* do not allow i to be paired with any other nucleotide (in context type) */
              if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
                for (j = 1; j < i; j++) {
                  hc->mx[n * i + j] &= ~type;
                  hc->mx[n * j + i] &= ~type;
                }
                for (j = i + 1; j <= n; j++) {
                  hc->mx[n * i + j] &= ~type;
                  hc->mx[n * j + i] &= ~type;
                }
              }

              hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            }
          }
        }
      }
    }
  }
}

PRIVATE void
default_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned int  i, j, n;
  vrna_hc_t     *hc;

  hc = fc->hc;

  {
    n = fc->length;

    for (j = n; j > 1; j--) {
      for (i = 1; i < j; i++) {
        hc->mx[n * i + j] = default_pair_constraint(fc, i, j);
        hc->mx[n * j + i] = hc->mx[n * i + j];
      }
    }
  }
}

PRIVATE void
prepare_hc_bp(vrna_fold_compound_t  *fc,
              unsigned int          options)
{
  unsigned char   option;
  unsigned int    i, j, k, p, q, n, actual_i, actual_j, strand_j, s, *ss;
  int             *idx;
  vrna_hc_t       *hc;
  vrna_hc_depot_t *depot;

  hc    = fc->hc;
  depot = hc->depot;
  ss    = fc->strand_start;

  if ((!depot) || (!depot->bp))
    return;

  {
    n   = fc->length;
    idx = fc->jindx;

    /* 2. apply constraints stored in depot */
    for (s = 0; s < depot->strands; s++) {
      for (actual_i = 1; actual_i <= depot->bp_size[s]; actual_i++) {
        for (k = 0; k < depot->bp[s][actual_i].list_size; k++) {
          option    = depot->bp[s][actual_i].context[k];
          strand_j  = depot->bp[s][actual_i].strand_j[k];
          actual_j  = depot->bp[s][actual_i].j[k];
          i         = ss[s] + actual_i - 1;         /* constraint position in current strand order */
          j         = ss[strand_j] + actual_j - 1;  /* constraint position in current strand order */

          if (i < j) {
            /* apply the constraint */
            hc->mx[n * i + j] = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;
            hc->mx[n * j + i] = option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS;

            /* is the ptype reset actually required??? */
            if ((fc->type == VRNA_FC_TYPE_SINGLE) &&
                (option & VRNA_CONSTRAINT_CONTEXT_ALL_LOOPS)) {
              /* reset ptype in case (i,j) is a non-canonical pair */
              if (fc->ptype[idx[j] + i] == 0)
                fc->ptype[idx[j] + i] = 7;
            }

            if (!(option & VRNA_CONSTRAINT_CONTEXT_NO_REMOVE)) {
              /*
               * remove all conflicting base pairs, i.e. do not allow i,j to pair
               * with any other nucleotide k
               */
              for (p = 1; p < i; p++) {
                hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

                for (q = i + 1; q < j; q++) {
                  hc->mx[n * p + q] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * q + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }
              for (p = i + 1; p < j; p++) {
                hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;

                for (q = j + 1; q <= n; q++) {
                  hc->mx[n * p + q] = VRNA_CONSTRAINT_CONTEXT_NONE;
                  hc->mx[n * q + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                }
              }
              for (p = j + 1; p <= n; p++) {
                hc->mx[n * i + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * j + p] = VRNA_CONSTRAINT_CONTEXT_NONE;
                hc->mx[n * p + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
              }
            }

            if (option & VRNA_CONSTRAINT_CONTEXT_ENFORCE) {
              /* do not allow i,j to be unpaired */
              hc->mx[n * i + i] = VRNA_CONSTRAINT_CONTEXT_NONE;
              hc->mx[n * j + j] = VRNA_CONSTRAINT_CONTEXT_NONE;
            }
          }
        }
      }
    }
  }
}

struct hc_bp {
  int           i;
  int           j;
  unsigned char options;
};

PRIVATE void
hc_reset_to_default(vrna_fold_compound_t *vc)
{
  vrna_hc_t *hc = vc->hc;

  /*
   * #########################
   * fill with default values
   * #########################
   */

  default_hc_up(vc, 0);

  default_hc_bp(vc, 0);

  /* should we reset the generalized hard constraint feature here? */
  if (hc->f || hc->data) {
    if (hc->free_data)
      hc->free_data(hc->data);

    hc->f         = NULL;
    hc->data      = NULL;
    hc->free_data = NULL;
  }

  hc_prepare_eval_f(hc);
}

PRIVATE void
hc_update_up(vrna_fold_compound_t *vc)
{
  unsigned int  i, n;
  vrna_hc_t     *hc;

  n   = vc->length;
  hc  = vc->hc;

  if (hc->type == VRNA_HC_WINDOW) {
    /* do nothing for now! */
  } else {
    for (hc->up_ext[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in exterior loop */
      hc->up_ext[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) ? 1 +
                      hc->up_ext[i + 1] : 0;

    for (hc->up_hp[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in hairpin loop */
      hc->up_hp[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) ? 1 +
                     hc->up_hp[i + 1] : 0;

    for (hc->up_int[n + 1] = 0, i = n; i > 0; i--) /* unpaired stretch in internal loop */
      hc->up_int[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) ? 1 +
                      hc->up_int[i + 1] : 0;

    for (hc->up_ml[n + 1] = 0, i = n; i > 0; i--)  /* unpaired stretch in multibranch loop */
      hc->up_ml[i] = (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) ? 1 +
                     hc->up_ml[i + 1] : 0;

    /*
     *  loop arround once more until we find a nucleotide that mustn't
     *  be unpaired (needed for circular folding)
     *  Note, circular fold is only possible for single strand predictions
     */
    if (vc->strands < 2) {
      if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
        hc->up_ext[n + 1] = hc->up_ext[1];
        for (i = n; i > 0; i--) {
          if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP)
            hc->up_ext[i] = MIN2(n, 1 + hc->up_ext[i + 1]);
          else
            break;
        }
      }

      if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP) {
        hc->up_hp[n + 1] = hc->up_hp[1];
        for (i = n; i > 0; i--) {
          if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_HP_LOOP)
            hc->up_hp[i] = MIN2(n, 1 + hc->up_hp[i + 1]);
          else
            break;
        }
      }

      if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP) {
        hc->up_int[n + 1] = hc->up_int[1];
        for (i = n; i > 0; i--) {
          if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_INT_LOOP)
            hc->up_int[i] = MIN2(n, 1 + hc->up_int[i + 1]);
          else
            break;
        }
      }

      if (hc->mx[n + 1] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP) {
        hc->up_ml[n + 1] = hc->up_ml[1];
        for (i = n; i > 0; i--) {
          if (hc->mx[n * i + i] & VRNA_CONSTRAINT_CONTEXT_MB_LOOP)
            hc->up_ml[i] = MIN2(n, 1 + hc->up_ml[i + 1]);
          else
            break;
        }
      }
    }
  }
}

PRIVATE void
hc_prepare_eval_f(vrna_hc_t *hc)
{
  hc->eval_ext  = dispatch_wrap_ext_f;
  hc->eval_hp   = dispatch_wrap_hp_f;
  hc->eval_int  = dispatch_wrap_int_f;
  hc->eval_mb   = dispatch_wrap_mb_f;
}

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

#endif
