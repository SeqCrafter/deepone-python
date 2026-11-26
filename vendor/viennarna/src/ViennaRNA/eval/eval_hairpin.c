#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/sequences/alphabet.h"
#include "ViennaRNA/eval/hairpin.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/hairpin_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE INLINE int
eval_hp_loop(vrna_fold_compound_t *fc,
             unsigned int         i,
             unsigned int         j,
             unsigned int         options);

PRIVATE INLINE int
eval_ext_hp_loop(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 unsigned int         options);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_hairpin(unsigned int size,
               unsigned int type,
               int          si1,
               int          sj1,
               const char   *sequence,
               vrna_param_t *P)
{
  int energy, salt_correction;

  energy = INF;

  if (P) {
    salt_correction = 0;

    if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
      if (size <= MAXLOOP)
        salt_correction = P->SaltLoop[size + 1];
      else
        salt_correction = vrna_salt_loop_int(size + 1,
                                             P->model_details.salt,
                                             P->temperature + K0,
                                             P->model_details.backbone_length);
    }

    if (size <= 30)
      energy = P->hairpin[size];
    else
      energy = P->hairpin[30] + (int)(P->lxc * log((size) / 30.));

    energy += salt_correction;

    if (size < 3)
      return energy;            /* should only be the case when folding alignments */

    if ((sequence) &&
        (P->model_details.special_hp)) {
      if (size == 4) {
        /* check for tetraloop bonus */
        char tl[7] = {
          0
        }, *ts;
        memcpy(tl, sequence, sizeof(char) * 6);
        tl[6] = '\0';
        if ((ts = strstr(P->Tetraloops, tl)))
          return P->Tetraloop_E[(ts - P->Tetraloops) / 7] + salt_correction;
      } else if (size == 6) {
        char tl[9] = {
          0
        }, *ts;
        memcpy(tl, sequence, sizeof(char) * 8);
        tl[8] = '\0';
        if ((ts = strstr(P->Hexaloops, tl)))
          return P->Hexaloop_E[(ts - P->Hexaloops) / 9] + salt_correction;
      } else if (size == 3) {
        char tl[6] = {
          0
        }, *ts;
        memcpy(tl, sequence, sizeof(char) * 5);
        tl[5] = '\0';
        if ((ts = strstr(P->Triloops, tl)))
          return P->Triloop_E[(ts - P->Triloops) / 6] + salt_correction;

        return energy + (type > 2 ? P->TerminalAU : 0);
      }
    }

    if ((si1 >= 0) &&
        (sj1 >= 0))
      energy += P->mismatchH[type][si1][sj1];
  }

  return energy;
}

/* evaluate a hairpin loop with all constraints applied */
PUBLIC int
vrna_eval_hairpin(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j,
                  unsigned int          options)
{
  if ((fc) &&
      (i > 0) &&
      (j > 0)) {
    /* is this base pair allowed to close a hairpin (like) loop ? */
    if ((options & VRNA_EVAL_LOOP_NO_HC) ||
        (fc->hc == NULL) ||
        (fc->hc->eval_hp(i, j, fc->hc)))
      return (i > j) ? eval_ext_hp_loop(fc, j, i, options) : eval_hp_loop(fc, i, j, options);
  }

  return INF;
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE int
eval_hp_loop(vrna_fold_compound_t *fc,
             unsigned int         i,
             unsigned int         j,
             unsigned int         options)
{
  char              **Ss;
  unsigned int      *sn, **a2s, u, s, n_seq, type, noGUclosure;
  short             *S, *S2, **SS, **S5, **S3;
  int               e, en;
  vrna_param_t      *P;
  vrna_md_t         *md;
  vrna_ud_t         *domains_up;
  struct sc_hp_dat  sc_wrapper;

  sn          = fc->strand_number;
  P           = fc->params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  domains_up  = fc->domains_up;
  e           = INF;

  if (sn[j] != sn[i])
    return e;

  /* regular hairpin loop */
  switch (fc->type) {
    /* sequence alignments */
    default:
      S     = fc->sequence_encoding;
      S2    = fc->sequence_encoding2;
      u     = j - i - 1;
      type  = vrna_get_ptype_md(S2[i], S2[j], md);

      if ((noGUclosure) &&
          ((type == 3) || (type == 4)) &&
          (!(options & VRNA_EVAL_LOOP_NO_HC)))
        break;

      e = vrna_E_hairpin(u, type, S[i + 1], S[j - 1], fc->sequence + i - 1, P);

      break;
  }

  if (e != INF) {
    if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
      init_sc_hp(fc, &sc_wrapper);

      if (sc_wrapper.pair)
        e += sc_wrapper.pair(i, j, &sc_wrapper);

      free_sc_hp(&sc_wrapper);
    }

    /* consider possible ligand binding */
    if (domains_up && domains_up->energy_cb) {
      en = domains_up->energy_cb(fc,
                                 i + 1, j - 1,
                                 VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
                                 domains_up->data);
      if (en != INF)
        en += e;

      e = MIN2(e, en);
    }
  }

  return e;
}

PRIVATE int
eval_ext_hp_loop(vrna_fold_compound_t *fc,
                 unsigned int         i,
                 unsigned int         j,
                 unsigned int         options)
{
  char              **Ss, loopseq[10] = {
    0
  };
  short             *S, *S2, **SS, **S5, **S3;
  unsigned int      **a2s, u1, u2, s, type, n_seq, length, noGUclosure, circ;
  int               e;
  vrna_param_t      *P;
  vrna_md_t         *md;
  struct sc_hp_dat  sc_wrapper;

  length      = fc->length;
  P           = fc->params;
  md          = &(P->model_details);
  noGUclosure = md->noGUclosure;
  circ        = md->circ;
  e           = INF;

  if (!circ)
    return e;

  u1  = length - j;
  u2  = i - 1;

  if ((u1 + u2) < 3)
    return e;

  switch (fc->type) {
    /* sequence alignments */
    default:
      S     = fc->sequence_encoding;
      S2    = fc->sequence_encoding2;
      type  = vrna_get_ptype_md(S2[j], S2[i], md);

      if ((noGUclosure) &&
          ((type == 3) || (type == 4)) &&
          (!(options & VRNA_EVAL_LOOP_NO_HC)))
        break;

      /* maximum special hp loop size: 6 */
      if ((u1 + u2) < 7) {
        memcpy(loopseq, fc->sequence + j - 1, sizeof(char) * (u1 + 1));
        memcpy(loopseq + u1 + 1, fc->sequence, sizeof(char) * (u2 + 1));
        loopseq[u1 + u2 + 2] = '\0';
      }

      e = vrna_E_hairpin(u1 + u2, type, S[j + 1], S[i - 1], loopseq, P);

      break;
  }

  if (e != INF) {
    if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
      init_sc_hp(fc, &sc_wrapper);

      if (sc_wrapper.pair_ext)
        e += sc_wrapper.pair_ext(i, j, &sc_wrapper);

      free_sc_hp(&sc_wrapper);
    }

    /* unstructured domains for circular case?! */
  }

  return e;
}

/*
 #####################################
 # DEPRECATED functions below        #
 #####################################
 */
#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#endif
