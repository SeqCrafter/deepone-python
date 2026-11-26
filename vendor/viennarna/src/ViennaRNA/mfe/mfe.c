/*
 *                minimum free energy
 *                RNA secondary structure prediction
 *
 *                c Ivo Hofacker, Chrisoph Flamm
 *                original implementation by
 *                Walter Fontana
 *                g-quadruplex support and threadsafety
 *                by Ronny Lorenz
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/eval/exterior.h"
#include "ViennaRNA/eval/hairpin.h"
#include "ViennaRNA/eval/internal.h"
#include "ViennaRNA/mfe/exterior.h"
#include "ViennaRNA/mfe/internal.h"
#include "ViennaRNA/mfe/multibranch.h"
#include "ViennaRNA/sequences/alphabet.h"

#include "ViennaRNA/intern/grammar_dat.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#define MAXSECTORS        500     /* dimension for a backtrack array */

#define   ADD_OR_INF(a, b)     (((a) != INF) && ((b) != INF) ?  (a) + (b) : INF)

#define M2_FORWARD

struct aux_arrays {
  int                   *cc;  /* auxilary arrays for canonical structures     */
  int                   *cc1; /* auxilary arrays for canonical structures     */
  vrna_mx_mfe_aux_ml_t  ml_helpers;
};

#include "ViennaRNA/constraints/exterior_sc.inc"

struct ms_helpers {
  struct sc_f5_dat      sc_wrapper;
};

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

PRIVATE int
fill_arrays(vrna_fold_compound_t  *fc,
            struct ms_helpers     *ms_dat);

PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               struct aux_arrays    *aux,
               struct ms_helpers    *ms_dat);

PRIVATE INLINE struct aux_arrays *
get_aux_arrays(unsigned int length);

PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux,
                  unsigned int      length);

PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux);

PRIVATE struct ms_helpers *
get_ms_helpers(vrna_fold_compound_t *fc);

PRIVATE void
free_ms_helpers(struct ms_helpers *ms_dat,
                size_t            strands);

PRIVATE void
update_fms5_arrays(vrna_fold_compound_t *fc,
                   int                  i,
                   struct ms_helpers    *ms_dat);

PRIVATE void
update_fms3_arrays(vrna_fold_compound_t *fc,
                   unsigned int         s,
                   struct ms_helpers    *ms_dat);

PRIVATE int
pair_multi_strand(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  struct ms_helpers     *ms_dat);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_mfe(vrna_fold_compound_t *fc,
         char                 *structure)
{
  unsigned int      length;
  int               energy;
  float             mfe;
  struct ms_helpers *ms_dat;

  mfe = (float)(INF / 100.);

  if (fc) {
    length  = fc->length;
    ms_dat  = NULL;

    if (!vrna_fold_compound_prepare(fc, VRNA_OPTION_MFE)) {
      vrna_log_warning("vrna_mfe@mfe.c: Failed to prepare vrna_fold_compound");
      return mfe;
    }

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(fc, VRNA_STATUS_MFE_PRE, fc->auxdata);

    /* call user-defined grammar pre-condition callback function */
    if (fc->aux_grammar) {
      for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->cbs_status); i++)
        if (fc->aux_grammar->cbs_status[i])
          fc->aux_grammar->cbs_status[i](fc, VRNA_STATUS_MFE_PRE, fc->aux_grammar->datas[i]);
    }

    if (fc->strands > 1)
      ms_dat = get_ms_helpers(fc);

    energy = fill_arrays(fc, ms_dat);

    /* DeepOne folds linear RNA using fixed default model details. */

    /* DeepOne requests energy only: structure traceback is omitted. */

    /* call user-defined recursion status callback function */
    if (fc->stat_cb)
      fc->stat_cb(fc, VRNA_STATUS_MFE_POST, fc->auxdata);

    /* call user-defined grammar post-condition callback function */
    if (fc->aux_grammar) {
      for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->cbs_status); i++)
        if (fc->aux_grammar->cbs_status[i])
          fc->aux_grammar->cbs_status[i](fc, VRNA_STATUS_MFE_POST, fc->aux_grammar->datas[i]);
    }

    switch (fc->params->model_details.backtrack_type) {
      case 'C':
        mfe = (float)fc->matrices->c[fc->jindx[length] + 1] / 100.;
        break;

      case 'M':
        mfe = (float)fc->matrices->fML[fc->jindx[length] + 1] / 100.;
        break;

      default:
        mfe = (float)energy / 100.;;

        break;
    }

    free_ms_helpers(ms_dat, fc->strands);
  }

  return mfe;
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

/* fill DP matrices */
PRIVATE int
fill_arrays(vrna_fold_compound_t  *fc,
            struct ms_helpers     *ms_dat)
{
  unsigned int      *sn;
  unsigned int      i, j, length, uniq_ML;
  int               ij, *indx, *f5, *c, *fML, *fM1, *fM2_real;
  vrna_param_t      *P;
  vrna_md_t         *md;
  vrna_mx_mfe_t     *matrices;
  vrna_ud_t         *domains_up;
  struct aux_arrays *helper_arrays;

  length      = (int)fc->length;
  indx        = fc->jindx;
  P           = fc->params;
  md          = &(P->model_details);
  uniq_ML     = md->uniq_ML;
  matrices    = fc->matrices;
  f5          = matrices->f5;
  c           = matrices->c;
  fML         = matrices->fML;
  fM1         = matrices->fM1;
  fM2_real    = matrices->fM2_real;
  domains_up  = fc->domains_up;
  sn          = fc->strand_number;

  /* allocate memory for all helper arrays */
  helper_arrays = get_aux_arrays(length);

  /* pre-processing ligand binding production rule(s) */
  if (domains_up && domains_up->prod_cb)
    domains_up->prod_cb(fc, domains_up->data);

  /* prefill matrices with init contributions */
  for (i = 1; i <= length; i++) {
    c[indx[i] + i] = fML[indx[i] + i] = INF;

    if (fM1)
      fM1[indx[i] + i] = INF;

    if (fM2_real)
      for (j = i; j <= length; j++)
        fM2_real[indx[j] + i] = INF;
  }

  /* start recursion */
  if (length <= ((fc->strands > 1) ? fc->strands : (unsigned int)md->min_loop_size)) {
    /* clean up memory */
    free_aux_arrays(helper_arrays);

    /* return free energy of unfolded chain */
    return 0;
  }

  for (i = length - 1; i >= 1; i--) {
    if ((fc->strands > 1) &&
        (sn[i] != sn[i + 1]))
      update_fms3_arrays(fc, sn[i + 1], ms_dat);

    for (j = i + 1; j <= length; j++) {
      ij = indx[j] + i;

      /* decompose subsegment [i, j] with pair (i, j) */
      c[ij] = decompose_pair(fc, i, j, helper_arrays, ms_dat);

      /* decompose subsegment [i, j] that is multibranch loop part with at least two branches */
      if (fM2_real)
        fM2_real[ij] = vrna_mfe_multibranch_m2_fast(fc,
                                                    i,
                                                    j,
                                                    helper_arrays->ml_helpers);

      /* decompose subsegment [i, j] that is multibranch loop part with at least one branch */
      fML[ij] = vrna_mfe_multibranch_stems_fast(fc, i, j, helper_arrays->ml_helpers);

      /* decompose subsegment [i, j] that is multibranch loop part with exactly one branch */
      if (fM1)
        fM1[ij] = vrna_mfe_multibranch_m1(fc, i, j);

      if (fc->aux_grammar)
        /* call auxiliary grammar rules */
        for (size_t i = 0; i < vrna_array_size(fc->aux_grammar->aux); i++)
          if (fc->aux_grammar->aux[i].cb)
            (void)fc->aux_grammar->aux[i].cb(fc, i, j, fc->aux_grammar->aux[i].data);
    } /* end of j-loop */

    rotate_aux_arrays(helper_arrays, length);

    if (fc->strands > 1)
      update_fms5_arrays(fc, i, ms_dat);
  } /*
     * end of i-loop
     * calculate energies of 5' fragments
     */
  (void)vrna_mfe_exterior_f5(fc);

  /* clean up memory */
  free_aux_arrays(helper_arrays);

  return f5[length];
}

/* post-processing step for circular RNAs */

/*
 **************************
 *  Construct fM_d5 array
 **************************
 */

/*
 **************************
 *  Construct fM_d3 array
 **************************
 */

#if 1

PRIVATE struct ms_helpers *
get_ms_helpers(vrna_fold_compound_t *fc)
{
  struct ms_helpers *dat = (struct ms_helpers *)vrna_alloc(sizeof(struct ms_helpers));

  init_sc_f5(fc, &(dat->sc_wrapper));

  return dat;
}

PRIVATE void
free_ms_helpers(struct ms_helpers *ms_dat,
                size_t            strands)
{
  if (strands > 1)
    free_sc_f5(&(ms_dat->sc_wrapper));

  free(ms_dat);
}

PRIVATE void
update_fms5_arrays(vrna_fold_compound_t *fc,
                   int                  i,
                   struct ms_helpers    *ms_dat)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *se, type;
  int                   e, tmp, **fms5, *c, *idx, n, end, dangle_model, base;
  vrna_param_t          *params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  struct sc_f5_dat      *sc_wrapper;
  sc_ext_red_cb         sc_spl;
  sc_ext_red_cb         sc_red_stem;
  sc_ext_red_cb         sc_red_ext;

  n             = fc->length;
  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  idx           = fc->jindx;
  c             = fc->matrices->c;
  fms5          = fc->matrices->fms5;
  sn            = fc->strand_number;
  se            = fc->strand_end;
  params        = fc->params;
  md            = &(params->model_details);
  dangle_model  = md->dangles;
  hc            = fc->hc;
  sc_wrapper    = &(ms_dat->sc_wrapper);
  sc_spl        = sc_wrapper->decomp;
  sc_red_stem   = sc_wrapper->red_stem;
  sc_red_ext    = sc_wrapper->red_ext;

  for (size_t s = 0; s < fc->strands; s++) {
    e   = tmp = INF;
    end = se[s];

    if (i < end) {
      if (hc->eval_ext(i, end, i + 1, end, VRNA_DECOMP_EXT_EXT, hc)) {
        tmp = fms5[s][i + 1];
        if ((sc_red_ext) && (tmp != INF))
          tmp += sc_red_ext(i, end, i + 1, end, sc_wrapper);
      }
    } else {
      tmp = 0;
    }

    e = MIN2(e, tmp);

    for (int k = i + 1; k < end; k++) {
      if ((hc->eval_ext(i, end, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT, hc)) &&
          (fms5[s][k + 1] != INF) &&
          (c[idx[k] + i] != INF)) {
        type = vrna_get_ptype_md(S2[i], S2[k], md);

        switch (dangle_model) {
          case 2:
            s5  = ((i > 1) && (sn[i - 1] == sn[i])) ? S1[i - 1] : -1;
            s3  = ((k < n) && (sn[k] == sn[k + 1])) ? S1[k + 1] : -1;
            break;

          default:
            s5 = s3 = -1;
            break;
        }

        base = vrna_E_exterior_stem(type, s5, s3, params);

        tmp = base +
              fms5[s][k + 1] +
              c[idx[k] + i];

        if (sc_red_stem)
          tmp += sc_red_stem(i, k, i, k, sc_wrapper);

        if (sc_spl)
          tmp += sc_spl(i, end, k, k + 1, sc_wrapper);

        e = MIN2(e, tmp);
      }
    }

    if ((hc->eval_ext(i, end, i, end, VRNA_DECOMP_EXT_STEM, hc)) &&
        (c[idx[end] + i] != INF)) {
      type = vrna_get_ptype_md(S2[i], S2[end], md);
      switch (dangle_model) {
        case 2:
          s5  = ((i > 1) && (sn[i - 1] == sn[i])) ? S1[i - 1] : -1;
          s3  = -1;
          break;

        default:
          s5 = s3 = -1;
          break;
      }

      base = vrna_E_exterior_stem(type, s5, s3, params);

      tmp = base +
            c[idx[end] + i];

      if (sc_red_stem)
        tmp += sc_red_stem(i, end, i, end, sc_wrapper);

      e = MIN2(e, tmp);
    }

    /* damn odd dangles below */
    if (dangle_model % 2) {
      s5 = S1[i];

      for (int k = i + 2; k + 1 < end; k++) {
        if ((hc->eval_ext(i, end, k, k + 1, VRNA_DECOMP_EXT_STEM_EXT1, hc)) &&
            (fms5[s][k + 1] != INF) &&
            (c[idx[k] + i + 1] != INF)) {
          type  = vrna_get_ptype_md(S2[i + 1], S2[k], md);
          base  = vrna_E_exterior_stem(type, s5, -1, params);
          tmp   = base +
                  fms5[s][k + 1] +
                  c[idx[k] + i + 1];

          if (sc_red_stem)
            tmp += sc_red_stem(i, k, i + 1, k, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(i, end, k, k + 1, sc_wrapper);

          e = MIN2(e, tmp);
        }
      }

      for (int k = i + 1; k + 1 < end; k++) {
        s3 = S1[k + 1];

        if ((hc->eval_ext(i, end, k, k + 2, VRNA_DECOMP_EXT_STEM_EXT, hc)) &&
            (fms5[s][k + 2] != INF) &&
            (c[idx[k] + i] != INF)) {
          type  = vrna_get_ptype_md(S2[i], S2[k], md);
          base  = vrna_E_exterior_stem(type, -1, s3, params);
          tmp   = base +
                  fms5[s][k + 2] +
                  c[idx[k] + i];

          if (sc_red_stem)
            tmp += sc_red_stem(i, k + 1, i, k, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(i, end, k + 1, k + 2, sc_wrapper);

          e = MIN2(e, tmp);
        }

        if ((hc->eval_ext(i, end, k, k + 2, VRNA_DECOMP_EXT_STEM_EXT1, hc)) &&
            (fms5[s][k + 2] != INF) &&
            (c[idx[k] + i + 1] != INF) &&
            (i + 1 < k)) {
          type  = vrna_get_ptype_md(S2[i + 1], S2[k], md);
          base  = vrna_E_exterior_stem(type, s5, s3, params);
          tmp   = base +
                  fms5[s][k + 2] +
                  c[idx[k] + i + 1];

          if (sc_red_stem)
            tmp += sc_red_stem(i, k + 1, i + 1, k, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(i, end, k + 1, k + 2, sc_wrapper);

          e = MIN2(e, tmp);
        }
      } /* end for (k... */

      s5  = S1[i];
      s3  = S1[end];

      if ((hc->eval_ext(i, end, i, end - 1, VRNA_DECOMP_EXT_STEM, hc)) &&
          (c[idx[end - 1] + i] != INF) &&
          (i + 1 < end)) {
        type  = vrna_get_ptype_md(S2[i], S2[end - 1], md);
        base  = vrna_E_exterior_stem(type, -1, s3, params);
        tmp   = base +
                c[idx[end - 1] + i];

        if (sc_red_stem)
          tmp += sc_red_stem(i, end, i, end - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }

      if ((hc->eval_ext(i, end, i + 1, end, VRNA_DECOMP_EXT_STEM, hc)) &&
          (c[idx[end] + i + 1] != INF) &&
          (i + 1 < end)) {
        type  = vrna_get_ptype_md(S2[i + 1], S2[end], md);
        base  = vrna_E_exterior_stem(type, s5, -1, params);
        tmp   = base +
                c[idx[end] + i + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(i, end, i + 1, end, sc_wrapper);

        e = MIN2(e, tmp);
      }

      if ((hc->eval_ext(i, end, i + 1, end - 1, VRNA_DECOMP_EXT_STEM, hc)) &&
          (c[idx[end - 1] + i + 1] != INF) &&
          (i + 2 < end)) {
        type  = vrna_get_ptype_md(S2[i + 1], S2[end - 1], md);
        base  = vrna_E_exterior_stem(type, s5, s3, params);
        tmp   = base +
                c[idx[end - 1] + i + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(i, end, i + 1, end - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }
    }

    fms5[s][i] = e;
  }
}

PRIVATE void
update_fms3_arrays(vrna_fold_compound_t *fc,
                   unsigned int         s,
                   struct ms_helpers    *ms_dat)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *ss, type;
  int                   *c, **fms3, base, e, tmp, j, k, start, n, *idx, dangle_model;
  vrna_param_t          *params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  struct sc_f5_dat      *sc_wrapper;
  sc_ext_red_cb         sc_spl;
  sc_ext_red_cb         sc_red_stem;
  sc_ext_red_cb         sc_red_ext;

  n             = fc->length;
  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  idx           = fc->jindx;
  c             = fc->matrices->c;
  fms3          = fc->matrices->fms3;
  sn            = fc->strand_number;
  ss            = fc->strand_start;
  params        = fc->params;
  md            = &(params->model_details);
  dangle_model  = md->dangles;
  hc            = fc->hc;
  sc_wrapper    = &(ms_dat->sc_wrapper);
  sc_spl        = sc_wrapper->decomp;
  sc_red_stem   = sc_wrapper->red_stem;
  sc_red_ext    = sc_wrapper->red_ext;

  start = ss[s];

  for (j = start; j <= n; j++) {
    e = tmp = INF;

    if (start < j) {
      if (hc->eval_ext(start, j, start, j - 1, VRNA_DECOMP_EXT_EXT, hc)) {
        tmp = fms3[s][j - 1];

        if (sc_red_ext)
          tmp += sc_red_ext(start, j, start, j - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }
    } else {
      e = 0;
    }

    if (hc->eval_ext(start, j, start, j, VRNA_DECOMP_EXT_STEM, hc) &&
        (c[idx[j] + start] != INF)) {
      type = vrna_get_ptype_md(S2[start], S2[j], md);

      switch (dangle_model) {
        case 2:
          s5  = -1;
          s3  = (sn[j] == sn[j + 1]) ? S1[j + 1] : -1;
          break;

        default:
          s5 = s3 = -1;
          break;
      }

      base = vrna_E_exterior_stem(type, s5, s3, params);

      tmp = base +
            c[idx[j] + start];

      if (sc_red_stem)
        tmp += sc_red_stem(start, j, start, j, sc_wrapper);

      e = MIN2(e, tmp);
    }

    for (k = start; k < j; k++) {
      if (hc->eval_ext(start, j, k, k + 1, VRNA_DECOMP_EXT_EXT_STEM, hc) &&
          (fms3[s][k] != INF) &&
          (c[idx[j] + k + 1] != INF)) {
        type = vrna_get_ptype_md(S2[k + 1], S2[j], md);

        switch (dangle_model) {
          case 2:
            s5  = (sn[k] == sn[k + 1]) ? S1[k] : -1;
            s3  = (sn[j] == sn[j + 1]) ? S1[j + 1] : -1;
            break;

          default:
            s5 = s3 = -1;
            break;
        }

        base = vrna_E_exterior_stem(type, s5, s3, params);

        tmp = base +
              fms3[s][k] +
              c[idx[j] + k + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(k + 1, j, k + 1, j, sc_wrapper);

        if (sc_spl)
          tmp += sc_spl(start, j, k, k + 1, sc_wrapper);

        e = MIN2(e, tmp);
      }
    }

    if (dangle_model % 2) {
      s5  = S1[start];
      s3  = S1[j];

      if (hc->eval_ext(start, j, start, j - 1, VRNA_DECOMP_EXT_STEM, hc) &&
          (c[idx[j - 1] + start] != INF) &&
          (start + 1 < j)) {
        type  = vrna_get_ptype_md(S2[start], S2[j - 1], md);
        base  = vrna_E_exterior_stem(type, -1, s3, params);
        tmp   = base +
                c[idx[j - 1] + start];

        if (sc_red_stem)
          tmp += sc_red_stem(start, j, start, j - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }

      if (hc->eval_ext(start, j, start + 1, j, VRNA_DECOMP_EXT_STEM, hc) &&
          (c[idx[j] + start + 1] != INF) &&
          (start + 1 < j)) {
        type  = vrna_get_ptype_md(S2[start + 1], S2[j], md);
        base  = vrna_E_exterior_stem(type, s5, -1, params);
        tmp   = base +
                c[idx[j] + start + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(start, j, start + 1, j, sc_wrapper);

        e = MIN2(e, tmp);
      }

      if (hc->eval_ext(start, j, start + 1, j - 1, VRNA_DECOMP_EXT_STEM, hc) &&
          (c[idx[j - 1] + start + 1] != INF) &&
          (start + 2 < j)) {
        type  = vrna_get_ptype_md(S2[start + 1], S2[j - 1], md);
        base  = vrna_E_exterior_stem(type, s5, s3, params);
        tmp   = base +
                c[idx[j - 1] + start + 1];

        if (sc_red_stem)
          tmp += sc_red_stem(start, j, start + 1, j - 1, sc_wrapper);

        e = MIN2(e, tmp);
      }

      for (k = start; k < j; k++) {
        s5 = S1[k + 1];

        if (hc->eval_ext(start, j, k, k + 2, VRNA_DECOMP_EXT_EXT_STEM, hc) &&
            (fms3[s][k] != INF) &&
            (c[idx[j] + k + 2] != INF) &&
            k + 2 < j) {
          type  = vrna_get_ptype_md(S2[k + 2], S2[j], md);
          base  = vrna_E_exterior_stem(type, s5, -1, params);
          tmp   = base +
                  fms3[s][k] +
                  c[idx[j] + k + 2];

          if (sc_red_stem)
            tmp += sc_red_stem(k + 1, j, k + 2, j, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(start, j, k, k + 1, sc_wrapper);

          e = MIN2(e, tmp);
        }

        if (hc->eval_ext(start, j, k, k + 1, VRNA_DECOMP_EXT_EXT_STEM1, hc) &&
            (fms3[s][k] != INF) &&
            (c[idx[j - 1] + k + 1] != INF) &&
            (k + 2 < j)) {
          type  = vrna_get_ptype_md(S2[k + 1], S2[j - 1], md);
          base  = vrna_E_exterior_stem(type, -1, s3, params);
          tmp   = base +
                  fms3[s][k] +
                  c[idx[j - 1] + k + 1];

          if (sc_red_stem)
            tmp += sc_red_stem(k + 1, j, k + 1, j - 1, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(start, j, k, k + 1, sc_wrapper);

          e = MIN2(e, tmp);
        }

        if (hc->eval_ext(start, j, k, k + 2, VRNA_DECOMP_EXT_EXT_STEM1, hc) &&
            (fms3[s][k] != INF) &&
            (c[idx[j - 1] + k + 2] != INF) &&
            (k + 3 < j)) {
          type  = vrna_get_ptype_md(S2[k + 2], S2[j - 1], md);
          base  = vrna_E_exterior_stem(type, s5, s3, params);
          tmp   = base +
                  fms3[s][k] +
                  c[idx[j - 1] + k + 2];

          if (sc_red_stem)
            tmp += sc_red_stem(k + 1, j, k + 2, j - 1, sc_wrapper);

          if (sc_spl)
            tmp += sc_spl(start, j, k, k + 1, sc_wrapper);

          e = MIN2(e, tmp);
        }
      }
    }

    fms3[s][j] = e;
  }
}

PRIVATE int
pair_multi_strand(vrna_fold_compound_t  *fc,
                  int                   i,
                  int                   j,
                  struct ms_helpers     *ms_dat)
{
  short                 *S1, *S2, s5, s3;
  unsigned int          *sn, *ends, type, nick;
  int                   start, contribution, **fms5, **fms3, base, tmp, tmp2, dangle_model;
  vrna_param_t          *params;
  vrna_md_t             *md;
  vrna_hc_t             *hc;
  struct sc_f5_dat      *sc_wrapper;
  sc_ext_red_cb         sc_red_stem;

  contribution  = INF;
  S1            = fc->sequence_encoding;
  S2            = fc->sequence_encoding2;
  params        = fc->params;
  md            = &(params->model_details);
  dangle_model  = md->dangles;
  sn            = fc->strand_number;
  ends          = fc->strand_end;
  fms5          = fc->matrices->fms5;
  fms3          = fc->matrices->fms3;
  hc            = fc->hc;
  sc_wrapper    = &(ms_dat->sc_wrapper);
  sc_red_stem   = sc_wrapper->red_stem;

  if ((sn[i] != sn[j]) &&
      (hc->eval_ext(i, j, i, j, VRNA_DECOMP_EXT_STEM, hc))) {
    /* most obious strand nick is at end of sn[i] and start of sn[j] */
    type = vrna_get_ptype_md(S2[j], S2[i], md);

    switch (dangle_model) {
      case 2:
        s5  = (sn[j] == sn[j - 1]) ? S1[j - 1] : -1;
        s3  = (sn[i] == sn[i + 1]) ? S1[i + 1] : -1;
        break;

      default:
        s5 = s3 = -1;
        break;
    }

    base = vrna_E_exterior_stem(type, s5, s3, params) +
           params->DuplexInit;

    if (sc_red_stem)
      base += sc_red_stem(j, i, j, i, sc_wrapper);

    tmp = INF;

    /*
     *  if (evaluate(i + 1,
     *               j - 1,
     *               ends[sn[i]],
     *               ends[sn[i]] + 1,
     *               VRNA_DECOMP_EXT_EXT_EXT,
     *               &hc_dat_local))
     */
    if (sn[i] != sn[i + 1]) {
      if ((sn[j - 1] != sn[j]) &&
          (i + 1 == j)) {
        tmp2  = 0;
        tmp   = MIN2(tmp, tmp2);
      } else if (sn[j - 1] == sn[j]) {
        tmp2  = fms3[sn[i + 1]][j - 1];
        tmp   = MIN2(tmp, tmp2);
      }
    } else if (sn[j - 1] != sn[j]) {
      tmp2  = fms5[sn[j - 1]][i + 1];
      tmp   = MIN2(tmp, tmp2);
    } else {
      if ((fms5[sn[i]][i + 1] != INF) &&
          (fms3[sn[ends[sn[i]] + 1]][j - 1] != INF)) {
        tmp2 = 0;

        if (ends[sn[i]] > (unsigned int)i)
          tmp2 += fms5[sn[i]][i + 1];

        if ((unsigned int)(j - 1) > ends[sn[i]])
          tmp2 += fms3[sn[ends[sn[i]] + 1]][j - 1];

        tmp = MIN2(tmp, tmp2);
      }

      /* check whether we find more strand nicks between i and j */
      nick = ends[sn[i]] + 1;
      while (sn[nick] != sn[j]) {
        /*
         *    if (evaluate(i + 1,
         *                 j - 1,
         *                 ends[sn[nick]],
         *                 ends[sn[nick]] + 1,
         *                 VRNA_DECOMP_EXT_EXT_EXT,
         *                 &hc_dat_local))
         */
        if ((fms5[sn[nick]][i + 1] != INF) &&
            (fms3[sn[ends[sn[nick]] + 1]][j - 1] != INF)) {
          tmp2 = 0;
          if ((unsigned int)(i + 1) <= ends[sn[nick]])
            tmp2 += fms5[sn[nick]][i + 1];

          if (ends[sn[nick]] + 1 <= (unsigned int)(j - 1))
            tmp2 += fms3[sn[ends[sn[nick]] + 1]][j - 1];

          tmp = MIN2(tmp, tmp2);
        }

        nick = ends[sn[nick]] + 1;
      }
    }

    if (tmp != INF)
      contribution = base + tmp;

    /* odd dangles below */
    if (dangle_model % 2) {
      s5  = (sn[j] == sn[j - 1]) ? S1[j - 1] : -1;
      s3  = (sn[i] == sn[i + 1]) ? S1[i + 1] : -1;

      if ((j > i + 1) &&
          (sn[i] != sn[i + 1]) &&
          (sn[j - 1] == sn[j])) {
        tmp = vrna_E_exterior_stem(type, s5, -1, params) +
              params->DuplexInit;

        if (sc_red_stem)
          tmp += sc_red_stem(j - 1, i, j, i, sc_wrapper);

        if (sn[j - 2] == sn[j]) {
          if (fms3[sn[i + 1]][j - 2] != INF) {
            tmp           += fms3[sn[i + 1]][j - 2];
            contribution  = MIN2(contribution, tmp);
          }
        } else {
          contribution = MIN2(contribution, tmp);
        }
      } else if ((i + 1 < j) &&
                 (sn[j - 1] != sn[j]) &&
                 (sn[i] == sn[i + 1])) {
        tmp = vrna_E_exterior_stem(type, -1, s3, params) +
              params->DuplexInit;

        if (sc_red_stem)
          tmp += sc_red_stem(j, i + 1, j, i, sc_wrapper);

        if (sn[i] == sn[i + 2]) {
          if (fms5[sn[j - 1]][i + 2] != INF) {
            tmp           += fms5[sn[j - 1]][i + 2];
            contribution  = MIN2(contribution, tmp);
          }
        } else {
          contribution = MIN2(contribution, tmp);
        }
      } else if ((sn[i] == sn[i + 1]) &&
                 (sn[j - 1] == sn[j])) {
        /* 1st case, i + 1 and j - 1 dangle on enclosing pair */
        tmp   = INF;
        base  = vrna_E_exterior_stem(type, s5, s3, params) +
                params->DuplexInit;

        if (sc_red_stem)
          base += sc_red_stem(j - 1, i + 1, j, i, sc_wrapper);

        start = i;
        nick  = ends[sn[start]] + 1;
        do {
          if ((fms5[sn[start]][i + 2] != INF) &&
              (fms3[sn[nick]][j - 2] != INF)) {
            tmp2 = 0;

            if (nick > (unsigned int)(i + 2))
              tmp2 += fms5[sn[start]][i + 2];

            if ((unsigned int)j > nick + 1)
              tmp2 += fms3[sn[nick]][j - 2];

            tmp = MIN2(tmp, tmp2);
          }

          start = nick;
          nick  = ends[sn[nick]] + 1;
        } while (sn[nick] != sn[j]);

        if (tmp != INF) {
          tmp           += base;
          contribution  = MIN2(contribution, tmp);
        }

        /* 2nd case, only i + 1 may dangle on enclosing pair */
        tmp   = INF;
        base  = vrna_E_exterior_stem(type, -1, s3, params) +
                params->DuplexInit;

        if (sc_red_stem)
          base += sc_red_stem(j, i + 1, j, i, sc_wrapper);

        start = i;
        nick  = ends[sn[start]] + 1;
        do {
          if ((fms5[sn[start]][i + 2] != INF) &&
              (fms3[sn[nick]][j - 1] != INF)) {
            tmp2 = 0;

            if (nick > (unsigned int)(i + 2))
              tmp2 += fms5[sn[start]][i + 2];

            if ((unsigned int)j > nick)
              tmp2 += fms3[sn[nick]][j - 1];

            tmp = MIN2(tmp, tmp2);
          }

          start = nick;
          nick  = ends[sn[nick]] + 1;
        } while (sn[nick] != sn[j]);

        if (tmp != INF) {
          tmp           += base;
          contribution  = MIN2(contribution, tmp);
        }

        /* 3rd case, only j - 1 may dangle on enclosing pair */
        tmp   = INF;
        base  = vrna_E_exterior_stem(type, s5, -1, params) +
                params->DuplexInit;

        if (sc_red_stem)
          base += sc_red_stem(j - 1, i, j, i, sc_wrapper);

        start = i;
        nick  = ends[sn[start]] + 1;
        do {
          if ((fms5[sn[start]][i + 1] != INF) &&
              (fms3[sn[nick]][j - 2] != INF)) {
            tmp2 = 0;

            if (nick > (unsigned int)(i + 1))
              tmp2 += fms5[sn[start]][i + 1];

            if ((unsigned int)j > nick + 1)
              tmp2 += fms3[sn[nick]][j - 2];

            tmp = MIN2(tmp, tmp2);
          }

          start = nick;
          nick  = ends[sn[nick]] + 1;
        } while (sn[nick] != sn[j]);

        if (tmp != INF) {
          tmp           += base;
          contribution  = MIN2(contribution, tmp);
        }
      }
    }
  }

  return contribution;
}

#endif

/**
*** trace back through the "c", "f5" and "fML" arrays to get the
*** base pairing list. No search for equivalent structures is done.
*** This is fast, since only few structure elements are recalculated.
***
*** normally s=0.
*** If s>0 then s items have been already pushed onto the bt_stack
**/

PRIVATE INLINE int
decompose_pair(vrna_fold_compound_t *fc,
               int                  i,
               int                  j,
               struct aux_arrays    *aux,
               struct ms_helpers    *ms_dat)
{
  unsigned char hc_decompose;
  unsigned int  n;
  int           e, new_c, energy, stackEnergy, ij, dangle_model, noLP, *cc, *cc1;

  n             = fc->length;
  ij            = fc->jindx[j] + i;
  dangle_model  = fc->params->model_details.dangles;
  noLP          = fc->params->model_details.noLP;
  hc_decompose  = fc->hc->mx[n * i + j];
  cc            = aux->cc;
  cc1           = aux->cc1;
  e             = INF;

  /* do we evaluate this pair? */
  if (hc_decompose) {
    new_c = INF;

    /* check for hairpin loop */
    energy  = vrna_eval_hairpin(fc, i, j, VRNA_EVAL_LOOP_DEFAULT);
    new_c   = MIN2(new_c, energy);

    /* check for multibranch loops */
    energy  = vrna_mfe_multibranch_loop_fast(fc, i, j, aux->ml_helpers);
    new_c   = MIN2(new_c, energy);

    if (dangle_model == 3) {
      /* coaxial stacking */
      energy  = vrna_mfe_multibranch_loop_stack(fc, i, j);
      new_c   = MIN2(new_c, energy);
    }

    /* check for internal loops */
    energy  = vrna_mfe_internal(fc, i, j);
    new_c   = MIN2(new_c, energy);

    /* multi-strand decomposition */
    if (fc->strands > 1) {
      energy  = pair_multi_strand(fc, i, j, ms_dat);
      new_c   = MIN2(new_c, energy);
    }

    /* remember stack energy for --noLP option */
    if (noLP) {
      stackEnergy = vrna_eval_stack(fc, i, j, VRNA_EVAL_LOOP_DEFAULT);
      new_c       = MIN2(new_c, cc1[j - 1] + stackEnergy);
      cc[j]       = new_c;
      if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
          (cc[j] != INF))
        cc[j] -= fc->pscore[ij];

      e = cc1[j - 1] + stackEnergy;
    } else {
      e = new_c;
    }

    /* finally, check for auxiliary grammar rule(s) */
    if (fc->aux_grammar) {
      for (size_t c = 0; c < vrna_array_size(fc->aux_grammar->c); c++) {
        if (fc->aux_grammar->c[c].cb) {
          energy  = fc->aux_grammar->c[c].cb(fc, i, j, fc->aux_grammar->c[c].data);
          e       = MIN2(e, energy);
        }
      }
    }

    if ((fc->type == VRNA_FC_TYPE_COMPARATIVE) &&
        (e != INF))
      e -= fc->pscore[ij];
  } /* end >> if (pair) << */

  return e;
}

PRIVATE INLINE struct aux_arrays *
get_aux_arrays(unsigned int length)
{
  struct aux_arrays *aux = (struct aux_arrays *)vrna_alloc(sizeof(struct aux_arrays));

  aux->cc   = (int *)vrna_alloc(sizeof(int) * (length + 2));    /* auxilary arrays for canonical structures     */
  aux->cc1  = (int *)vrna_alloc(sizeof(int) * (length + 2));    /* auxilary arrays for canonical structures     */

  /* get fast multiloop decomposition helpers */
  aux->ml_helpers = vrna_mfe_multibranch_fast_init(length);

  return aux;
}

PRIVATE INLINE void
rotate_aux_arrays(struct aux_arrays *aux,
                  unsigned int      length)
{
  unsigned int  j;
  int           *FF;

  FF        = aux->cc1;
  aux->cc1  = aux->cc;
  aux->cc   = FF;
  for (j = 1; j <= length; j++)
    aux->cc[j] = INF;

  vrna_mfe_multibranch_fast_rotate(aux->ml_helpers);
}

PRIVATE INLINE void
free_aux_arrays(struct aux_arrays *aux)
{
  free(aux->cc);
  free(aux->cc1);
  vrna_mfe_multibranch_fast_free(aux->ml_helpers);

  free(aux);
}
