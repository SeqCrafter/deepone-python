/*
 *                Model Details structure creation/modification/destruction
 *
 *                This file contains everything which is necessary to
 *                obtain, modify, and destroy the model_details datastructure
 *                used in the folding recurrences throughout the ViennaRNA
 *                Package
 *
 *                c Ronny Lorenx
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>

#include "ViennaRNA/utils/basic.h"

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*  below are the evil global variables that will vanish
 *  as soon as we drop backward compatibility in ViennaRNA
 *  Package v3
 */

/* below are some more deprecated global symbols we need to get rid off */

    /* internal loops of size 2 get energy 0.8Kcal and
                                           * no mismatches (no longer used) */
 /* TODO: compile ribosums into program
                                           * Warning: this variable will vanish */
    /*generate comma seperated output*/

 /* base pairing prob. matrix */
 /* pr[i,j] -> pr[iindx[i]-j] */
    /* fold with constraints */

#endif

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

#define   BP_REV_DEFAULT        { 0, 2, 1, 4, 3, 6, 5, 7 }

#define   BP_ALIAS_DEFAULT      { 0, 1, 2, 3, 4, 3, 2, 0 }

#define   BP_ENCODING_DEFAULT \
  /*  _  A  C  G  U  X  K  I */ \
  { { 0, 0, 0, 0, 0, 0, 0, 0 }, \
    { 0, 0, 0, 0, 5, 0, 0, 5 }, \
    { 0, 0, 0, 1, 0, 0, 0, 0 }, \
    { 0, 0, 2, 0, 3, 0, 0, 0 }, \
    { 0, 6, 0, 4, 0, 0, 0, 6 }, \
    { 0, 0, 0, 0, 0, 0, 2, 0 }, \
    { 0, 0, 0, 0, 0, 1, 0, 0 }, \
    { 0, 6, 0, 0, 5, 0, 0, 0 } }

#define   DM_DEFAULT \
  { { 0, 0, 0, 0, 0, 0, 0 }, /* hamming distance between pairs */ \
    { 0, 0, 2, 2, 1, 2, 2 } /* CG */, \
    { 0, 2, 0, 1, 2, 2, 2 } /* GC */, \
    { 0, 2, 1, 0, 2, 1, 2 } /* GU */, \
    { 0, 1, 2, 2, 0, 2, 1 } /* UG */, \
    { 0, 2, 2, 1, 2, 0, 2 } /* AU */, \
    { 0, 2, 2, 2, 1, 2, 0 } /* UA */ }

PRIVATE vrna_md_t defaults = {
  .temperature = VRNA_MODEL_DEFAULT_TEMPERATURE,
  .betaScale = 1.,
  .pf_smooth = VRNA_MODEL_DEFAULT_PF_SMOOTH,
  .dangles  = VRNA_MODEL_DEFAULT_DANGLES,
  .special_hp  = VRNA_MODEL_DEFAULT_SPECIAL_HP,
  .noLP  = VRNA_MODEL_DEFAULT_NO_LP,
  .noGU  = VRNA_MODEL_DEFAULT_NO_GU,
  .noGUclosure  = VRNA_MODEL_DEFAULT_NO_GU_CLOSURE,
  .logML  = VRNA_MODEL_DEFAULT_LOG_ML,
  .circ  = VRNA_MODEL_DEFAULT_CIRC,
  .circ_penalty  = VRNA_MODEL_DEFAULT_CIRC_PENALTY,
  .gquad  = VRNA_MODEL_DEFAULT_GQUAD,
  .uniq_ML  = VRNA_MODEL_DEFAULT_UNIQ_ML,
  .energy_set  = VRNA_MODEL_DEFAULT_ENERGY_SET,
  .backtrack  = VRNA_MODEL_DEFAULT_BACKTRACK,
  .backtrack_type  = VRNA_MODEL_DEFAULT_BACKTRACK_TYPE,
  .compute_bpp  = VRNA_MODEL_DEFAULT_COMPUTE_BPP,
  .nonstandards  = { 0 },
  .max_bp_span  = VRNA_MODEL_DEFAULT_MAX_BP_SPAN,
  .min_loop_size  = TURN,
  .window_size  = VRNA_MODEL_DEFAULT_WINDOW_SIZE,
  .oldAliEn  = VRNA_MODEL_DEFAULT_ALI_OLD_EN,
  .ribo  = VRNA_MODEL_DEFAULT_ALI_RIBO,
  .cv_fact  = VRNA_MODEL_DEFAULT_ALI_CV_FACT,
  .nc_fact  = VRNA_MODEL_DEFAULT_ALI_NC_FACT,
  .sfact  = 1.07,
  .rtype  = BP_REV_DEFAULT,
  .alias  = BP_ALIAS_DEFAULT,
  .pair  = BP_ENCODING_DEFAULT,
  .pair_dist  = DM_DEFAULT,
  .salt  = VRNA_MODEL_DEFAULT_SALT,
  .saltMLLower  = VRNA_MODEL_DEFAULT_SALT_MLLOWER,
  .saltMLUpper  = VRNA_MODEL_DEFAULT_SALT_MLUPPER,
  .saltDPXInit  = VRNA_MODEL_DEFAULT_SALT_DPXINIT,
  .saltDPXInitFact  = VRNA_MODEL_DEFAULT_SALT_DPXINIT_FACT,
  .helical_rise  = VRNA_MODEL_DEFAULT_HELICAL_RISE,
  .backbone_length  = VRNA_MODEL_DEFAULT_BACKBONE_LENGTH,
  .circ_alpha0  = VRNA_MODEL_DEFAULT_CIRC_ALPHA0
};

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

/* Fill the base pair type encodings according to the model details */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_md_t *
vrna_md_copy(vrna_md_t        *md_to,
             const vrna_md_t  *md_from)
{
  int       i;
  vrna_md_t *md;

  md = NULL;

  /* only process if md_from is non-NULL */
  if (md_from) {
    if (!md_to)
      /* create container to be filled */
      md = (vrna_md_t *)vrna_alloc(sizeof(vrna_md_t));
    else
      /* or directly write to target */
      md = md_to;

    /* check if not the same object */
    if (md_to != md_from) {
      /* copy simple members */
      memcpy(md, md_from, sizeof(vrna_md_t));
      /* copy arrays */
      memcpy(md->rtype, &(md_from->rtype[0]), 8 * sizeof(int));
      memcpy(md->alias, &(md_from->alias[0]), (MAXALPHA + 1) * sizeof(short));
      memcpy(md->nonstandards, &(md_from->nonstandards[0]), 64 * sizeof(char));
      /* copy matrices */
      for (i = 0; i <= MAXALPHA; i++)
        memcpy(md->pair[i], (md_from->pair[i]), (MAXALPHA + 1) * sizeof(int));
      /* copy pair dists */
      for (i = 0; i <= 6; i++)
        memcpy(md->pair_dist[i], (md_from->pair_dist[i]), 7 * sizeof(float));
    }
  }

  return md;
}

PUBLIC void
vrna_md_set_default(vrna_md_t *md)
{
  if (md) /* copy defaults */
    vrna_md_copy(md, &defaults);
}

/* below are the setter functions for global default settings */

/*
 *  This function updates the pair/alias/rtype arrays according to model settings.
 *  It should be called whenever there is a change in the following model settings:
 *  - energy_set
 *  - noGU
 *  - nonstandards
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 * ###########################################
 * # deprecated functions below              #
 *###########################################
 */

#endif
