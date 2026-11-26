/*
 *
 *                c Ivo Hofacker
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/salt.h"

#define SPEEDUP_PARAMS    1

/**
 *** \file ViennaRNA/params/basic.c
 *** <P>
 *** This file provides functions that return temperature scaled energy parameters and
 *** Boltzmann weights packed in datastructures
 *** </P>
 ***/

/*------------------------------------------------------------------------*/
#define SCALE 10
/**
 *** dangling ends should never be destabilizing, i.e. expdangle>=1<BR>
 *** specific heat needs smooth function (2nd derivative)<BR>
 *** we use a*(sin(x+b)+1)^2, with a=2/(3*sqrt(3)), b=Pi/6-sqrt(3)/2,
 *** in the interval b<x<sqrt(3)/2
 */
#define CLIP_NEGATIVE(X) ((X) < 0 ? 0 : (X))
#define SMOOTH(X) ((!pf_smooth)               ?   CLIP_NEGATIVE(X) : \
                   ((X) / SCALE < -1.2283697) ?                 0  : \
                   ((X) / SCALE > 0.8660254) ?                (X) : \
                   SCALE *0.38490018   \
                   * (sin((X) / SCALE - 0.34242663) + 1) \
                   * (sin((X) / SCALE - 0.34242663) + 1) \
                   )

/* #define SMOOTH(X) ((X)<0 ? 0 : (X)) */

/*
 * If the global use_mfelike_energies flag is set, truncate doubles to int
 * values and cast back to double. This makes the energy parameters of the
 * partition (folding get_scaled_exp_params()) compatible with the mfe folding
 * parameters (get_scaled_exp_params()), e.g. for explicit partition function
 * computations.
 */
#define TRUNC_MAYBE(X) ((!pf_smooth) ? (double)((int)(X)) : (X))

/* Rescale Free energy contribution according to deviation of temperature from measurement conditions */
#define RESCALE_dG(dG, dH, dT)   ((dH) - ((dH) - (dG)) * dT)

/*
 * Rescale Free energy contribution according to deviation of temperature from measurement conditions
 * and convert it to Boltzmann Factor for specific kT
 */
#define RESCALE_BF(dG, dH, dT, kT)          ( \
    exp( \
      -TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) \
      * 10. \
      / kT \
      ) \
    )

#define RESCALE_BF_SMOOTH(dG, dH, dT, kT)   ( \
    exp(  \
      SMOOTH( \
        -TRUNC_MAYBE((double)RESCALE_dG((dG), (dH), (dT))) \
        ) \
      * 10. \
      / kT \
      ) \
    )

#define saltT md->temperature+K0

/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */
PRIVATE vrna_param_t p;
PRIVATE int               id = -1;
/* variables for partition function */
PRIVATE vrna_exp_param_t  pf;

#if SPEEDUP_PARAMS

PRIVATE vrna_param_t      p_pre       = { 0 };
PRIVATE unsigned char     p_pre_init  = 0;

#endif

#ifdef _OPENMP
#pragma omp threadprivate(id, pf_id)
#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE vrna_param_t *
get_scaled_params(vrna_md_t *md);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC vrna_param_t *
vrna_params(vrna_md_t *md_p)
{
  vrna_param_t  *cp = NULL;
  vrna_md_t     md;

  if (md_p == NULL) {
    vrna_md_set_default(&md);
    md_p = &md;
  }

#if SPEEDUP_PARAMS
  if (p_pre_init) {
    /* ignore window_size, max_bp_span, and min_loop_size here */
    p_pre.model_details.window_size   = md_p->window_size;
    p_pre.model_details.min_loop_size = md_p->min_loop_size;
    p_pre.model_details.max_bp_span   = md_p->max_bp_span;

    if (memcmp(md_p, &(p_pre.model_details), sizeof(vrna_md_t)) == 0) {
      return vrna_params_copy(&p_pre);
    }
  }

  cp = get_scaled_params(md_p);

  /* store with current settings */
  memcpy(&p_pre, cp, sizeof(vrna_param_t));

  p_pre_init = (unsigned char)1;

  return cp;
#else
  return get_scaled_params(md_p);
#endif
}

PUBLIC vrna_param_t *
vrna_params_copy(vrna_param_t *par)
{
  vrna_param_t *copy = NULL;

  if (par) {
    copy = (vrna_param_t *)vrna_alloc(sizeof(vrna_param_t));
    memcpy(copy, par, sizeof(vrna_param_t));
  }

  return copy;
}

PUBLIC void
vrna_params_prepare(vrna_fold_compound_t  *fc,
                    unsigned int          options)
{
  if (fc) {
    vrna_md_t *md_p;

    /*
     *  every vrna_fold_compound_t must have a vrna_paramt_t structure attached
     *  to it that holds the current model details. So we just use this here as
     *  the reference model
     */
    md_p = &(fc->params->model_details);

    ;
  }
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE vrna_param_t *
get_scaled_params(vrna_md_t *md)
{
  unsigned int  i, j, k, l;
  double        tempf;
  vrna_param_t  *params;
  double salt = md->salt;
  double saltStandard = VRNA_MODEL_DEFAULT_SALT;

  params = (vrna_param_t *)vrna_alloc(sizeof(vrna_param_t));

  memset(params->param_file, '\0', 256);

  params->model_details = *md;  /* copy over the model details */
  params->temperature   = md->temperature;
  tempf                 = ((params->temperature + K0) / Tmeasure);

  params->ninio[2]              = RESCALE_dG(ninio37, niniodH, tempf);
  params->lxc                   = lxc37 * tempf;
  params->TripleC               = RESCALE_dG(TripleC37, TripleCdH, tempf);
  params->MultipleCA            = RESCALE_dG(MultipleCA37, MultipleCAdH, tempf);
  params->MultipleCB            = RESCALE_dG(MultipleCB37, MultipleCBdH, tempf);
  params->TerminalAU            = RESCALE_dG(TerminalAU37, TerminalAUdH, tempf);
  params->DuplexInit            = RESCALE_dG(DuplexInit37, DuplexInitdH, tempf);
  params->MLbase                = RESCALE_dG(ML_BASE37, ML_BASEdH, tempf);
  params->MLclosing             = RESCALE_dG(ML_closing37, ML_closingdH, tempf);
  params->gquadLayerMismatch    = RESCALE_dG(GQuadLayerMismatch37, GQuadLayerMismatchH, tempf);
  params->gquadLayerMismatchMax = GQuadLayerMismatchMax;

  for (i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++)
    for (j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3 * VRNA_GQUAD_MAX_LINKER_LENGTH; j++) {
      double  GQuadAlpha_T  = RESCALE_dG(GQuadAlpha37, GQuadAlphadH, tempf);
      double  GQuadBeta_T   = RESCALE_dG(GQuadBeta37, GQuadBetadH, tempf);
      double  GQuadEnergy   = (double)GQuadAlpha_T * (double)(i - 1) +
                              (double)GQuadBeta_T * log((double)j - 2.);
      params->gquad[i][j]   = (int)GQuadEnergy;
    }

  for (i = 0; i < 31; i++)
    params->hairpin[i] = RESCALE_dG(hairpin37[i], hairpindH[i], tempf);

  for (i = 0; i <= MIN2(30, MAXLOOP); i++) {
    params->bulge[i]          = RESCALE_dG(bulge37[i], bulgedH[i], tempf);
    params->internal_loop[i]  = RESCALE_dG(internal_loop37[i], internal_loopdH[i], tempf);
  }

  for (; i <= MAXLOOP; i++) {
    params->bulge[i] = params->bulge[30] +
                       (int)(params->lxc * log((double)(i) / 30.));
    params->internal_loop[i] = params->internal_loop[30] +
                               (int)(params->lxc * log((double)(i) / 30.));
    params->SaltLoopDbl[i]   = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i]      = (int) (params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i]<0));
  }

  for (i = 0; i <= MIN2(31, MAXLOOP+1); i++) {
    params->SaltLoopDbl[i]    = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i]       = (int) (params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i]<0));
  }

  for (; i <= MAXLOOP; i++) {
    params->SaltLoopDbl[i]   = (salt==saltStandard) ? 0. : vrna_salt_loop(i, salt, saltT, md->backbone_length);
    params->SaltLoop[i]      = (int) (params->SaltLoopDbl[i] + 0.5 - (params->SaltLoopDbl[i]<0));
  }

  for (i = 0; (i * 7) < strlen(Tetraloops); i++)
    params->Tetraloop_E[i] = RESCALE_dG(Tetraloop37[i], TetraloopdH[i], tempf);

  for (i = 0; (i * 5) < strlen(Triloops); i++)
    params->Triloop_E[i] = RESCALE_dG(Triloop37[i], TriloopdH[i], tempf);

  for (i = 0; (i * 9) < strlen(Hexaloops); i++)
    params->Hexaloop_E[i] = RESCALE_dG(Hexaloop37[i], HexaloopdH[i], tempf);

  for (i = 0; i <= NBPAIRS; i++)
    params->MLintern[i] = RESCALE_dG(ML_intern37, ML_interndH, tempf);

  /* stacks    G(T) = H - [H - G(T0)]*T/T0 */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      params->stack[i][j] = RESCALE_dG(stack37[i][j],
                                       stackdH[i][j],
                                       tempf);

  /* mismatches */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++)
      for (k = 0; k < 5; k++) {
        int mm;
        params->mismatchI[i][j][k] = RESCALE_dG(mismatchI37[i][j][k],
                                                mismatchIdH[i][j][k],
                                                tempf);
        params->mismatchH[i][j][k] = RESCALE_dG(mismatchH37[i][j][k],
                                                mismatchHdH[i][j][k],
                                                tempf);
        params->mismatch1nI[i][j][k] = RESCALE_dG(mismatch1nI37[i][j][k],
                                                  mismatch1nIdH[i][j][k],
                                                  tempf);
        params->mismatch23I[i][j][k] = RESCALE_dG(mismatch23I37[i][j][k],
                                                  mismatch23IdH[i][j][k],
                                                  tempf);
        if (md->dangles) {
          mm = RESCALE_dG(mismatchM37[i][j][k],
                          mismatchMdH[i][j][k],
                          tempf);
          params->mismatchM[i][j][k]  = (mm > 0) ? 0 : mm;
          mm                          = RESCALE_dG(mismatchExt37[i][j][k],
                                                   mismatchExtdH[i][j][k],
                                                   tempf);
          params->mismatchExt[i][j][k] = (mm > 0) ? 0 : mm;
        } else {
          params->mismatchM[i][j][k] = params->mismatchExt[i][j][k] = 0;
        }
      }

  /* dangles */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j < 5; j++) {
      int dd;
      dd = RESCALE_dG(dangle5_37[i][j],
                      dangle5_dH[i][j],
                      tempf);
      params->dangle5[i][j] = (dd > 0) ? 0 : dd;  /* must be <= 0 */
      dd                    = RESCALE_dG(dangle3_37[i][j],
                                         dangle3_dH[i][j],
                                         tempf);
      params->dangle3[i][j] = (dd > 0) ? 0 : dd;  /* must be <= 0 */
    }

  /* internal 1x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++)
          params->int11[i][j][k][l] = RESCALE_dG(int11_37[i][j][k][l],
                                                 int11_dH[i][j][k][l],
                                                 tempf);

  /* internal 2x1 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m;
          for (m = 0; m < 5; m++)
            params->int21[i][j][k][l][m] = RESCALE_dG(int21_37[i][j][k][l][m],
                                                      int21_dH[i][j][k][l][m],
                                                      tempf);
        }

  /* internal 2x2 loops */
  for (i = 0; i <= NBPAIRS; i++)
    for (j = 0; j <= NBPAIRS; j++)
      for (k = 0; k < 5; k++)
        for (l = 0; l < 5; l++) {
          int m, n;
          for (m = 0; m < 5; m++)
            for (n = 0; n < 5; n++)
              params->int22[i][j][k][l][m][n] = RESCALE_dG(int22_37[i][j][k][l][m][n],
                                                           int22_dH[i][j][k][l][m][n],
                                                           tempf);
        }

  strncpy(params->Tetraloops, Tetraloops, 281);
  strncpy(params->Triloops, Triloops, 241);
  strncpy(params->Hexaloops, Hexaloops, 361);

  /* Salt correction for stack and multiloop */
  params->SaltStack = (salt==saltStandard) ? 0 : vrna_salt_stack(salt, saltT, md->helical_rise);
  if (salt == saltStandard)
    params->SaltMLbase = params->SaltMLclosing = 0;
  else
    vrna_salt_ml(params->SaltLoopDbl, md->saltMLLower, md->saltMLUpper, &params->SaltMLbase, &params->SaltMLclosing);
  
  params->MLclosing += params->SaltMLbase;
  params->MLclosing += params->SaltMLclosing;
  params->MLbase += params->SaltMLbase;
  for (i = 0; i <= NBPAIRS; i++)
    params->MLintern[i] += params->SaltMLbase;

  params->SaltDPXInit = 0;
  if (salt != saltStandard)
  {
    if (md->saltDPXInit != VRNA_MODEL_DEFAULT_SALT_DPXINIT)
      params->SaltDPXInit = md->saltDPXInit;
    else if (md->saltDPXInit)
      params->SaltDPXInit = vrna_salt_duplex_init(md);
  }
  params->DuplexInit += params->SaltDPXInit;

  params->id = ++id;
  return params;
}

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */

#endif
