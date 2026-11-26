/** \file dp_matrices.c **/

/*
 *                Dynamic Programming matrix related functions
 *
 *                This file contains everything necessary to
 *                obtain and destroy data structures representing
 *                dynamic programming (DP) matrices used in the folding
 *                recurrences throughout the VienneRNA paclage
 *
 *                c Ronny Lorenz
 *
 *                ViennaRNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/mfe/gquad.h"

/*
 #################################
 # PRIVATE MACROS                #
 #################################
 */

/* the definitions below indicate which arrays should be allocated upon retrieval of a matrices data structure */
#define ALLOC_NOTHING     0
#define ALLOC_F           1
#define ALLOC_F5          2
#define ALLOC_F3          4
#define ALLOC_FC          8
#define ALLOC_C           16
#define ALLOC_FML         32
#define ALLOC_PROBS       256
#define ALLOC_AUX         512

#define ALLOC_CIRC        1024
#define ALLOC_MULTISTRAND 2048
#define ALLOC_UNIQ        4096

#define ALLOC_MFE_DEFAULT         (ALLOC_F5 | ALLOC_C | ALLOC_FML)
#define ALLOC_MFE_LOCAL           (ALLOC_F3 | ALLOC_C | ALLOC_FML)

#define ALLOC_PF_WO_PROBS         (ALLOC_F | ALLOC_C | ALLOC_FML)
#define ALLOC_PF_DEFAULT          (ALLOC_PF_WO_PROBS | ALLOC_PROBS | ALLOC_AUX)

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
PRIVATE unsigned int
get_mx_alloc_vector(vrna_fold_compound_t  *fc,
                    vrna_mx_type_e        type,
                    unsigned int          options);

PRIVATE unsigned int
get_mx_mfe_alloc_vector_current(vrna_mx_mfe_t   *mx,
                                vrna_mx_type_e  mx_type);

PRIVATE void
mfe_matrices_free_default(vrna_mx_mfe_t *self);

PRIVATE int
add_mfe_matrices(vrna_fold_compound_t *vc,
                 vrna_mx_type_e       type,
                 unsigned int         alloc_vector);

PRIVATE vrna_mx_mfe_t *
init_mx_mfe_default(vrna_fold_compound_t  *fc,
                    unsigned int          alloc_vector);

PRIVATE INLINE void
nullify_mfe(vrna_mx_mfe_t *mx);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC void
vrna_mx_mfe_free(vrna_fold_compound_t *vc)
{
  if (vc) {
    vrna_mx_mfe_t *self = vc->matrices;
    if (self) {
      switch (self->type) {
        case VRNA_MX_DEFAULT:
          mfe_matrices_free_default(self);
          break;

        default:                /* do nothing */
          break;
      }
      free(self);
      vc->matrices = NULL;
    }
  }
}

PUBLIC int
vrna_mx_add(vrna_fold_compound_t  *vc,
            vrna_mx_type_e        mx_type,
            unsigned int          options)
{
  int ret;

  ret = 1;

  if (options & VRNA_OPTION_MFE)
    ret &= vrna_mx_mfe_add(vc, mx_type, options);

  ;;

  return ret;
}

PUBLIC int
vrna_mx_mfe_add(vrna_fold_compound_t  *vc,
                vrna_mx_type_e        mx_type,
                unsigned int          options)
{
  unsigned int mx_alloc_vector;

  if (vc->params) {
    options |= VRNA_OPTION_MFE;

    mx_alloc_vector = get_mx_alloc_vector(vc,
                                          mx_type,
                                          options);
    vrna_mx_mfe_free(vc);
    return add_mfe_matrices(vc, mx_type, mx_alloc_vector);
  }

  return 0;
}

PUBLIC int
vrna_mx_prepare(vrna_fold_compound_t  *vc,
                unsigned int          options)
{
  int             ret, realloc;
  unsigned int    mx_alloc_vector, mx_alloc_vector_current;
  vrna_mx_type_e  mx_type;

  ret = 1;

  if (vc) {
    /*  check whether we have the correct DP matrices attached, and if there is
     *  enough memory allocated
     */
    if (options & VRNA_OPTION_MFE) {
      /* prepare for MFE computation */
      mx_type = VRNA_MX_DEFAULT;;

      if (vc->strands > 1)
        options |= VRNA_OPTION_HYBRID;

      realloc = 0;

      if (!vc->matrices || (vc->matrices->type != mx_type) || (vc->matrices->length < vc->length)) {
        realloc = 1;
      } else {
        mx_alloc_vector =
          get_mx_alloc_vector(vc, mx_type, options);
        mx_alloc_vector_current = get_mx_mfe_alloc_vector_current(vc->matrices, mx_type);
        if ((mx_alloc_vector & mx_alloc_vector_current) != mx_alloc_vector)
          realloc = 1;
      }

      if (realloc) /* Add DP matrices, if not they are not present */
        ret &= vrna_mx_mfe_add(vc, mx_type, options);
    }

    ;
  } else {
    ret = 0;
  }

  return ret;
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */
PRIVATE unsigned int
get_mx_mfe_alloc_vector_current(vrna_mx_mfe_t   *mx,
                                vrna_mx_type_e  mx_type)
{
  unsigned int mx_alloc_vector = ALLOC_NOTHING;

  if (mx) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        if (mx->f5)
          mx_alloc_vector |= ALLOC_F5;

        if (mx->f3)
          mx_alloc_vector |= ALLOC_F3;

        if ((mx->fms5) ||
            (mx->fms3))
          mx_alloc_vector |= ALLOC_MULTISTRAND;

        if (mx->c)
          mx_alloc_vector |= ALLOC_C;

        if (mx->fML)
          mx_alloc_vector |= ALLOC_FML;

        if ((mx->fM1) ||
            (mx->fM2_real))
          mx_alloc_vector |= ALLOC_UNIQ;

        if ((mx->fM1_new) &&
            (mx->fM2_real))
          mx_alloc_vector |= ALLOC_CIRC;

        break;

      default:
        break;
    }
  }

  return mx_alloc_vector;
}

PRIVATE int
add_mfe_matrices(vrna_fold_compound_t *vc,
                 vrna_mx_type_e       mx_type,
                 unsigned int         alloc_vector)
{
  if (vc) {
    switch (mx_type) {
      case VRNA_MX_DEFAULT:
        vc->matrices = init_mx_mfe_default(vc, alloc_vector);
        break;

      default:
        return 0;
    }

    if (!vc->matrices)
      return 0;

    ;
  }

  return 1;
}

PRIVATE unsigned int
get_mx_alloc_vector(vrna_fold_compound_t  *fc,
                    vrna_mx_type_e        mx_type,
                    unsigned int          options)
{
  unsigned int  v;
  vrna_md_t     *md_p;

  md_p = &(fc->params->model_details);

  v = ALLOC_NOTHING;

  /* default MFE matrices ? */
  if (options & VRNA_OPTION_MFE)
    v |= (mx_type == VRNA_MX_WINDOW) ? ALLOC_MFE_LOCAL : ALLOC_MFE_DEFAULT;

  /* default PF matrices ? */
  ;;

  if ((fc->strands > 1) || (options & VRNA_OPTION_HYBRID))
    v |= ALLOC_MULTISTRAND;

  /* matrices for circular folding ? */
  ;;

  /* unique ML decomposition ? */
  if (md_p->uniq_ML)
    v |= ALLOC_UNIQ;

  return v;
}

PRIVATE void
mfe_matrices_free_default(vrna_mx_mfe_t *self)
{
  free(self->f5);
  free(self->f3);

  if (self->fms5)
    for (unsigned int s = 0; s < self->strands; s++)
      free(self->fms5[s]);

  free(self->fms5);

  if (self->fms3)
    for (unsigned int s = 0; s < self->strands; s++)
      free(self->fms3[s]);

  free(self->fms3);

  free(self->c);
  free(self->fML);
  free(self->fM1);
  free(self->fM2);
  free(self->fM1_new);
  free(self->fM2_real);
#ifndef VRNA_DISABLE_C11_FEATURES
#else
#endif
}

PRIVATE vrna_mx_mfe_t *
init_mx_mfe_default(vrna_fold_compound_t  *fc,
                    unsigned int          alloc_vector)
{
  unsigned int  n, size, lin_size, s, strands;
  vrna_mx_mfe_t *mx;
  vrna_mx_mfe_t init = {
    .type = VRNA_MX_DEFAULT
  };

  n = fc->length;

  if ((int)(n * n) >= (int)INT_MAX) {
    vrna_log_warning("init_mx_mfe_default(): "
                         "sequence length %d exceeds addressable range",
                         n);
    return NULL;
  }

  mx = vrna_alloc(sizeof(vrna_mx_mfe_t));

  if (mx) {
    memcpy(mx, &init, sizeof(vrna_mx_mfe_t));
    nullify_mfe(mx);

    strands     = fc->strands;
    size        = ((n + 1) * (n + 2)) / 2;
    lin_size    = n + 2;
    mx->length  = n;
    mx->strands = strands;

    if (alloc_vector & ALLOC_F5)
      mx->f5 = (int *)vrna_alloc(sizeof(int) * lin_size);

    if (alloc_vector & ALLOC_F3)
      mx->f3 = (int *)vrna_alloc(sizeof(int) * lin_size);

    if (alloc_vector & ALLOC_MULTISTRAND) {
      mx->fms5  = (int **)vrna_alloc(sizeof(int *) * strands);
      mx->fms3  = (int **)vrna_alloc(sizeof(int *) * strands);

      for (s = 0; s < strands; s++) {
        mx->fms5[s] = (int *)vrna_alloc(sizeof(int) * (n + 1));
        mx->fms3[s] = (int *)vrna_alloc(sizeof(int) * (n + 1));
      }
    }

    if (alloc_vector & ALLOC_C)
      mx->c = (int *)vrna_alloc(sizeof(int) * size);

    if (alloc_vector & ALLOC_FML)
      mx->fML = (int *)vrna_alloc(sizeof(int) * size);

    if ((alloc_vector & ALLOC_UNIQ) ||
        (alloc_vector & ALLOC_CIRC)) {
      if (!(alloc_vector & ALLOC_CIRC))
        mx->fM1 = (int *)vrna_alloc(sizeof(int) * size);

      mx->fM2_real = (int *)vrna_alloc(sizeof(int) * size);
    }

    if (alloc_vector & ALLOC_CIRC)
      mx->fM1_new  = (int *)vrna_alloc(sizeof(int) * lin_size);
  }

  return mx;
}

PRIVATE INLINE void
nullify_mfe(vrna_mx_mfe_t *mx)
{
  if (mx) {
    mx->length  = 0;
    mx->strands = 0;

    switch (mx->type) {
      case VRNA_MX_DEFAULT:
        mx->c     = NULL;
        mx->f5    = NULL;
        mx->f3    = NULL;
        mx->fms5  = NULL;
        mx->fms3  = NULL;
        mx->fML   = NULL;
        mx->fM1   = NULL;
        mx->fM2   = NULL;
        mx->fM1_new   = NULL;
        mx->fM2_real  = NULL;
        mx->Fc    = INF;
        mx->FcH   = INF;
        mx->FcI   = INF;
        mx->FcM   = INF;
        mx->c_gq  = NULL;
        break;

      }
  }
}

