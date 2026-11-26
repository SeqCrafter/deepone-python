#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/eval/gquad.h"
#include "ViennaRNA/eval/internal.h"

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

#include "ViennaRNA/constraints/internal_sc.inc"

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
eval_stack(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           unsigned int         options);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC int
vrna_E_internal(unsigned int  n1,
                unsigned int  n2,
                unsigned int  type,
                unsigned int  type_2,
                int           si1,
                int           sj1,
                int           sp1,
                int           sq1,
                vrna_param_t  *P)
{
  /* compute energy of degree 2 loop (stack bulge or internal) */
  unsigned int  nl, ns, u, backbones, no_close;
  int           energy, salt_stack_correction, salt_loop_correction;

  energy = INF;

  if (P) {
    no_close              = 0;
    salt_stack_correction = P->SaltStack;
    salt_loop_correction  = 0;

    if ((P->model_details.noGUclosure) &&
        ((type_2 == 3) || (type_2 == 4) || (type == 3) || (type == 4)))
      no_close = 1;

    if (n1 > n2) {
      nl  = n1;
      ns  = n2;
    } else {
      nl  = n2;
      ns  = n1;
    }

    if (nl == 0)  /* stack */
      return P->stack[type][type_2] +
             salt_stack_correction;

    if (no_close)
      return INF;

    if (P->model_details.salt != VRNA_MODEL_DEFAULT_SALT) {
      /* salt correction for loop */
      backbones = nl + ns + 2;
      if (backbones <= MAXLOOP + 1)
        salt_loop_correction = P->SaltLoop[backbones];
      else
        salt_loop_correction = vrna_salt_loop_int(backbones,
                                                  P->model_details.salt,
                                                  P->temperature + K0,
                                                  P->model_details.backbone_length);
    }

    energy = 0;

    switch (ns) {
      case 0:
        /* bulge */
        energy = (nl <= MAXLOOP) ?
                 P->bulge[nl] :
                 (P->bulge[30] + (int)(P->lxc * log(nl / 30.)));
        if (nl == 1) {
          /* add stacking energy for 1-bulges */
          energy += P->stack[type][type_2];
        } else {
          if (type > 2)
            energy += P->TerminalAU;

          if (type_2 > 2)
            energy += P->TerminalAU;
        }

        break;

      case 1:
        if (nl == 1) {
          /* 1x1 loop */
          energy = P->int11[type][type_2][si1][sj1];
        } else if (nl == 2) {
          /* 2x1 loop */
          if (n1 == 1)
            energy = P->int21[type][type_2][si1][sq1][sj1];
          else
            energy = P->int21[type_2][type][sq1][si1][sp1];
        } else {
          /* 1xn loop */
          energy = (nl + 1 <= MAXLOOP) ?
                   (P->internal_loop[nl + 1]) :
                   (P->internal_loop[30] + (int)(P->lxc * log((nl + 1) / 30.)));
          energy  += MIN2(MAX_NINIO, (int)(nl - ns) * P->ninio[2]);
          energy  += P->mismatch1nI[type][si1][sj1] +
                     P->mismatch1nI[type_2][sq1][sp1];
        }

        break;

      case 2:
        if (nl == 2) {
          /* 2x2 loop */
          energy = P->int22[type][type_2][si1][sp1][sq1][sj1];
          break;
        } else if (nl == 3) {
          /* 2x3 loop */
          energy = P->internal_loop[5] +
                   P->ninio[2];
          energy += P->mismatch23I[type][si1][sj1] +
                    P->mismatch23I[type_2][sq1][sp1];
          break;
        }

      /* fall through */

      default:
        /* generic internal loop (no else here!)*/
        u       = nl + ns;
        energy  += (u <= MAXLOOP) ?
                   (P->internal_loop[u]) :
                   (P->internal_loop[30] + (int)(P->lxc * log((u) / 30.)));

        energy += MIN2(MAX_NINIO, (int)(nl - ns) * P->ninio[2]);

        energy += P->mismatchI[type][si1][sj1] +
                  P->mismatchI[type_2][sq1][sp1];

        break;
    }

    energy += salt_loop_correction;
  }

  return energy;
}

PUBLIC int
vrna_eval_stack(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j,
                unsigned int          options)
{
  if ((fc) &&
      (i > 0) &&
      (i + 2 < j)) {
    /* is this base pair allowed to close a hairpin (like) loop ? */
    if ((options & VRNA_EVAL_LOOP_NO_HC) ||
        (fc->hc == NULL) ||
        (fc->hc->eval_int(i, j, i + 1, j - 1, fc->hc)))
      return eval_stack(fc, i, j, options);
  }

  return INF;
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

PRIVATE int
eval_stack(vrna_fold_compound_t *fc,
           unsigned int         i,
           unsigned int         j,
           unsigned int         options)
{
  unsigned char     sliding_window;
  char              *ptype, **ptype_local;
  short             **SS;
  unsigned int      p, q, s, n_seq, *sn, type, type_2;
  int               e, ij, pq, *rtype, *indx;
  vrna_param_t      *P;
  vrna_md_t         *md;
  struct sc_int_dat sc_wrapper;

  e               = INF;
  sliding_window  = (fc->hc->type == VRNA_HC_WINDOW) ? 1 : 0;
  p               = i + 1;
  q               = j - 1;
  n_seq           = (fc->type == VRNA_FC_TYPE_SINGLE) ? 1 : fc->n_seq;
  sn              = fc->strand_number;
  SS              = (fc->type == VRNA_FC_TYPE_SINGLE) ? NULL : fc->S;
  ptype           = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                    (sliding_window ? NULL : fc->ptype) :
                    NULL;
  ptype_local = (fc->type == VRNA_FC_TYPE_SINGLE) ?
                (sliding_window ? fc->ptype_local : NULL) :
                NULL;
  P     = fc->params;
  md    = &(P->model_details);
  rtype = &(md->rtype[0]);
  indx  = (sliding_window) ? NULL : fc->jindx;
  ij    = (sliding_window) ? 0 : indx[j] + i;
  pq    = (sliding_window) ? 0 : indx[q] + p;

  if ((sn[p] != sn[i]) ||
      (sn[j] != sn[q]))
    return e;

  switch (fc->type) {
    case VRNA_FC_TYPE_SINGLE:
      type = sliding_window ?
             vrna_get_ptype_window(i, j, ptype_local) :
             vrna_get_ptype(ij, ptype);
      type_2 = sliding_window ?
               rtype[vrna_get_ptype_window(p, q, ptype_local)] :
               rtype[vrna_get_ptype(pq, ptype)];

      e = P->stack[type][type_2];

      break;

    }

  if (e != INF) {
    if (!(options & VRNA_EVAL_LOOP_NO_SC)) {
      init_sc_int(fc, &sc_wrapper);

      /* add soft constraints */
      if (sc_wrapper.pair)
        e += sc_wrapper.pair(i, j, p, q, &sc_wrapper);

      free_sc_int(&sc_wrapper);
    }
  }

  return e;
}

/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#endif
