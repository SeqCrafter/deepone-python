/* constraints handling */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <config.h>

#include "ViennaRNA/utils/basic.h"

#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

#define STATE_CLEAN         (unsigned char)0
#define STATE_DIRTY_UP_MFE  (unsigned char)1
#define STATE_DIRTY_UP_PF   (unsigned char)2
#define STATE_DIRTY_BP_MFE  (unsigned char)4
#define STATE_DIRTY_BP_PF   (unsigned char)8

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
prepare_sc_up_mfe(vrna_sc_t     *sc,
                  unsigned int  n,
                  unsigned int  options);

PRIVATE void
prepare_sc_bp_mfe(vrna_sc_t     *sc,
                  unsigned int  n,
                  int           *idx,
                  unsigned int  options);

PRIVATE int
prepare_sc_user_cb(vrna_fold_compound_t *fc,
                   vrna_sc_t            *sc,
                   unsigned int         options);

PRIVATE INLINE void
populate_sc_up_mfe(vrna_sc_t    *sc,
                   unsigned int i,
                   unsigned int n);

PRIVATE INLINE int
get_stored_bp_contributions(vrna_sc_bp_storage_t  *container,
                            unsigned int          j);

PRIVATE INLINE void
populate_sc_bp_mfe(vrna_sc_t    *sc,
                   unsigned int i,
                   unsigned int maxdist,
                   unsigned int n,
                   int          *idx);

PRIVATE INLINE void
free_sc_up(vrna_sc_t *sc);

PRIVATE INLINE void
free_sc_bp(vrna_sc_t *sc);

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC int
vrna_sc_prepare(vrna_fold_compound_t  *fc,
                unsigned int          options)
{
  int ret = 0;

  if (fc) {
    unsigned int  n, s;

    n = fc->length;
    switch (fc->type) {
      case VRNA_FC_TYPE_SINGLE:
        if (options & VRNA_OPTION_MFE) {
          prepare_sc_up_mfe(fc->sc, n, options);
          prepare_sc_bp_mfe(fc->sc, n, fc->jindx, options);
        }

        ;

        ret &= prepare_sc_user_cb(fc, fc->sc, options);

        break;

      }

  }

  return ret;
}

PUBLIC void
vrna_sc_free(vrna_sc_t *sc)
{
  if (sc) {
    free_sc_up(sc);
    free_sc_bp(sc);

    free(sc->energy_stack);
    free(sc->exp_energy_stack);

    if (sc->free_data)
      sc->free_data(sc->data);

    free(sc);
  }
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

/* pupulate sc->energy_up array at position i from sc->up_storage data */
PRIVATE INLINE void
populate_sc_up_mfe(vrna_sc_t    *sc,
                   unsigned int i,
                   unsigned int n)
{
  unsigned int  j;

  sc->energy_up[i][0] = 0;
  for (j = 1; j <= n; j++)
    sc->energy_up[i][j] = sc->energy_up[i][j - 1] +
                          sc->up_storage[i + j - 1];
}

PRIVATE INLINE int
get_stored_bp_contributions(vrna_sc_bp_storage_t  *container,
                            unsigned int          j)
{
  unsigned int  cnt;
  int           e;

  e = 0;

  /* go through list of constraints for current position i */
  for (cnt = 0; container[cnt].interval_start != 0; cnt++) {
    if (container[cnt].interval_start > j)
      break; /* only constraints for pairs (i,q) with q > j left */

    if (container[cnt].interval_end < j)
      continue; /* constraint for pairs (i,q) with q < j */

    /* constraint has interval [p,q] with p <= j <= q */
    e += container[cnt].e;
  }

  return e;
}

PRIVATE INLINE void
populate_sc_bp_mfe(vrna_sc_t    *sc,
                   unsigned int i,
                   unsigned int maxdist,
                   unsigned int n,
                   int          *idx)
{
  unsigned int  j, k;
  int           e;

  if (sc->bp_storage[i]) {
    for (k = 1; k < maxdist; k++) {
      j = i + k;

      if (j > n)
        break;

      e = get_stored_bp_contributions(sc->bp_storage[i], j);

      switch (sc->type) {
        case VRNA_SC_DEFAULT:
          sc->energy_bp[idx[j] + i] = e;
          break;

        case VRNA_SC_WINDOW:
          sc->energy_bp_local[i][j - i] = e;
          break;
      }
    }
  } else {
    for (k = 1; k < maxdist; k++) {
      j = i + k;
      if (j > n)
        break;

      switch (sc->type) {
        case VRNA_SC_DEFAULT:
          sc->energy_bp[idx[j] + i] = 0;
          break;

        case VRNA_SC_WINDOW:
          sc->energy_bp_local[i][j - i] = 0;
          break;
      }
    }
  }
}

PRIVATE INLINE void
free_sc_up(vrna_sc_t *sc)
{
  unsigned int i;

  if (sc) {
    free(sc->up_storage);

    sc->up_storage = NULL;

    if (sc->type == VRNA_SC_DEFAULT) {
      if (sc->energy_up)
        for (i = 0; i <= sc->n + 1; i++)
          free(sc->energy_up[i]);

      if (sc->exp_energy_up)
        for (i = 0; i <= sc->n + 1; i++)
          free(sc->exp_energy_up[i]);
    }

    free(sc->energy_up);
    sc->energy_up = NULL;

    free(sc->exp_energy_up);
    sc->exp_energy_up = NULL;

    sc->state &= ~(STATE_DIRTY_UP_MFE | STATE_DIRTY_UP_PF);
  }
}

PRIVATE INLINE void
free_sc_bp(vrna_sc_t *sc)
{
  unsigned int i;

  if (sc) {
    if (sc->bp_storage) {
      for (i = 1; i <= sc->n; i++)
        free(sc->bp_storage[i]);
      free(sc->bp_storage);
      sc->bp_storage = NULL;
    }

    switch (sc->type) {
      case VRNA_SC_DEFAULT:
        free(sc->energy_bp);
        sc->energy_bp = NULL;

        free(sc->exp_energy_bp);
        sc->energy_bp = NULL;

        break;

      case VRNA_SC_WINDOW:
        free(sc->energy_bp_local);
        sc->energy_bp_local = NULL;

        free(sc->exp_energy_bp_local);
        sc->exp_energy_bp_local = NULL;

        break;
    }
    sc->state &= ~(STATE_DIRTY_BP_MFE | STATE_DIRTY_BP_PF);
  }
}

/* populate sc->energy_up arrays for usage in MFE computations */
PRIVATE void
prepare_sc_up_mfe(vrna_sc_t     *sc,
                  unsigned int  n,
                  unsigned int  options)
{
  unsigned int i;

  /* prepare sc for unpaired nucleotides only if we actually have some to apply */
  if (sc) {
    if (sc->up_storage) {
      if (sc->state & STATE_DIRTY_UP_MFE) {
        /*  allocate memory such that we can access the soft constraint
         *  energies of a subsequence of length j starting at position i
         *  via sc->energy_up[i][j]
         */
        sc->energy_up = (int **)vrna_realloc(sc->energy_up, sizeof(int *) * (n + 2));

        {
          for (i = 1; i <= n; i++)
            sc->energy_up[i] = (int *)vrna_realloc(sc->energy_up[i], sizeof(int) * (n - i + 2));

          sc->energy_up[0]      = (int *)vrna_realloc(sc->energy_up[0], sizeof(int));
          sc->energy_up[n + 1]  = (int *)vrna_realloc(sc->energy_up[n + 1], sizeof(int));

          /* now add soft constraints as stored in container for unpaired sc */
          for (i = 1; i <= n; i++)
            populate_sc_up_mfe(sc, i, n - i + 1);

          sc->energy_up[0][0]     = 0;
          sc->energy_up[n + 1][0] = 0;
        }

        sc->state &= ~STATE_DIRTY_UP_MFE;
      }
    } else if (sc->energy_up) {
      /* remove any unpaired sc if storage container is empty */
      free_sc_up(sc);
    }
  }
}

/* populate sc->exp_energy_up arrays for usage in partition function computations */

/* populate sc->energy_bp arrays for usage in MFE computations */
PRIVATE void
prepare_sc_bp_mfe(vrna_sc_t     *sc,
                  unsigned int  n,
                  int           *idx,
                  unsigned int  options)
{
  unsigned int  i;

  if (sc) {
    /* prepare sc for base paired positions only if we actually have some to apply */
    if (sc->bp_storage) {
      if (sc->state & STATE_DIRTY_BP_MFE) {
        {
          sc->energy_bp =
            (int *)vrna_realloc(sc->energy_bp, sizeof(int) * (((n + 1) * (n + 2)) / 2));

          for (i = 1; i < n; i++)
            populate_sc_bp_mfe(sc, i, n, n, idx);
        }

        sc->state &= ~STATE_DIRTY_BP_MFE;
      }
    } else {
      free_sc_bp(sc);
    }
  }
}

/* populate sc->exp_energy_bp arrays for usage in partition function computations */

PRIVATE int
prepare_sc_user_cb(vrna_fold_compound_t *fc,
                   vrna_sc_t            *sc,
                   unsigned int         options)
{
  int ret = 0;

  if ((sc) &&
      (sc->data) &&
      (sc->prepare_data) &&
      (sc->type == VRNA_SC_DEFAULT))
    ret = sc->prepare_data(fc, sc->data, options, NULL);

  return ret;
}

