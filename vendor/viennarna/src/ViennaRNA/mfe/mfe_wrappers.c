/*
 * Various wrappers to create simplified interfaces for Minimum Free Energy
 * prediction
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/mfe/global.h"

#ifndef INLINE
# ifdef __GNUC__
#   define INLINE inline
# else
#   define INLINE
# endif
#endif

/* wrappers for single sequences */
PUBLIC float
vrna_fold(const char  *string,
          char        *structure)
{
  float                 mfe;
  vrna_fold_compound_t  *vc;
  vrna_md_t             md;

  vrna_md_set_default(&md);
  vc  = vrna_fold_compound(string, &md, 0);
  mfe = vrna_mfe(vc, NULL);

  vrna_fold_compound_free(vc);

  return mfe;
}

/* wrappers for multiple sequence alignments */

/* wrappers for RNA-RNA cofolding interaction */

