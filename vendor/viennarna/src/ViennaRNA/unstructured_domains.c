/** \file unstructured_domains.c **/

/*
 *                Unstructured domains
 *
 *                This file contains everything necessary to
 *                deal with the default implementation for unstructured
 *                domains in secondary structures. This feature enables,
 *                for instance, ligand binding to unpaired stretches of
 *                an RNA secondary structure.
 *
 *                c 2016 Ronny Lorenz
 *
 *                ViennaRNA package
 */

#include <config.h>

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

/*
 #################################
 # PRIVATE VARIABLES/STRUCTS     #
 #################################
 */

struct binding_segment {
  unsigned int  start;
  unsigned int  end;
  unsigned int  type;
};

struct ud_bt_stack {
  unsigned int    from;
  vrna_ud_motif_t *motifs;
  unsigned int    motif_cnt;
  unsigned int    motif_size;
};

struct default_outside {
  int         motif_num;
  FLT_OR_DBL  exp_energy;
};

/*
 *  Default data structure for ligand binding to unpaired stretches
 */
struct ligands_up_data_default {
  /*
   **********************************
   * pre-computed position-wise
   * motif list
   **********************************
   */
  int         n;
  int         **motif_list_ext;
  int         **motif_list_hp;
  int         **motif_list_int;
  int         **motif_list_mb;

  int         *dG;
  FLT_OR_DBL  *exp_dG;
  int         *len;

  /*
   **********************************
   * below are DP matrices to store
   * the production rule results
   **********************************
   */
  int         *energies_ext;
  int         *energies_hp;
  int         *energies_int;
  int         *energies_mb;
  FLT_OR_DBL  *exp_energies_ext;
  FLT_OR_DBL  *exp_energies_hp;
  FLT_OR_DBL  *exp_energies_int;
  FLT_OR_DBL  *exp_energies_mb;

  /*
   **********************************
   * below are lists to store the
   * outside partition function for
   * each motif starting at each
   * position
   **********************************
   */
  unsigned int            *outside_ext_count;
  struct default_outside  **outside_ext;
  unsigned int            *outside_hp_count;
  struct default_outside  **outside_hp;
  unsigned int            *outside_int_count;
  struct default_outside  **outside_int;
  unsigned int            *outside_mb_count;
  struct default_outside  **outside_mb;

  FLT_OR_DBL              (*default_cb[32])(int,
                                            int,
                                            struct ligands_up_data_default *);
  FLT_OR_DBL              *exp_e_mx[32];
};

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE void
remove_ligands_up(vrna_fold_compound_t *vc);

/* default implementations for unstructured domains feature */

/* helper functions for default implementatations of unstructured domains feature */

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */

PUBLIC void
vrna_ud_remove(vrna_fold_compound_t *vc)
{
  if (vc && vc->domains_up)
    remove_ligands_up(vc);
}

/*
 #####################################
 # BEGIN OF STATIC HELPER FUNCTIONS  #
 #####################################
 */

PRIVATE void
remove_ligands_up(vrna_fold_compound_t *vc)
{
  int i;

  /* free auxiliary data */
  if (vc->domains_up->free_data)
    vc->domains_up->free_data(vc->domains_up->data);

  /* free motif sequences */
  for (i = 0; i < vc->domains_up->motif_count; i++)
    free(vc->domains_up->motif[i]);

  /* free motif names */
  for (i = 0; i < vc->domains_up->motif_count; i++)
    free(vc->domains_up->motif_name[i]);

  free(vc->domains_up->motif);
  free(vc->domains_up->motif_name);
  free(vc->domains_up->motif_size);
  free(vc->domains_up->motif_en);
  free(vc->domains_up->motif_type);

  free(vc->domains_up->uniq_motif_size);

  free(vc->domains_up);

  vc->domains_up = NULL;
}

/*
 *  Given a secondary structure in dot-bracket notation,
 *  extract all consecutive unpaired nucleotides as individual
 *  segments.
 *  After successful execution, this function returns a list of
 *  segments and the number of list elements is stored in the
 *  variable passed as segments_num
 */

/*
 **********************************
 * Default implementation for
 * ligand binding to unpaired
 * stretches follows below
 **********************************
 */

