#ifndef VIENNA_RNA_PACKAGE_LOOPS_HAIRPIN_H
#define VIENNA_RNA_PACKAGE_LOOPS_HAIRPIN_H

#include <math.h>
#include <string.h>
#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/params/salt.h>
#include <ViennaRNA/eval/basic.h>

#ifdef VRNA_WARN_DEPRECATED
# if defined(DEPRECATED)
#   undef DEPRECATED
# endif
# if defined(__clang__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated("", msg)))
# elif defined(__GNUC__)
#  define DEPRECATED(func, msg) func __attribute__ ((deprecated(msg)))
# else
#  define DEPRECATED(func, msg) func
# endif
#else
# define DEPRECATED(func, msg) func
#endif

#ifdef __GNUC__
# define INLINE inline
#else
# define INLINE
#endif

/**
 *
 *  @file     ViennaRNA/loops/hairpin.h
 *  @ingroup  eval, eval_loops, eval_loops_hp
 *  @brief    Energy evaluation of hairpin loops for MFE and partition function calculations
 */

/**
 *
 *  @addtogroup   eval_loops_hp
 *  @{
 */

/**
 *  @name Basic free energy interface
 *  @{
 */

/**
 *  @brief Retrieve the energy of a hairpin-loop
 *
 *  To evaluate the free energy of a hairpin-loop, several parameters have to be known.
 *  A general hairpin-loop has this structure:<BR>
 *  <PRE>
 *        a3 a4
 *      a2     a5
 *      a1     a6
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as there are
 *  six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
 *  a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
 *
 *  @warning  This function @b only evaluates the free energy of a hairpin loop according to the current
 *            Turner energy parameter set! No additional hard- or soft constraints are applied. See vrna_eval_hp_loop()
 *            for a function that also takes into account any user-supplied constraints!
 *
 *  @note     Whenever one of the mismatch base encodings @p si1 or @p sj1 is negative, terminal mismatch energies are not applied!
 *
 *  @note     The parameter @p sequence is a 0-terminated string of size @p size + 2 that contain the nucleic acid
 *            sequence of the loop in upper-case letters. This parameter is only required for loops of @p size below
 *            7, since it is used for look-up of unusually stable tri-, tetra- and hexa-loops, such as GNRA tetra loops.
 *            Those may have additional sequence-dependent tabulated free energies available.
 *
 *  @see      #vrna_param_t, vrna_eval_hp_loop()
 *
 *  @param  size      The size of the loop (number of unpaired nucleotides)
 *  @param  type      The pair type of the base pair closing the hairpin
 *  @param  si1       The 5'-mismatching nucleotide
 *  @param  sj1       The 3'-mismatching nucleotide
 *  @param  sequence  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P         The datastructure containing scaled energy parameters
 *  @return           The Free energy of the Hairpin-loop in dcal/mol
 */
int
vrna_E_hairpin(unsigned int size,
               unsigned int type,
               int          si1,
               int          sj1,
               const char   *sequence,
               vrna_param_t *P);

/**
 *  @brief Evaluate free energy of a hairpin loop
 *
 *  This function evaluates the free energy of a hairpin loop closed by a base pair
 *  (i,j). By default (@p options = #VRNA_EVAL_LOOP_DEFAULT), *all* user-supplied
 *  constraints will be taken into consideration. This means that any hard constraints
 *  that prohibit the formation of this loop will result in an energy contribution
 *  of @b #INF. On the other hand, if, given the set of constraints, the loop is
 *  allowed then its free energy is evaluated according to the Nearest Neighbor
 *  energy parameter set. On top of that, any user-supplied soft-constraints will
 *  be added, if applicable.
 *
 *  The @p options argument allows for (de-)activating certain aspects of the evaluation,
 *  e.g. hard constraints, soft constraints, etc.
 *
 *  @note   If sequence position @p i is larger than @p j, the function assumes a
 *          hairpin loop formed by a circular RNA, where the unpaired loop sequence
 *          spans the n,1-junction.
 *
 *  @note   By default, all user-supplied hard- and soft constraints will be taken
 *          into account! Use the #VRNA_EVAL_LOOP_NO_HC and #VRNA_EVAL_LOOP_NO_SC
 *          bit flags as input for @p options to change the default behavior if necessary.
 *
 *  @note   This function is polymorphic! The provided #vrna_fold_compound_t may be
 *          of type #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @see    vrna_E_hairpin(), vrna_exp_eval_hairpin(),
 *          #VRNA_EVAL_LOOP_NO_HC, #VRNA_EVAL_LOOP_NO_SC, #VRNA_EVAL_LOOP_NO_CONSTRAINTS
 *
 *  @param  fc        The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i         5'-position of the base pair
 *  @param  j         3'-position of the base pair
 *  @param  options   A bit-field that specifies which aspects (not) to consider during evaluation
 *  @returns          Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol or #INF if the loop is forbidden
 */
int
vrna_eval_hairpin(vrna_fold_compound_t  *fc,
                  unsigned int          i,
                  unsigned int          j,
                  unsigned int          options);

/* End basic interface */
/** @} */

/**
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */

/**
 *  @brief Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of a hairpin loop
 *
 *  This is the partition function variant of vrna_E_hp() that returns the Boltzmann
 *  weight @f$e^{-\Delta E/kT} @f$ instead of the energy @f$ E @f$.
 *
 *  @note   Whenever one of the mismatch base encodings @p si1 or @p sj1 is negative, terminal mismatch energies are not applied!
 *
 *  @note   Do not forget to scale this Bolztmann factor properly, e.g. by multiplying with scale[u+2]
 *
 *  @see vrna_exp_eval_hp_loop(), #vrna_exp_param_t, vrna_E_hp()
 *
 *  @param  size      The size of the loop (number of unpaired nucleotides)
 *  @param  type      The pair type of the base pair closing the hairpin
 *  @param  si1       The 5'-mismatching nucleotide
 *  @param  sj1       The 3'-mismatching nucleotide
 *  @param  sequence  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P         The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Hairpin-loop
 */

/**
 *  @brief High-Level function for hairpin loop energy evaluation (partition function variant)
 *
 *  This is the partition function variant of vrna_eval_hp_loop() that returns the Boltzmann
 *  weight @f$e^{-\Delta E/kT} @f$ instead of the energy @f$ E @f$. On top of all constraints
 *  application, this function already scales the Boltzmann factor, i.e. it multiplies the
 *  result with scale[u + 2]
 *
 *  The @p options argument allows for (de-)activating certain aspects of the evaluation,
 *  e.g. hard constraints, soft constraints, etc.
 *
 *  @note   If sequence position @p i is larger than @p j, the function assumes a
 *          hairpin loop formed by a circular RNA, where the unpaired loop sequence
 *          spans the n,1-junction.
 *
 *  @note   By default, all user-supplied hard- and soft constraints will be taken
 *          into account! Use the #VRNA_EVAL_LOOP_NO_HC and #VRNA_EVAL_LOOP_NO_SC
 *          bit flags to change the default behavior if necessary.
 *
 *  @note   This function is polymorphic! The provided #vrna_fold_compound_t may be
 *          of type #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @see    vrna_eval_hairpin(), vrna_exp_E_hairpin(),
 *          #VRNA_EVAL_LOOP_NO_HC, #VRNA_EVAL_LOOP_NO_SC, #VRNA_EVAL_LOOP_NO_CONSTRAINTS
 *
 *  @param  fc        The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i         5'-position of the base pair
 *  @param  j         3'-position of the base pair
 *  @param  options   A bit-field that specifies which aspects (not) to consider during evaluation
 *  @returns          Boltzmann factor of the free energy of the hairpin loop closed by @f$ (i,j) @f$ or 0. if the loop is forbidden
 */

/* End partition function interface */
/** @} */

/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup   eval_deprecated
 *  @{
 */

/**
 *  @brief Compute the Energy of a hairpin-loop
 *
 *  To evaluate the free energy of a hairpin-loop, several parameters have to be known.
 *  A general hairpin-loop has this structure:<BR>
 *  <PRE>
 *        a3 a4
 *      a2     a5
 *      a1     a6
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as there are
 *  six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
 *  a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
 *
 *  @note The parameter sequence should contain the sequence of the loop in capital letters of the nucleic acid
 *        alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
 *        which are treated differently (based on experimental data) if they are tabulated.
 *
 *  @see scale_parameters(), vrna_param_t
 *
 *  @warning  Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *            Energy evaluation may change due to updates in global variable "tetra_loop"
 *
 *  @param  size  The size of the loop (number of unpaired nucleotides)
 *  @param  type  The pair type of the base pair closing the hairpin
 *  @param  si1   The 5'-mismatching nucleotide
 *  @param  sj1   The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P     The datastructure containing scaled energy parameters
 *  @return The Free energy of the Hairpin-loop in dcal/mol
 */
DEPRECATED(PRIVATE INLINE int
           E_Hairpin(int          size,
                     int          type,
                     int          si1,
                     int          sj1,
                     const char   *string,
                     vrna_param_t *P),
           "Use vrna_E_hairpin() instead!");

DEPRECATED(int
           vrna_E_hp_loop(vrna_fold_compound_t  *fc,
                          int                   i,
                          int                   j),
           "Use vrna_eval_hairpin() instead!");

DEPRECATED(int
           vrna_E_ext_hp_loop(vrna_fold_compound_t  *fc,
                              int                   i,
                              int                   j),
           "Use vrna_eval_hairpin() instead!");

DEPRECATED(int
           vrna_eval_hp_loop(vrna_fold_compound_t *fc,
                             int                  i,
                             int                  j),
           "Use vrna_eval_hairpin() instead!");

DEPRECATED(int
           vrna_eval_ext_hp_loop(vrna_fold_compound_t *fc,
                                 int                  i,
                                 int                  j),
           "Use vrna_eval_hairpin() instead!");

/**
 *  @brief Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of a hairpin loop
 *
 *  @note multiply by scale[u+2]
 *
 *  @see get_scaled_pf_parameters(), vrna_exp_param_t, E_Hairpin()
 *
 *  @warning  Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *            Energy evaluation may change due to updates in global variable "tetra_loop"
 *
 *  @param  u       The size of the loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the hairpin
 *  @param  si1     The 5'-mismatching nucleotide
 *  @param  sj1     The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Hairpin-loop
 */
DEPRECATED(PRIVATE INLINE FLT_OR_DBL
           exp_E_Hairpin(int              u,
                         int              type,
                         short            si1,
                         short            sj1,
                         const char       *string,
                         vrna_exp_param_t *P),
           "Use vrna_exp_E_hairpin() instead!");

DEPRECATED(FLT_OR_DBL
           vrna_exp_E_hp_loop(vrna_fold_compound_t  *fc,
                              int                   i,
                              int                   j),
           "Use vrna_exp_eval_hairpin() instead!");

/** @} */

#endif

#endif
