#ifndef VIENNA_RNA_PACKAGE_EVAL_INTERNAL_H
#define VIENNA_RNA_PACKAGE_EVAL_INTERNAL_H

#include <math.h>

#include <ViennaRNA/utils/basic.h>
#include <ViennaRNA/params/default.h>
#include <ViennaRNA/datastructures/basic.h>
#include <ViennaRNA/fold_compound.h>
#include <ViennaRNA/params/basic.h>
#include <ViennaRNA/constraints/hard.h>
#include <ViennaRNA/constraints/soft.h>
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
 *  @file     ViennaRNA/eval/internal.h
 *  @ingroup  eval, eval_loops, eval_loops_int
 *  @brief    Energy evaluation of internal loops
 */

/**
 *  @addtogroup   eval_loops_int
 *  @{
 */

/**
 *  @name Basic free energy interface
 *  @{
 */

/**
 *  @brief  Compute the Energy of an internal loop
 *
 *  This function computes the free energy @f$ E @f$ of an internal-loop with the
 *  following structure: <BR>
 *  <PRE>
 *        3'  5'
 *        |   |
 *        U - V
 *    a_n       b_1
 *     .        .
 *     .        .
 *     .        .
 *    a_1       b_m
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  This general structure depicts an internal-loop that is closed by the base pair (X,Y).
 *  The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
 *  that constitute the loop. In this example, the length of the internal-loop is @f$(n+m)@f$
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:<BR>
 *  5'-mismatch: a_1<BR>
 *  3'-mismatch: b_m<BR>
 *  and for the enclosed base pair (V,U):<BR>
 *  5'-mismatch: b_1<BR>
 *  3'-mismatch: a_n<BR>
 *
 *  @note Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair
 *        must be 'turned arround' when evaluating the free energy of the internal-loop<br>
 *        This function is threadsafe
 *
 *  @see vrna_exp_E_internal()
 *
 *  @param  n1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  n2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the internal loop
 *  @param  type_2  The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled energy parameters
 *  @return The Free energy of the internal loop in dcal/mol
 */
int
vrna_E_internal(unsigned int  n1,
                unsigned int  n2,
                unsigned int  type,
                unsigned int  type_2,
                int           si1,
                int           sj1,
                int           sp1,
                int           sq1,
                vrna_param_t  *P);

/**
 *  @brief Evaluate the free energy contribution of an internal loop with delimiting
 *  base pairs @f$(i,j)@f$ and @f$(k,l)@f$
 *
 *  @note This function is polymorphic, i.e. it accepts #vrna_fold_compound_t of type
 *        #VRNA_FC_TYPE_SINGLE as well as #VRNA_FC_TYPE_COMPARATIVE
 */

int
vrna_eval_stack(vrna_fold_compound_t  *fc,
                unsigned int          i,
                unsigned int          j,
                unsigned int          options);

/* End basic interface */
/**@}*/

/**
 *  @name Boltzmann weight (partition function) interface
 *  @{
 */

/* End partition function interface */
/**@}*/

/**
 * @}
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/**
 *  @addtogroup   eval_deprecated
 *  @{
 */

#ifdef ON_SAME_STRAND
#undef ON_SAME_STRAND
#endif

#define ON_SAME_STRAND(I, J, C)  (((I) >= (C)) || ((J) < (C)))

/**
 *  <H2>Compute the Energy of an internal-loop</H2>
 *  This function computes the free energy @f$\Delta G@f$ of an internal-loop with the
 *  following structure: <BR>
 *  <PRE>
 *        3'  5'
 *        |   |
 *        U - V
 *    a_n       b_1
 *     .        .
 *     .        .
 *     .        .
 *    a_1       b_m
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  This general structure depicts an internal-loop that is closed by the base pair (X,Y).
 *  The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
 *  that constitute the loop. In this example, the length of the internal-loop is @f$(n+m)@f$
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:<BR>
 *  5'-mismatch: a_1<BR>
 *  3'-mismatch: b_m<BR>
 *  and for the enclosed base pair (V,U):<BR>
 *  5'-mismatch: b_1<BR>
 *  3'-mismatch: a_n<BR>
 *
 *  @note Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair
 *        must be 'turned arround' when evaluating the free energy of the internal-loop<br>
 *        This function is threadsafe
 *
 *  @see scale_parameters(), vrna_param_t
 *
 *  @param  n1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  n2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the internal loop
 *  @param  type_2  The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled energy parameters
 *  @return The Free energy of the Interior-loop in dcal/mol
 */
DEPRECATED(PRIVATE INLINE int
           E_IntLoop(int          n1,
                     int          n2,
                     int          type,
                     int          type_2,
                     int          si1,
                     int          sj1,
                     int          sp1,
                     int          sq1,
                     vrna_param_t *P),
           "Use vrna_E_internal() instead!");

/**
 *  <H2>Compute Boltzmann weight @f$e^{-\Delta G/kT} @f$ of internal loop</H2>
 *  multiply by scale[u1+u2+2] for scaling
 *
 *  @note This function is threadsafe
 *
 *  @see get_scaled_pf_parameters(), vrna_exp_param_t, E_IntLoop()
 *
 *  @param  u1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  u2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the internal loop
 *  @param  type2   The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled Boltzmann weights of the energy parameters
 *  @return The Boltzmann weight of the Interior-loop
 */
DEPRECATED(PRIVATE INLINE FLT_OR_DBL
           exp_E_IntLoop(int              u1,
                         int              u2,
                         int              type,
                         int              type2,
                         short            si1,
                         short            sj1,
                         short            sp1,
                         short            sq1,
                         vrna_exp_param_t *P),
           "Use vrna_exp_E_internal() instead!");

DEPRECATED(PRIVATE INLINE int
           E_IntLoop_Co(int           type,
                        int           type_2,
                        int           i,
                        int           j,
                        int           p,
                        int           q,
                        int           cutpoint,
                        short         si1,
                        short         sj1,
                        short         sp1,
                        short         sq1,
                        int           dangles,
                        vrna_param_t  *P),
           "This function is obsolete");

/*
 *  ugly but fast internal loop evaluation
 *
 *  Avoid including this function in your own code. It only serves
 *  as a fast inline block internally re-used throughout the RNAlib. It
 *  evalutes the free energy of internal loops in single sequences or sequence
 *  hybrids. Soft constraints are also applied if available.
 *
 *  NOTE: do not include into doxygen reference manual!
 */
DEPRECATED(PRIVATE INLINE int
           ubf_eval_int_loop(int            i,
                             int            j,
                             int            p,
                             int            q,
                             int            i1,
                             int            j1,
                             int            p1,
                             int            q1,
                             short          si,
                             short          sj,
                             short          sp,
                             short          sq,
                             unsigned char  type,
                             unsigned char  type_2,
                             int            *rtype,
                             int            ij,
                             int            cp,
                             vrna_param_t   *P,
                             vrna_sc_t      *sc),
           "This function is obsolete");

DEPRECATED(PRIVATE INLINE int
           ubf_eval_int_loop2(int           i,
                              int           j,
                              int           p,
                              int           q,
                              int           i1,
                              int           j1,
                              int           p1,
                              int           q1,
                              short         si,
                              short         sj,
                              short         sp,
                              short         sq,
                              unsigned char type,
                              unsigned char type_2,
                              int           *rtype,
                              int           ij,
                              unsigned int  *sn,
                              unsigned int  *ss,
                              vrna_param_t  *P,
                              vrna_sc_t     *sc),
           "This function is obsolete");

/*
 *  ugly but fast exterior internal loop evaluation
 *
 *  Avoid including this function in your own code. It only serves
 *  as a fast inline block internally re-used throughout the RNAlib. It
 *  evalutes the free energy of internal loops in single sequences or sequence
 *  hybrids. Soft constraints are also applied if available.
 *
 *  NOTE: do not include into doxygen reference manual!
 */
DEPRECATED(PRIVATE INLINE int
           ubf_eval_ext_int_loop(int            i,
                                 int            j,
                                 int            p,
                                 int            q,
                                 int            i1,
                                 int            j1,
                                 int            p1,
                                 int            q1,
                                 short          si,
                                 short          sj,
                                 short          sp,
                                 short          sq,
                                 unsigned char  type,
                                 unsigned char  type_2,
                                 int            length,
                                 vrna_param_t   *P,
                                 vrna_sc_t      *sc),
           "This function is obsolete");

DEPRECATED(int
           vrna_eval_int_loop(vrna_fold_compound_t  *fc,
                              int                   i,
                              int                   j,
                              int                   k,
                              int                   l),
           "Use vrna_eval_internal() instead!");

DEPRECATED(int
           vrna_E_stack(vrna_fold_compound_t  *fc,
                        int                   i,
                        int                   j),
           "Use vrna_eval_stack() instead!");

DEPRECATED(FLT_OR_DBL
           vrna_exp_E_interior_loop(vrna_fold_compound_t *fc,
                                    int                  i,
                                    int                  j,
                                    int                  k,
                                    int                  l),
           "Use vrna_exp_eval_internal() instead!");

/**
 * @}
 */

#endif

#endif
