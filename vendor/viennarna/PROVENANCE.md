# ViennaRNA provenance

Version: 2.7.2, upstream ViennaRNA/ViennaRNA.
Source: https://pypi.org/project/ViennaRNA/2.7.2/
Archive: https://files.pythonhosted.org/packages/d0/7a/2165d7d54234efc63a929a33d8a746cc23c3834747d195bd4d1b56cd1124/viennarna-2.7.2.tar.gz

COPYING is preserved verbatim. This is a private, modified energy-only subset,
not the full upstream Python or RNAlib distribution.

Modifications:

- Keep only the MFE scalar for single linear RNA with default model settings.
- Remove dot-bracket traceback and circular postprocessing from vrna_mfe.
- Prune comparative-alignment and G-quadruplex option branches, which are disabled
  by the fixed defaults used by the only entry point. Partition-function and
  window/2D-fold option branches are removed as well.
- Remove functions/variables not reachable from vrna_fold, and unused declarations.
- Remove parameter-file loading and filename bookkeeping; retain compiled default
  Turner 2004 energy tables and their required enthalpy data only.
- Remove unused includes and every source/header/data file outside the resulting
  compiler dependency closure. sources.txt is the exact C compilation allowlist.
- Keep the upstream config.h with HAVE_MALLOC_H disabled for macOS portability.
- Build no SWIG, C++, SVM/dlib, SIMD, OpenMP, plotting or external executables.

The Cython extension calls only vrna_fold and hides all internal C symbols.
Folding values are checked against independent RNAfold 2.7.2 CLI and Python
references, including randomized sequences. Builds need no pruning scripts,
Clang Python bindings, downloads or the upstream source tree.

Archive SHA-256: 44682b8a5f077faf32bf5d0353a232c74091ea47255a76a56495f9fb02e1e872

The final subset contains 30 C translation units. The default energy definitions
include intl11/intl21/intl22 and their enthalpy headers: these define globals used
across translation units and are required even when syntax-only checks pass
without them. Cephes const.c retains only MACHEP, MAXLOG, MAXNUM and PI. G-quadruplex
sparse-matrix cleanup is removed because those matrices are never allocated by
this fixed-default implementation. A clean linked-extension import and independent
reference comparisons verify these dependencies, in addition to syntax checks.
