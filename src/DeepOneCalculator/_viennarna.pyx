"""Minimal Cython interface to the vendored ViennaRNA 2.7.2 RNAlib."""
cdef extern from "ViennaRNA/mfe/global.h":
    float vrna_fold(const char *sequence, char *structure)


def fold_energy(str sequence):
    """Compute RNA minimum free energy using default Turner 2004 parameters.

    Parameters
    ----------
    sequence : str
        Nonempty RNA or DNA sequence containing only A, C, G, U or T.

    Returns
    -------
    float
        Minimum free energy in kcal/mol at 37 degrees Celsius.

    Raises
    ------
    ValueError
        If the sequence is empty or contains invalid bases.
    """
    normalized = sequence.upper().replace("T", "U")
    if not normalized or set(normalized) - set("ACGU"):
        raise ValueError("sequence must be nonempty and contain only A, C, G, U or T")
    cdef bytes encoded = normalized.encode("ascii")
    # Keep the GIL: RNAlib defaults and initialization include shared state.
    # The vendored energy-only implementation has no structure traceback.
    return float(vrna_fold(encoded, NULL))
