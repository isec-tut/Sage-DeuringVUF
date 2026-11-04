"""Derived protocol sizes without curve or torsion precomputation."""

from functools import lru_cache

from src.config.parameters import DEFAULT_PARAMETERS

SECURITY_BITS = DEFAULT_PARAMETERS.security_bits
TORSION_EXPONENT = DEFAULT_PARAMETERS.f
REPRESENTATION_EXPONENT = DEFAULT_PARAMETERS.representation_exponent
REPRESENTATION_COFACTOR = DEFAULT_PARAMETERS.representation_cofactor
AUXILIARY_NORM = DEFAULT_PARAMETERS.auxiliary_norm


@lru_cache(maxsize=None)
def secret_ideal_norm(security_bits=None):
    """Fixed public sampling norm; compute its next prime only once per size."""
    from sage.all import ZZ, next_prime
    bound = (DEFAULT_PARAMETERS.secret_ideal_bound if security_bits is None
             else ZZ(2)**(4*security_bits))
    return next_prime(ZZ(bound))
