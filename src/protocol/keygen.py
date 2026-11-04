"""DeuringVUF key generation (paper Algorithm 1)."""

if __package__ in (None, ""):
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import random
from sage.all import ZZ, gcd, inverse_mod, matrix

from src.config import context
from src.arithmetic import torsion
from src.config.constants import secret_ideal_norm
from src.encoding.encode import basis_matrix, canonical_basis
from src.encoding.serialization import encode_public_key, encode_secret_key
from src.arithmetic.projective import odd_linear_combination
from src.arithmetic.backend import (
    IdealToIsogenyKeygen as IdealToIsogeny,
    RandomEquivalentPrimeIdeal,
    RandomIdealGivenNorm,
)


def _mix_odd_basis(P, Q, PmQ, mixing, order):
    """Apply the column-based mixing matrix to an odd-prime differential basis.

    Return x(R), x(S), x(R-S) without affine normalization or y recovery.
    Reduce each coefficient modulo the known order to avoid negative ladders.
    """
    def combine(a, b):
        return odd_linear_combination(P, Q, PmQ, a, b, order)

    a, b = mixing[0, 0], mixing[1, 0]
    c, d = mixing[0, 1], mixing[1, 1]
    return combine(a, b), combine(c, d), combine(a-c, b-d)


class Keygen:
    """Create byte-encoded pk and sk with deterministic-basis matrix encodings."""

    def __init__(self):
        N, f = context.N, context.f
        secret_norm = secret_ideal_norm()
        while True:
            I_sk = RandomIdealGivenNorm(secret_norm, prime=True)
            success, I_sk = RandomEquivalentPrimeIdeal(I_sk)
            if not success:
                continue
            success, curve, P, Q, PmQ, U, V, UmV = IdealToIsogeny(I_sk)
            if success:
                break

        # Sample a change of basis in GL(2, Z/NZ).
        while True:
            mixing = matrix(ZZ, 2, 2, [random.randint(0, N-1) for _ in range(4)])
            if gcd(mixing.det(), N) == 1:
                break
        R, S, RmS = _mix_odd_basis(P, Q, PmQ, mixing, N)

        N_sk_inv = ZZ(inverse_mod(I_sk.norm(), 2**f))
        U, V = (N_sk_inv*U).y_recovery(N_sk_inv*V, N_sk_inv*UmV)
        canonical, hint = canonical_basis(curve)
        C0 = basis_matrix(torsion.ChangeOfBasis((U, V), canonical))
        U_can, V_can = canonical
        pk = (curve, R, S, RmS, hint)
        sk = (curve, I_sk, mixing, C0, U_can, V_can)
        self.pk = encode_public_key(pk)
        self.sk = encode_secret_key(pk, sk)


if __name__ == "__main__":
    from examples.keygen import main
    main()
