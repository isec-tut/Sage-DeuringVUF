"""Shared construction of the two-dimensional proof (paper Algorithm 5)."""

import random
from sage.all import ZZ, inverse_mod

from src.config.parameters import N, f
from src.config.constants import REPRESENTATION_COFACTOR, REPRESENTATION_EXPONENT
from src.encoding.encode import encode_proof
from src.arithmetic.backend import CouplePoint, IdealToIsogenyFromMatrix, Isogeny22ChainSqrt


def build_representation(I, I_x, U, V, source_matrix, auxiliary_ideal, *, kernel=None,
                         auxiliary_curve=None, rng_state=None, auxiliary_intersection=None):
    """Return (Ex, proof, Fx) using a kernel check or a cached complementary curve.

    Normal evaluation identifies Ex by killing its known kernel. Amortized
    evaluation identifies the other factor using the cached j(Fx). Both paths
    use the same ideal conversion, normalization and product-isogeny kernel.
    The proof is (E3, C3, hint3), where C3 expresses (U3, V3) as rows
    relative to the deterministic full 2^f basis on E3.
    """
    I_aux = (I.intersection(auxiliary_ideal) if auxiliary_intersection is None
             else auxiliary_intersection)
    L = I_aux.intersection(I_x)
    if rng_state is not None:
        random.setstate(rng_state)
    success, E3, U3, V3, difference = IdealToIsogenyFromMatrix(
        L, source_matrix, REPRESENTATION_EXPONENT
    )
    assert success

    inverse = ZZ(inverse_mod(I.norm(), 2**f))
    U3, V3 = (inverse*U3).y_recovery(inverse*V3, inverse*difference)
    P, Q = CouplePoint(U, U3), CouplePoint(V, V3)
    P = REPRESENTATION_COFACTOR*CouplePoint(N*P[0], P[1])
    Q = REPRESENTATION_COFACTOR*CouplePoint(N*Q[0], Q[1])
    points = [(kernel, E3(0))] if auxiliary_curve is None else []
    curves, images = Isogeny22ChainSqrt(P, Q, REPRESENTATION_EXPONENT, points)

    if auxiliary_curve is None:
        if _is_kummer_zero(images[0][0]):
            idx = 0
        elif _is_kummer_zero(images[0][1]):
            idx = 1
        else:
            raise ValueError("Neither codomain factor kills the input kernel")
    else:
        auxiliary_j = auxiliary_curve.j_invariant()
        if curves[0].j_invariant() == auxiliary_j:
            idx = 1
        elif curves[1].j_invariant() == auxiliary_j:
            idx = 0
        else:
            raise ValueError("Neither codomain factor matches the auxiliary curve")

    # return curves[idx], (E3, U3, V3), curves[1-idx]
    return curves[idx], encode_proof(E3, U3, V3), curves[1-idx]


def _is_kummer_zero(point):
    return point.is_zero() if hasattr(point, "is_zero") else point[1] == 0
