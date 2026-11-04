"""Torsion-pair encodings relative to deterministic full 2^f bases.

Matrices act on a column of points: rows are coordinates of the output pair.
This module represents Sage objects, not a byte-level serialization format.
"""

from functools import lru_cache
from sage.all import ZZ, Zmod, matrix

from src.config import context
from src.arithmetic import torsion
from src.config.parameters import f


@lru_cache(maxsize=1)
def standard_basis():
    """The deterministic source basis used by the upstream parameter adapter."""
    return tuple(context.iso_E0_Fp4_to_Et_Fp4(P)
                 for P in (context.U0_Fp4, context.V0_Fp4))


def canonical_basis(curve, hint=None):
    """Return (basis, hint), with canonical signs and no caller RNG use.

    A=0 uses the deterministic point sampler, since the Montgomery hint
    construction requires A!=0. None is its only accepted hint.
    """
    if curve.base_ring() != context.Fp4:
        raise ValueError('expected the configured quartic field')
    if any(curve.a_invariants()[i] != expected
           for i, expected in ((0, 0), (2, 0), (3, 1), (4, 0))):
        raise ValueError('expected a Montgomery model y^2=x^3+A*x^2+x')
    if curve.a2() == 0:
        if hint is not None:
            raise ValueError('the A=0 basis has no hint')
        basis = torsion.TorsionBasis_2e(curve, context.Fp2d, False, f)
    else:
        if hint is None:
            _, hint = torsion.TorsionBasisToHint(curve, f)
        else:
            hA, h = hint
            hA, h = ZZ(hA), ZZ(h)
            if hA not in (0, 1) or not 0 <= h < 128:
                raise ValueError('invalid torsion basis hint')
            hint = (hA, h)
        basis = torsion.TorsionBasisFromHintFull(curve, hint, f)
    if not torsion.is_full_two_power_basis(*basis, f):
        raise ValueError('hint did not produce a full torsion basis')
    return basis, hint


def basis_matrix(M):
    """Validate a GL(2,Z/2^f Z) basis matrix and normalize its coefficient ring."""
    if M.nrows() != 2 or M.ncols() != 2:
        raise ValueError('expected a 2 by 2 basis matrix')
    # Do not silently coerce matrices over a different finite ring.
    if M.base_ring() not in (ZZ, Zmod(2**f)):
        raise ValueError('incorrect basis matrix coefficient ring')
    M = matrix(Zmod(2**f), M)
    if ZZ(M.det()) % 2 != 1:
        raise ValueError('basis matrix must have odd determinant')
    return M


def encode_proof(curve, P, Q):
    """Encode a full pair relative to the deterministic basis on its curve."""
    basis, hint = canonical_basis(curve)
    M = basis_matrix(torsion.ChangeOfBasis(basis, (P, Q)))
    return curve, M, hint


def decode_proof(proof):
    """Recover a full 2^f pair; callers must reject malformed public data."""
    curve, M, hint = proof
    M = basis_matrix(M)
    basis, _ = canonical_basis(curve, hint)
    return curve, *torsion.EvalMatrix(M, basis)


def decode_secret_key(pk, sk):
    """Expand (E,I,M_N,C0,U_can,V_can) to the historical point tuple."""
    E, I, mixing, C0, U, V = sk
    if E != pk[0]:
        raise ValueError('secret and public curves differ')
    U0, V0 = torsion.EvalMatrix(basis_matrix(C0), standard_basis())
    return E, I, mixing, U0, V0, U, V
