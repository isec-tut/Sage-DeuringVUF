"""Arithmetic backend conventions and IdealToIsogeny images.

Run: sage -python -m tests.test_backend
"""

import random
from unittest.mock import patch

from sage.all import ZZ, inverse_mod, matrix, set_random_seed
from src.config import context
from src.config.constants import secret_ideal_norm
from src.protocol.keygen import _mix_odd_basis
from src.arithmetic.backend import (
    upstream_ec, REPRESENTATION_EXPONENT, _odd_basis_images,
    IdealToIsogenyKeygen, RandomEquivalentPrimeIdeal, RandomIdealGivenNorm,
)
from src.arithmetic.kummer import KummerPoint
from src.arithmetic.projective import JacobianPoint
from src.encoding.encode import standard_basis
from precompute.endomorphism import action, quaternion_change_basis


def test_action_conventions():
    assert context.f == context.DEFAULT_PARAMETERS.f
    assert REPRESENTATION_EXPONENT == ZZ(context.N).nbits()
    assert REPRESENTATION_EXPONENT <= context.f
    i, j, k = context.B.gens()
    for theta in (context.B(1), i, (i+j)/2, (1+k)/2):
        a, b, c, d = theta
        coefficients = tuple(map(ZZ, (a-d, b-c, 2*c, 2*d)))
        expected = tuple(
            context.iso_E0_Fp4_to_Et_Fp4(action(coefficients, P, context.zeta4, context.Fp4))
            for P in (context.U0_Fp4, context.V0_Fp4)
        )
        assert upstream_ec.EvalBasisPrecomp(theta) == expected
        actual = upstream_ec.EvalBasisPrecomp(theta, tuple(JacobianPoint(P) for P in standard_basis()))
        assert tuple(P.to_affine() for P in actual) == expected



def test_odd_basis_images():
    """Compare x-only images on the full non-smooth N-torsion with affine ones."""
    i, j, k = context.B.gens()
    corrections = (
        context.B(1), i, (i+j)/2, (1+k)/2,
        2-3*i+(i+j)/2+4*(1+k)/2,
        (ZZ(2)**300+7)+11*i-3*(i+j)/2+5*(1+k)/2,
    )
    for correction in corrections:
        coefficients = tuple(ZZ(v) for v in quaternion_change_basis(correction, context.O0))
        for degree in (ZZ(1), ZZ(17)):
            inverse = ZZ(inverse_mod(degree, context.N))
            expected = tuple(
                inverse*context.iso_E0_Fp4_to_Et_Fp4(
                    action(coefficients, point, context.zeta4, context.Fp4)
                )
                for point in (context.P0_fp4, context.Q0_fp4)
            )
            actual = _odd_basis_images(correction, degree)
            assert actual == expected or actual == tuple(-point for point in expected)
    print('odd N-basis x-only images: passed', flush=True)



def test_keygen_kummer_images():
    """Check both returned triples and odd mixing against full-point arithmetic."""
    random.seed(11)
    set_random_seed(11)
    success, ideal = RandomEquivalentPrimeIdeal(
        RandomIdealGivenNorm(secret_ideal_norm(), prime=True)
    )
    assert success
    success, E, P, Q, PmQ, U, V, UmV = IdealToIsogenyKeygen(ideal)
    assert success
    images = (P, Q, PmQ, U, V, UmV)
    assert all(isinstance(X, KummerPoint) for X in images)
    line = P.parent()
    assert line.curve() == E
    assert all(X.parent() == line for X in images)
    assert all(not X.is_zero() and (context.N*X).is_zero() for X in (P, Q, PmQ))
    assert all(not X.is_zero() and (2**context.f*X).is_zero() for X in (U, V, UmV))
    Pa, Qa = P.y_recovery(Q, PmQ)
    assert KummerPoint(line, Pa-Qa) == PmQ
    odd_pairing = context.P0_fp4.weil_pairing(context.Q0_fp4, context.N)
    assert Pa.weil_pairing(Qa, context.N) == odd_pairing**ZZ(ideal.norm())

    N = context.N
    matrices = (
        [[1, 0], [0, 1]], [[0, 1], [1, 0]], [[0, 1], [1, 1]],
        [[N-1, 0], [0, 1]], [[2, 5], [3, 7]], [[-3, 4], [5, -7]],
        [[1, 1], [0, 0]], [[0, 0], [0, 0]],
        [[2*N+3, -N-7], [N+5, 3*N+2]],
    )
    scaled = tuple(line((scale*X, scale*Z))
                   for point, scale in zip((P, Q, PmQ), (2, 3, 5))
                   for X, Z in (point.XZ(),))
    for entries in matrices:
        M = matrix(ZZ, entries)
        R = M[0, 0]*Pa+M[1, 0]*Qa
        S = M[0, 1]*Pa+M[1, 1]*Qa
        for basis in ((P, Q, PmQ), scaled):
            with (
                patch.object(KummerPoint, 'x', side_effect=AssertionError('affine normalization')),
                patch.object(KummerPoint, 'y_recovery', side_effect=AssertionError('y recovery')),
            ):
                actual = _mix_odd_basis(*basis, M, N)
            for X, expected in zip(actual, (R, S, R-S)):
                if expected.is_zero():
                    assert X.is_zero()
                else:
                    assert not X.is_zero() and X == KummerPoint(line, expected)
    print('keygen Kummer triples and projective odd-basis mixing: passed', flush=True)


def main():
    test_action_conventions()
    test_odd_basis_images()
    test_keygen_kummer_images()
    print('backend: all tests passed', flush=True)


if __name__ == '__main__':
    main()
