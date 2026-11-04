"""Projective odd-order chains, including exceptional multiples and signs."""

if __package__ in (None, ""):
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


from unittest.mock import patch

from sage.all import EllipticCurve, GF
from sage.misc.randstate import seed
from src.arithmetic.endomorphism_odd import OddTorsionBasis
from src.arithmetic.kummer import KummerPoint
from precompute.endomorphism import Frobenius, action, i_action


def _basis():
    field = GF(11**4, name='z')
    curve = EllipticCurve(field, [1, 0])
    zeta = min(field(-1).sqrt(all=True))
    with seed(0):
        for _ in range(100):
            R = 24*curve.random_point()
            P = Frobenius(R)+2*R
            if not P.is_zero():
                Q = i_action(P, zeta)
                return P, Q, zeta, OddTorsionBasis(P, Q, 5, 2, zeta)
    raise AssertionError('failed to find test eigenbasis')


def test_projective_chains():
    P, Q, zeta, prepared = _basis()
    for a in range(-6, 7):
        for b in range(-6, 7):
            assert prepared._image(a, b) == prepared.parent(a*P+b*Q)

    coefficients = (
        (0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0),
        (0, 0, 1, 0), (0, 0, 0, 1), (2, -3, 1, 4),
        (-5, 2, -3, 1), (2**300+7, 11, -3, 5),
    )
    for alpha in coefficients:
        A, B = (action(alpha, point, zeta, P.base_ring()) for point in (P, Q))
        for scalar in (0, 1, 2, -1, 5, 17):
            expected = tuple(prepared.parent(scalar*point) for point in (A, B, A-B))
            with patch.object(KummerPoint, 'x', side_effect=AssertionError('affine normalization')):
                images = prepared.images(alpha, scalar)
            assert images == expected, (alpha, scalar)

    field = P.base_ring()
    for name, factor in (('P', field.gen()+1), ('Q', field.gen()+2),
                         ('difference', field.gen()+3), ('sum', field.gen()+4)):
        X, Z = getattr(prepared, name)
        setattr(prepared, name, (factor*X, factor*Z))
    for a in range(5):
        for b in range(5):
            assert prepared._image(a, b) == prepared.parent(a*P+b*Q)


def test_rejects_incorrect_basis():
    P, Q, zeta, _ = _basis()
    for wrong_Q, wrong_order, wrong_lambda in ((P, 5, 2), (Q, 4, 2), (Q, 5, 1)):
        try:
            OddTorsionBasis(P, wrong_Q, wrong_order, wrong_lambda, zeta)
        except ValueError:
            pass
        else:
            raise AssertionError('invalid eigenbasis accepted')


def test_general_alpha_and_kappa():
    P, Q, zeta, _ = _basis()
    for base, lam, mu, t in ((P, 2, 0, 1), (P+Q, 0, 2, 0), (P+Q, 0, 2, 1)):
        second = i_action(base, zeta)+t*Frobenius(base)
        prepared = OddTorsionBasis(base, second, 5, lam, zeta,
                                  frobenius_i=mu, kappa_frobenius=t)
        for coefficients in ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0),
                             (0, 0, 0, 1), (2, -3, 1, 4)):
            A, B = (action(coefficients, R, zeta, R.base_ring()) for R in (base, second))
            for scalar in (1, -1, 3):
                expected = tuple(prepared.parent(scalar*R) for R in (A, B, A-B))
                with patch.object(KummerPoint, 'x', side_effect=AssertionError('affine normalization')):
                    assert prepared.images(coefficients, scalar) == expected


def main():
    test_projective_chains()
    test_rejects_incorrect_basis()
    test_general_alpha_and_kappa()
    print('odd endomorphism: all tests passed')


if __name__ == '__main__':
    main()
