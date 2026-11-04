"""Current Kummer extensions and projective odd-basis validation.

Run: sage -python -m tests.test_kummer
"""

from unittest.mock import patch

from sage.all import EllipticCurve, GF

from precompute.torsion_points import point_ord
from src.arithmetic.kummer import KummerLine, KummerPoint
from src.arithmetic.projective import odd_linear_combination, valid_odd_basis


def basis():
    field = GF(19**2, name='i', modulus=[1, 0, 1])
    curve = EllipticCurve(field, [1, 0])
    P, counter = point_ord(curve, field, False, 5, 1, return_next=True)
    for _ in range(100):
        Q, counter = point_ord(curve, field, False, 5, 1,
                               start=counter, return_next=True)
        if P.weil_pairing(Q, 5) != 1:
            return curve, P, Q
    raise AssertionError('failed to find an independent basis')


def test_kummer_extensions():
    curve, P, Q = basis()
    line = KummerLine(curve)
    for A, B in ((P, Q), (P+Q, Q), (-P, Q), (P, -Q)):
        triple = tuple(line(T) for T in (A, B, A-B))
        for scales in ((1, 1, 1), (2, 3, 4)):
            R, S, D = tuple(line(tuple(scale*c for c in point.XZ()))
                            for point, scale in zip(triple, scales))
            recovered = R.y_recovery(S, D)
            assert recovered in ((A, B), (-A, -B))
            assert line(R.lift()) == R
            for a, b in ((0, 0), (0, 1), (1, 0), (1, 1), (1, -1),
                         (-7, 13), (5, 10), (2**80+7, -2**90-3)):
                expected = line(a*A+b*B)
                with (
                    patch.object(KummerPoint, 'lift', side_effect=AssertionError('lifting')),
                    patch.object(KummerPoint, 'x', side_effect=AssertionError('normalization')),
                ):
                    assert odd_linear_combination(R, S, D, a, b, 5) == expected


def test_public_basis_validation():
    curve, P, Q = basis()
    line = KummerLine(curve)
    triple = tuple(line(T) for T in (P, Q, P-Q))
    assert valid_odd_basis(curve, triple, 5)
    assert valid_odd_basis(curve, (line(P), line(Q), line(P+Q)), 5)
    assert not valid_odd_basis(curve, triple, 3)
    R, S, D = triple
    for malformed in ((R, S), (R, S, D, R),
                      (line.zero(), S, D), (R, line((0, 0)), D),
                      (R, R, D), (R, S, R)):
        assert not valid_odd_basis(curve, malformed, 5)
    # Rational two-torsion lies on the curve but is not killed by odd order 5.
    assert not valid_odd_basis(curve, (line(curve(0, 0)), S, D), 5)
    outsider = next(x for x in curve.base_ring() if not curve.is_x_coord(x))
    assert not valid_odd_basis(curve, (line((outsider, 1)), S, D), 5)
    scaled = tuple(line(tuple(scale*c for c in point.XZ()))
                   for point, scale in zip(triple, (2, 3, 4)))
    with patch.object(KummerPoint, 'lift', side_effect=AssertionError('lifting')):
        assert valid_odd_basis(curve, scaled, 5)


def main():
    test_kummer_extensions()
    test_public_basis_validation()
    print('Kummer recovery, signed linear combinations, and public basis validation: passed')


if __name__ == '__main__':
    main()
