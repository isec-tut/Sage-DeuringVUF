"""Exact-order PointCoords and the legacy bi_dlp coefficient convention."""


from sage.all import EllipticCurve, GF, ZZ, Zmod, discrete_log, matrix

from precompute.endomorphism import action, action_matrices
from precompute.torsion_points import TorsionBasis_2e, point_ord
from src.arithmetic.torsion import ChangeOfBasis, PointCoords, bi_dlp, point_coords_pairing
from src.arithmetic.projective import JacobianPoint


def legacy_coordinates(P, Q, R, S, n):
    """Independent Sage Weil/discrete_log oracle, using the former convention."""
    zeta = R.weil_pairing(S, n)
    return tuple(ZZ(discrete_log(value, zeta, ord=n)) for value in (
        P.weil_pairing(S, n), R.weil_pairing(P, n),
        Q.weil_pairing(S, n), R.weil_pairing(Q, n),
    ))


def check_exact_order(basis, n):
    P, Q = basis
    R, S = 3*P+2*Q, 5*P+7*Q  # determinant 11, coprime to all tested orders.
    zeta = point_coords_pairing(R, S, order=n)
    coefficients = ((0, 0, 1, 0), (-17, 19, n-1, n-2), (n, n+1, -3, -5))
    for a, b, c, d in coefficients:
        X, Y = a*R+b*S, c*R+d*S
        expected = tuple(ZZ(v) % n for v in (a, b, c, d))
        assert PointCoords(X, R, S, zeta, order=n) == expected[:2]
        assert PointCoords(Y, R, S, zeta, order=n) == expected[2:]
        assert bi_dlp(X, Y, R, S, n) == legacy_coordinates(X, Y, R, S, n) == expected
        assert ChangeOfBasis((R, S), (X, Y), order=n) == matrix(Zmod(n), 2, 2, expected)


def main():
    F = GF(31**2, name='i')
    E = EllipticCurve(F, [1, 0])
    full = TorsionBasis_2e(E, F, False, 5)
    for n in (ZZ(4), ZZ(16), ZZ(32)):
        basis = tuple((32//n)*P for P in full)
        check_exact_order(basis, n)
    JP, JQ = 3*JacobianPoint(full[0]), 5*JacobianPoint(full[1])
    JX = -17*JP+19*JQ
    zeta = point_coords_pairing(JP, JQ, order=32)
    assert PointCoords(JX, JP, JQ, zeta, order=32) == (15, 19)

    # The endomorphism precomputation also works for a basis smaller than the
    # default 2^240 basis: order=16 must not apply the default f-e scaling.
    basis = tuple(2*P for P in full)
    i = F(-1).sqrt()
    matrices = action_matrices(basis, ZZ(16), i, F)
    for coefficients, M in zip(((1, 0, 0, 0), (0, 1, 0, 0),
                                (0, 0, 1, 0), (0, 0, 0, 1)), matrices):
        for col, P in enumerate(basis):
            assert M[0, col]*basis[0]+M[1, col]*basis[1] == action(coefficients, P, i, F)

    # Preserve the old helper's ability to handle small odd/composite orders.
    Fodd = GF(19**2, name='j')
    Eodd = EllipticCurve(Fodd, [1, 0])
    P, counter = point_ord(Eodd, Fodd, False, 5, 1, return_next=True)
    while True:
        Q, counter = point_ord(Eodd, Fodd, False, 5, 1, start=counter, return_next=True)
        if P.weil_pairing(Q, 5) != 1:
            break
    check_exact_order((P, Q), ZZ(5))
    P4, Q4 = TorsionBasis_2e(Eodd, Fodd, False, 2)
    check_exact_order((P4+P, Q4+Q), ZZ(20))

    for kwargs in ({'order': 0}, {'order': 1}, {'order': -5}, {'e': 4, 'order': 16}):
        for call in (
            lambda: PointCoords(basis[0], *basis, **kwargs),
            lambda: point_coords_pairing(*basis, **kwargs),
            lambda: ChangeOfBasis(basis, basis, **kwargs),
        ):
            try:
                call()
            except ValueError:
                pass
            else:
                raise AssertionError(f'invalid order options accepted: {kwargs}')
    print('PointCoords exact orders, legacy coordinates and precomputation matrices: passed', flush=True)


if __name__ == '__main__':
    main()
