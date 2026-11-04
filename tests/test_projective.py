"""Projective arithmetic versus independent Sage elliptic-curve operations."""


from unittest.mock import patch
from sage.all import EllipticCurve, GF, ZZ
from src.arithmetic.projective import JacobianPoint


def test_gluing_coordinate_view():
    """The unchanged gluing API sees the same coordinates as Sage points."""
    from third_party.sqisign_prism_v2.theta.theta_structures.couple_point import CouplePoint
    from third_party.sqisign_prism_v2.theta.theta_isogenies.gluing_isogeny import GluingThetaIsogeny
    E = EllipticCurve(GF(11**2, name='i'), [1, 0])
    points = [P for P in E if (4*P).is_zero() and not (2*P).is_zero()]
    U = next(P for P in points if 2*P == E(0, 0))
    V = next(P for P in points if 2*P != 2*U)

    def scaled(P, scale):
        J = JacobianPoint(P)
        return J._new(scale**2*J._X, scale**3*J._Y, scale*J._Z)

    T1, T2 = CouplePoint(U, V), CouplePoint(V, U)
    J1 = CouplePoint(scaled(U, 2), scaled(V, 3))
    J2 = CouplePoint(scaled(V, 5), scaled(U, 7))
    expected = GluingThetaIsogeny.get_base_change_matrix(T1, T2)
    assert GluingThetaIsogeny.get_base_change_matrix(J1, J2) == expected
    gluing = object.__new__(GluingThetaIsogeny)
    gluing._base_change_matrix = expected
    assert gluing.base_change(J1) == gluing.base_change(T1)
    assert gluing.base_change(J2) == gluing.base_change(T2)
    for P, J in ((U, J1[0]), (V, J1[1])):
        assert tuple(J[i] for i in range(3)) == tuple(P)
        assert J[2] == 1
    print('theta gluing coordinate view matches Sage affine inputs: passed', flush=True)


def main():
    F = GF(11)
    curves = [EllipticCurve(F, coefficients) for coefficients in
              ([0, 0, 0, 1, 0], [0, 1, 0, 1, 0], [0, 3, 0, 1, 0], [4, 3])]
    for E in curves:
        for P in E:
            JP = JacobianPoint(P)
            scaled = JP._new(4*JP._X, 8*JP._Y, 2*JP._Z)
            assert scaled == JP and scaled.to_affine() == P
            for Q in E:
                JQ = JacobianPoint(Q)
                with patch.object(JacobianPoint, 'to_affine', side_effect=AssertionError('normalization')):
                    total = scaled+JQ
                    difference = scaled-JQ
                    doubled = scaled.double()
                assert total.to_affine() == P+Q
                assert difference.to_affine() == P-Q
                assert doubled.to_affine() == 2*P
            for scalar in (*range(-13, 14), 2**80+12345, -(2**80+12345), 2**240-1, 2**239):
                with patch.object(JacobianPoint, 'to_affine', side_effect=AssertionError('normalization')):
                    actual = ZZ(scalar)*scaled
                assert actual.to_affine() == scalar*P
    print('Jacobian addition, doubling, signed scalars and projective rescaling: passed', flush=True)
    test_gluing_coordinate_view()


if __name__ == '__main__':
    main()
