# SageMath imports
from sage.all import (
    EllipticCurve,
    PolynomialRing,
    matrix,
    identity_matrix
)
from src.arithmetic.torsion import PointCoords, point_coords_pairing
# Preserve the old import location with the shared PointCoords compatibility wrapper.
from src.arithmetic.torsion import bi_dlp as bi_dlp


def quaternion_change_basis(gamma, order):
    """Coordinates of gamma in the Z-basis of a quaternion order."""
    from sage.all import Matrix, QQ, vector

    basis = order.basis()
    change = Matrix(QQ, [b.coefficient_tuple() for b in basis]).inverse()
    coefficients = vector(QQ, gamma.coefficient_tuple())*change
    assert gamma == sum(c*b for c, b in zip(coefficients, basis))
    return coefficients

# action of square root of -1
def i_action(P, zeta4):
    F = P.base_ring()
    E = P.curve()
    assert zeta4 in F
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    X, Y, Z = P
    return E([-X, zeta4*Y, Z])

# Frobenius endomorphism
def Frobenius(P):
    p = P.base_ring().characteristic()
    E = P.curve()
    X, Y, Z = P
    return E([X**p, Y**p, Z**p])

# return retP s.t. 2*retP = P. Note that retP is over an extension field.
def half_point(P, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    assert F.is_subring(F2)

    E2 = EllipticCurve(F2, [1, 0])
    R = PolynomialRing(F2, name="X")
    X = R.gens()[0]
    if P.is_zero():
        f = X**3 + X
    else:
        f = P[2]*(X**2 - 1)**2 - P[0]*4*(X**3 + X)
    xs = f.roots(multiplicities=False)
    assert len(xs) > 0
    x = xs[0]
    y = (x**3 + x).sqrt()
    retP = E2([x, y])
    if not E(2*retP) == P:
        retP = -retP
    assert E(2*retP) == P
    return retP

# the action of (i + j)/2
def i_j_2_action(P, zeta4, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    halfP = half_point(P, F2)
    return E(i_action(halfP, zeta4) + Frobenius(halfP))

# the action of (1 + ij)/2
def one_ij_2_action(P, zeta4, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    halfP = half_point(P, F2)
    return E(halfP + i_action(Frobenius(halfP), zeta4))

# the action of a + bi + c(i + j)/2 + d(1 + ij)/2
def action(alpha, P, zeta4, F2):
    F = P.base_ring()
    E = P.curve()
    assert E == EllipticCurve(F, [1,0]) # P should be on the curve y^2 = x^3 + x
    a, b, c, d = alpha
    ret = a*P
    ret += b*i_action(P, zeta4)
    ret += c*i_j_2_action(P, zeta4, F2)
    ret += d*one_ij_2_action(P, zeta4, F2)
    return ret

# matrix of the multiplication by alpha w.r.t. basis of N-torsion subgroup
def action_matrix(alpha, basis, N, zeta4, F2):
    P, Q = basis
    aP, aQ = [action(alpha, R, zeta4, F2) for R in [P, Q]]
    zeta = point_coords_pairing(P, Q, order=N)
    a, b = PointCoords(aP, P, Q, zeta, order=N)
    c, d = PointCoords(aQ, P, Q, zeta, order=N)
    # Precomputation stores image coordinates in columns (upstream uses rows).
    return matrix([[a, c], [b, d]])

# matrices of the multiplications by 1, i, (i + j)/2, (1 + ij)/2
def action_matrices(basis, N, zeta4, F2):
    one = [1,0,0,0]
    Ms = [identity_matrix(2)]
    for alpha in [one[-i:]+one[:-i] for i in range(1,4)]:
        Ms.append(action_matrix(alpha, basis, N, zeta4, F2))
    return Ms

