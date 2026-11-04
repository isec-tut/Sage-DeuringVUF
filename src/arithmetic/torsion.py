"""Canonical torsion bases, hints and Weil-pairing basis coordinates."""

from sage.all import ZZ, Zmod, kronecker, matrix
from src.arithmetic.montgomery import Ladder, ProjectiveDifference, lift_x, order_twof
from precompute.torsion_points import (
    deterministic_field_element, point_ord, TorsionBasis_2e,
)


def sqrt_minus_one(F):
    q = ZZ(F.cardinality())
    i = F.multiplicative_generator()**((q - 1)//4)
    if i**2 != -1:
        i = F(-1).sqrt(extend=False)
    assert i**2 == -1
    return i


def is_montgomery_x(A, x):
    return (x**3 + A*x**2 + x).is_square()


def torsion_hint_coefficient(A):
    if hasattr(A, "a2"):
        A = A.a2()
    F = A.parent()
    p = ZZ(F.characteristic())
    if ZZ(F.degree()) > 2:
        if A**(p**2) != A:
            raise ValueError("Montgomery coefficient is not defined over the quadratic subfield")
        F = F.subfield(2)
        A = F(A)
    return A


def torsion_basis_from_hint_x(A, hint, e):
    A = torsion_hint_coefficient(A)
    F = A.parent()
    hA, h = hint
    hA = ZZ(hA)
    h = ZZ(h)

    if h == 0:
        basis, _ = TorsionBasisToHint(A, e)
        return basis

    if hA == 1:
        i = sqrt_minus_one(F)
        xR = -A/(F(1) + i*F(h))
    else:
        xR = F(h)*A
    xRS = -xR - A

    scalar = (ZZ(F.characteristic()) + 1)//(2**e)
    E24 = (A + F(2), F(4))
    xR = Ladder((xR, F(1)), E24, scalar)
    xRS = Ladder((xRS, F(1)), E24, scalar)
    xS = ProjectiveDifference(xR, xRS, (A, F(1)))
    return xR, xS, xRS


def TorsionBasisToHint(A, e):
    A = torsion_hint_coefficient(A)
    F = A.parent()
    assert A != 0

    hA = ZZ(1) if A.is_square() else ZZ(0)
    h = ZZ(0)
    i = sqrt_minus_one(F) if hA == 1 else None
    while True:
        h += 1
        if hA == 1:
            # The field characteristic is already prime; Kronecker avoids
            # repeating the primality checks performed by legendre_symbol.
            if kronecker(ZZ(1 + h**2), ZZ(F.characteristic())) != -1:
                continue
            xR = -A/(F(1) + i*F(h))
        else:
            xR = F(h)*A
        if is_montgomery_x(A, xR):
            break

    basis = torsion_basis_from_hint_x(A, (hA, h), e)
    if h >= 128:
        h = ZZ(0)
    return basis, (hA, h)


def TorsionBasisFromHint(A, hint, e):
    return torsion_basis_from_hint_x(A, hint, e)


def is_full_two_power_basis(P, Q, e):
    from src.arithmetic.projective import kummer_point
    if e < 1:
        return False
    P, Q = kummer_point(P), kummer_point(Q)
    half = ZZ(2)**(e-1)
    P2, Q2 = (P, Q) if e == 1 else (half*P, half*Q)
    return (not P2.is_zero() and not Q2.is_zero()
            and P2.double().is_zero() and Q2.double().is_zero() and P2 != Q2)


def LiftXOnlyTorsionBasis(E, basis, e):
    F = E.base_ring()
    basis = tuple((F(X), F(Z)) for X, Z in basis)
    xR, xS, xRS = basis
    candidates = [(xR, xS, xRS), (xR, xRS, xS)]
    for P, Q, PmQ in candidates:
        try:
            U, V = lift_x(E, P, Q, PmQ)
        except (AssertionError, ArithmeticError, ValueError, ZeroDivisionError):
            continue
        if is_full_two_power_basis(U, V, e):
            return U, V
    raise ValueError("x-only torsion basis could not be lifted")


def TorsionBasisFromHintFull(E, hint, e):
    basis = TorsionBasisFromHint(E.a2(), hint, e)
    P, Q, PmQ = basis
    return lift_x(E, P, Q, PmQ)


def bi_dlp(P, Q, R, S, n):
    """Compatibility wrapper: P=aR+bS, Q=cR+dS modulo n, via PointCoords."""
    zeta = point_coords_pairing(R, S, order=n)
    a, b = PointCoords(P, R, S, zeta, order=n)
    c, d = PointCoords(Q, R, S, zeta, order=n)
    return a, b, c, d


def _coordinate_order(e, order):
    """Resolve (target order, source-basis scaling) without context initialization."""
    if order is not None:
        if e is not None:
            raise ValueError('provide either e or order, not both')
        order = ZZ(order)
        if order <= 1:
            raise ValueError('expected a torsion order greater than one')
        return order, ZZ(1)
    from src.config.parameters import f
    exponent = f if e is None else ZZ(e)
    if not 1 <= exponent <= f:
        raise ValueError('expected a torsion exponent between 1 and f')
    return ZZ(2)**exponent, ZZ(2)**(f-exponent)


def point_coords_pairing(P, Q, e=None, *, order=None):
    """Shared Weil base; order=n uses an already n-torsion basis without scaling."""
    from third_party.sqisign_prism_v2.theta.utilities.discrete_log import weil_pairing_pari
    n, scale = _coordinate_order(e, order)
    from src.arithmetic.projective import JacobianPoint
    P, Q = ((scale*JacobianPoint(T)).to_affine() for T in (P, Q))
    return weil_pairing_pari(P, Q, n)


def PointCoords(X, P, Q, zeta=None, e=None, *, order=None):
    """PointCoords semantics using upstream BiDLP's Weil/PARI implementation.

    With e, return a,b such that X=2^(f-e)*(aP+bQ) in a full 2^f basis.
    With order=n, P,Q are already an n-torsion basis and X=aP+bQ modulo n.
    These modes are mutually exclusive; omitting both means e=f.
    The signature implementation's
    Fp2 Tate reduction is not valid as a drop-in on the quartic twist; see
    docs/basis_encoding.md for a reproduced top-bit error. zeta is the Weil base from
    point_coords_pairing(P,Q,e), rather than an upstream Tate pairing value.
    """
    from third_party.sqisign_prism_v2.theta.utilities.discrete_log import BiDLP
    n, scale = _coordinate_order(e, order)
    from src.arithmetic.projective import JacobianPoint
    JX, JP, JQ = JacobianPoint(X), JacobianPoint(P), JacobianPoint(Q)
    scaled = tuple((scale*T).to_affine() for T in (JP, JQ))
    a, b = BiDLP(JX.to_affine(), *scaled, n, ePQ=zeta)
    if scale*(a*JP+b*JQ) != JX:
        raise ArithmeticError('point coordinates did not reconstruct the input')
    return a, b


def ChangeOfBasis(B1, B2, e=None, *, order=None):
    """Return M with B2 = 2^(f-e) M B1 (column of points, matrix rows).

    B1 is a full 2^f basis (None selects the standard source basis);
    B2 is a 2^e basis. With e=f this is the same
    matrix convention as upstream ec.ChangeOfBasis, using Weil/PARI BiDLP.
    With order=n, both bases have order n and B2=M*B1 without scaling.
    """
    n, _ = _coordinate_order(e, order)
    if B1 is None:
        if order is not None:
            raise ValueError('order mode requires an explicit source basis')
        from src.encoding.encode import standard_basis
        B1 = standard_basis()
    P, Q = B1
    zeta = point_coords_pairing(P, Q, e, order=order)
    rows = [PointCoords(X, P, Q, zeta, e=e, order=order) for X in B2]
    if len(rows) != 2:
        raise ValueError('expected a pair of target points')
    return matrix(Zmod(n), rows)


def EvalMatrix(M, basis):
    """Use upstream matrix semantics, Jacobian arithmetic, then return full points."""
    from third_party.sqisign_prism_v2.ec import EvalMatrix as upstream_eval
    from src.arithmetic.projective import JacobianPoint
    images = upstream_eval(M, tuple(JacobianPoint(P) for P in basis))
    return tuple(P.to_affine() for P in images)
