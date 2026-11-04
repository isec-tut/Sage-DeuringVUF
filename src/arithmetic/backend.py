"""DeuringVUF adapters for the pinned SQIsign / PRISM arithmetic.

The local fork adds Kummer helpers and a parameter adapter. Keep twist/field
and API conventions here; see third_party/sqisign_prism_v2/UPSTREAM.md.
"""

from functools import lru_cache

from sage.all import ZZ, Matrix, gcd, inverse_mod
from src.config import context
from third_party.sqisign_prism_v2 import ec as upstream_ec
from third_party.sqisign_prism_v2 import hd, qlapoti as upstream_qlapoti
from third_party.sqisign_prism_v2 import quaternions as qt
from third_party.sqisign_prism_v2.theta.theta_structures.couple_point import CouplePoint as UpstreamCouplePoint
from third_party.sqisign_prism_v2.theta.theta_isogenies.product_isogeny import EllipticProductIsogeny
# Preserve the local Kummer API (y_recovery and LadderBiscalar extensions).
from src.arithmetic.kummer import KummerLine, KummerPoint
from src.arithmetic.endomorphism_odd import OddTorsionBasis
from precompute.endomorphism import quaternion_change_basis
from src.arithmetic.torsion import PointCoords, point_coords_pairing
from src.encoding.encode import basis_matrix as validate_basis_matrix
from src.arithmetic.projective import JacobianPoint, kummer_point, lift_kummer
from src.arithmetic.montgomery import order_twof

from src.config.constants import (
    TORSION_EXPONENT, REPRESENTATION_EXPONENT, REPRESENTATION_COFACTOR,
)

def CouplePoint(P, Q):
    """Use signed Jacobian points internally; lift Kummer inputs only here."""
    P = lift_kummer(P) if isinstance(P, KummerPoint) else P
    Q = lift_kummer(Q) if isinstance(Q, KummerPoint) else Q
    return UpstreamCouplePoint(JacobianPoint(P), JacobianPoint(Q))


def RandomIdealGivenNorm(norm, prime):
    norm = ZZ(norm)
    for _ in range(100):
        result = qt.RandomIdealGivenNorm(norm, prime)
        if result is not False:
            if result.norm() != norm or result.left_order() != context.O0:
                raise ArithmeticError("upstream ideal has incorrect norm or order")
            return result
    raise RuntimeError("RandomIdealGivenNorm failed after 100 attempts")


def RandomEquivalentPrimeIdeal(ideal):
    # Preserve the existing (success, ideal) interface.
    return True, qt.RandomEquivalentPrimeIdeal(ideal)


def _kernel_tuple(P, Q):
    return ((P[0], P[1]), (Q[0], Q[1]))


def Isogeny22Chain(P, Q, exponent, points):
    phi = EllipticProductIsogeny((P, Q), exponent)
    return phi.codomain(), [phi(CouplePoint(*pt), lift=False) for pt in points]


def Isogeny22ChainSqrt(P, Q, exponent, points):
    phi = hd.Dim2Iso(_kernel_tuple(P, Q), exponent)
    return phi.codomain(), [phi(CouplePoint(*pt), lift=False) for pt in points]


def _embedded_images(phi, degree, basis):
    """EmbeddedIsogeny with Weil pairing over the quartic twist field.

    The upstream reduced Tate exponent assumes Fp2. Here both the 2-power
    points and odd-order points live on Et over Fp4, so use the full-order
    Weil pairing to identify the degree-d factor instead.
    """
    P, Q = basis
    images = _product_pair_images(phi, P, Q)
    target = P.weil_pairing(Q, 2**TORSION_EXPONENT)**degree
    for idx in (0, 1):
        R, S = _recover_pair(phi.codomain()[idx], images, idx)
        if R.weil_pairing(S, 2**TORSION_EXPONENT) == target:
            return idx, R, S
    raise ArithmeticError("neither product factor has the required pairing")


def _product_pair_images(phi, P, Q):
    """Evaluate a differential triple without lifting either product factor."""
    P, Q = JacobianPoint(P), JacobianPoint(Q)
    zero = phi.E2(0)
    return tuple(phi(CouplePoint(point, zero), lift=False) for point in (P, Q, P-Q))


def _recover_pair(curve, images, idx):
    """Recover one consistent pair with a single square root on this factor."""
    P, Q, difference = _kummer_component(curve, images, idx)
    return P.y_recovery(Q, difference)


def _kummer_component(curve, images, idx):
    """Select a product factor, retaining the projective differential triple."""
    line = KummerLine(curve)
    return tuple(KummerPoint(line, image[idx]) for image in images)


def _kummer_images(curve, P, Q):
    line = KummerLine(curve)
    P, Q = JacobianPoint(P), JacobianPoint(Q)
    return tuple(kummer_point(T, line) for T in (P, Q, P-Q))


@lru_cache(maxsize=1)
def _odd_action_basis():
    """Derive the odd action from alpha and kappa without an N-torsion DLP."""
    parameters = context.DEFAULT_PARAMETERS
    eigenvalue = ZZ(parameters.alpha_scalar*inverse_mod(parameters.alpha_frobenius, context.N))
    return OddTorsionBasis(
        context.P0_fp4, context.Q0_fp4, context.N, eigenvalue,
        context.zeta4, field=context.Fp2d,
        frobenius_i=ZZ(parameters.alpha_i*inverse_mod(parameters.alpha_frobenius, context.N)),
        kappa_frobenius=parameters.kappa_frobenius,
    )


@lru_cache(maxsize=1)
def _standard_basis_pairing():
    P0 = context.iso_E0_Fp4_to_Et_Fp4(context.U0_Fp4)
    Q0 = context.iso_E0_Fp4_to_Et_Fp4(context.V0_Fp4)
    return point_coords_pairing(P0, Q0, TORSION_EXPONENT)


def _odd_basis_images(correction, degree):
    """Evaluate degree^-1 * correction on the N-basis using x-only arithmetic.

    Keep P, Q and P-Q together to preserve relative signs. Affine points are
    recovered only at the boundary with the existing twist/product-isogeny API;
    the recovered pair may differ by a simultaneous sign from affine evaluation.
    The prepared eigenbasis turns the action into two-scalar projective chains
    in the quadratic subfield, without field division or half-point searches.
    """
    N = context.N
    coefficients = tuple(
        ZZ(value) % N for value in quaternion_change_basis(correction, context.O0)
    )
    inverse_N = ZZ(inverse_mod(degree, N))
    # Fold inverse-degree multiplication into the action matrix modulo N.
    P, Q, difference = _odd_action_basis().images(coefficients, scalar=inverse_N)
    P, Q = P.y_recovery(Q, difference)
    return tuple(context.iso_E0_Fp4_to_Et_Fp4(point) for point in (P, Q))


def _ideal_to_isogeny(ideal, basis=None, odd_basis=False, *, basis_matrix=None):
    """Upstream IdealToIsogeny extended to input bases and odd-order points.

    Qlapoti always uses the available 2^f torsion. The smaller exponent used
    in the VUF proof is independent of this internal ideal conversion.
    """
    if basis is not None and basis_matrix is not None:
        raise ValueError('provide points or a basis matrix, not both')
    if basis_matrix is not None:
        basis_matrix = validate_basis_matrix(basis_matrix)
    for _ in range(100):
        mu1, mu2, theta, norm, reduced, beta = upstream_qlapoti.qlapoti(ideal)
        degree = ZZ(mu1.reduced_norm()/norm)
        # The odd images require division modulo N as well as modulo 2^f.
        if not odd_basis or gcd(degree, context.N) == 1:
            break
    else:
        raise RuntimeError("could not find a degree coprime to the odd torsion")

    inverse = ZZ(inverse_mod(degree, 2**TORSION_EXPONENT))
    P0 = context.iso_E0_Fp4_to_Et_Fp4(context.U0_Fp4)
    Q0 = context.iso_E0_Fp4_to_Et_Fp4(context.V0_Fp4)
    base = JacobianPoint(P0), JacobianPoint(Q0)
    theta_P, theta_Q = upstream_ec.EvalBasisPrecomp(inverse*theta, base)
    phi = hd.Dim2Iso(((base[0], theta_P), (base[1], theta_Q)), TORSION_EXPONENT)
    idx, P, Q = _embedded_images(phi, degree, (P0, Q0))
    curve = phi.codomain()[idx]

    correction = mu1*beta/norm
    P, Q = upstream_ec.EvalBasisPrecomp(
        inverse*correction, (JacobianPoint(P), JacobianPoint(Q))
    )
    if basis is not None:
        pairing = _standard_basis_pairing()
        a, b = PointCoords(basis[0], P0, Q0, pairing, e=TORSION_EXPONENT)
        c, d = PointCoords(basis[1], P0, Q0, pairing, e=TORSION_EXPONENT)
        basis_matrix = validate_basis_matrix(Matrix(ZZ, [[a, b], [c, d]]))
    if basis_matrix is not None:
        P, Q = upstream_ec.EvalMatrix(basis_matrix, (P, Q))

    # Check the degree of the recovered map, including its basis correction.
    # e(M B) = e(B)^det(M): matrix input avoids decoding then doing two DLPs.
    determinant = 1 if basis_matrix is None else ZZ(basis_matrix.det())
    expected = _standard_basis_pairing()**ZZ(ideal.norm()*determinant)
    if P.weil_pairing(Q, 2**TORSION_EXPONENT) != expected:
        raise ArithmeticError("IdealToIsogeny images have incorrect pairing")
    kummer = _kummer_images(curve, P, Q)
    if not odd_basis:
        return True, curve, *kummer

    # Endomorphism evaluation and inverse-degree scaling use Kummer coordinates.
    odd_points = _odd_basis_images(correction, degree)
    images = _product_pair_images(phi, *odd_points)
    odd_kummer = _kummer_component(curve, images, idx)
    return True, curve, *odd_kummer, *kummer


def IdealToIsogenyKeygen(ideal):
    """Return success, E, x(P), x(Q), x(P-Q), x(U), x(V), x(U-V).

    All six images are projective KummerPoint objects on the same line.
    P,Q are N-torsion images; U,V are the full 2^f-torsion images.
    """
    return _ideal_to_isogeny(ideal, odd_basis=True)


def IdealToIsogeny(ideal, U0, V0, UmV0, representation_exponent):
    if not 3 <= representation_exponent <= TORSION_EXPONENT:
        raise ValueError("invalid representation exponent")
    return _ideal_to_isogeny(ideal, basis=(U0, V0))


def IdealToIsogenyFromMatrix(ideal, change_of_basis, representation_exponent):
    """Convert using images of M*(P0,Q0), without decoding the source pair."""
    if not 3 <= representation_exponent <= TORSION_EXPONENT:
        raise ValueError('invalid representation exponent')
    return _ideal_to_isogeny(ideal, basis_matrix=change_of_basis)


def IsogVerif2dim(data, j):
    E1, P1, Q1, E3, P3, Q3 = data
    P = CouplePoint(context.N*P1, P3)
    Q = CouplePoint(context.N*Q1, Q3)
    exponent = REPRESENTATION_EXPONENT
    if not all(order_twof(kummer_point(T), exponent) for T in (*P.points(), *Q.points())):
        return False
    try:
        phi = hd.Dim2Iso(_kernel_tuple(P, Q), exponent)
    except (AssertionError, ArithmeticError, ValueError, ZeroDivisionError):
        return False
    for idx, curve in enumerate(phi.codomain()):
        if curve.j_invariant() == j:
            return phi, idx
    return False
