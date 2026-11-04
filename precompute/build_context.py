"""Construct and cache the curve/quaternion context from scalar parameters."""

from functools import lru_cache
from types import SimpleNamespace

from sage.all import EllipticCurve, GF, QuaternionAlgebra, ZZ
from sage.misc.randstate import seed

from third_party.sqisign_prism_v2.ideals import O0_leftideal_create
from src.arithmetic.kummer import KummerLine, KummerPoint
from src.config.parameters import DEFAULT_PARAMETERS
from precompute.endomorphism import action, action_matrices, quaternion_change_basis
from precompute.torsion_points import (
    calcFields, Fp2ToFp2d, Fp2dToFp2, TorsionBasis_2e,
)


def points_order_N(curve, N):
    """Deterministically find an N-torsion point on the quadratic twist.

    Its exponent is p-1. Avoid factoring the full curve order just to find a
    point of prime order N; exact order is checked by scalar multiplication.
    """
    field = curve.base_ring()
    cofactor = (field.characteristic()-1)//N
    counter = 0
    while True:
        counter += 1
        x = field([1, counter])
        if not curve.is_x_coord(x):
            continue
        point = min(curve.lift_x(x, all=True), key=lambda pt: list(pt.y()))
        point = cofactor*point
        if not point.is_zero():
            if not (N*point).is_zero():
                raise ArithmeticError("point is not in the requested N-torsion")
            point.set_order(N)
            return point


def _build_quaternions(parameters):
    B = QuaternionAlgebra(-1, -parameters.p)
    i, j, k = B.gens()
    O0 = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])
    alpha = (parameters.alpha_scalar + parameters.alpha_i*i
             + parameters.alpha_frobenius*j)
    if alpha not in O0 or alpha.reduced_norm() % parameters.N:
        raise ValueError("alpha must lie in O0 and have norm divisible by N")
    I_P0 = O0_leftideal_create(alpha.conjugate(), parameters.N, O0)
    kappa = i + parameters.kappa_frobenius*j
    if parameters.require_nonsquare_kappa:
        discriminant = ZZ(kappa.reduced_trace()**2-4*kappa.reduced_norm())
        if discriminant.kronecker(parameters.N) != -1:
            raise ValueError("kappa discriminant must be a nonsquare modulo N")
    return SimpleNamespace(B=B, i=i, j=j, k=k, O0=O0, alpha=alpha,
                           I_P0=I_P0, kappa=kappa)


def _build_base_curves(parameters):
    Fp2 = GF(parameters.p**2, modulus=[1, 0, 1], names="i")
    E0_fp2 = EllipticCurve(Fp2, [1, 0])
    Et_fp2 = E0_fp2.quadratic_twist()
    if Et_fp2.order() % parameters.N:
        raise ValueError("the twist does not contain the requested odd torsion")
    return SimpleNamespace(Fp2=Fp2, E0_fp2=E0_fp2, Et_fp2=Et_fp2)


def _extend_curves(parameters, curves):
    Fp4, Fp2d, zeta4 = calcFields(parameters.p)
    E0_Fp4 = EllipticCurve(Fp4, [1, 0])
    Et_Fp4 = EllipticCurve(Fp4, [
        Fp2dToFp2(curves.Et_fp2.a4(), zeta4),
        Fp2dToFp2(curves.Et_fp2.a6(), zeta4),
    ])
    return SimpleNamespace(
        **vars(curves), Fp4=Fp4, Fp2d=Fp2d, zeta4=zeta4,
        E0_Fp4=E0_Fp4, Et_Fp4=Et_Fp4,
        iso_Et_Fp4_to_E0_Fp4=Et_Fp4.isomorphisms(E0_Fp4)[0],
        iso_E0_Fp4_to_Et_Fp4=E0_Fp4.isomorphisms(Et_Fp4)[0],
    )


def _build_odd_torsion(parameters, quaternions, curves, R0):
    N, Fp4, zeta4 = parameters.N, curves.Fp4, curves.zeta4
    R0_Fp4 = curves.Et_Fp4(*(Fp2dToFp2(v, zeta4) for v in R0.xy()))
    R0_E0 = curves.iso_Et_Fp4_to_E0_Fp4(R0_Fp4)
    alpha_coeffs = [ZZ(v) for v in quaternion_change_basis(quaternions.alpha, quaternions.O0)]
    kappa_coeffs = [ZZ(v) for v in quaternion_change_basis(quaternions.kappa, quaternions.O0)]
    P0_fp4 = action(alpha_coeffs, R0_E0, zeta4, Fp4)
    Q0_fp4 = action(kappa_coeffs, P0_fp4, zeta4, Fp4)
    P0_Et_Fp4 = curves.iso_E0_Fp4_to_Et_Fp4(P0_fp4)
    Q0_Et_Fp4 = curves.iso_E0_Fp4_to_Et_Fp4(Q0_fp4)
    P0_Et, Q0_Et = (
        curves.Et_fp2(*(Fp2ToFp2d(v, zeta4, curves.Fp2.gen()) for v in point.xy()))
        for point in (P0_Et_Fp4, Q0_Et_Fp4)
    )
    for point in (P0_Et, Q0_Et):
        if point.is_zero() or not (N*point).is_zero():
            raise ArithmeticError("computed odd-torsion basis has incorrect order")
        point.set_order(N)
    return SimpleNamespace(
        R0=R0, R0_Fp4=R0_Fp4, R0_E0=R0_E0,
        alpha_coeffs=alpha_coeffs, kappa_coeffs=kappa_coeffs,
        P0_fp4=P0_fp4, Q0_fp4=Q0_fp4,
        P0_Et_Fp4=P0_Et_Fp4, Q0_Et_Fp4=Q0_Et_Fp4,
        P0_Et=P0_Et, Q0_Et=Q0_Et,
    )


def _build_even_torsion(parameters, curves):
    E0_fp2d = EllipticCurve(curves.Fp2d, [1, 0])
    U0_Fp4, V0_Fp4 = TorsionBasis_2e(
        curves.E0_Fp4, curves.Fp2d, False, parameters.f
    )
    Ms = action_matrices((U0_Fp4, V0_Fp4), parameters.torsion_order,
                         curves.zeta4, curves.Fp4)
    U0_Fp2, V0_Fp2 = (
        curves.E0_fp2(*(Fp2ToFp2d(v, curves.zeta4, curves.Fp2.gen()) for v in point.xy()))
        for point in (U0_Fp4, V0_Fp4)
    )
    line = KummerLine(curves.E0_fp2)
    difference = U0_Fp2-V0_Fp2
    return SimpleNamespace(
        E0_fp2d=E0_fp2d, U0_Fp4=U0_Fp4, V0_Fp4=V0_Fp4, Ms=Ms,
        U0_Fp2=U0_Fp2, V0_Fp2=V0_Fp2, E0_Fp2_kummer=line, UmV0_Fp2=difference,
        U0_Kummer=KummerPoint(line, U0_Fp2), V0_Kummer=KummerPoint(line, V0_Fp2),
        UmV0_Kummer=KummerPoint(line, difference),
    )


def get_context(parameters=DEFAULT_PARAMETERS):
    """Build once per parameter set, preserving the caller's random states.

    The signature arithmetic adapter uses DEFAULT_PARAMETERS. Other parameter
    sets can be constructed here for experiments without rebinding that adapter.
    """
    return _build_context(parameters)


@lru_cache(maxsize=None)
def _build_context(parameters):
    with seed(parameters.precomputation_seed):
        if not ZZ(parameters.p).is_prime() or not ZZ(parameters.N).is_prime():
            raise ValueError("p and N must both be prime")
        quaternions = _build_quaternions(parameters)
        curves = _build_base_curves(parameters)
        R0 = points_order_N(curves.Et_fp2, parameters.N)
        curves = _extend_curves(parameters, curves)
        odd = _build_odd_torsion(parameters, quaternions, curves, R0)
        even = _build_even_torsion(parameters, curves)
    return SimpleNamespace(
        parameters=parameters, p=parameters.p, f=parameters.f, c=parameters.c, N=parameters.N,
        **vars(quaternions), **vars(curves), **vars(odd), **vars(even),
    )
