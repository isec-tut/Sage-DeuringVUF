"""Compatibility access to the single default parameter/context definition.

New code should use src.config.parameters for integers and src.config.context for arithmetic
objects. Importing this facade does not initialize curves or reseed randomness.
"""

from src.config import context as _context
from src.config.parameters import DEFAULT_PARAMETERS, c, f, N, p

__all__ = [
    *_context.__all__, "points_order_N", "quaternion_change_basis", "Fp2ToFp4",
]


def __getattr__(name):
    return getattr(_context, name)


def __dir__():
    return sorted(set(globals()) | set(__all__))


def points_order_N(E, N):
    from precompute.build_context import points_order_N as build_point
    return build_point(E, N)


def quaternion_change_basis(gamma, order):
    from precompute.endomorphism import quaternion_change_basis as coordinates
    return coordinates(gamma, order)


def Fp2ToFp4(x, zeta):
    from precompute.torsion_points import Fp2dToFp2
    return Fp2dToFp2(x, zeta)


def main():
    _context.main()


if __name__ == "__main__":
    main()
