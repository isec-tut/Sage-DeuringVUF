"""DeuringVUF verification (paper Algorithm 3)."""

if __package__ in (None, ""):
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from src.config import context
from src.config.constants import REPRESENTATION_COFACTOR
from src.encoding.encode import canonical_basis, decode_proof
from src.protocol.hashing import f_in
from src.arithmetic.backend import CouplePoint, IsogVerif2dim
from src.arithmetic.projective import JacobianPoint, odd_linear_combination, valid_odd_basis
from src.encoding import serialization as wire


def _verification_setup(pk):
    """Validate public points and prepare the canonical proof-torsion basis."""
    try:
        curve, R, S, RmS, hint = pk
        if not valid_odd_basis(curve, (R, S, RmS), context.N):
            return None
        (U, V), _ = canonical_basis(curve, hint)
    except (AssertionError, ArithmeticError, AttributeError, TypeError, ValueError):
        return None
    return (REPRESENTATION_COFACTOR*JacobianPoint(U),
            REPRESENTATION_COFACTOR*JacobianPoint(V))


def _verify_proof(pk, pi, x, v):
    """Validate the public key and verify one proof."""
    basis = _verification_setup(pk)
    if basis is None:
        return False
    U, V = basis
    curve, R, S, RmS, _ = pk
    try:
        E3, U3, V3 = decode_proof(pi)
    except (AssertionError, ArithmeticError, AttributeError, TypeError, ValueError):
        return False
    r, s = f_in(x)
    kernel = odd_linear_combination(R, S, RmS, r, s, context.N)
    if kernel.is_zero():
        return False
    data = (curve, U, V, E3, REPRESENTATION_COFACTOR*JacobianPoint(U3),
            REPRESENTATION_COFACTOR*JacobianPoint(V3))
    result = IsogVerif2dim(data, v)
    if not result:
        return False
    phi, idx = result
    reference = R if r == 0 else S
    reference_image = phi(CouplePoint(reference, E3(0)), lift=False)[idx]
    kernel_image = phi(CouplePoint(kernel, E3(0)), lift=False)[idx]
    return not reference_image[1].is_zero() and kernel_image[1].is_zero()


class Verify:
    """Decode public data and verify, validating the key and preparing bases."""

    def __init__(self, pk, pi, x, v):
        self.pk, self.pi, self.x, self.v = pk, pi, x, v
        try:
            pk = wire.decode_public_key(pk) if wire.is_bytes(pk) else pk
            pi = wire.decode_proof(pi) if wire.is_bytes(pi) else pi
            v = wire.decode_output(v) if wire.is_bytes(v) else v
        except (AssertionError, ArithmeticError, AttributeError, TypeError, ValueError):
            self.verify = False
            return
        self.verify = _verify_proof(pk, pi, x, v)


if __name__ == "__main__":
    from examples.verify import main
    main()
