"""DeuringVUF wire format, following SQIsign v2 section 4.6.

Parameters are agreed out of band. See docs/serialization.md for the exact
layout and the DeuringVUF-specific extensions. No pickle or Sage payloads.
"""

from sage.all import ZZ, Zmod, EllipticCurve, gcd, matrix

from src.config import context
from src.encoding.encode import basis_matrix, canonical_basis
from src.arithmetic.kummer import KummerLine, KummerPoint

FP_BYTES = (int(context.p).bit_length() + 7) // 8
TORSION_BYTES = (context.f + 7) // 8
ODD_BYTES = (int(context.N).bit_length() + 7) // 8
PUBLIC_KEY_BYTES = 6 * FP_BYTES + 2
PROOF_BYTES = 2 * FP_BYTES + 4 * TORSION_BYTES + 1
OUTPUT_BYTES = 2 * FP_BYTES
SECRET_KEY_BYTES = PUBLIC_KEY_BYTES + 5 * FP_BYTES + 4 * ODD_BYTES + 4 * TORSION_BYTES


class _Reader:
    def __init__(self, data, size):
        if not isinstance(data, (bytes, bytearray, memoryview)):
            raise TypeError('expected a byte string')
        self.data = bytes(data)
        if len(self.data) != size:
            raise ValueError(f'expected exactly {size} bytes')
        self.offset = 0

    def take(self, size):
        part = self.data[self.offset:self.offset + size]
        if len(part) != size:
            raise ValueError('truncated encoding')
        self.offset += size
        return part

    def integer(self, size, bound=None, signed=False):
        value = ZZ(int.from_bytes(self.take(size), 'little', signed=signed))
        if bound is not None and not 0 <= value < bound:
            raise ValueError('noncanonical integer')
        return value


def _integer(value, size, bound=None, signed=False):
    value = ZZ(value)
    if bound is not None and not 0 <= value < bound:
        raise ValueError('integer out of range')
    try:
        return int(value).to_bytes(size, 'little', signed=signed)
    except OverflowError as exc:
        raise ValueError('integer does not fit the wire format') from exc


def encode_fp2(value):
    """Encode a+bi as LE(a)||LE(b), using the configured i^2=-1."""
    value = context.Fp4(value)
    if value**context.p**2 != value:
        raise ValueError('value is outside Fp2')
    imag_unit = context.Fp4(context.zeta4)
    real = ZZ((value + value**context.p) / 2)
    imag = ZZ((value - value**context.p) / (2 * imag_unit))
    return _integer(real, FP_BYTES, context.p) + _integer(imag, FP_BYTES, context.p)


def _read_fp2(reader):
    real = reader.integer(FP_BYTES, context.p)
    imag = reader.integer(FP_BYTES, context.p)
    return context.Fp4(real) + imag * context.Fp4(context.zeta4)


def decode_fp2(data):
    return _read_fp2(_Reader(data, OUTPUT_BYTES))


def _curve_bytes(curve):
    if curve.base_ring() != context.Fp4 or tuple(curve.a_invariants()) != (
        0, curve.a2(), 0, 1, 0
    ):
        raise ValueError('expected a configured Montgomery curve')
    return encode_fp2(curve.a2())


def _read_curve(reader):
    return EllipticCurve(context.Fp4, [0, _read_fp2(reader), 0, 1, 0])


def _hint_bytes(curve, hint):
    if curve.a2() == 0:
        if hint is not None:
            raise ValueError('A=0 must use the implicit basis')
        return b'\0'
    hA, h = map(ZZ, hint)
    if hA not in (0, 1) or not 0 <= h < 128:
        raise ValueError('invalid basis hint')
    return bytes([int(hA + 2*h)])


def _read_hint(reader, curve):
    value = int(reader.integer(1))
    if curve.a2() == 0:
        if value:
            raise ValueError('noncanonical A=0 hint')
        return None
    return value & 1, value >> 1


def _matrix_bytes(value, width, modulus):
    if value.nrows() != 2 or value.ncols() != 2:
        raise ValueError('expected a 2 by 2 matrix')
    if gcd(ZZ(value.det()), modulus) != 1:
        raise ValueError('noninvertible basis matrix')
    return b''.join(_integer(c, width, modulus) for c in value.list())


def _read_matrix(reader, width, modulus):
    result = matrix(Zmod(modulus), 2, 2,
                    [reader.integer(width, modulus) for _ in range(4)])
    if gcd(ZZ(result.det()), modulus) != 1:
        raise ValueError('noninvertible basis matrix')
    return result


def _difference_candidates(curve, r, s):
    """Return the two affine differential roots in a fixed wire order.

    The square root is taken in Fp4, then checked to lie in Fp2. This also
    supports points on the quadratic twist without recovering either y.
    """
    denominator = (r-s)**2
    if not denominator:
        raise ValueError('public basis must have distinct x coordinates')
    middle = (r*s+1)*(r+s)+2*curve.a2()*r*s
    discriminant = middle**2-(r*s-1)**2*denominator
    if not discriminant:
        raise ValueError('public basis has a repeated differential root')
    try:
        root = discriminant.sqrt(extend=False)
    except (ValueError, ArithmeticError) as exc:
        raise ValueError('invalid public differential roots') from exc
    # Fix the sign using the wire embedding rather than Sage's root ordering.
    if encode_fp2(-root) < encode_fp2(root):
        root = -root
    return (middle+root)/denominator, (middle-root)/denominator


def encode_public_key(pk):
    curve, R, S, RmS, hint = pk
    encoded = _curve_bytes(curve) + _hint_bytes(curve, hint)
    coordinates = []
    for point in (R, S, RmS):
        if point.parent().curve() != curve:
            raise ValueError('public point has the wrong curve')
        X, Z = point.XZ()
        if not Z:
            raise ValueError('cannot encode a public point at infinity')
        x = X / Z
        encode_fp2(x)  # Require every input coordinate to lie in Fp2.
        coordinates.append(x)
    r, s, difference = coordinates
    candidates = _difference_candidates(curve, r, s)
    if difference not in candidates:
        raise ValueError('inconsistent public differential triple')
    branch = candidates.index(difference)
    return encoded + encode_fp2(r) + encode_fp2(s) + bytes([branch])


def decode_public_key(data):
    reader = _Reader(data, PUBLIC_KEY_BYTES)
    curve = _read_curve(reader)
    hint = _read_hint(reader, curve)
    line = KummerLine(curve)
    r, s = _read_fp2(reader), _read_fp2(reader)
    branch = int(reader.integer(1, bound=2))
    difference = _difference_candidates(curve, r, s)[branch]
    points = tuple(KummerPoint(line, (x, 1)) for x in (r, s, difference))
    return curve, *points, hint


def _ideal_bytes(ideal):
    from third_party.sqisign_prism_v2.ideals import IdealGenerator
    if ideal.left_order() != context.O0:
        raise ValueError('expected a left O0 ideal')
    norm = ZZ(ideal.norm())
    if norm <= 0 or norm % 2 == 0:
        raise ValueError('expected a positive odd ideal norm')
    generator = IdealGenerator(ideal)
    # SQIsign section 4.6 footnote 1: clear the invertible denominator.
    denominator = ZZ(generator.denominator())
    if gcd(denominator, norm) != 1:
        raise ValueError('ideal generator denominator is not invertible')
    generator *= denominator
    # Center coefficients modulo norm to keep the fixed signed width.
    coefficients = [((ZZ(c) + norm//2) % norm) - norm//2 for c in generator]
    generator = context.B(coefficients)
    recovered = context.O0.left_ideal(
        [b*generator for b in context.O0.basis()] + [b*norm for b in context.O0.basis()]
    )
    if recovered != ideal:
        raise ValueError('generator does not recover the secret ideal')
    return _integer(norm, FP_BYTES, signed=True) + b''.join(
        _integer(c, FP_BYTES, signed=True) for c in coefficients
    )


def _read_ideal(reader):
    norm = reader.integer(FP_BYTES, signed=True)
    if norm <= 0 or norm % 2 == 0:
        raise ValueError('invalid ideal norm')
    coefficients = [reader.integer(FP_BYTES, signed=True) for _ in range(4)]
    if any(not -(norm//2) <= c <= norm//2 for c in coefficients):
        raise ValueError('noncanonical ideal generator')
    generator = context.B(coefficients)
    ideal = context.O0.left_ideal(
        [b*generator for b in context.O0.basis()] + [b*norm for b in context.O0.basis()]
    )
    if ideal.norm() != norm:
        raise ValueError('incorrect encoded ideal norm')
    return ideal


def encode_secret_key(pk, sk):
    curve, ideal, mixing, C0, U, V = sk
    if curve != pk[0] or (U, V) != tuple(canonical_basis(curve, pk[-1])[0]):
        raise ValueError('secret key does not match the public basis')
    return (encode_public_key(pk) + _ideal_bytes(ideal)
            + _matrix_bytes(mixing, ODD_BYTES, context.N)
            + _matrix_bytes(basis_matrix(C0), TORSION_BYTES, 2**context.f))


def decode_secret_key(data, pk=None):
    reader = _Reader(data, SECRET_KEY_BYTES)
    embedded_pk = reader.take(PUBLIC_KEY_BYTES)
    if pk is not None:
        expected = bytes(pk) if isinstance(pk, (bytes, bytearray, memoryview)) else encode_public_key(pk)
        if embedded_pk != expected:
            raise ValueError('secret and public keys differ')
    public = decode_public_key(embedded_pk)
    curve = public[0]
    ideal = _read_ideal(reader)
    mixing = matrix(ZZ, _read_matrix(reader, ODD_BYTES, context.N))
    C0 = _read_matrix(reader, TORSION_BYTES, 2**context.f)
    (U, V), _ = canonical_basis(curve, public[-1])
    return curve, ideal, mixing, C0, U, V


def encode_proof(proof):
    curve, C3, hint = proof
    return (_curve_bytes(curve)
            + _matrix_bytes(basis_matrix(C3), TORSION_BYTES, 2**context.f)
            + _hint_bytes(curve, hint))


def decode_proof(data):
    reader = _Reader(data, PROOF_BYTES)
    curve = _read_curve(reader)
    C3 = _read_matrix(reader, TORSION_BYTES, 2**context.f)
    return curve, C3, _read_hint(reader, curve)


encode_output = encode_fp2
decode_output = decode_fp2


def is_bytes(value):
    return isinstance(value, (bytes, bytearray, memoryview))
