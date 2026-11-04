"""Wire-format regressions. Run: sage -python -m tests.test_serialization."""


import random
from unittest.mock import patch
from sage.all import set_random_seed
from sage.misc.randstate import current_randstate

from src.config import context
from src.encoding import serialization as wire
from src.encoding.encode import decode_proof as expand_proof, encode_proof as compress_proof
from src.protocol.eval import Eval, Amortized_Eval
from src.protocol.keygen import Keygen
from src.protocol.verify import Verify


def rejects(function, *args):
    try:
        function(*args)
    except (TypeError, ValueError, ArithmeticError):
        return
    raise AssertionError(f'{function.__name__} accepted malformed data')


def main():
    random.seed(0)
    set_random_seed(0)
    # Known byte vectors fix coefficient order, endianness and the i embedding.
    value = context.Fp4(258) + 3*context.Fp4(context.zeta4)
    expected = (258).to_bytes(wire.FP_BYTES, 'little') + (3).to_bytes(wire.FP_BYTES, 'little')
    assert wire.encode_fp2(value) == expected
    assert wire.decode_fp2(expected) == value
    for value in (0, 1, context.p-1, context.Fp4(context.zeta4)):
        assert wire.decode_fp2(wire.encode_fp2(value)) == value
    rejects(wire.encode_fp2, context.Fp4.gen())
    rejects(wire.decode_fp2, int(context.p).to_bytes(wire.FP_BYTES, 'little') + bytes(wire.FP_BYTES))

    key = Keygen()
    assert isinstance(key.pk, bytes) and isinstance(key.sk, bytes)
    assert len(key.pk) == wire.PUBLIC_KEY_BYTES
    assert len(key.sk) == wire.SECRET_KEY_BYTES
    pk = wire.decode_public_key(key.pk)
    sk = wire.decode_secret_key(key.sk, key.pk)
    assert wire.encode_public_key(pk) == key.pk
    python_state = random.getstate()
    sage_state = current_randstate()
    assert wire.encode_secret_key(pk, sk) == key.sk
    assert random.getstate() == python_state
    assert current_randstate() is sage_state
    # Independent projective rescaling must not affect the affine x encoding.
    E, *rest = pk
    R, S, RmS, hint = rest
    scaled = tuple(P.parent()(tuple(c*t for c in P.XZ()))
                   for P, t in zip((R, S, RmS), (2, 3, 5)))
    assert wire.encode_public_key((E, *scaled, hint)) == key.pk
    assert wire.PUBLIC_KEY_BYTES == 194
    assert wire.SECRET_KEY_BYTES == 594
    # Both relative signs must reconstruct the corresponding differential point.
    Ra, Sa = R.y_recovery(S, RmS)
    for branch in (0, 1):
        encoded = key.pk[:-1] + bytes([branch])
        with (
            patch.object(wire.KummerPoint, 'lift', side_effect=AssertionError('lifting')),
            patch.object(wire.KummerPoint, 'y_recovery', side_effect=AssertionError('y recovery')),
        ):
            decoded = wire.decode_public_key(encoded)
            assert wire.encode_public_key(decoded) == encoded
        expected_point = Ra-Sa if branch == key.pk[-1] else Ra+Sa
        assert decoded[3] == R.parent()(expected_point)
    rejects(wire.decode_public_key, key.pk[:-1]+b'\2')
    rejects(wire.decode_public_key, key.pk[:-1]+b'\xff')
    offset = 2*wire.FP_BYTES+1
    width = 2*wire.FP_BYTES
    repeated = (key.pk[:offset+width] + key.pk[offset:offset+width] + key.pk[-1:])
    rejects(wire.decode_public_key, repeated)
    repeated_root = key.pk[:offset] + bytes(width) + key.pk[offset+width:]
    rejects(wire.decode_public_key, repeated_root)
    rejects(wire.encode_public_key, (E, R, S, R, hint))
    state = None
    for x in (b'\0\0', b'wire regression'):
        evaluation = Eval(key.pk, key.sk, x)
        cached = Amortized_Eval(key.pk, key.sk, x, evaluation.amortized,
                                rng_state=evaluation.rng_state)
        if state is None:
            amortized = Amortized_Eval(key.pk, key.sk, x)
        else:
            # Reuse must not decode keys or regenerate their canonical basis.
            with (
                patch.object(wire, 'decode_public_key', side_effect=AssertionError('public key decoded again')),
                patch.object(wire, 'decode_secret_key', side_effect=AssertionError('secret key decoded again')),
            ):
                amortized = Amortized_Eval(key.pk, key.sk, x, state=state)
        state = amortized.state
        for result in (evaluation, cached, amortized):
            assert isinstance(result.pi, bytes) and isinstance(result.v, bytes)
            assert len(result.pi) == wire.PROOF_BYTES
            assert len(result.v) == wire.OUTPUT_BYTES
            assert wire.decode_output(result.v) == wire.decode_output(evaluation.v)
            proof = wire.decode_proof(result.pi)
            assert wire.encode_proof(proof) == result.pi
            assert compress_proof(*expand_proof(proof)) == proof
            assert Verify(key.pk, result.pi, x, result.v).verify
        # Raw decoded objects remain accepted by the arithmetic compatibility API.
        assert Verify(pk, wire.decode_proof(evaluation.pi), x,
                      wire.decode_output(evaluation.v)).verify
        for bad in (key.pk[:-1], key.pk+b'\0'):
            assert not Verify(bad, evaluation.pi, x, evaluation.v).verify
        for bad in (evaluation.pi[:-1], evaluation.pi+b'\0'):
            assert not Verify(key.pk, bad, x, evaluation.v).verify
        for bad in (evaluation.v[:-1], evaluation.v+b'\0'):
            assert not Verify(key.pk, evaluation.pi, x, bad).verify
        zero_matrix = (evaluation.pi[:2*wire.FP_BYTES]
                       + bytes(4*wire.TORSION_BYTES) + evaluation.pi[-1:])
        assert not Verify(key.pk, zero_matrix, x, evaluation.v).verify
        assert not Verify(key.pk, evaluation.pi, x+b'changed', evaluation.v).verify
    for decode, data in ((wire.decode_public_key, key.pk),
                         (wire.decode_secret_key, key.sk),
                         (wire.decode_proof, evaluation.pi),
                         (wire.decode_output, evaluation.v)):
        rejects(decode, data[:-1])
        rejects(decode, data+b'\0')
    bad_pk = bytearray(key.pk)
    bad_pk[0] ^= 1
    rejects(wire.decode_secret_key, key.sk, bad_pk)
    rejects(lambda: Amortized_Eval(bytes(bad_pk), key.sk, b'wrong key', state=state))
    rejects(lambda: Amortized_Eval(key.pk, key.sk+b'\0', b'wrong key', state=state))
    bad_sk = bytearray(key.sk)
    bad_sk[wire.PUBLIC_KEY_BYTES:wire.PUBLIC_KEY_BYTES+wire.FP_BYTES] = bytes(wire.FP_BYTES)
    rejects(wire.decode_secret_key, bad_sk)
    # A=0 has a unique zero hint encoding and implicit deterministic basis.
    from src.encoding.encode import canonical_basis
    basis, hint = canonical_basis(context.E0_Fp4)
    proof = compress_proof(context.E0_Fp4, *basis)
    data = wire.encode_proof(proof)
    assert wire.decode_proof(data) == proof
    rejects(wire.decode_proof, data[:-1]+b'\1')
    print(f'serialization ({context.DEFAULT_PARAMETERS.name}): all tests passed', flush=True)


if __name__ == '__main__':
    main()
