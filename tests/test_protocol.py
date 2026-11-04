"""End-to-end DeuringVUF evaluation, verification, and basis encodings.

Run: sage -python -m tests.test_protocol
"""

import random
from unittest.mock import patch

from sage.all import ZZ, Zmod, identity_matrix, matrix, set_random_seed
from src.config import context
from src.protocol.keygen import Keygen
from src.protocol.eval import Eval, Amortized_Eval
from src.protocol.verify import Verify
from src.encoding import serialization as wire
from src.encoding.encode import (
    canonical_basis, decode_proof, decode_secret_key, encode_proof, standard_basis,
)
from src.arithmetic.backend import qt
from src.arithmetic.kummer import KummerPoint
from src.arithmetic.projective import odd_linear_combination, valid_odd_basis
from src.arithmetic.torsion import ChangeOfBasis, EvalMatrix, PointCoords, point_coords_pairing


def test_protocol(seed):
    random.seed(seed)
    set_random_seed(seed)
    key = Keygen()
    key.sk = wire.decode_secret_key(key.sk, key.pk)
    key.pk = wire.decode_public_key(key.pk)
    E, R, S, RmS, hint = key.pk
    assert all(isinstance(P, KummerPoint) for P in (R, S, RmS))
    assert valid_odd_basis(E, (R, S, RmS), context.N)
    with (
        patch.object(KummerPoint, 'x', side_effect=AssertionError('affine normalization')),
        patch.object(KummerPoint, 'lift', side_effect=AssertionError('lifting')),
        patch.object(KummerPoint, 'y_recovery', side_effect=AssertionError('y recovery')),
    ):
        assert valid_odd_basis(E, (R, S, RmS), context.N)
        projective_kernel = odd_linear_combination(R, S, RmS, 0, 1, context.N)
        assert projective_kernel == S
    Ra, Sa = R.y_recovery(S, RmS)
    assert key.sk[2].det() % context.N != 0
    odd_pairing = context.P0_fp4.weil_pairing(context.Q0_fp4, context.N)
    assert Ra.weil_pairing(Sa, context.N) == odd_pairing**ZZ(key.sk[1].norm()*key.sk[2].det())
    reduced = qt.LLLReducedBasis(key.sk[1])
    assert all(b in key.sk[1] for b in reduced)
    expanded = decode_secret_key(key.pk, key.sk)
    canonical, _ = canonical_basis(E, hint)
    assert len(key.sk) == 6
    assert key.sk[3].base_ring() == Zmod(2**context.f)
    assert key.sk[4:] == canonical
    assert expanded[3:5] == EvalMatrix(key.sk[3], standard_basis())
    assert ChangeOfBasis(None, expanded[3:5]) == key.sk[3]
    assert expanded[-2:] == canonical
    _check_pointcoords_basis(*canonical)
    initial = Amortized_Eval(key.pk, key.sk, b'initial amortized evaluation')
    initial_eval = Eval(key.pk, key.sk, initial.x)
    # v is the j-invariant: compare outputs for the same key and input.
    assert initial.v == initial_eval.v, (
        seed, initial.x, 'initial', initial_eval.v, initial.v
    )
    assert Verify(key.pk, initial.pi, initial.x, initial.v).verify
    state = initial.state
    for x in (b'\0\0', b'\0\1', b'protocol regression'):
        evaluation = Eval(key.pk, key.sk, x)
        check = Verify(key.pk, evaluation.pi, x, evaluation.v)
        assert check.verify, (seed, x)
        if x == b'\0\0':
            line = R.parent()
            scaled = tuple(line((scale*X, scale*Z))
                           for P, scale in zip((R, S, RmS), (2, 3, 5))
                           for X, Z in (P.XZ(),))
            scaled_pk = (E, *scaled, hint)
            assert Verify(scaled_pk, evaluation.pi, x, evaluation.v).verify
            with patch('src.protocol.verify.f_in', return_value=(0, 0)):
                assert not Verify(key.pk, evaluation.pi, x, evaluation.v).verify
            if seed == 0:
                # Exercise the two verifier reference choices, including r=0.
                for scalars in ((0, 1), (1, 0)):
                    with (patch('src.protocol.eval.f_in', return_value=scalars),
                          patch('src.protocol.verify.f_in', return_value=scalars)):
                        edge = Eval(key.pk, key.sk, x)
                        assert Verify(key.pk, edge.pi, x, edge.v).verify
        cached = Amortized_Eval(key.pk, key.sk, x, evaluation.amortized)
        assert cached.v == evaluation.v, (seed, x, 'cached', evaluation.v, cached.v)
        assert Verify(key.pk, cached.pi, x, cached.v).verify
        with (
            patch('src.protocol.eval.RandomIdealGivenNorm', side_effect=AssertionError('resampling K')),
            patch('src.protocol.eval.odd_linear_combination', side_effect=AssertionError('computing Rx')),
        ):
            reused = Amortized_Eval(key.pk, key.sk, x, state=state)
        assert reused.state[0] is state[0]
        assert reused.state[1] is state[1]
        assert reused.v == evaluation.v, (seed, x, 'reused', evaluation.v, reused.v)
        assert Verify(key.pk, reused.pi, x, reused.v).verify
        state = reused.state
        assert not Verify(key.pk, evaluation.pi, x+b'changed', evaluation.v).verify
        assert not Verify(key.pk, evaluation.pi, x, wire.decode_output(evaluation.v)+1).verify
        proof = wire.decode_proof(evaluation.pi)
        E3, C3, hint3 = proof
        _, U3, V3 = decode_proof(proof)
        assert encode_proof(E3, U3, V3) == proof
        assert C3.base_ring() == Zmod(2**context.f) and ZZ(C3.det()) % 2 == 1
        for bad_pi in (
            (E3, matrix(Zmod(2**context.f), 2, 2), hint3),
            (E3, C3, (2, 1)),
            (E3, C3, (0, -1)),
            (E3, C3, (0, 128)),
            (E3, matrix(ZZ, 3, 3), hint3),
            (E3, matrix(Zmod(17), 2, 2, [ZZ(c) for c in C3.list()]), hint3),
            (E3, matrix(Zmod(2**context.f), [[3, 0], [0, 1]])*C3, hint3),
            (E3, E3(0), V3),  # Historical point tuple is no longer a proof encoding.
        ):
            assert not Verify(key.pk, bad_pi, x, evaluation.v).verify
        zero, invalid = R.parent().zero(), R.parent()((0, 0))
        for bad_pk in (
            (E, zero, S, RmS, hint), (E, R, zero, RmS, hint),
            (E, R, S, zero, hint), (E, R, S, R, hint),
            (E, invalid, S, RmS, hint), (E, R, S, invalid, hint),
            (E, E(0), S, RmS, hint),
        ):
            rejected = Verify(bad_pk, evaluation.pi, x, evaluation.v)
            assert not rejected.verify
    print(f'protocol seed {seed}: passed', flush=True)



def _check_pointcoords_basis(P, Q):
    """Exercise full/reduced torsion and a nontrivial change of basis."""
    R, S = 3*P+2*Q, 5*P+7*Q
    for e in (1, 4, context.f):
        n = ZZ(2)**e
        pairing = point_coords_pairing(R, S, e)
        for a, b in ((0, 0), (-17, n//2+7), (n-1, n-2)):
            X = 2**(context.f-e)*(a*R+b*S)
            assert PointCoords(X, R, S, pairing, e=e) == (a % n, b % n)
        expected = matrix(Zmod(n), [[3, 2], [5, 7]])
        reduced = tuple(2**(context.f-e)*X for X in (R, S))
        assert ChangeOfBasis((P, Q), reduced, e=e) == expected
        assert tuple(2**(context.f-e)*X for X in EvalMatrix(expected, (P, Q))) == reduced



def test_canonical_zero_coefficient():
    basis, hint = canonical_basis(context.E0_Fp4)
    assert hint is None
    assert canonical_basis(context.E0_Fp4) == (basis, hint)
    proof = encode_proof(context.E0_Fp4, *basis)
    assert proof[1] == identity_matrix(Zmod(2**context.f), 2)
    assert decode_proof(proof) == (context.E0_Fp4, *basis)


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--seeds', type=int, nargs='+', default=[0, 1, 7])
    args = parser.parse_args()
    print(f'params: {context.DEFAULT_PARAMETERS.name}', flush=True)
    test_canonical_zero_coefficient()
    for seed in args.seeds:
        test_protocol(seed)
    print('protocol: all tests passed', flush=True)


if __name__ == '__main__':
    main()
