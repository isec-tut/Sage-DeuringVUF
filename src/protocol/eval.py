"""Ordinary and amortized DeuringVUF evaluation (paper Algorithms 2 and 5)."""

if __package__ in (None, ""):
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import random

from src.config import context
from src.config.constants import AUXILIARY_NORM
from src.protocol.hashing import f_in
from third_party.sqisign_prism_v2.ideals import O0_leftideal_create
from src.arithmetic.backend import RandomIdealGivenNorm
from src.protocol.representation import build_representation
from src.arithmetic.projective import odd_linear_combination
from src.encoding import serialization as wire

ALPHA_CONJUGATE = context.alpha.conjugate()


def _input_ideal(pk, sk, x, *, compute_kernel=True):
    """Translate the input to its public kernel point and secret O0-ideal."""
    _, R, S, RmS, _ = pk
    mixing = sk[2]
    r, s = f_in(x)
    gamma = ((r*mixing[0, 0] + s*mixing[0, 1])
             + (r*mixing[1, 0] + s*mixing[1, 1])*context.kappa)
    ideal = O0_leftideal_create(
        ALPHA_CONJUGATE*gamma.conjugate(), context.N, context.O0
    )
    kernel = (odd_linear_combination(R, S, RmS, r, s, context.N)
              if compute_kernel else None)
    return ideal, kernel


class Eval:
    """Decode keys and evaluate once, storing byte-encoded v and pi plus setup."""

    def __init__(self, pk, sk, x):
        self.pk, self.sk, self.x = pk, sk, x
        pk = wire.decode_public_key(pk) if wire.is_bytes(pk) else pk
        sk = wire.decode_secret_key(sk, pk) if wire.is_bytes(sk) else sk
        _, ideal, _, C0, U, V = sk
        I_x, kernel = _input_ideal(pk, sk, x)
        # Ex, E3, U3, V3 = self.IsogenyRepresentation(
        Ex, E3, C3, hint3 = self.IsogenyRepresentation(
            ideal, I_x, U, V, C0, kernel
        )
        # self.pi = encode_proof(E3, U3, V3)
        self.pi = wire.encode_proof((E3, C3, hint3))
        self.v = wire.encode_output(Ex.j_invariant())
        self.amortized = (self.K, self.E_aux)

    def IsogenyRepresentation(self, I, I_x, U, V, C0, kernel):
        """Return (Ex, E3, C3, hint3) with matrix-encoded source and proof bases."""
        self.K = RandomIdealGivenNorm(AUXILIARY_NORM, prime=False)
        self.rng_state = random.getstate()
        curve, proof, self.E_aux = build_representation(
            I, I_x, U, V, C0, self.K, kernel=kernel
        )
        return curve, *proof


class AmortizedState(tuple):
    """Pair (I_aux, E_aux) with private, key-bound decoded-key caches.

    This is an in-memory state, not a wire format. Indexing and unpacking
    retain the historical two-element state interface.
    """

    def __new__(cls, I_aux, E_aux, pk, sk, decoded_pk, decoded_sk):
        instance = super().__new__(cls, (I_aux, E_aux))
        instance._pk = bytes(pk) if wire.is_bytes(pk) else pk
        instance._sk = bytes(sk) if wire.is_bytes(sk) else sk
        instance._decoded_pk = decoded_pk
        instance._decoded_sk = decoded_sk
        return instance


class Amortized_Eval:
    """Algorithm 1: initialize or reuse state (I_aux, E_aux) under one key.

    Pass ``state=None`` for the first evaluation and pass the returned ``state``
    to subsequent evaluations, including those on different inputs. The old
    positional ``amortized=(K, E_aux)`` argument remains supported.
    Optional ``rng_state`` restores Python's RNG before ideal conversion;
    it does not restore SageMath's RNG.
    The returned state caches decoded keys; keep the state object intact to
    skip key decoding on subsequent evaluations under the same key.
    """

    # def __init__(self, pk, sk, x, amortized=None, rng_state=None, *, state=None):
    def __init__(self, pk, sk, x, amortized=None, *, state=None, rng_state=None):
        self.pk, self.sk, self.x = pk, sk, x
        input_pk, input_sk = pk, sk
        if isinstance(state, AmortizedState):
            if pk != state._pk or sk != state._sk:
                raise ValueError('amortized state belongs to a different key pair')
            pk, sk = state._decoded_pk, state._decoded_sk
        else:
            pk = wire.decode_public_key(pk) if wire.is_bytes(pk) else pk
            sk = wire.decode_secret_key(sk, pk) if wire.is_bytes(sk) else sk
        _, ideal, _, C0, U, V = sk
        if state is not None and amortized is not None:
            raise ValueError("Pass either state or legacy amortized, not both")
        # Lines 8-11: sample K only once; store I_aux rather than K.
        if state is not None:
            I_aux, E_aux = state
        elif amortized is not None:
            K, E_aux = amortized
            I_aux = ideal.intersection(K)
        else:
            K = RandomIdealGivenNorm(AUXILIARY_NORM, prime=False)
            I_aux, E_aux = ideal.intersection(K), None
        I_x, kernel = _input_ideal(pk, sk, x, compute_kernel=E_aux is None)
        # Lines 12-36: first identify by the kernel, then by cached j(E_aux).
        Ex, proof, complementary_curve = build_representation(
            ideal, I_x, U, V, C0, None, kernel=kernel,
            # auxiliary_curve=E_aux, rng_state=rng_state, auxiliary_intersection=I_aux,
            auxiliary_curve=E_aux, auxiliary_intersection=I_aux, rng_state=rng_state,
        )
        # self.pi = encode_proof(*proof)
        self.pi = wire.encode_proof(proof)
        self.v = wire.encode_output(Ex.j_invariant())
        self.state = AmortizedState(
            I_aux, complementary_curve if E_aux is None else E_aux,
            input_pk, input_sk, pk, sk,
        )


if __name__ == "__main__":
    from benchmarks.benchmark_vuf import main
    main()
