"""Local DeuringVUF ideal helpers for the SQIsign fork.

Keep construction independent of params: public precomputation uses it before
the default arithmetic context is ready.
"""


def O0_leftideal_create(alpha, N, O):
    """Construct and check the left ideal O<alpha,N> of norm N."""
    if N <= 0 or alpha not in O or alpha.reduced_norm() % N:
        raise ValueError('expected alpha in O with norm divisible by positive N')
    ideal = O.left_ideal([b*alpha for b in O.basis()] + [b*N for b in O.basis()])
    if ideal.left_order() != O or ideal.norm() != N:
        raise ArithmeticError('constructed ideal has incorrect order or norm')
    return ideal


def IdealGenerator(ideal):
    """Use SQIsign's reduced basis and SmallishGenerator for an odd-norm ideal."""
    from .quaternions import LLLReducedBasis, SmallishGenerator
    from sage.misc.randstate import seed

    norm = ideal.norm()
    if norm <= 0 or norm % 2 == 0:
        raise ValueError('expected a positive odd ideal norm')
    # The upstream sampler is randomized. Fix its local seed for reproducible
    # key encoding while preserving the caller's Sage RNG state.
    with seed(0):
        generator = SmallishGenerator(ideal, norm, LLLReducedBasis(ideal))
    if O0_leftideal_create(generator, norm, ideal.left_order()) != ideal:
        raise ArithmeticError('generator does not recover the supplied ideal')
    return generator
