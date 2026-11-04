"""Ideal construction and the SQIsign-backed serialization generator.

Run: sage -python -m tests.test_ideals
"""

import random

from sage.all import ZZ, gcd
from sage.misc.randstate import current_randstate, seed

from src.config import context
from third_party.sqisign_prism_v2.ideals import IdealGenerator, O0_leftideal_create


def rejects(function, error, *args):
    try:
        function(*args)
    except error:
        return
    raise AssertionError(f'{function.__name__} accepted invalid input')


def test_generators():
    order = context.O0
    i = context.B.gen(0)
    for alpha, norm in ((2+i, 5), (3+2*i, 13), ((2+i)*(3+2*i), 65),
                        (context.alpha.conjugate(), context.N)):
        ideal = O0_leftideal_create(alpha, ZZ(norm), order)
        assert ideal.norm() == norm and ideal.left_order() == order
        assert alpha in ideal
        with seed(17):
            python_state = random.getstate()
            sage_state = current_randstate()
            generator = IdealGenerator(ideal)
            assert current_randstate() is sage_state
            assert random.getstate() == python_state
        with seed(91):
            assert IdealGenerator(ideal) == generator
        assert generator in ideal
        assert gcd(generator.reduced_norm(), norm**2) == norm
        # Independent reconstruction through order multiplication and addition.
        assert order*generator + order*norm == ideal
        basis = ideal.basis()
        equivalent_basis = (basis[0]+basis[1], basis[1], -basis[2], basis[3])
        same_ideal = order.left_ideal(equivalent_basis)
        assert same_ideal == ideal
        assert IdealGenerator(same_ideal) == generator


def test_invalid_construction():
    order, i = context.O0, context.B.gen(0)
    for alpha, norm in ((2+i, 0), (2+i, -5), (context.B(1), 5), (i/3, 5)):
        rejects(O0_leftideal_create, ValueError, alpha, norm, order)
    # Divisibility alone does not imply that O<alpha,N> has norm N.
    rejects(O0_leftideal_create, ArithmeticError, context.B(5), 5, order)
    rejects(IdealGenerator, ValueError, order*context.B(2))


def main():
    test_generators()
    test_invalid_construction()
    print('ideal construction, generator reconstruction, and RNG preservation: passed')


if __name__ == '__main__':
    main()
