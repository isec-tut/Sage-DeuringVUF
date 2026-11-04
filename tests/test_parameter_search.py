"""Mathematical regression checks for the parameter search (SageMath)."""

from math import gcd
from sage.all import GF, QuaternionAlgebra

from precompute.search_parameters import find_alpha, find_general_alpha, find_kappa


def main():
    # p=31, N=5: the old a+b*j form exists, but kappa=i has square discriminant.
    p, N = 31, 5
    a, b = find_alpha(p, N)
    assert b % N and gcd(a*a+p*b*b, N*N) == N
    assert GF(N)(-4).is_square()
    t = find_kappa(p, N)
    assert gcd(1+p*t*t, N) == 1
    assert not GF(N)(-4*(1+p*t*t)).is_square()

    # p=71, N=7: a+b*j is unavailable, while a+u*i+j retains kappa=i.
    p, N = 71, 7
    assert not GF(N)(-p).is_square()
    a, u, b = find_general_alpha(p, N)
    B = QuaternionAlgebra(-1, -p)
    i, j, k = B.gens()
    alpha, kappa = a+u*i+b*j, i
    assert alpha.reduced_norm() == a*a+u*u+p*b*b
    assert gcd(int(alpha.reduced_norm()), N*N) == N
    assert not GF(N)(kappa.reduced_trace()**2-4*kappa.reduced_norm()).is_square()
    print("parameter search: both alpha/kappa families passed", flush=True)


if __name__ == "__main__":
    main()
