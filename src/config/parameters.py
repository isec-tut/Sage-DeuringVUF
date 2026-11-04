"""Scalar DeuringVUF parameters; importing this module needs no SageMath.

The protocol uses one fixed parameter set. Curves, quaternion objects and
torsion bases are constructed separately on first use.
"""

from dataclasses import dataclass


@dataclass(frozen=True)
class DeuringVUFParameters:
    """Fixed integers and the deterministic seed for public precomputations."""

    name: str = "original"
    security_bits: int = 128
    c: int = 305
    f: int = 240
    odd_torsion_cofactor: int = 1182
    alpha_scalar: int = 597586948685359749757390890506743252
    alpha_frobenius: int = 314328096788165707865159343442888735
    precomputation_seed: int = 0
    alpha_i: int = 0
    kappa_frobenius: int = 0
    secret_norm_strategy: str = "legacy"
    require_nonsquare_kappa: bool = False

    def __post_init__(self):
        if self.security_bits <= 0 or self.c <= 0 or self.c % 2 != 1 or self.f < 3:
            raise ValueError("security_bits and odd c must be positive; f must be at least 3")
        if self.odd_torsion_cofactor <= 0 or (self.p-1) % self.odd_torsion_cofactor:
            raise ValueError("odd_torsion_cofactor must divide p-1")
        if self.N <= 1 or self.N % 2 != 1:
            raise ValueError("N must be an odd integer greater than one")
        if not 3 <= self.representation_exponent <= self.f:
            raise ValueError("insufficient 2-power torsion for the VUF representation")
        if self.secret_norm_strategy not in ("legacy", "p-squared"):
            raise ValueError("unknown secret norm strategy")

    @property
    def p(self):
        return self.c*2**self.f-1

    @property
    def N(self):
        return (self.p-1)//self.odd_torsion_cofactor

    @property
    def torsion_order(self):
        return 2**self.f

    @property
    def representation_exponent(self):
        return self.N.bit_length()

    @property
    def representation_cofactor(self):
        return 2**(self.f-self.representation_exponent)

    @property
    def auxiliary_norm(self):
        return 2**self.representation_exponent-self.N

    @property
    def secret_ideal_bound(self):
        return self.p**2 if self.secret_norm_strategy == "p-squared" else 2**(4*self.security_bits)


DEFAULT_PARAMETERS = DeuringVUFParameters()

# Mathematical names retained for callers which only need scalar parameters.
p, f, c, N = (
    DEFAULT_PARAMETERS.p, DEFAULT_PARAMETERS.f, DEFAULT_PARAMETERS.c,
    DEFAULT_PARAMETERS.N,
)


def main():
    parameters = DEFAULT_PARAMETERS
    for name in (
        "name", "security_bits", "c", "f", "p", "odd_torsion_cofactor", "N",
        "representation_exponent", "representation_cofactor", "auxiliary_norm",
        "alpha_i", "kappa_frobenius", "secret_norm_strategy", "precomputation_seed",
    ):
        print(f"{name}: {getattr(parameters, name)}")


if __name__ == "__main__":
    main()
