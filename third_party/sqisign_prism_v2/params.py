"""DeuringVUF parameters, replacing the upstream signature parameter presets.

Uses the cached default context from src/config/context.py. Its first quaternion access
initializes the arithmetic objects. Upstream ec matrices act on rows, whereas
the local precomputations store columns.
"""

from sage.all import ZZ, Zmod, matrix, next_prime
from src.config import context as _setup

p = ZZ(_setup.p)
f = ZZ(_setup.f)
B, O0 = _setup.B, _setup.O0
QUAT_equiv_bound_coeff = 64
QUAT_prime_cofactor = next_prime(2**p.nbits())
Fp2 = _setup.Fp2
Fp = Fp2.prime_subfield()
Fp2_i = Fp2.gen()
E0 = _setup.Et_Fp4
P0 = _setup.iso_E0_Fp4_to_Et_Fp4(_setup.U0_Fp4)
Q0 = _setup.iso_E0_Fp4_to_Et_Fp4(_setup.V0_Fp4)
mat_1, mat_i, mat_ij2, mat_1k2 = (
    matrix(Zmod(2**f), M.transpose()) for M in _setup.Ms
)
