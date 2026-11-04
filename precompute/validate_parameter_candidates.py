"""Build actual public curve/torsion data for search_parameters candidates.

This is isolated from the fixed production backend and does not rebind defaults.
"""

import argparse
import json
from pathlib import Path
from types import SimpleNamespace

from sage.all import QuaternionAlgebra, ZZ
from sage.misc.randstate import seed

from third_party.sqisign_prism_v2.ideals import O0_leftideal_create
from src.config.parameters import DeuringVUFParameters
from precompute.build_context import (
    _build_base_curves, _extend_curves, _build_odd_torsion,
    _build_even_torsion, points_order_N,
)
from precompute.search_parameters import validate_candidate


def validate_public_context(candidate):
    validate_candidate(candidate)
    params = DeuringVUFParameters(**candidate["parameters"])
    with seed(params.precomputation_seed):
        B = QuaternionAlgebra(-1, -params.p)
        i, j, k = B.gens()
        O0 = B.quaternion_order([1, i, (i+j)/2, (1+k)/2])
        alpha = params.alpha_scalar + candidate["alpha_i"]*i + params.alpha_frobenius*j
        kappa = i + candidate["kappa_frobenius"]*j
        ideal = O0_leftideal_create(alpha.conjugate(), params.N, O0)
        quaternion = SimpleNamespace(alpha=alpha, kappa=kappa, O0=O0)
        curves = _build_base_curves(params)
        R0 = points_order_N(curves.Et_fp2, params.N)
        curves = _extend_curves(params, curves)
        odd = _build_odd_torsion(params, quaternion, curves, R0)
        even = _build_even_torsion(params, curves)
        N, T = params.N, 2**params.f
        odd_pairing = odd.P0_fp4.weil_pairing(odd.Q0_fp4, N)
        even_pairing = even.U0_Fp4.weil_pairing(even.V0_Fp4, T)
        checks = {
            "alpha_in_maximal_order": alpha in O0,
            "kappa_in_maximal_order": kappa in O0,
            "kernel_ideal_norm_N": ideal.norm() == N,
            "odd_P_nonzero_order_N": not odd.P0_fp4.is_zero() and (N*odd.P0_fp4).is_zero(),
            "odd_Q_nonzero_order_N": not odd.Q0_fp4.is_zero() and (N*odd.Q0_fp4).is_zero(),
            "odd_basis_pairing_nontrivial": odd_pairing != 1 and odd_pairing**N == 1,
            "even_U_exact_order": (T*even.U0_Fp4).is_zero() and not ((T//2)*even.U0_Fp4).is_zero(),
            "even_V_exact_order": (T*even.V0_Fp4).is_zero() and not ((T//2)*even.V0_Fp4).is_zero(),
            "even_basis_pairing_primitive": even_pairing**T == 1 and even_pairing**(T//2) != 1,
            "discriminant_nonsquare": ZZ(kappa.reduced_trace()**2-4*kappa.reduced_norm()).kronecker(N) == -1,
        }
        if not all(checks.values()):
            raise ArithmeticError(f"public context failed: {checks}")
        return {name: bool(value) for name, value in checks.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--all", action="store_true", help="validate every candidate, rather than first per level")
    args = parser.parse_args()
    result = json.loads(args.input.read_text(encoding="utf-8"))
    seen = set()
    for candidate in result["candidates"]:
        if not args.all and candidate["level"] in seen:
            continue
        seen.add(candidate["level"])
        print(f"Building {candidate['parameters']['name']}", flush=True)
        candidate["public_context_checks"] = validate_public_context(candidate)
        args.input.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
        print("public context: all checks passed", flush=True)


if __name__ == "__main__":
    main()
