"""Search DeuringVUF parameters at SQIsign-v3 field sizes (SageMath).

Run: sage -python -m precompute.search_parameters --output parameter_search.json
This writes research candidates, without changing DEFAULT_PARAMETERS.
"""

import argparse
from dataclasses import asdict
import json
from math import gcd
from pathlib import Path

from sage.all import GF, ZZ, matrix, next_prime

from src.config.parameters import DeuringVUFParameters


# Exact sizes computed from the primes in SQIsign v3, Section 5.2.
TARGETS = {"I": (128, 3, 324), "III": (192, 27, 500), "V": (256, 17, 664)}


def find_alpha(p, N):
    """Find small a,b with v_N(a^2+p*b^2)=1 and b invertible mod N."""
    root = ZZ(GF(N)(-p).sqrt())
    lattice = matrix(ZZ, [[N, 0], [root, 1]]).LLL()
    # A primitive short vector normally works; include a second vector and sums.
    vectors = [lattice[0], lattice[1], lattice[0]+lattice[1], lattice[0]-lattice[1]]
    for a, b in vectors:
        a, b = int(a), int(b)
        norm = a*a + p*b*b
        if b % N and gcd(norm, N*N) == N:
            return a, b
    raise ArithmeticError("could not construct a rank-one alpha modulo N")


def find_kappa(p, N):
    """kappa=i+t*j: norm=1+p*t^2, discriminant=-4*(1+p*t^2)."""
    field = GF(N)
    for t in range(1, 10000):
        norm = 1+p*t*t
        discriminant = -4*norm
        if norm % N and not field(discriminant).is_square():
            return t
    raise ArithmeticError("could not find a nonsquare kappa discriminant")


def find_general_alpha(p, N):
    """For N=3 mod 4, retain kappa=i and use alpha=a+u*i+j."""
    field = GF(N)
    for u in range(1, 10000):
        value = field(-p-u*u)
        if not value or not value.is_square():
            continue
        a = int(value.sqrt())
        a = min(a, N-a)
        if gcd(a*a+u*u+p, N*N) != N:
            a += N
        if gcd(a*a+u*u+p, N*N) == N:
            return a, u, 1
    raise ArithmeticError("could not construct general alpha")


def validate_candidate(candidate):
    """Prove scalar conditions; do not confuse these with a security proof."""
    params = DeuringVUFParameters(**candidate["parameters"])
    p, N, f = params.p, params.N, params.f
    a, b = params.alpha_scalar, params.alpha_frobenius
    u = candidate["alpha_i"]
    t = candidate["kappa_frobenius"]
    checks = {
        "p_prime_proven": bool(ZZ(p).is_prime(proof=True)),
        "N_prime_proven": bool(ZZ(N).is_prime(proof=True)),
        "p_mod_4": p % 4 == 3,
        "N_divides_half_p_minus_1": ((p-1)//2) % N == 0,
        "N_less_than_2f": N < 2**f,
        "N_at_least_2lambda_bits": N.bit_length() >= 2*params.security_bits,
        "field_at_least_reference_prime": p >= candidate["reference_prime"],
        "alpha_norm_valuation_one": gcd(a*a+u*u+p*b*b, N*N) == N,
        "alpha_frobenius_invertible": gcd(b, N) == 1,
        "kappa_norm_coprime_N": gcd(1+p*t*t, N) == 1,
        "kappa_discriminant_nonsquare": not GF(N)(-4*(1+p*t*t)).is_square(),
        "auxiliary_positive_odd_coprime_N": (
            params.auxiliary_norm > 0 and params.auxiliary_norm % 2 == 1
            and gcd(params.auxiliary_norm, N) == 1
        ),
    }
    if not all(checks.values()):
        raise ArithmeticError(f"candidate failed checks: {checks}")
    return checks


def make_candidate(level, c, f, cofactor, origin):
    security, ref_c, ref_f = TARGETS[level]
    p = c*2**f-1
    N = (p-1)//cofactor
    if N % 4 == 1:
        a, b = find_alpha(p, N)
        u, t = 0, find_kappa(p, N)
        adaptation = "kappa differs from fixed i; adapt odd-action basis conventions"
    else:
        a, u, b = find_general_alpha(p, N)
        t = 0
        adaptation = "alpha has an i component; generalize alpha representation and odd-action optimization"
    params = DeuringVUFParameters(
        name=f"deuringvuf-v3-size-{level.lower()}-c{c}-f{f}",
        security_bits=security, c=c, f=f, odd_torsion_cofactor=cofactor,
        alpha_scalar=a, alpha_frobenius=b,
        alpha_i=u, kappa_frobenius=t, secret_norm_strategy="p-squared",
        require_nonsquare_kappa=True,
    )
    result = {
        "level": level, "origin": origin,
        "reference_prime": ref_c*2**ref_f-1,
        "parameters": asdict(params),
        "p": p, "p_bits": p.bit_length(), "N": N, "N_bits": N.bit_length(),
        "representation_exponent": params.representation_exponent,
        "auxiliary_norm": params.auxiliary_norm,
        "alpha_i": u,
        "kappa_frobenius": t,
        "kappa": f"i + {t}*j",
        "recommended_secret_ideal_norm": int(next_prime(ZZ(p)**2, proof=True)),
        "production_backend_compatible": False,
        "note": adaptation + "; scalar parameters alone are not a drop-in preset; no security certification",
    }
    result["checks"] = validate_candidate(result)
    return result


def search_prime(level, c, f, max_cofactor, origin):
    security, ref_c, ref_f = TARGETS[level]
    p = c*2**f-1
    if p < ref_c*2**ref_f-1 or not ZZ(p).is_pseudoprime():
        return None
    # p-1 has exactly one factor of two. Search both residues N=1,3 mod 4.
    for cofactor in range(2, max_cofactor+1, 4):
        if (p-1) % cofactor:
            continue
        N = (p-1)//cofactor
        if N >= 2**f or N.bit_length() < 2*security:
            continue
        if ZZ(N).is_pseudoprime():
            return make_candidate(level, c, f, cofactor, origin)
    return None


def search(level, count, max_c, max_cofactor, extra_bits=0):
    _, ref_c, ref_f = TARGETS[level]
    reference = ref_c*2**ref_f-1
    bits = reference.bit_length()
    found = []
    candidate = search_prime(level, ref_c, ref_f, max_cofactor, "exact SQIsign-v3 prime")
    if candidate:
        found.append(candidate)
        print(f"{level}: exact reference prime has a scalar candidate", flush=True)
    else:
        print(f"{level}: no exact-prime candidate within cofactor bound {max_cofactor}", flush=True)
    # Prefer small c; permit an explicit, bounded increase in field size.
    for c in range(1, max_c+1, 2):
        if len(found) >= count:
            break
        for increase in range(extra_bits+1):
            if len(found) >= count:
                break
            f = bits+increase-c.bit_length()
            if f < 3 or (c, f) == (ref_c, ref_f):
                continue
            candidate = search_prime(level, c, f, max_cofactor, "p >= reference, bounded bit length")
            if candidate:
                found.append(candidate)
                print(f"{level}: p={c}*2^{f}-1, h={candidate['parameters']['odd_torsion_cofactor']}, "
                      f"p_bits={candidate['p_bits']}, N_bits={candidate['N_bits']}", flush=True)
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--levels", nargs="+", choices=TARGETS, default=["I"])
    parser.add_argument("--count", type=int, default=3)
    parser.add_argument("--max-c", type=int, default=4095)
    parser.add_argument("--max-cofactor", type=int, default=65536)
    parser.add_argument("--extra-bits", type=int, default=8,
                        help="maximum field-bit increase over the reference (0 for exact size)")
    parser.add_argument("--output", type=Path, default=Path("parameter_search.json"))
    args = parser.parse_args()
    if min(args.count, args.max_c, args.max_cofactor) <= 0 or args.extra_bits < 0:
        parser.error("search limits must be positive")
    candidates = []
    for level in args.levels:
        candidates.extend(search(level, args.count, args.max_c, args.max_cofactor, args.extra_bits))
    args.output.write_text(json.dumps({
        "source": "https://sqisign.org/spec/sqisign-20260901.pdf",
        "scope": "Scalar research candidates with alpha=a+u*i+b*j, kappa=i+t*j; defaults unchanged",
        "limits": {"levels": args.levels, "count": args.count,
                   "max_c": args.max_c, "max_cofactor": args.max_cofactor,
                   "extra_bits": args.extra_bits},
        "candidates": candidates,
    }, indent=2)+"\n", encoding="utf-8")
    print(f"Wrote {len(candidates)} candidates to {args.output}", flush=True)
    if not candidates:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
