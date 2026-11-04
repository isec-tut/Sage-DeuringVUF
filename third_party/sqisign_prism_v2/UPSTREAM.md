# Local SQIsign-v2.0 fork

- Repository: https://github.com/KULeuven-COSIC/sqisign_prism_v2
- Commit: `297ddf1f46fa677e0a2164fbcf51d433e12ccf7d`
- License: MIT; retain both `LICENSE` and `theta/LICENSE`.
- Imported: `quaternions.py`, `qlapoti.py`, `sos.py`, `ec.py`, `hd.py`, and
  the complete Python `theta/` library.

`UPSTREAM.json` retains the original upstream Git blob hashes in `files`.
The `local_modifications` entries separately record the hashes and reasons for
modified imported files; `tests/test_dependencies.py` checks the resulting
fork against this manifest. The only modified imported file is currently
`theta/montgomery_isogenies/kummer_line.py`, which adds `LadderBiscalar`,
`lift`, and `y_recovery` while preserving existing upstream methods.
The package initializers and this documentation are local. `params.py` is a
local adapter replacing the upstream signature presets and hardcoded constants.

`ideals.py` is a local extension. Its checked constructor builds O<alpha,N>
without importing `params`, so public precomputation can use it without a
circular dependency. Serialization uses the unchanged upstream
`LLLReducedBasis` and `SmallishGenerator` functions, with a fixed local Sage
seed for reproducible encodings that preserve the caller's RNG state.
The former `src/arithmetic/ideal.py` and its unused Gram/LLL helpers were removed.

The bundled theta implementation originates from
[ThetaIsogenies/two-isogenies](https://github.com/ThetaIsogenies/two-isogenies),
the code accompanying [ePrint 2023/1747](https://eprint.iacr.org/2023/1747.pdf).
Both (2,2)-isogeny calculations and Deuring correspondence algorithms use this
SQIsign fork; no separate theta source tree is required.

## DeuringVUF integration

`src/arithmetic/backend.py` connects the upstream arithmetic to the existing
DeuringVUF parameters and interfaces:

1. Random ideals, equivalent prime ideals, and LLL reduction use
   `quaternions.py`, which uses the upstream `sos.py` implementation.
2. Ideal conversion uses upstream `qlapoti.qlapoti` and `ec.EvalBasisPrecomp`.
   The correction is applied to an optional input basis, as required by
   Algorithm 5 of ePrint 2023/1251. Key generation also evaluates the
   odd-order points through the same product isogeny.
3. Product isogenies use `hd.py` and its packaged `theta/` implementation.
   Local differential-triple recovery evaluates with `lift=False` and recovers
   only the needed product factor, reducing square-root computations while
   preserving relative signs.

`IdealToIsogenyKeygen` returns `(success,E,P,Q,PmQ,U,V,UmV)`, with both
the odd and 2-power images represented by projective Kummer differential
triples. The odd images are selected directly from `lift=False` output.
Keygen mixes the odd basis and computes its new difference projectively,
and keeps that differential triple in `pk=(E,R,S,RmS,hint)`. Eval and Verify
compute their odd input kernels projectively. The Eval conversion still
returns `(success,E,U,V,UmV)`.

Local `projective.JacobianPoint` adapts signed full points to the unchanged
theta API. Endomorphism/basis matrices, representation scalars, and elliptic
kernel additions/doublings use Jacobian arithmetic. Indexing supplies an
affine view only at the gluing coordinate boundary; PARI pairings also still
require normalized full points. Kummer kernel/reference points are lifted
only when passed to the signed theta API. Public-basis validation uses
projective membership, torsion and differential-relation checks; it does not
run a Weil independence test.

The `Keygen`, `Eval`, `Amortized_Eval`, `Verify` call interfaces
and output `j(Ex)` are retained. The secret source pair and proof pair are stored
as matrices: `sk=(E,I,M,C0,U_can,V_can)`, `pi=(E3,C3,hint3)`.
The canonical public-side pair remains stored as points. Evaluation caches are `(K,Fx)`.
See [basis encoding and Tate diagnostics](../../docs/basis_encoding.md).
Proof randomness and curve models need not match the previous implementation
for a given seed.

### Mathematical conventions

- Keep the existing prime `p = 305*2^240 - 1`, odd prime `N`, quaternion order,
  field embeddings, twist, and precomputed action matrices. Upstream security
  level presets and their `constants.py` are not used.
- `params.f = 240` is the available torsion exponent used internally by Qlapoti.
  `REPRESENTATION_EXPONENT = N.nbits()` is the separate proof exponent. The
  proof's auxiliary ideal has norm `2^REPRESENTATION_EXPONENT - N`.
- The precomputed local matrices store images in columns. Upstream `ec` uses
  rows, so `params.py` transposes these matrices once.
- The internal product curve is the existing twist over Fp4. Upstream
  `hd.EmbeddedIsogeny` and `ec.ChangeOfBasis` use a reduced Tate exponent for
  Fp2. The adapter instead selects the factor with full-order Weil pairings
  and uses the bundled `theta.utilities.discrete_log.BiDLP` implementation
  (PARI Weil pairings and `fflog`) for input basis changes. The local
  `torsion.PointCoords` preserves the coordinate/scaling convention of the
  signature API; its optional pairing base is a Weil value. The raw upstream
  Tate-based `ec.PointCoords` is not a drop-in on the quartic twist: its Fp2
  final exponent gives a non-bilinear map there, and reproduced coordinate errors
  of 2^(f-1). The proper quartic exponent loses one bit of pairing order on the
  standard source pair. Transporting the same failing input to E0 over the
  quadratic subfield makes the raw upstream function work. These cases are
  regression-tested; the VUF protocol uses Weil/PARI instead.
  Local `torsion.ChangeOfBasis` follows the upstream row convention, and matrix
  evaluation reuses upstream `ec.EvalMatrix`. Eval passes the stored source
  matrix directly to the backend, avoiding point decoding followed by two DLPs.
  Endomorphism precomputation also uses local PointCoords, with `order=n` for
  an already n-torsion basis. The `e` mode retains upstream's 2^(f-e) scaling.
  The old bi_dlp entry points delegate to PointCoords rather than maintaining
  a second Sage discrete-log implementation. Importing the coordinate helper
  does not access the cached context, so it can run during context construction.
- `hd.Dim2Iso` uses the square-root chain and needs exactly the kernel torsion,
  without the two additional levels needed by `EllipticProductIsogeny`.
- Odd-order correction divides modulo `N`, whereas 2-power correction divides
  modulo `2^f`. A Qlapoti result whose embedded degree is not coprime to `N` is
  retried before evaluating the odd-order points.
- Key generation's N-torsion basis satisfies `Q0=i(P0)`, `j(P0)=lambda*P0`,
  `j(Q0)=-lambda*Q0`, where `lambda=alpha_scalar/alpha_frobenius mod N`.
  These relations are verified once when preparing the public basis. They
  determine the action matrix without an N-torsion DLP. `endomorphism_odd.py`
  evaluates its columns and their difference with projective differential
  chains over the quadratic subfield. Inverse-degree scaling is folded into
  the matrix modulo N. The hot path has no finite-field division, normalization,
  root extraction or half-point search. A consistent affine pair is recovered
  and transported to the twist at the product-isogeny API boundary. A
  simultaneous sign change of the pair is allowed.

The extended Kummer API is provided by this fork's theta package. The three
additional methods were ported from the former local Kummer copy; both copies
carry the same theta MIT copyright notice. No historical top-level theta
modules or path injection are needed at the compatibility boundary.
The canonical torsion hints and the VUF protocol itself also remain local.
No upstream signature protocol is substituted for DeuringVUF.

## Validation

Run from the repository root, with assertions enabled:

```sh
sage -python -m tests.test_dependencies
sage -python -m tests.test_backend
sage -python -m tests.test_protocol
sage -python -m tests.test_projective
sage -python -m tests.test_pointcoords
sage -python -m benchmarks.benchmark_vuf --trials 3 --seed 0
sage -python -m tests.test_odd_endomorphism
```

The migration regressions check upstream hashes, action conventions, torsion
exponents, three seeds and three inputs per key,
ordinary/amortized evaluation, ordinary verification, and rejection of altered
inputs, outputs, proof matrices/hints, and zero public-key points. They also
check basis-encoding round trips, the stored canonical public-side points,
Kummer image types/orders and projective odd-basis mixing against affine
arithmetic (including zero coefficients, signed scalars and rescaled inputs),
the local Weil coordinate solver. The benchmark
script exits with an error if evaluation or verification fails.
