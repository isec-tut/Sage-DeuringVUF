# DeuringVUF SageMath Implementation

An experimental implementation of DeuringVUF, covering key generation,
evaluation, amortized evaluation, and verification. The protocol code is
separated from the adapters to the SQIsign / PRISM arithmetic libraries.

The construction and example parameters follow
[ePrint 2023/1251](https://eprint.iacr.org/2023/1251.pdf).

## Arithmetic dependencies and attribution

Deuring correspondence algorithms, including the Qlapoti-based IdealToIsogeny
machinery, and the (2,2)-isogeny computations use a locally maintained fork of
[SQIsign-v2.0 / PRISM](https://github.com/KULeuven-COSIC/sqisign_prism_v2),
vendored in `third_party/sqisign_prism_v2/` at upstream commit
`297ddf1f46fa677e0a2164fbcf51d433e12ccf7d`.
Both the ordinary and square-root variants of product-isogeny chains use
that fork's `theta/` package. Kummer points, pairings, and basis conversion
helpers also come from the same package.

The theta implementation originates from
[two-isogenies](https://github.com/ThetaIsogenies/two-isogenies), accompanying
[An Algorithmic Approach to (2,2)-isogenies in the Theta Model and Applications
to Isogeny-based Cryptography](https://eprint.iacr.org/2023/1747.pdf).
This project uses the version bundled with SQIsign-v2.0.

Local changes add the DeuringVUF parameter adapter and the Kummer methods
`LadderBiscalar`, `lift`, and `y_recovery`. A local `ideals.py` helper uses
SQIsign's `LLLReducedBasis` and `SmallishGenerator` for secret ideal serialization
and provides checked ideal construction without initializing the default context.
The original upstream hashes and
the modified file's hash are recorded separately. See
[upstream integration notes](third_party/sqisign_prism_v2/UPSTREAM.md) and retain
the upstream [MIT license](third_party/sqisign_prism_v2/LICENSE) and
[theta MIT license](third_party/sqisign_prism_v2/theta/LICENSE) when redistributing.

## Execution environment

The full test suite was verified with the following environment:

| Component | Tested version |
| --- | --- |
| SageMath | 10.6 |
| Python (bundled with SageMath) | 3.12.5 |
| Runtime | Linux container using the `sagemath/sagemath` Docker image on Windows |

Run the protocol, examples, and tests with SageMath's Python interpreter via
`sage -python`. A standalone Python installation does not provide the required
`sage` modules. The arithmetic dependencies are included in `third_party/`;
no separate checkout of SQIsign or the theta library is needed.

With SageMath installed, check the interpreter versions with:

```sh
sage --version
sage -python --version
```

Alternatively, use Docker. The image used for the full test run was
`sagemath/sagemath@sha256:19995db6194f4a4bab18ce9a88556fd15b9ed5e916b4504fefe618a7796ddbdb`.
The digest fixes the tested image rather than following the changing `latest` tag.
From the repository root, run the tests in PowerShell with:

```powershell
$sageImage = 'sagemath/sagemath@sha256:19995db6194f4a4bab18ce9a88556fd15b9ed5e916b4504fefe618a7796ddbdb'
docker run --rm --mount "type=bind,source=$($PWD.Path),target=/work,readonly" -w /work $sageImage sage -python -B -u -m tests.run_all
```

Docker must be running with Linux containers enabled. The source directory is
mounted read-only, and `-B` disables Python bytecode writes. To run an example
instead, replace `tests.run_all` with, for example, `examples.eval`.

## Running the code

Run the following commands from the repository root in the SageMath environment.

Run all tests with:

```sh
sage -python -m tests.run_all
```

The runner discovers `tests/test_*.py` and runs each module in a fresh Python
process using the same interpreter. It streams test output, continues after
failures, and reports each module's result and elapsed time. It exits with code
0 when all tests pass, 1 if any test fails, or 130 if interrupted. Separate
processes isolate imports, cached precomputation, and RNG state between modules.

The integration checks are grouped by responsibility: `test_dependencies.py`
checks the fork's integrity and module connections, `test_backend.py` checks
arithmetic conventions and IdealToIsogeny images, and `test_protocol.py`
checks key generation, ordinary/amortized evaluation, and verification.

### Examples

Run the key generation, evaluation, and verification demos with:

```sh
sage -python -m examples.keygen
sage -python -m examples.eval
sage -python -m examples.verify
```

Each demo generates its own key pair. The key generation demo prints serialized
keys in hexadecimal alongside the decoded curve, ideal, and basis information.
The evaluation demo prints the input, serialized output and proof in hexadecimal,
the decoded output j-invariant and proof curve, basis matrix, and hint, and
evaluation time excluding key generation and display decoding.

### Benchmarks

Run the protocol benchmark with a short trial count:

```sh
sage -python -m benchmarks.benchmark_vuf --trials 3 --seed 0
```

For a longer run, or to replay the saved Python RNG state during amortized
evaluation:

```sh
sage -python -m benchmarks.benchmark_vuf --trials 100 --seed 0
sage -python -m benchmarks.benchmark_vuf --trials 100 --seed 0 --rng-state
```

| Argument | Default | Meaning |
| --- | --- | --- |
| `--trials N` | `100` | Number of calls to each operation: Keygen, Eval, Amortized Eval, and Verify. Use a positive integer; the current two-byte input counter supports up to 65,536 trials. |
| `--seed N` | `0` | Initial seed for both Python's and SageMath's RNGs. The RNGs are seeded once at the start of the run. |
| `--rng-state` | Disabled | Pass each ordinary Eval's saved Python RNG state to the corresponding Amortized Eval. This does not restore SageMath's RNG or guarantee identical proofs. |
| `-h`, `--help` | -- | Show command-line help and exit. |

Each key-generation call creates a key pair. All evaluations and verifications
use the **first key pair** and inputs encoded as two-byte big-endian counters
starting at zero. Amortized evaluation starts with `state=None` and reuses the
returned state across subsequent inputs. Its reported statistics include the
first call, which initializes the state and caches decoded keys; the totals
also list first-call and subsequent-call times separately.

The benchmark reports the platform, Python/SageMath versions, parameters, and
seed, followed by total times and per-call mean, population variance, median,
minimum, and maximum. Times are in seconds; variance is in seconds squared.
Module imports and initial public precomputation happen before the timed calls.
Key/proof/output encoding and decoding performed inside those calls are included,
as are public key validation and canonical basis preparation during Verify.

Every input is checked for agreement between ordinary and amortized outputs.
Ordinary proof verification is timed; amortized proofs are additionally verified
outside the timed sections. A mismatch or failed verification terminates the
run with an error. Randomized searches and runtime caches can affect timings,
so compare runs under the same environment and options.

Running `src/protocol/keygen.py` or `src/protocol/verify.py` directly launches
the corresponding demo in `examples/`. Running `src/protocol/eval.py` directly
launches the protocol benchmark.

## API

Use the API from SageMath's Python environment:

```python
from src import Keygen, Eval, Amortized_Eval, Verify

key = Keygen()
x = b"example"

evaluation = Eval(key.pk, key.sk, x)
check = Verify(key.pk, evaluation.pi, x, evaluation.v)
assert check.verify

cached = Amortized_Eval(
    key.pk, key.sk, x, evaluation.amortized, rng_state=evaluation.rng_state
)
assert cached.v == evaluation.v
assert Verify(key.pk, cached.pi, x, cached.v).verify

first = Amortized_Eval(key.pk, key.sk, b"first")
second = Amortized_Eval(key.pk, key.sk, b"second", state=first.state)
assert Verify(key.pk, second.pi, b"second", second.v).verify
```

Keys, proofs, and outputs are byte strings. The source-side power-of-two torsion
basis in the secret key and the basis in the proof are stored as change-of-basis
matrices relative to deterministic bases. These canonical bases are reconstructed
during decoding. See the [binary format](docs/serialization.md) for byte order,
lengths, and the relationship to the SQIsign v2 format.

| Value | Contents |
| --- | --- |
| `key.pk` | Bytes encoding the curve, basis hint, `x(R), x(S)`, and differential branch bit |
| `key.sk` | Bytes encoding the public key, secret ideal, `M_N`, and `C0` |
| `evaluation.pi` | Bytes encoding `(E3, C3, hint3)` |
| `evaluation.v` | Bytes encoding `j(Ex)` as an element of Fp2 |
| `evaluation.amortized` | Arithmetic objects `(K, Fx)` |

Each row of a basis encoding matrix contains the coefficients of one point.
`Eval` decodes the keys and reconstructs `(U_can, V_can)`. `Amortized_Eval`
caches the decoded keys in its returned `state` and skips decoding when that
state is reused under the same key pair.

The key generation adapter returns
`(success, curve, P, Q, PmQ, U, V, UmV)`, with both torsion triples represented
by Kummer points. Odd-order basis mixing and the input kernel computations in
Eval and Verify use x-only arithmetic. Signed point operations at the theta
library boundary use a Jacobian coordinate adapter.

See [basis encoding](docs/basis_encoding.md) for reconstruction details and a
reproducible explanation of why the upstream Tate-based `PointCoords` cannot
be used directly with the current field and twist conventions.

Every verification includes public key validation and preparation of the
canonical basis. `rng_state` captures Python's RNG state. Arithmetic also uses
SageMath's RNG, so this state alone does not guarantee reproducible proofs.

### Amortized evaluation

`Amortized_Eval` evaluates multiple inputs under the same key pair while reusing
auxiliary arithmetic and decoded keys. Unlike calling `Eval` independently for
each input, it samples the auxiliary ideal once and retains its intersection
with the secret ideal (`I_aux`) and the complementary curve (`E_aux`). Later
calls use the cached curve to identify the output component without computing
the input kernel again. Each call still computes and returns a proof for its
input; this is sequential reuse of setup rather than a batch proof.

Start with `state=None` (the default), then pass the returned `state` to the next
call. No preceding ordinary `Eval` is required:

```python
key = Keygen()
state = None
for x in (b"first", b"second", b"third"):
    result = Amortized_Eval(key.pk, key.sk, x, state=state)
    assert Verify(key.pk, result.pi, x, result.v).verify
    state = result.state
```

| Result | Meaning |
| --- | --- |
| `result.v` | Serialized VUF output; equals ordinary `Eval` for the same key and input |
| `result.pi` | Serialized proof, verified with the ordinary `Verify` API |
| `result.state` | In-memory `(I_aux, E_aux)` state with decoded-key caches |

Keep the returned state object intact to retain its decoded-key caches. It is
bound to the original key pair, and reusing it with different keys raises
`ValueError`. The state contains private key material and must be kept private;
it is not a serialized key or proof format. Start a fresh state for a new key pair.
The first call includes setup costs; subsequent calls reuse that setup.

For compatibility, the fourth positional argument (`amortized`) also accepts
the `(K, E_aux)` pair from an ordinary `Eval` result's `amortized` attribute,
as shown above. This legacy pair differs from `result.state`, which stores the
intersection `I_aux` rather than `K`. Pass either `amortized` or `state`, not both.
The optional keyword `rng_state` restores only Python's RNG before ideal
conversion; proofs may differ from ordinary evaluation even when outputs agree.

## Key, proof, and output sizes

The fixed parameters match Section 4.2 of
[the paper](https://eprint.iacr.org/2023/1251.pdf): `p` has 249 bits, `N` has
239 bits, and `f = 240`. The sizes below are serialized byte lengths. They
exclude parameter identifiers, transport headers, and Python/Sage object overhead.

Let `L = ceil(bit_length(p)/8) = 32`, `D = ceil(bit_length(N)/8) = 30`, and
`T = ceil(f/8) = 30`. An Fp2 element occupies two Fp coefficients, or 64 bytes.

| Object | Breakdown in bytes | Total |
| --- | --- | ---: |
| Public key `key.pk` | Curve coefficient `A`: 64; basis hint: 1; `x(R), x(S)`: 128; differential branch: 1 | **194 bytes** |
| Secret key `key.sk` | Public key: 194; ideal norm: 32; four generator coefficients: 128; four `M_N` coefficients: 120; four `C0` coefficients: 120 | **594 bytes** |
| Proof `evaluation.pi` | Curve coefficient `A3`: 64; four `C3` coefficients: 120; basis hint: 1 | **185 bytes** |
| Output `evaluation.v` | Fp2 encoding of `j(Ex)`: 64 | **64 bytes** |
| Output and proof `(v, pi)` | 64 + 185 | **249 bytes** |

The general formulas are `6L+2` for the public key,
`len(pk)+5L+4D+4T` for the secret key, `2L+4T+1` for the proof, and `2L` for
the output. The secret key includes the public key. Ordinary and amortized
evaluation produce proofs and outputs of the same size. The reusable `state`
contains in-memory arithmetic objects and is excluded from these sizes.

The constants are defined in [serialization.py](src/encoding/serialization.py).
See the [binary format](docs/serialization.md) for field order and byte order.

### Relationship to Section 4.3, "Sizes"

Section 4.3 counts compressed curves and points. This implementation represents
the same mathematical data with the following encoding choices:

| Object | Paper, Section 4.3 | This implementation |
| --- | --- | --- |
| Public key | Curve j-invariant and the x-coordinates and signs of two N-torsion points: three Fp2 elements plus two bits in the text | Curve coefficient `A` and two x-coordinates: three Fp2 elements, with one byte each for the basis hint and differential branch |
| Proof | Curve j-invariant and two x-coordinates with a relative sign for a 2^f-torsion basis: three Fp2 elements plus one bit | Curve coefficient `A3`, a 2-by-2 matrix `C3` relative to a deterministic basis, and a basis hint |
| Output | One Fp2 element, `j(Ex)` | Same |
| Secret key | No size calculation in this section | Ideal generator and change-of-basis matrices, as listed above |

The public key stores `x(R), x(S)` and a branch bit. Decoding reconstructs
`x(R-S)` with one square root using the biquadratic relation, without recovering
y-coordinates or solving discrete logarithms. The branch uses one bit of a
separate byte. The proof matrix `C3` contains four 30-byte coefficients. The
additional compression to "3 scalars" mentioned at the end of Section 4.3 is
not implemented.

The paper's abstract formulas assume `log2(p) = 2*lambda`, giving `4*lambda`
bits for the output and `12*lambda+1` bits for the proof. The public key is
listed as `12*lambda+2` bits in the text and `12*lambda+3` bits in Table 3.
These formulas do not specify this implementation's fixed-width binary format.
For example, each 249-bit Fp coefficient is stored in a 256-bit slot here.

Table 1 reports a 192-byte public key and a 256-byte proof. In particular,
Table 1 compares VRFs: the VUF-to-VRF transformation in Section 5.1 includes
the VUF output in the proof. The corresponding payload here is **249 bytes
for `(v, pi)`**, rather than `pi` alone. Table 1's values also do not exactly
match the formulas including sign bits in Section 4.3; distinguish those
reported sizes from the exact byte lengths of this implementation. The returned
`v` is the VUF output itself; it is not the hashed VRF output of Section 5.1.

## Repository layout

| Directory | Purpose |
| --- | --- |
| `src/protocol/` | Key generation, ordinary/amortized evaluation, verification, and input hashing |
| `src/arithmetic/` | Point and torsion arithmetic, endomorphisms, and SQIsign backend adapters |
| `src/config/` | Parameters, derived constants, and the public arithmetic context |
| `src/encoding/` | Basis encoding, serialization, and display helpers |
| `precompute/` | Public precomputation, parameter search, and candidate validation |
| `third_party/` | The local SQIsign-v2.0 fork and upstream licenses |
| `tests/` | Regression tests and the test runner |
| `benchmarks/` | Protocol timing benchmark |
| `examples/` | Key generation, evaluation, and verification demos |
| `docs/` | Encoding specifications and implementation notes |

The internal Qlapoti torsion exponent `f` (240 for `original`) is distinct from
the proof exponent `REPRESENTATION_EXPONENT = N.nbits()`. The auxiliary ideal
norm is `2^REPRESENTATION_EXPONENT - N`. See
[UPSTREAM.md](third_party/sqisign_prism_v2/UPSTREAM.md) for the upstream arithmetic,
basis, and twist conventions.

## Parameters and precomputation

The protocol uses one fixed parameter set, with no runtime selection option:

| Name | c | f | p bits | N bits |
| --- | --- | --- | --- | --- |
| `original` | 305 | 240 | 249 | 239 |

Inspect the scalar parameters with ordinary Python:

```sh
python -m src.config.parameters
```

The definitions and parameter type are in `src/config/parameters.py`.

| Parameter | Meaning or derived value |
| --- | --- |
| `c`, `f` | `p = c*2^f - 1` |
| `odd_torsion_cofactor` | `N = (p-1)/odd_torsion_cofactor` |
| `N` | Determines the proof exponent `N.bit_length()` and auxiliary ideal norm |
| `alpha_scalar`, `alpha_i`, `alpha_frobenius` | `alpha = alpha_scalar + alpha_i*i + alpha_frobenius*j` |
| `kappa_frobenius` | `kappa = i + kappa_frobenius*j` |
| `secret_norm_strategy` | The fixed set samples the secret ideal using `next_prime(2^(4*lambda))` |
| `precomputation_seed` | Seed for reproducible public precomputation, including finite field models |

Importing `src.config.context` alone does not run precomputation. The first
access to an arithmetic object such as `context.O0` or `context.U0_Fp4` invokes
`precompute/build_context.py`, which constructs the quaternions, base curve,
extension fields, odd-order torsion, and power-of-two torsion. Results are cached
per parameter set. Precomputation uses a local random seed and preserves the
caller's Python and SageMath RNG states.

Parameter objects are immutable and validate integer relations when created.
Precomputation also checks primality of `p` and `N`, the membership and norm of
`alpha`, and point orders. Changing the prime requires recomputing a compatible
`alpha`. The protocol always uses the fixed default values, including when
called as a library; environment variables do not select parameters.

For experiments, `context.get_context(parameters)` can build precomputations
for another parameter set. It does not switch the active protocol parameters:
the protocol and upstream adapters use `DEFAULT_PARAMETERS`.

## Optimizing odd-order point images

The relation `alpha = a + u*i + b*j` gives
`j(P0) = (a/b)*P0 - (u/b)*i(P0)`. The change to the basis
`Q0 = (i+t*j)(P0)` is also known, so the action matrix can be computed without
a discrete logarithm. For `original`, `u = t = 0`, giving the special case of
a Frobenius eigenbasis. These basis relations are checked once during setup.
Subsequent endomorphism evaluations use a known 2-by-2 matrix modulo N and
projective double-scalar multiplication. The image of `P0-Q0` is evaluated too,
preserving the relative sign needed for reconstruction.

This specialized evaluation path avoids normalization, finite field division,
square roots, and point halving during odd-order endomorphism evaluation.
Final point reconstruction and the isogeny API boundary still use affine
coordinates.

Basis conversion is centralized in `torsion.PointCoords`, using the bundled
`BiDLP` and PARI discrete logarithms. Precomputation uses
`PointCoords(..., order=n)` directly on bases of order n. The legacy `bi_dlp`
is a compatibility wrapper around two calls to `PointCoords`.

For the current twist over a quartic extension, directly applying the upstream
`ec.PointCoords` Tate reduction exponent for quadratic fields breaks bilinearity
and introduces a coefficient discrepancy of `2^239`. The local implementation
uses Weil pairings instead. If supplying `zeta`, use the Weil pairing value
returned by `point_coords_pairing`.

Isogeny images are also evaluated as differential triples. Reconstructing only
the required product component reduces the number of square roots. The fixed
public norm used for secret ideal generation is computed once per size.

The protocol benchmark measures the current implementation. Timings depend on
the environment and Qlapoti's randomized search.
