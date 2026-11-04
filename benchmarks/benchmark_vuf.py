"""Benchmark keygen, evaluation, amortized evaluation and verification."""

if __package__ in (None, ""):
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import argparse
import platform
import random
import statistics
import time


from sage.all import set_random_seed
from sage.env import SAGE_VERSION

from src.config.parameters import DEFAULT_PARAMETERS


def average(values):
    return sum(values) / len(values)


def make_x(index):
    return index.to_bytes(2, byteorder="big")


def print_environment(trials, seed):
    print("----- environment -----")
    print(f"platform: {platform.platform()}")
    print(f"python:   {platform.python_version()}")
    print(f"sage:     {SAGE_VERSION}")
    print(f"trials:   {trials}")
    print(f"seed:     {seed}")
    parameters = DEFAULT_PARAMETERS
    print(f"params:   {parameters.name}")
    print(f"p:        {parameters.p.bit_length()} bits (c={parameters.c}, f={parameters.f})")
    print(f"N:        {parameters.N.bit_length()} bits")
    print("condition: all operations timed for trials calls; evaluations use the first key pair")
    print("validation: amortized proofs are additionally verified outside timing")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=100)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--rng-state", action="store_true",
        help="restore each Eval's saved Python RNG state in Amortized_Eval",
    )
    args = parser.parse_args()

    if args.trials <= 0:
        raise ValueError("--trials must be positive")

    from src.protocol.keygen import Keygen
    from src.protocol.eval import Eval, Amortized_Eval
    from src.protocol.verify import Verify
    from src.encoding.serialization import decode_output

    random.seed(args.seed)
    set_random_seed(args.seed)
    print_environment(args.trials, args.seed)
    print(f"rng_state replay: {args.rng_state}")

    keygen_times = []
    eval_times = []
    amortized_eval_times = []
    verify_times = []

    print("----- keygen start -----")
    for i in range(args.trials):
        start = time.perf_counter()
        key = Keygen()
        end = time.perf_counter()
        keygen_times.append(end - start)
        if i == 0:
            pk, sk = key.pk, key.sk
    print("----- keygen finish -----")

    # The first timed amortized evaluation initializes the reusable state.
    state = None

    print("----- eval and verify start -----")
    for i in range(args.trials):
        x = make_x(i)

        start = time.perf_counter()
        eval = Eval(pk, sk, x)
        end = time.perf_counter()
        eval_times.append(end - start)
        v = eval.v
        pi = eval.pi

        start = time.perf_counter()
        amortized_eval = Amortized_Eval(
            pk, sk, x, state=state,
            rng_state=eval.rng_state if args.rng_state else None,
        )
        end = time.perf_counter()
        amortized_eval_times.append(end - start)
        state = amortized_eval.state
        # Decode the serialized j-invariants before comparing field elements.
        eval_j = decode_output(eval.v)
        amortized_j = decode_output(amortized_eval.v)
        # if amortized_eval.v != v:
        if amortized_j != eval_j:
            raise AssertionError(
                f"{i}: j-invariant mismatch: x={x!r}, "
                f"Eval.j={eval_j}, Amortized_Eval.j={amortized_j}"
            )

        start = time.perf_counter()
        verify = Verify(pk, pi, x, v)
        end = time.perf_counter()
        verify_times.append(end - start)
        if verify.verify:
            print(f"{i}: ok (Eval.j == Amortized_Eval.j)")
        else:
            raise AssertionError(f"{i}: verification failed")
        if not Verify(pk, amortized_eval.pi, x, amortized_eval.v).verify:
            raise AssertionError(f"{i}: amortized verification failed")

    print("----- eval and verify finish -----")

    print("----- total -----")
    print(f"keygen:           {sum(keygen_times):.6f} seconds ({len(keygen_times)} calls)")
    print(f"eval:             {sum(eval_times):.6f} seconds ({len(eval_times)} calls)")
    print(f"amortized_eval:   {sum(amortized_eval_times):.6f} seconds "
          f"({len(amortized_eval_times)} calls, including first evaluation)")
    print(f"  first:          {amortized_eval_times[0]:.6f} seconds (1 call)")
    print(f"  subsequent:     {sum(amortized_eval_times[1:]):.6f} seconds "
          f"({len(amortized_eval_times) - 1} calls)")
    print(f"verify:           {sum(verify_times):.6f} seconds ({len(verify_times)} timed calls)")

    print("----- statistics per call -----")
    print("times in seconds; variance is population variance in seconds^2")
    print("amortized_eval includes the first evaluation")
    print(f"{'operation':<17} {'mean':>12} {'variance':>14} {'median':>12} {'min':>12} {'max':>12}")
    for name, times in (
        ("keygen", keygen_times),
        ("eval", eval_times),
        ("amortized_eval", amortized_eval_times),
        ("verify", verify_times),
    ):
        print(
            f"{name:<17} {average(times):>12.6f} "
            f"{statistics.pvariance(times):>14.6g} "
            f"{statistics.median(times):>12.6f} "
            f"{min(times):>12.6f} {max(times):>12.6f}"
        )


if __name__ == "__main__":
    main()
