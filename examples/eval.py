"""Run the evaluation demonstration."""

if __package__ in (None, ""):
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import random
import time
from src.protocol.keygen import Keygen
from src.protocol.eval import Eval
from src.encoding.serialization import decode_output, decode_proof
from examples.keygen import print_fp2


def main():
    random.seed(0)
    kg = Keygen()
    x = bytes.fromhex("0000")
    eval_s = time.perf_counter()
    evaluation = Eval(kg.pk, kg.sk, x)
    eval_e = time.perf_counter()
    print("----- eval -----")
    print(f"Input x: {x.hex()}")
    print(f"Output v ({len(evaluation.v)} bytes): {evaluation.v.hex()}")
    print_fp2("j(Ex)", decode_output(evaluation.v))
    print(f"Proof pi ({len(evaluation.pi)} bytes): {evaluation.pi.hex()}")
    E3, C3, hint3 = decode_proof(evaluation.pi)
    print_fp2("E3 coefficient A3 (y^2 = x^3 + A3*x^2 + x)", E3.a2())
    print("C3 (relative to the deterministic proof basis):")
    print(C3)
    print(f"hint3: {hint3}")
    print(f"Evaluation time: {eval_e - eval_s} seconds")


if __name__ == "__main__":
    main()
