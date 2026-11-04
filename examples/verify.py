"""Run the verify demonstration."""

if __package__ in (None, ""):
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import random
import time
from src.protocol.keygen import Keygen
from src.protocol.eval import Eval
from src.protocol.verify import Verify


def main():
    random.seed(0)
    kg = Keygen()
    pk = kg.pk
    sk = kg.sk
    x = bytes.fromhex("0000")
    eval = Eval(pk, sk, x)
    pi = eval.pi
    v = eval.v
    verify_s = time.time()
    verify = Verify(pk, pi, x, v)
    verify_e = time.time()
    print("----- verify -----")
    print(f"Verification result: {verify.verify}")
    print(f"Verification time: {verify_e - verify_s} seconds")


if __name__ == "__main__":
    main()
