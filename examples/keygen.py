"""Run the keygen demonstration."""

if __package__ in (None, ""):
    from pathlib import Path
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import random
import time
from src.protocol.keygen import Keygen
from src.config.context import Fp2, zeta4
from precompute.torsion_points import Fp2ToFp2d
from src.encoding import printer
from src.encoding.serialization import decode_public_key, decode_secret_key


def print_fp2(label, value):
    a, b = Fp2ToFp2d(value, zeta4, Fp2.gen())
    print(f"{label}: {hex(a)} + {hex(b)}*i")


def main():
    random.seed(0)
    keygen_s = time.time()
    kg = Keygen()
    keygen_e = time.time()
    print("----- pk -----")
    print(f"pk ({len(kg.pk)} bytes): {kg.pk.hex()}")
    E_pk, R, S, RmS, torsion_hint = decode_public_key(kg.pk)
    print_fp2("E_pk coefficient A (y^2 = x^3 + A*x^2 + x)", E_pk.a2())
    for label, point in (("xR", R), ("xS", S), ("x(R-S)", RmS)):
        print_fp2(label, point.x())
    print(f"torsion_hint: {torsion_hint}")
    print("----- sk -----")
    print(f"sk ({len(kg.sk)} bytes): {kg.sk.hex()}")
    _, I_sk, M, C0, U_can, V_can = decode_secret_key(kg.sk, kg.pk)
    print("I_sk:")
    printer.print_ideal(I_sk)
    print("M:")
    print(M)
    print("C0 (relative to the deterministic source basis):")
    print(C0)
    print_fp2("xU_can", U_can.x())
    print_fp2("xV_can", V_can.x())
    print("--------------")
    print(f"Key generation took {keygen_e - keygen_s} seconds")


if __name__ == "__main__":
    main()
