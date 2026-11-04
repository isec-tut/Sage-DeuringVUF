from sage.all import *
from sage.schemes.elliptic_curves.ell_point import EllipticCurvePoint_field
from src.arithmetic.kummer import KummerPoint

def print_ideal(I):
    print(f"norm: {I.norm()}")
    dem = I.basis_matrix().denominator()
    print(f"denominator: {dem}")
    print("basis:")
    gens = I.gens()
    gen1 = gens[0]*dem
    gen2 = gens[1]*dem
    gen3 = gens[2]*dem
    gen4 = gens[3]*dem
    for n in range(4):
        print(f"{gen1[n]} {gen2[n]} {gen3[n]} {gen4[n]}")
    return

def print_Montgomery_curve(E):
    A = E.a2()
    print(f"{hex(A[0])} + {hex(A[1])}*i")
    return

def print_point(P, Q, PmQ=None):
    if isinstance(P, EllipticCurvePoint_field) and isinstance(Q, EllipticCurvePoint_field):
        print(f"{hex(P[0][0])} + {hex(P[0][1])}*i")
        print(f"{hex(Q[0][0])} + {hex(Q[0][1])}*i")
        if PmQ is not None:
            print(f"{hex(PmQ[0][0])} + {hex(PmQ[0][1])}*i")
    elif isinstance(P, KummerPoint) and isinstance(Q, KummerPoint):
        xP = P.x()
        xQ = Q.x()
        print(f"{hex(xP[0])} + {hex(xP[1])}*i")
        print(f"{hex(xQ[0])} + {hex(xQ[1])}*i")
        if PmQ is not None:
            xPmQ = PmQ.x()
            print(f"{hex(xPmQ[0])} + {hex(xPmQ[1])}*i")
    return