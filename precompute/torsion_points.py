from sage.all import *
import hashlib

# return GF(p^4), GF(p^2), square root of -1 in GF(p^2)
def calcFields(p):
    assert p % 4 == 3
    R = PolynomialRing(GF(p), name="x")
    x = R.gens()[0]
    Fp2 = GF(p**2, modulus=x**2+1, name="i")
    z = Fp2.random_element()
    while z.is_square():
        z = Fp2.random_element()
    t = ZZ(z + z**p)
    n = ZZ(z**(p+1))
    R = PolynomialRing(ZZ, name="x")
    x = R.gens()[0]
    Fp4 = GF(p**4, modulus=x**4 - t*x**2 + n, name="z")
    Fp2 = Fp4.subfield(2)
    i = Fp2.gen()**((p**2 - 1)//4)
    assert i**2 == -1
    return Fp4, Fp2, i

# transform an element x in Fp(zeta4) to an element in Fp2(zeta4d)
def Fp2ToFp2d(x, zeta4, zeta4d):
    p = x.base_ring().characteristic()
    return ZZ((x + x**p)/2) + ZZ((x - x**p)/(2*zeta4)) * zeta4d

# transform an element x in Fp2(zeta4d) to an element in Fp(zeta4)
def Fp2dToFp2(x, zeta4):
    x0 = x[0]
    x1 = x[1]
    return x0 + x1*zeta4

def fp2_sqrt(a, zeta4=None):
    K = a.parent()
    if zeta4 != None:
        K = K.base_ring()
    p = K.characteristic()
    gen = K.gen()
    if zeta4 != None:
        R = PolynomialRing(GF(p), name="x")
        x = R.gens()[0]
        Fp2 = GF(p**2, modulus=x**2+1, name="i")
        gen = Fp2.gen()
        a = Fp2ToFp2d(a, zeta4, gen)
    
    i = gen
    assert i**2 == -1

    a0 = a[0]
    a1 = a[1]

    x0 = (a0**2 + a1**2)**((p+1)//4)
    if a1 == 0:
        x0 = a0
    x0 = x0 + a0
    t0 = x0 + x0
    x1 = t0**((p-3)//4)
    x0 = x0*x1
    x1 = x1*a1
    t1 = (2*x0)**2
    if t1 == t0:
        x = x0 + x1*i
    else:
        x = x1 - x0*i

    if ZZ(t0) % 2 == 1 or (t0 == 0 and ZZ(t1) % 2 == 1):
        x = -x

    if zeta4 != None:
        x = Fp2dToFp2(x, zeta4)
        a = Fp2dToFp2(a, zeta4)

    #print(x**2)
    #print(a)
    assert x**2 == a

    return x

def deterministic_field_element(F, counter):
    p = ZZ(F.characteristic())
    degree = ZZ(F.degree())
    nbytes = (p.nbits() + 7) // 8 + 8

    coeffs = []
    for i in range(degree):
        data = b"torsion-basis:" + str(counter).encode() + b":" + str(i).encode()
        buf = b""
        block = 0
        while len(buf) < nbytes:
            buf += hashlib.sha256(data + b":" + str(block).encode()).digest()
            block += 1
        coeffs.append(ZZ(int.from_bytes(buf, byteorder="big") % p))

    if degree == 1:
        return F(coeffs[0])
    return F(coeffs)

def point_ord(E, Fp2, is_twist, l, e, start=0, return_next=False):
    p = E.base_ring().characteristic()
    if is_twist:
        ord = p - 1
    else:
        ord = p + 1
    assert is_prime(l)
    assert ord % l**e == 0

    counter = ZZ(start)
    while True:
        xP = deterministic_field_element(Fp2, counter)
        counter += 1
        try:
            P = E.lift_x(xP)
        except (ArithmeticError, ValueError):
            continue
        y = P.xy()[1]
        if (is_twist and not y in Fp2) or (not is_twist and y in Fp2):
            P = ord//l**e * P
            if not (l**(e-1)*P).is_zero():
                break
    assert (l**e*P).is_zero()

    # Choose a canonical sign so the same x-coordinate always gives the same point.
    if P[1][0] >= (p + 1)//2 or (P[1][0] == 0 and P[1][1] >= (p + 1)//2):
        P = -P
    if return_next:
        return P, counter
    return P

def TorsionBasis_2e(E, Fp2, is_twist, e):
    P, counter = point_ord(E, Fp2, is_twist, 2, e, return_next=True)
    Q, counter = point_ord(E, Fp2, is_twist, 2, e, start=counter, return_next=True)
    while (2**(e-1)*P).weil_pairing(2**(e-1)*Q, 2) == 1:
        Q, counter = point_ord(E, Fp2, is_twist, 2, e, start=counter, return_next=True)
    return P, Q

if __name__ == "__main__":
    from src.config.parameters import p, f

    print(f"p = {ZZ(p+1).factor()} - 1")
    Fp2d = GF(p**2, modulus=[1, 0, 1], names="i")

    Fp4, Fp2, zeta4 = calcFields(p)
    E0_fp4 = EllipticCurve(Fp4, [1, 0])
    P0_fp4, Q0_fp4 = TorsionBasis_2e(E0_fp4, Fp2, False, f)

    E0_fp2 = EllipticCurve(Fp2d, [1, 0])
    P0x_fp2, P0y_fp2 = [Fp2ToFp2d(v, zeta4, Fp2d.gen()) for v in P0_fp4.xy()]
    P0_fp2 = E0_fp2(P0x_fp2, P0y_fp2)
    Q0x_fp2, Q0y_fp2 = [Fp2ToFp2d(v, zeta4, Fp2d.gen()) for v in Q0_fp4.xy()]
    Q0_fp2 = E0_fp2(Q0x_fp2, Q0y_fp2)
    print(f"{P0_fp2 = }")
    print(f"{Q0_fp2 = }")

    assert P0_fp2.order() == 1 << f
    assert Q0_fp2.order() == 1 << f
    e0 = P0_fp2.weil_pairing(Q0_fp2, 1 << f)
    assert e0 ** (1 << f - 1) == -1

    # assert P0 == P and Q0 == Q
