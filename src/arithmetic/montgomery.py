"""Projective Montgomery arithmetic and recovery of affine point pairs."""

from sage.all import Integer
from src.arithmetic.kummer import KummerPoint, sqrt_Fp2


def xDBL(P, E):
    """
    function for Montgomery doubling with projective curve constant
        
    Input:  projective point P = (X:Z), curve constants (A:C) 
    Output: projective point [2]P = (X2:Z2)

    Cost: 4M + 2S + 8a
    """
    X, Z = P
    A, C = E

    t0 = X - Z      
    t1 = X + Z
    t0 = t0**2
    t1 = t1**2
    Z2 = C * t0
    Z2 = Z2 + Z2
    Z2 = Z2 + Z2
    X2 = Z2 * t1
    t1 = t1 - t0
    t0 = C + C
    t0 = t0 + A
    t0 = t0 * t1
    Z2 = Z2 + t0
    Z2 = Z2 * t1
        
    return X2, Z2


def xADD(P, Q, PmQ):
    """
    function for Montgomery differential addition
        
    Input:  projective coordinates P = (XP : ZP), 
           Q=(XQ : ZQ), and their difference 
            x(P-Q) = (xPQ : zPQ) 
    Output: coordinates of sum P + Q = (XQP : ZQP)

    Cost: 4M + 2S + 6a
    """
    XP, ZP = P
    XQ, ZQ = Q
    xPQ, zPQ = PmQ

    t0 = XP + ZP     
    t1 = XP - ZP   
    XP = XQ - ZQ   
    ZP = XQ + ZQ   
    t0 = XP * t0   
    t1 = ZP * t1   
    ZP = t0 - t1   
    XP = t0 + t1   
    ZP = ZP**2     
    XQP = XP**2    
    ZQP = xPQ * ZP 
    XQP = XQP * zPQ
        
    return XQP, ZQP


def xDBLADD(P, Q, PmQ, E24):
    """
    function for step in Montgomery ladder
    simultaneous doubling and differential addition
        
    Input: projective coordinates P=(XP:ZP) and Q=(XQ:ZQ), 
            projective difference P-Q=(xPQ:zPQ) and 
            curve constant A24/C24=(A+2C)/4C.   
    Output: projective coordinates of 2P=(X2P:Z2P)
            and Q+P=(XQP:ZQP)

    Cost: 8M + 4S + 8A
    """
    XP, ZP = P
    XQ, ZQ = Q
    xPQ, zPQ = PmQ
    A24, C24 = E24
        
    t0 = XP + ZP                  
    t1 = XP - ZP 
    X2P = t0**2
    t2 = XQ - ZQ
    XQP = XQ + ZQ
    t0 = t0 * t2
    Z2P = t1**2
    t1 = t1 * XQP
    t2 = X2P - Z2P
    Z2P = Z2P * C24
    X2P = X2P * Z2P
    XQP = A24 * t2
    ZQP = t0 - t1
    Z2P = XQP + Z2P
    XQP = t0 + t1
    Z2P = Z2P * t2
    ZQP = ZQP**2
    XQP = XQP**2
    ZQP = xPQ * ZQP
    XQP = XQP * zPQ

    return X2P, Z2P, XQP, ZQP


def Ladder(P, E, m):
    m = Integer(m)
    if not isinstance(m, (int, Integer)):
        try:
            m = Integer(m)
        except:
            raise TypeError(f"Cannot coerce input scalar {m = } to an integer")
    if m == 0:
        return (1, 0)
    m = abs(m)

    X0, Z0 = 1, 0
    X1, Z1 = P
    XP, ZP = P
    A24, C24 = E
    for i in reversed(range(m.nbits())):
        mi = (m >> i) & 1
        if mi == 1:
            X1, Z1, X0, Z0 = xDBLADD([X1, Z1], [X0, Z0], [XP, ZP], [A24, C24])
        else:
            X0, Z0, X1, Z1 = xDBLADD([X0, Z0], [X1, Z1], [XP, ZP], [A24, C24])
    return (X0, Z0)


def LadderBiscalar(P, Q, PmQ, E24, m, n):
    """
    Output: [m]P + [n]Q
    """
    if not isinstance(m, (int, Integer)):
        try:
            m = Integer(m)
        except:
            raise TypeError(f"Cannot coerce input scalar {m = } to an integer")
    if not isinstance(n, (int, Integer)):
        try:
            n = Integer(n)
        except:
            raise TypeError(f"Cannot coerce input scalar {n = } to an integer")
    if m == 0 and n == 0:
        return (1, 0)
    m = abs(m)
    n = abs(n)

    # k is the common bitlength (as in the spec)
    k = max(m.bit_length(), n.bit_length())

    # Step 2–5: choose initial sigma
    if (m % 2 == 0) and (n % 2 == 1):
        sigma0, sigma1 = 1, 0
    else:
        sigma0, sigma1 = 0, 1

    # Step 6: m', n'
    mp = m
    np_ = n

    # Step 7–8: make m' odd if needed
    if m % 2 == 0:
        mp -= 1

    # Step 9–10: make n' odd if needed
    if n % 2 == 0:
        np_ -= 1

    # Step 11: build b0, b1 with a leading 0
    # b0 = (0, m'_{k-1}, ..., m'_0)
    # b1 = (0, n'_{k-1}, ..., n'_0)
    # indices: b?[0] is the leading 0, b?[i+1] is bit i
    b0 = [0] * (k + 1)
    b1 = [0] * (k + 1)
    for i in range(k):
        b0[i + 1] = (mp >> i) & 1
        b1[i + 1] = (np_ >> i) & 1

    # r will have length 2k: r_0,...,r_{2k-1}
    r = [0] * (2 * k)

    # We'll mutate sigma0,sigma1 during the loop
    for i in range(k):
        # pick the right b-array by sigma
        arr0 = b0 if sigma0 == 0 else b1
        arr1 = b0 if sigma1 == 0 else b1

        # Step 13–14
        r[2 * i] = arr0[i] ^ arr0[i + 1]
        r[2 * i + 1] = arr1[i] ^ arr1[i + 1]

        # Step 15–18: swap sigma0,sigma1 if r_{2i+1} = 1
        if r[2 * i + 1] == 1:
            sigma0, sigma1 = sigma1, sigma0

    # At this point sigma0,sigma1 are their final values
    # --------------------------
    # 2. Evaluation stage
    # --------------------------

    # Step 20: R0 = (1:0)  (point at infinity in x-only coords)
    R0 = (1, 0)

    # Step 21: T = (T0,T1) = (P,Q)
    T0 = P
    T1 = Q

    T = [T0, T1]  # just for clarity

    # Step 22–23: R1, R2 from final sigma
    T_sigma0 = T[sigma0]
    T_sigma1 = T[(sigma0 + 1) % 2]

    R1 = T_sigma0
    R2 = T_sigma1

    # Step 24–25: D1, D2
    D1 = R1
    D2 = R2

    # Step 26: R2 = xADD(R1, R2, P-Q)
    R2 = xADD(R1, R2, PmQ)

    # Step 27–28: F1, F2
    F1 = R2
    F2 = PmQ

    # Main loop: i from k-1 down to 0
    for i in range(k - 1, -1, -1):
        # Step 30: h = r_{2i} + r_{2i+1}
        r0 = r[2 * i]
        r1 = r[2 * i + 1]
        h = r0 + r1

        # R-array as a list for indexing R_h
        R = [R0, R1, R2]

        # Step 31: T0 <- R_{h (mod 2)}
        T0 = R[h % 2]

        # Step 32: T = (T0, T1) <- (T0, R2)
        T1 = R2
        T = [T0, T1]

        # Step 33: T0 <- xDBL(T_{floor(h/2)})
        # floor(h/2) is 0 for h=0,1 and 1 for h=2
        idx_dbl = h // 2
        T0 = xDBL(T[idx_dbl], E24)
        T[0] = T0

        # Step 34–35: T1, T2 from R
        T1 = R[r1]          # r1 is 0 or 1
        T2 = R[r1 + 1]      # 1 or 2

        # Step 36–39: maybe swap D1, D2
        if r1 == 1:
            D1, D2 = D2, D1

        # Step 40: T1 <- xADD(T1, T2, D1)
        T1 = xADD(T1, T2, D1)

        # Step 41: T2 <- xADD(R0, R2, F1)
        T2 = xADD(R0, R2, F1)

        # Step 42–45: maybe swap F1, F2 if h is odd
        if (h % 2) == 1:
            F1, F2 = F2, F1

        # Step 46–48: update R0,R1,R2
        R0 = T0
        R1 = T1
        R2 = T2

    # After the loop, the result is R0 (as in Ladder/Ladder3pt style)
    idx = ((m & 1) ^ 1) + ((n & 1) ^ 1)
    X_res, Z_res = R[idx]
    return (X_res, Z_res)


def ec_is_zero(P):
    X, Z = P
    return Z == 0


def order_twof(P, t):
    #print(P)
    if P.is_zero():
        return False
    #for i in range(t-1):
    #    print((P*2**(i+1)).x())
    test = P*2**(t-1)
    #print(test)
    if test.is_zero():
        return False
    test = test*2
    return test.is_zero()


def ProjectiveDifference(P, Q, E):
    XP, ZP = P
    XQ, ZQ = Q
    A, C = E
    BXX = C*(XP*XQ-ZP*ZQ)**2
    BXZ = C*(XP*XQ+ZP*ZQ)*(XP*ZQ+ZP*XQ)+2*A*XP*XQ*ZP*ZQ
    BZZ = C*(XP*ZQ-ZP*XQ)**2
    gamma = C*(C.conjugate() * ZP.conjugate() * ZQ.conjugate())**2
    BXX = gamma*BXX
    BXZ = gamma*BXZ
    BZZ = gamma*BZZ
    delta = sqrt_Fp2(BXZ**2-BXX*BZZ, canonical=True)
    XPQ = delta+BXZ
    ZPQ = BZZ
    return (XPQ, ZPQ)


def lift_x(E, P, Q, PmQ):
    if isinstance(P, KummerPoint) and isinstance(Q, KummerPoint) and isinstance(PmQ, KummerPoint):
        XP, ZP = P.XZ()
        XQ, ZQ = Q.XZ()
        XPQ, ZPQ = PmQ.XZ()
    else:
        XP, ZP = P
        XQ, ZQ = Q
        XPQ, ZPQ = PmQ

    A = E.a2()
    A24 = A + 2
    C24 = 4
    B = 1
    xP = XP / ZP
    yP2 = xP**3 + E.a2()*xP**2 + xP
    yP = sqrt_Fp2(yP2, canonical=True)
    assert yP**2 == yP2
    P_full = E([xP, yP])

    _, _, XPpQ, ZPpQ = xDBLADD((XP, ZP), (XQ, ZQ), (XPQ, ZPQ), (A24, C24))
    v1 = xP*ZQ
    v2 = XQ + v1
    v3 = XQ - v1
    v3 = v3*v3
    v3 = v3*XPpQ
    v1 = 2*A*ZQ
    v2 = v2 + v1
    v4 = xP*XQ
    v4 = v4 + ZQ
    v2 = v2*v4
    v1 = v1*ZQ
    v2 = v2 - v1
    v2 = v2*ZPpQ
    Y_ = v2 - v3
    v1 = 2*B*yP
    v1 = v1*ZQ
    v1 = v1*ZPpQ
    X_ = v1*XQ
    Z_ = v1*ZQ

    xQ = X_/Z_
    yQ = Y_/Z_
    Q_full = E([xQ, yQ])

    xPpQ = XPpQ/ZPpQ
    xPmQ = XPQ/ZPQ
    assert (P_full + Q_full)[0] == xPpQ
    assert (P_full - Q_full)[0] == xPmQ
    return P_full, Q_full
