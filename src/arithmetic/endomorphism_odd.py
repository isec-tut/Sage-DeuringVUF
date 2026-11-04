"""Division-free endomorphism evaluation on a prepared odd-torsion basis.

For alpha=a+u*i+b*j and Q=(i+t*j)(P), the known relation
j(P)=(a/b)*P-(u/b)*i(P) determines the action matrix without a DLP.
Its columns and their difference are evaluated
by projective differential chains. No halving, roots, DLP or affine x() is
needed in the evaluation loop.
"""

from operator import index

from src.arithmetic.kummer import KummerLine


def _double(P):
    """Projective doubling for y^2=x^3+x, with the common factor 2 retained."""
    X, Z = P
    plus = (X+Z)**2
    minus = (X-Z)**2
    product = plus*minus
    return product+product, (plus-minus)*(plus+minus)


def _add(P, Q, difference):
    XP, ZP = P
    XQ, ZQ = Q
    XD, ZD = difference
    if not ZP:
        return Q
    if not ZQ:
        return P
    if not ZD:
        return _double(P)
    # There is no nonzero two-torsion in an odd-order subgroup, so the
    # usual differential formula has no exceptional x(P-Q)=0 case here.
    plus = (XP+ZP)*(XQ-ZQ)
    minus = (XP-ZP)*(XQ+ZQ)
    return ZD*(plus+minus)**2, XD*(plus-minus)**2


def _multiply(n, P):
    one = P[0].parent().one()
    zero = P[0].parent().zero()
    result, next_point = (one, zero), P
    for bit in bin(n)[2:]:
        if bit == '0':
            next_point = _add(result, next_point, P)
            result = _double(result)
        else:
            result = _add(result, next_point, P)
            next_point = _double(next_point)
    return result


def _biscalar(P, Q, difference, a, b):
    """Return [a]P+[b]Q for nonnegative coefficients using a Euclidean chain."""
    if not a:
        return _multiply(b, Q)
    if not b:
        return _multiply(a, P)
    x0, x1 = Q, P
    s0, s1 = b, a
    while s0:
        if s1 < s0:
            s0, s1 = s1, s0
            x0, x1 = x1, x0
        if not s0:
            break
        if s1 <= 4*s0:
            s1 -= s0
            x0, difference = _add(x1, x0, difference), x0
        elif s0 % 2 == s1 % 2:
            s1 = (s1-s0)//2
            x0, x1 = _add(x1, x0, difference), _double(x1)
        elif s1 % 2 == 0:
            s1 //= 2
            x1, difference = _double(x1), _add(x1, difference, x0)
        else:
            s0 //= 2
            x0, difference = _double(x0), _add(x0, difference, x1)
    return x1 if s1 == 1 else _multiply(s1, x1)


class OddTorsionBasis:
    """Prepared public basis for one odd prime order N.

    The caller supplies a verified prime N and an optional quadratic subfield
    for the differential arithmetic. Basis relationships are checked once.
    """

    def __init__(self, P, Q, order, eigenvalue, zeta4, field=None, *,
                 frobenius_i=0, kappa_frobenius=0):
        self.order = index(order)
        if self.order <= 2 or self.order % 2 == 0:
            raise ValueError('order must be an odd prime')
        self.eigenvalue = index(eigenvalue) % self.order
        self.frobenius_i = index(frobenius_i) % self.order
        t = index(kappa_frobenius) % self.order
        curve = P.curve()
        if Q.curve() != curve or curve.a_invariants() != (0, 0, 0, 1, 0):
            raise ValueError('expected a basis on y^2=x^3+x')
        if P.is_zero() or Q.is_zero() or order*P != curve(0) or order*Q != curve(0):
            raise ValueError('basis points must have the supplied order')
        p = curve.base_ring().characteristic()
        iP = curve(-P[0], zeta4*P[1])
        jP = curve(*(v**p for v in P.xy()))
        lam, mu = self.eigenvalue, self.frobenius_i
        if Q != iP+t*jP or jP != lam*P-mu*iP:
            raise ValueError('points do not match the supplied alpha and kappa')
        if (lam**2+mu**2+int(p)) % self.order:
            raise ValueError('invalid Frobenius action')
        # Columns of T express (P, kappa(P)) in the reference basis (P,i(P)).
        self.basis_shear = t*lam % self.order
        self.basis_scale = (1-t*mu) % self.order
        if self.basis_scale == 0:
            raise ValueError('dependent odd-torsion basis')
        self.inverse_basis_scale = pow(self.basis_scale, -1, self.order)

        self.parent = KummerLine(curve)
        self.field = field if field is not None else curve.base_ring()
        self.P = (self.field(P[0]), self.field.one())
        self.Q = (self.field(Q[0]), self.field.one())
        self.difference = (self.field((P-Q)[0]), self.field.one())
        self.sum = _add(self.P, self.Q, self.difference)
        self.half = (self.order+1)//2

    def matrix(self, coefficients, scalar=1):
        """Row-major action matrix; columns are the images of P and Q."""
        if len(coefficients) != 4:
            raise ValueError('expected four integer order-basis coefficients')
        a, b, c, d = map(index, coefficients)
        scalar = index(scalar)
        N, half, lam = self.order, self.half, self.eigenvalue
        mu = self.frobenius_i
        common = a+d*half
        c_half, d_lam_half = c*half, d*lam*half
        # Reference matrices: i=[[0,-1],[1,0]],
        # j=[[lam,-mu],[-mu,-lam]], k=i*j. Return T^-1 M T.
        A = common+lam*c_half+mu*d*half
        B = -b-c_half-mu*c_half+d_lam_half
        C = b+c_half-mu*c_half+d_lam_half
        D = common-lam*c_half-mu*d*half
        h, z, inverse_z = self.basis_shear, self.basis_scale, self.inverse_basis_scale
        lower_left = C*inverse_z % N
        lower_right = (C*h+D*z)*inverse_z % N
        return tuple(scalar*v % N for v in (
            A-h*lower_left, A*h+B*z-h*lower_right,
            lower_left, lower_right,
        ))

    def _image(self, a, b):
        N = self.order
        a, b = a % N, b % N
        if a > N//2:
            a -= N
        if b > N//2:
            b -= N
        difference = self.sum if (a < 0) != (b < 0) else self.difference
        X, Z = _biscalar(self.P, self.Q, difference, abs(a), abs(b))
        # Promotion to the quartic field does not normalize (X:Z).
        return self.parent((X, Z))

    def images(self, coefficients, scalar=1):
        """Return x(tau(P)), x(tau(Q)), x(tau(P-Q)), scaled by scalar."""
        a, c, b, d = self.matrix(coefficients, scalar)
        return self._image(a, b), self._image(c, d), self._image(a-c, b-d)
