"""Projective point arithmetic at the boundary with the unchanged theta library.

JacobianPoint keeps a signed full point in Jacobian coordinates on the short
Weierstrass model x_short=x_Montgomery+A/3. The theta code's indexing API gets
an affine view only when it needs coordinates for the gluing matrices.
"""

from functools import lru_cache
from sage.all import ZZ
from src.arithmetic.kummer import KummerLine, KummerPoint


@lru_cache(maxsize=32)
def _short_coefficients(curve):
    a1, A, a3, a4, a6 = curve.a_invariants()
    if a1 != 0 or a3 != 0:
        raise ValueError('expected a curve with a1=a3=0')
    shift = A/3
    return shift, a4-A*shift


class JacobianPoint:
    """A full signed point with division-free addition and scalar multiplication."""

    __slots__ = ('_curve', '_X', '_Y', '_Z', '_shift', '_a', '_affine')

    def __init__(self, point):
        if isinstance(point, JacobianPoint):
            for name in self.__slots__:
                setattr(self, name, getattr(point, name))
            return
        self._curve = point.curve()
        self._shift, self._a = _short_coefficients(self._curve)
        F = self.base_ring()
        self._affine = point
        if point.is_zero():
            self._X, self._Y, self._Z = F.one(), F.one(), F.zero()
        else:
            self._X, self._Y, self._Z = point[0]+self._shift, point[1], F.one()

    def _new(self, X, Y, Z):
        result = object.__new__(type(self))
        result._curve = self._curve
        result._shift, result._a = self._shift, self._a
        result._X, result._Y, result._Z = X, Y, Z
        result._affine = None
        return result

    def zero(self):
        F = self.base_ring()
        return self._new(F.one(), F.one(), F.zero())

    def curve(self):
        return self._curve

    def base_ring(self):
        return self._curve.base_ring()

    def is_zero(self):
        return self._Z == 0

    def __bool__(self):
        return not self.is_zero()

    def to_affine(self):
        """Normalize only for APIs that need affine coordinates or PARI points."""
        if self._affine is None:
            if self.is_zero():
                self._affine = self._curve(0)
            else:
                inverse = 1/self._Z
                inverse2 = inverse*inverse
                self._affine = self._curve(
                    self._X*inverse2-self._shift, self._Y*inverse2*inverse
                )
        return self._affine

    def __getitem__(self, index):
        """Sage-compatible coordinates (x,y,1), not the internal Jacobian XYZ.

        Upstream gluing reads T[0] and assumes T[2]==1 for a nonzero point.
        Normalize here so its affine coordinate formulas remain unchanged.
        Addition, doubling and multiplication use _X,_Y,_Z directly.
        """
        return self.to_affine()[index]

    def xy(self):
        return self.to_affine().xy()

    def kummer(self, line=None):
        """Project to (X:Z) without normalizing the full point."""
        line = KummerLine(self._curve) if line is None else line
        if self.is_zero():
            return line.zero()
        Z2 = self._Z*self._Z
        return line((self._X-self._shift*Z2, Z2))

    def __eq__(self, other):
        if not isinstance(other, JacobianPoint):
            try:
                other = JacobianPoint(other)
            except (AttributeError, TypeError, ValueError):
                return False
        if self._curve != other._curve:
            return False
        if self.is_zero() or other.is_zero():
            return self.is_zero() and other.is_zero()
        z1, z2 = self._Z**2, other._Z**2
        return (self._X*z2 == other._X*z1
                and self._Y*z2*other._Z == other._Y*z1*self._Z)

    def __neg__(self):
        return self._new(self._X, -self._Y, self._Z)

    def double(self):
        if self.is_zero() or self._Y == 0:
            return self.zero()
        X, Y, Z = self._X, self._Y, self._Z
        XX, YY, ZZ2 = X*X, Y*Y, Z*Z
        YYYY = YY*YY
        S = 2*((X+YY)**2-XX-YYYY)
        M = 3*XX+self._a*ZZ2**2
        X3 = M*M-2*S
        return self._new(X3, M*(S-X3)-8*YYYY, 2*Y*Z)

    def __add__(self, other):
        if not isinstance(other, JacobianPoint):
            other = JacobianPoint(other)
        if self._curve != other._curve:
            raise ValueError('points are on different curves')
        if self.is_zero():
            return other
        if other.is_zero():
            return self
        z1 = self._Z**2
        if other._Z == 1:
            U1, U2 = self._X, other._X*z1
            S1, S2 = self._Y, other._Y*z1*self._Z
            z_product = self._Z
        else:
            z2 = other._Z**2
            U1, U2 = self._X*z2, other._X*z1
            S1, S2 = self._Y*z2*other._Z, other._Y*z1*self._Z
            z_product = self._Z*other._Z
        H, r = U2-U1, S2-S1
        if H == 0:
            return self.double() if r == 0 else self.zero()
        HH = H*H
        HHH, V = H*HH, U1*HH
        X3 = r*r-HHH-2*V
        return self._new(X3, r*(V-X3)-S1*HHH, H*z_product)

    def __sub__(self, other):
        return self + -JacobianPoint(other)

    def __mul__(self, scalar):
        scalar = ZZ(scalar)
        if scalar < 0:
            return (-self)*(-scalar)
        if scalar == 0:
            return self.zero()
        if scalar & (scalar-1) == 0:
            result = self
            for _ in range(scalar.nbits()-1):
                result = result.double()
            return result
        if scalar.nbits() >= 16:
            # Width-4 signed windows reduce additions for the long scalars used
            # by endomorphism matrices and the representation's N multiplier.
            digits = []
            while scalar:
                digit = ZZ(0)
                if scalar % 2:
                    digit = scalar % 16
                    if digit >= 8:
                        digit -= 16
                    scalar -= digit
                digits.append(digit)
                scalar //= 2
            twice = self.double()
            odd = {1: self}
            for digit in (3, 5, 7):
                odd[digit] = odd[digit-2]+twice
            result = self.zero()
            for digit in reversed(digits):
                result = result.double()
                if digit:
                    result = result+(odd[digit] if digit > 0 else -odd[-digit])
            return result
        result = self.zero()
        for bit in bin(scalar)[2:]:
            result = result.double()
            if bit == '1':
                result = result+self
        return result

    def __rmul__(self, scalar):
        return self*scalar

    def order(self):
        return self.to_affine().order()

    def weil_pairing(self, other, order):
        return self.to_affine().weil_pairing(JacobianPoint(other).to_affine(), order)


def kummer_point(point, line=None):
    if isinstance(point, KummerPoint):
        return point
    return JacobianPoint(point).kummer(line)


def lift_kummer(point):
    """Lift one Kummer point at the signed theta/PARI API boundary."""
    if point.is_zero():
        return point.parent().curve()(0)
    return point.lift()


def odd_linear_combination(P, Q, difference, a, b, order):
    """x(aP+bQ) in an independent odd-prime basis, entirely projectively."""
    a, b = ZZ(a) % order, ZZ(b) % order
    if not a:
        return b*Q
    if not b:
        return a*P
    # The historical local ladder's m,n arguments compute [n]P+[m]Q.
    return P.LadderBiscalar(Q, difference, b, a)


def valid_odd_basis(curve, basis, order):
    """Check curve membership, N-torsion and the projective differential relation.

    This checks a consistent nonzero differential triple, not Weil independence.
    """
    if len(basis) != 3:
        return False
    A = curve.a2()
    for P in basis:
        if not isinstance(P, KummerPoint) or P.parent().curve() != curve:
            return False
        X, Z = P.XZ()
        if Z == 0 or not (X*Z*(X*X+A*X*Z+Z*Z)).is_square():
            return False
        if not (order*P).is_zero():
            return False
    (X, Z), (Y, T), (D, W) = (P.XZ() for P in basis)
    delta = X*T-Y*Z
    if delta == 0:
        return False
    fX = X*(X*X+A*X*Z+Z*Z)
    fY = Y*(Y*Y+A*Y*T+T*T)
    numerator = (D*Z*T+W*(A*Z*T+X*T+Y*Z))*delta**2
    numerator -= W*(fX*T**3+fY*Z**3)
    return numerator**2 == 4*fX*fY*W**2*Z**3*T**3
