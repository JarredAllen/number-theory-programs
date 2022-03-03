import math
import euclid


def ModularIntegers(n):
    """Returns the integers modulo n"""
    assert n > 0, "Must be modulo a positive integer"

    class Modular:
        """A class representing arithmetic modulo n"""

        def __init__(self, a):
            self.a = a % n

        @classmethod
        def __iter__(self):
            def iterate():
                for i in range(n):
                    yield Modular(i)

            return iterate()

        def __radd__(self, other):
            if type(other) == int:
                other = Modular(other)
            return Modular(self.a + other.a)

        def __add__(self, other):
            if type(other) == int:
                other = Modular(other)
            return Modular(self.a + other.a)

        def __rsub__(self, other):
            if type(other) == int:
                other = Modular(other)
            return Modular(self.a - other.a)

        def __sub__(self, other):
            if type(other) == int:
                other = Modular(other)
            return Modular(self.a - other.a)

        def __neg__(self):
            return Modular(n - self.a)

        def __rmul__(self, other):
            if type(other) == int:
                other = Modular(other)
            return Modular(self.a * other.a)

        def __mul__(self, other):
            if type(other) == int:
                other = Modular(other)
            return Modular(self.a * other.a)

        def __truediv__(self, other):
            return self * other.multiplicative_inverse()

        def __pow__(self, p):
            return Modular(pow(self.a, p, n))

        def __eq__(self, other):
            if type(other) is int:
                self.a == other % n
            else:
                return self.a == other.a and n == other.get_modular_base()

        def to_nth_root_of_unity(self):
            """Map to the nth root of unity by the projection homomorphism which maps 1 mod n to exp(2 pi i / n)"""
            re = math.cos(2 * math.pi * self.a / n)
            if abs(re) < 1e-12:
                re = 0
            im = math.sin(2 * math.pi * self.a / n)
            if abs(im) < 1e-12:
                im = 0
            return complex(re, im)

        @classmethod
        def get_modular_base(cls):
            """Get the modulus for this ring"""
            return n

        def multiplicative_inverse(self):
            """Compute the multiplicative inverse, raising an error if this value is not a unit"""
            if math.gcd(self.a, n) != 1:
                raise ValueError(f"{self.a} does not have an inverse modulo {n}")
            x, y = euclid.extended_euclidean(self.a, n)
            if self.a * x + n * y == 1:
                return Modular(x)
            else:
                return -Modular(x)

        @classmethod
        def natural_project_from(cls, value):
            """Perform the natural projection into this additive group from the additive group the value is in"""
            m = n // value.get_modular_base()
            assert m * value.get_modular_base() == n, "Value is not from a subring"
            return Modular(value.a * m)

        def __str__(self):
            return f"{self.a} mod {n}"

        def __repr__(self):
            return f"ModularIntegers({n})({self.a})"

    return Modular


def DihedralGroup(n):
    """Get the dihedral group of order n"""
    assert n % 2 == 0 and n > 0, "Dihedral groups must have positive, even order"

    class Dihedral:
        f"The dihedral group of order {n}"

        e = Dihedral(0, 0)
        s = Dihedral(1, 0)
        r = Dihedral(0, 1)

        def __init__(self, s, r):
            self._s = s % 2
            self._r = r % (n // 2)

        @classmethod
        def __iter__(self):
            def iterate():
                for g in range(n):
                    yield Dihedral(g % 2, g // 2)

            return iterate()

        def __mul__(self, other):
            if other._s == 0:
                return Dihedral(self._s, self._r + other._r)
            else:
                return Dihedral(self._s + 1, -self._r + other._r)

        def __truediv__(self, other):
            return self * other.multiplicative_inverse()

        def __pow__(self, p):
            if self._s == 0:
                return Dihedral(0, self._r * p)
            else:
                if p % 2 == 0:
                    return Dihedral(0, 0)
                else:
                    return self

        def __eq__(self, other):
            return (
                self._s == other._s and self._r == other._r and n == other.get_order()
            )

        @classmethod
        def get_order(cls):
            """Get the modulus for this ring"""
            return n

        def multiplicative_inverse(self):
            """Compute the multiplicative inverse, raising an error if this value is not a unit"""
            if math.gcd(self.a, n) != 1:
                raise ValueError(f"{self.a} does not have an inverse modulo {n}")
            x, y = euclid.extended_euclidean(self.a, n)
            if self.a * x + n * y == 1:
                return Modular(x)
            else:
                return -Modular(x)

        @classmethod
        def natural_project_from(cls, value):
            """Perform the natural projection into this dihedral group from the dihedral group the value is in"""
            m = n // value.get_order()
            assert m * value.get_order() == n, "Value is not from a subring"
            return Modular(value.a * m)

        def __str__(self):
            if self._s == 0:
                return f"r^{self._r}"
            else:
                return f"sr^{self._s}"

        def __repr__(self):
            return f"DihedralGroup({n})({self._s}, {self._r})"

    return Dihedral
