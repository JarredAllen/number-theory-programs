import math
import euclid


def ModularIntegers(n):
    """Returns the integers modulo n"""

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

        def get_modular_base(self):
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

        def __str__(self):
            return f"{self.a} mod {n}"

        def __repr__(self):
            return f"ModularIntegers({n})({self.a})"

    return Modular
