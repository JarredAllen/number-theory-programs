# Group interface: methods provided on groups in here
# and <group>.e should be identity


def CyclicGroup(n):
    """Produce the cyclic group of order n"""
    assert n > 0, "Must have positive order"

    class Cyclic:
        f"The cyclic group of order {n}"

        # Functions to make one

        def __init__(self, a):
            self._a = r % n

        @classmethod
        def natural_project_from(cls, value):
            """Perform the natural projection into this dihedral group from the dihedral group the value is in"""
            m = n // value.get_order()
            if m * value.get_order() != n:
                raise TypeError(f"{value} is not from a subgroup of {cls.__doc__}")
            return Modular(value.a * m)

        @classmethod
        def __iter__(self):
            def iterate():
                for a in range(n):
                    yield Cyclic(a)

            return iterate()

        # Group-theoretic manipulations

        def __mul__(self, other):
            try:
                return Cyclic(self._a + other._a)
            except AttributeError:
                return NotImplemented

        def __truediv__(self, other):
            try:
                return self * other.multiplicative_inverse()
            except AttributeError:
                return NotImplemented

        def __pow__(self, p):
            try:
                return Cyclic(self._a * p)
            except AttributeError:
                return NotImplemented

        def __eq__(self, other):
            return self._a == other._a and n == other.get_order()

        def multiplicative_inverse(self):
            """Compute the multiplicative inverse"""
            return Cyclic(-self._a)

        # Display control

        def __str__(self):
            if self._a == 0:
                return "1"
            else:
                return "x^{{%s}}" % self._a

        def __repr__(self):
            return f"CyclicGroup({n})({self._a})"

        # Type-level manipulations

        @classmethod
        def get_order(cls):
            """Get the order of this group"""
            return n

    Cyclic.e = Cyclic(0)
    Cyclic.x = Cyclic(1)

    return Cyclic


def DihedralGroup(n):
    """Get the dihedral group of order n"""
    assert n % 2 == 0 and n > 0, "Dihedral groups must have positive, even order"

    class Dihedral:
        f"The dihedral group of order {n}"

        # Functions to make one

        def __init__(self, s, r):
            self._s = s % 2
            self._r = r % (n // 2)

        @classmethod
        def natural_project_from(cls, value):
            """Perform the natural projection into this dihedral group from the dihedral group the value is in"""
            m = n // value.get_order()
            if m * value.get_order() != n:
                raise TypeError(f"{value} is not from a subgroup of {cls.__doc__}")
            return Dihedral(value._s, value._r * m)

        @classmethod
        def __iter__(self):
            def iterate():
                for g in range(n):
                    yield Dihedral(g // (n // 2), g % (n // 2))

            return iterate()

        # Group-theoretic manipulations

        def __mul__(self, other):
            try:
                if other._s == 0:
                    return Dihedral(self._s, self._r + other._r)
                else:
                    return Dihedral(self._s + 1, -self._r + other._r)
            except AttributeError:
                return NotImplemented

        def __truediv__(self, other):
            try:
                return self * other.multiplicative_inverse()
            except AttributeError:
                return NotImplemented

        def __pow__(self, p):
            try:
                if self._s == 0:
                    return Dihedral(0, self._r * p)
                else:
                    if p % 2 == 0:
                        return Dihedral(0, 0)
                    else:
                        return self
            except AttributeError:
                return NotImplemented

        def __eq__(self, other):
            return (
                self._s == other._s and self._r == other._r and n == other.get_order()
            )

        def multiplicative_inverse(self):
            """Compute the multiplicative inverse, raising an error if this value is not a unit"""
            if self._s == 0:
                return Dihedral(0, -self._r)
            else:
                return self

        # Display control

        def __str__(self):
            s = "s" if self._s else ""
            if self._r == 0:
                r = ""
            elif self._r == 1:
                r = "r"
            else:
                r = f"r^{self._r}"
            return (s + r) or "e"

        def __repr__(self):
            return f"DihedralGroup({n})({self._s}, {self._r})"

        # Type-level manipulations

        @classmethod
        def get_order(cls):
            """Get the modulus for this ring"""
            return n

    Dihedral.e = Dihedral(0, 0)
    Dihedral.s = Dihedral(1, 0)
    Dihedral.r = Dihedral(0, 1)

    return Dihedral
