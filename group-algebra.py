def GroupRing(R, G):
    """Returns the group ring associated with R and G

    R must be some ring with identity, and it must also be a function
    which maps values into its type, with 0 mapped to the additive
    identity and 1 mapped to the multiplicative identity.

    G must be some group.
    """

    class Ring:
        f"The group ring {R}{G}"

        # Construct values in the ring

        def __init__(self, value):
            """Create a new group ring element.

            If value is a list, then it is taken as a list in the same
            order as the group elements.

            If value is a ring element, then it is taken as the the
            coefficient on the group identity (such that 0 maps to the
            additive identity and 1 maps to the multiplicative
            identity).
            """
            if type(value) == list:
                self._values = [R(v) for v in value]
            else:
                r = R(value)
                self._values = [r if g == G.e else 0 for g in G.__iter__()]

        @classmethod
        def from_function(cls, func):
            return cls([func(g) for g in G.__iter__()])

        @classmethod
        def natural_project_from(cls, value):
            """Attempt a natural projection into the group ring."""
            if type(value) == cls:
                # If projecting an element of the group ring, return it directly
                return value
            elif type(value) == R:
                # If projecting an element of the base ring, project
                return cls(value)
            try:
                # Try projecting from an element of a subgroup
                g = G.natural_project_from(value)
                return cls.from_function(lambda h: R(1 if g == h else 0))
            except TypeError:
                pass
            try:
                elements = list(G.__iter__())
                values = [R(0) for _ in elements]
                for g in iter(value.get_base_group()):
                    index = elements.index(G.natural_project_from(g))
                    assert values[index] == R(0), "Group projection is not injective"
                    values[index] += value(g)
                return cls(values)
            except TypeError:
                pass
            raise TypeError(f"{value} is not from a subring of {R}{G}")

        # Manipulate values in the ring

        def __call__(self, arg):
            """Get the coefficient on the given group element"""
            elements = list(G.__iter__())
            return self._values[elements.index(arg)]

        def __add__(self, other):
            """Add the two functions"""
            return Ring([a + b for (a, b) in zip(self._values, other._values)])

        def __sub__(self, other):
            """Subtract the two functions"""
            return Ring([a - b for (a, b) in zip(self._values, other._values)])

        def __neg__(self):
            """Additive inverse"""
            return Ring([-a for a in self._values])

        def __mul__(self, other):
            """Convolve the two functions"""
            self = Ring.natural_project_from(self)
            other = Ring.natural_project_from(other)
            elements = list(G.__iter__())
            answers = [R(0) for _ in elements]
            for e1 in elements:
                for e2 in elements:
                    i = elements.index(e1 * e2)
                    answers[i] += self(e1) * other(e2)
            return Ring(answers)

        def __rmul__(self, other):
            return Ring.__mul__(other, self)

        def __eq__(self, other):
            return all(
                self(g) == Ring.natural_project_from(other)(g) for g in G.__iter__()
            )

        # Represent as strings

        def __str__(self):
            return " + ".join(f"{self(g)}{g}" for g in G.__iter__())

        def __repr__(self):
            return f"GroupRing({R}, {G})({self._values})"

        # Class-level manipulation

        @classmethod
        def get_base_ring(cls):
            """Return the underlying ring"""
            return R

        @classmethod
        def get_base_group(cls):
            """Return the underlying group"""
            return G

    return Ring


# CyclicGroupRing isn't needed anymore, consult GroupRing instead
# This class is kept as a reminder of functionality which hasn't been ported over yet
def CyclicGroupRing(n):
    """Returns the cyclic group ring associated with Z_n and the complex numbers"""
    base_group = ModularIntegers(n)

    class GroupRing:
        f"A class representing the group ring of the complex numbers over Z_{n}"

        @classmethod
        def get_characters(cls):
            """Returns the characters of the group ring"""
            return [
                cls.from_function(lambda i: (x * i).to_nth_root_of_unity())
                for x in base_group.__iter__()
            ]

        @classmethod
        def get_primitive_idempotents(cls):
            """Returns the primitive idempotents of the group ring

            Any function in this list, convolved with itself, gives itself back.
            All such functions can be produced by adding some set of functions
            in this list."""
            return [
                cls.from_function(lambda i: (x * i).to_nth_root_of_unity() / n)
                for x in base_group.__iter__()
            ]

        @classmethod
        def from_dft(cls, dft):
            """Get the function with the given DFT"""
            return cls.from_function(
                lambda i: sum(d * c(i) for (c, d) in zip(cls.get_characters(), dft)) / n
            )

        def inner_product(self, other):
            """Take the inner product of the two group ring elements"""
            return sum(self(i) * other(i).conjugate() for i in base_group.__iter__())

        def get_dft(self):
            """Get the Discrete Fourier transform of this group ring element"""
            return [
                complex_zero_approx(self.inner_product(character))
                for character in GroupRing.get_characters()
            ]
