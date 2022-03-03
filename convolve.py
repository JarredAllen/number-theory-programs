from groups import ModularIntegers


def complex_zero_approx(c, epsilon=1e-15):
    """Convert the input to a complex number. If either the real or imaginary parts are less than epsilon, set them to zero"""
    re = complex(c).real
    if abs(re) < epsilon:
        re = 0
    im = complex(c).imag
    if abs(im) < epsilon:
        im = 0
    return complex(re, im)


def CyclicGroupRing(n):
    """Returns the cyclic group ring associated with Z_n and the complex numbers"""
    base_group = ModularIntegers(n)

    class GroupRing:
        """A class representing the group ring of the complex numbers over Z_""" + str(
            n
        )

        def __init__(self, values):
            "Create a new function for the given values. The argument must have length " + str(
                n
            ) + "and represents the values on 1,x,..."
            self._values = [complex_zero_approx(v) for v in values]

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
        def from_function(cls, func):
            """Create a new function according the the function passed as an argument"""
            return cls([func(i) for i in base_group.__iter__()])

        @classmethod
        def natural_project_from(cls, value):
            """Perform the natural projection into this group ring from the subring associated with a subgroup"""
            outputs = [0 for _ in range(n)]
            for input in value.get_base_group().__iter__():
                outputs[base_group.natural_project_from(input).a] = value(input)
            return cls(outputs)

        @classmethod
        def get_base_group(cls):
            """Get the base group that this group ring is over"""
            return base_group

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

        def __call__(self, arg):
            """Evaluate the given character at the given input"""
            return self._values[arg.a]

        def __add__(self, other):
            """Add the two functions"""
            return GroupRing(self._values + other._values)

        def __mul__(self, other):
            """Convolve the two functions"""
            return GroupRing(
                [
                    sum(
                        self._values[(x + i) % n] * other._values[(-i) % n]
                        for i in range(n)
                    )
                    for x in range(n)
                ]
            )

        def __str__(self):
            return str(self._values)

        def __repr__(self):
            return f"CyclicGroupRing({n})({self._values})"

    return GroupRing
