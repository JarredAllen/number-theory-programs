from fractions import Fraction


class EllipticCurve:
    """A class which implements elliptic curves"""

    def __init__(self, a, b, c):
        self.a = a
        self.b = b
        self.c = c

        class Point:
            def __init__(self, x, y):
                """Create a new point at (x, y). If the point at infinity
                is wanted, pass x=None, y=None
                """
                if x is not None:
                    assert y is not None
                    assert y ** 2 == x ** 3 + a * x ** 2 + b * x + c
                # Convert to fractions on input to avoid floating-point errors
                self.x = Fraction(x)
                self.y = Fraction(y)

            def __eq__(self, other):
                return self.x == other.x and self.y == other.y

            def __neq__(self, other):
                return not self.__eq__(other)

            def __neg__(self):
                return Point(self.x, -self.y)

            def __add__(self, other):
                # Identity point: get other point
                if self.x is None and self.y is None:
                    return other
                elif other.x is None and other.y is None:
                    return self
                # Vertical lines: Get the identity point
                elif self.x == other.x and self.y != other.y:
                    return Point(None, None)
                elif self.x == other.x and self.y == 0 and other.y == 0:
                    return Point(None, None)
                # Otherwise, get tangent or secant line, and use it to
                # find the third point of intersection on the line
                else:
                    if self.x == other.x:
                        m = (3 * self.x ** 2 + 2 * a * self.x + b) / (2 * self.y)
                    else:
                        m = (self.y - other.y) / (self.x - other.x)
                    line_b = -m * self.x + self.y
                    # (mx + ...)^2 = x^3 + ax^2 + ...
                    # So we re-arrange and find:
                    # 0 = x^3 - (m^2 - a)x^2 + ...
                    # And so the three x values sum to m^2-a
                    x = m ** 2 - a - self.x - other.x
                    y = m * x + line_b
                    return Point(x, -y)

            def __repr__(self):
                return "<point>(%s, %s)" % (repr(self.x), repr(self.y))

            def __str__(self):
                return "(%s, %s)" % (self.x, self.y)

        self.point = Point

    def count_solutions_mod_p(self, p):
        """Count the number of solutions to this elliptic curve modulo p"""
        from jacobi import jacobi

        count = 0
        for x in range(p):
            l = jacobi(x ** 3 + self.a * x ** 2 + self.b * x + self.c, p)
            count += l + 1
        return count


def output_nonzero_p_defects(curve, stop):
    from prime import eratosthenes

    primes = eratosthenes(stop)
    print("\\begin{ttabular}")
    for p in primes:
        defect = curve.count_solutions_mod_p(p) - p
        if defect != 0:
            print(f"\\bfseries {p} & {defect} \\\\")
    print("\\end{ttabular}")
