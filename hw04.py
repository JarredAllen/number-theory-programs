from functools import reduce


def extended_euclidean(a, b):
    """Given integers a,b, find integers x,y such that a*x+b*y = gcd(a,b)"""
    coeffs = [(0, 1), (1, 0)]
    while b != 0:
        q, r = divmod(a, b)
        x = coeffs[-2][1] - q * coeffs[-1][1]
        y = coeffs[-2][0] - q * coeffs[-1][0]
        a, b = b, r
        coeffs += [(y, x)]
    return coeffs[-2][::-1]


def crt(b, m, c, n):
    """Find an integer x such that:
    - $$x \\equiv b (mod m)$$
    - $$x \\equiv c (mod n)$$.
    """
    x, y = extended_euclidean(m, n)
    return x * m * c + y * n * b


def euclid_list(p, maximum=1000):
    """Use the algorithm to produce a list of prime numbers, terminating
    when the produced number has no factors below the given maximum and
    is greater than the maximum.
    """
    while True:
        a = reduce(lambda x, y: x * y, p) + 1
        for i in range(2, maximum):
            if a % i == 0:
                p += [i]
                break
        else:
            return p


def pow(a, k, m=None):
    """Compute $$a^k$$ by the method of successive squaring"""
    if m is None:
        b = 1
        while k >= 1:
            if k & 1:
                b *= a
            a **= 2
            k //= 2
        return b
    else:
        b = 1
        while k >= 1:
            if k & 1:
                b = (b * a) % m
            a = (a * a) % m
            k //= 2
        return b
