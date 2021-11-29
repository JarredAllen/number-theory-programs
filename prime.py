from math import gcd, isqrt


def eratosthenes(n):
    """Perform the Sieve of Eratosthenes to identify all prime numbers
    up to n.

    Returns all such primes as a list.
    """
    sieve = [False, False] + [True] * (n - 1)
    for i in range(isqrt(n)):
        if not sieve[i]:
            continue
        for j in range(i * i, n + 1, i):
            sieve[j] = False
    return [n for (n, p) in enumerate(sieve) if p]


def is_prime(n):
    """Returns true iff n is prime"""
    n = abs(n)
    if n == 0 or n == 1:
        return False
    for p in eratosthenes(isqrt(n)):
        if n % p == 0:
            return False
    return True


def prime_factor(n):
    """Returns the prime factorization of n"""
    factors = []
    for p in eratosthenes(int(n ** 0.5 + 1)):
        while n % p == 0:
            n //= p
            factors += [p]
    if n > 1:
        factors += [n]
    return factors


def rabin_miller(n, iters=100):
    """Perform `iter` iterations on the Rabin Miller test, returning
    True iff n can still be a prime"""
    if iters > (n + 3) // 4:
        iters = (n + 3) // 4
    q = n - 1
    k = 0
    while q % 2 == 0:
        q //= 2
        k += 1
    for i in range(1, iters + 1):
        if i % n == 0:
            continue
        if not rabin_miller_test_case(n, q, k, i):
            return False
    return True


def rabin_miller_test_case(n, q, k, a):
    """Check the base a in the Rabin-Miller primality test, returning
    True iff a does not prevent $n = 1+q2^k$ from being prime"""
    base = pow(a, q, n)
    if base == 1 or base == n - 1:
        return True
    for _ in range(k):
        base = pow(base, 2, n)
        if base == n - 1:
            return True
        if base == 1:
            return False
    return False


def korselts(n, display=False):
    """Return true iff n passes Korselt's Criterion"""

    def output(s, *args, **kwargs):
        if display:
            print(s.format(*args, **kwargs))

    if n % 2 == 0:
        output("{} is even, so it is not Carmichael.", n)
        return False
    factors = prime_factor(n)
    if len(factors) == 1:
        output("{} is prime, so it is not Carmichael.", n)
        return False
    output("{} has prime factors: ${}$\n", n, factors)
    for factor in factors:
        if n % (factor ** 2) == 0:
            output(
                "${p}$ is a prime factor of ${n}$, and ${p}^2|n$, "
                "so $n$ cannot be Carmichael.",
                p=factor,
                n=n,
            )
            return False
        if (n - 1) % (factor - 1) != 0:
            output(
                "${p}$ is a prime factor of ${n}$, and ${p}-1\\nmid{n}-1$, "
                "so $n$ cannot be Carmichael.",
                p=factor,
                n=n,
            )
            return False
        output(
            "${p}$ is a prime factor of ${n}$, and both ${p}\\parallel{n}$ "
            "and ${p}-1\\mid{n}-1$.\n",
            p=factor,
            n=n,
        )
    output("All prime factors pass, so ${n}$ is Carmichael.", n=n)
    return True


def write_korselts_list(nums):
    """Given a list of numbers, output the LaTeX for finding whether or
    not they're Carmichael"""
    for n in nums:
        print("\\item ", end="")
        korselts(n, True)
